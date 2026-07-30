#!/usr/bin/env python
# coding: utf-8
from __future__ import annotations

# Limit BLAS/OpenMP threads inside each process to avoid outer x inner oversubscription.
# Set these before numpy/scipy/sklearn imports when possible.
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")


# # Run scCST on parameter-sweep replicated simulations
# 
# This notebook is adapted from `sccst_sim_multiple_v3.ipynb` and is configured for the datasets created by `create_sim_parameter_sweeps_replicates.ipynb`.
# 
# Expected simulation directory:
# 
# ```python
# /home/anilprakash/labs/Mei/projects/anil/srda/notebooks/data/scrna_seq/simulation
# ```
# 
# It reads `simulation_manifest_de_parameter_sweeps.csv` when available, runs scCST on each generated `.h5ad`, and writes one result object plus per-cell summaries for each dataset and neighborhood-size setting. The output filenames include `sweep_name` and `replicate` so the repeated default-parameter datasets from different sweeps do not overwrite each other.
# 

# In[3]:


#sccst performance on simulated data.

import scanpy as sc
import scipy.sparse as sp
import time
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad
import scipy.io
import os
from scipy.sparse import coo_matrix, csr_matrix
from pathlib import Path


# In[1]:


"""
Latent-neighborhood matched biological-replicate partial-list RRA
with optional evenly spaced anchor selection.

Main behavior
-------------
1. If n_anchors_per_celltype is None and anchor_indices is None:
       Run focal-cell RRA over range(adata.n_obs), i.e. every cell.

2. If n_anchors_per_celltype is an int/dict or anchor_indices is provided:
       Run RRA only on selected anchor cells.

"""


from typing import Any, Dict, Optional, Tuple, Union, List

import numpy as np
import scanpy as sc
import scipy.stats as stats

from scipy.sparse import issparse
from joblib import Parallel, delayed
from scipy.optimize import linear_sum_assignment
from sklearn.neighbors import NearestNeighbors


# =============================================================================
# Multiple testing
# =============================================================================

def apply_bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """
    Applies Benjamini-Hochberg FDR correction.
    """
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)

    if n == 0:
        return p_values

    idx = np.argsort(p_values)
    sorted_p = p_values[idx]

    q_vals = sorted_p * (n / np.arange(1, n + 1))
    q_vals = np.minimum.accumulate(q_vals[::-1])[::-1]

    adj_p_values = np.zeros_like(q_vals)
    adj_p_values[idx] = np.minimum(q_vals, 1.0)

    return adj_p_values


# =============================================================================
# Bonferroni-corrected RRA rho p-values
# =============================================================================

def rra_bonferroni_pvalues(normalized_ranks: np.ndarray) -> np.ndarray:
    """
    Computes Bonferroni-corrected RRA rho p-values.

    Parameters
    ----------
    normalized_ranks:
        Array of shape (n_lists, n_genes), values in [0, 1].
        Smaller values indicate stronger ranking.

    Returns
    -------
    p_rra:
        Bonferroni-corrected RRA rho p-values for each gene.

    Notes
    -----
    For each gene, RRA computes beta order-statistic p-values:
        beta_{m,k}(r_(k))
    where m is the number of ranked lists and r_(k) is the kth smallest
    normalized rank.

    The rho score is:
        rho = min_k beta_{m,k}(r_(k))

    This function converts rho to a conservative p-value using:
        p = min(rho * m, 1)

    This is not the exact RobustRankAggreg p-value. It is a conservative
    RRA-style rho-score p-value.
    """
    normalized_ranks = np.asarray(normalized_ranks, dtype=float)

    if normalized_ranks.ndim != 2:
        raise ValueError("normalized_ranks must have shape (n_lists, n_genes).")

    n_lists, _ = normalized_ranks.shape

    if n_lists < 1:
        raise ValueError("RRA requires at least one ranked list.")

    ranks = np.clip(normalized_ranks, 0.0, 1.0)
    sorted_ranks = np.sort(ranks, axis=0)

    beta_scores = np.zeros_like(sorted_ranks)

    for k in range(1, n_lists + 1):
        beta_scores[k - 1, :] = stats.beta.cdf(
            sorted_ranks[k - 1, :],
            k,
            n_lists - k + 1,
        )

    rho = np.min(beta_scores, axis=0)
    p_rra = np.minimum(rho * n_lists, 1.0)

    return p_rra


# =============================================================================
# Local latent calipers
# =============================================================================

def compute_celltype_local_calipers(
    latent_matrix: np.ndarray,
    y_ctypes: np.ndarray,
    k_ref: int = 20,
    percentile: float = 75.0,
    min_cells: int = 30,
) -> Dict[int, float]:
    """
    Computes a cell-type-specific local latent-distance caliper.

    For each cell type:
        1. Find each cell's k_ref nearest same-cell-type neighbors.
        2. Compute the mean distance to those neighbors.
        3. Use the requested percentile of those mean local distances as the caliper.

    Returns
    -------
    calipers:
        Dict mapping numeric cell-type id -> latent-distance caliper.
    """
    calipers: Dict[int, float] = {}

    for ctype in np.unique(y_ctypes):
        idx = np.where(y_ctypes == ctype)[0]

        if len(idx) < max(min_cells, k_ref + 1):
            continue

        X = latent_matrix[idx]

        nn = NearestNeighbors(n_neighbors=k_ref + 1, metric="euclidean")
        nn.fit(X)

        distances, _ = nn.kneighbors(X)
        neighbor_distances = distances[:, 1:]
        mean_local_distance = neighbor_distances.mean(axis=1)

        calipers[int(ctype)] = float(np.percentile(mean_local_distance, percentile))

    return calipers


# =============================================================================
# Anchor selection
# =============================================================================

def select_evenly_spaced_anchors(
    latent_matrix: np.ndarray,
    y_ctypes: np.ndarray,
    n_anchors_per_celltype: Union[int, Dict[Union[int, str], int]],
    ctype_map: Optional[Dict[Any, int]] = None,
    random_state: int = 0,
    min_cells_per_celltype: int = 30,
) -> np.ndarray:
    """
    Selects evenly spaced anchors within each cell type using farthest-point
    sampling in latent space.

    Parameters
    ----------
    latent_matrix:
        Latent embedding, shape (n_cells, n_latent_dims).

    y_ctypes:
        Numeric cell-type labels, shape (n_cells,).

    n_anchors_per_celltype:
        Either:
            - int: same number of anchors per cell type.
            - dict: mapping numeric cell-type id or original cell-type label
              to number of anchors.

    ctype_map:
        Optional mapping from original cell-type label -> numeric id.
        Needed only if n_anchors_per_celltype uses original cell-type labels.

    random_state:
        Random seed for the first anchor in each cell type.

    min_cells_per_celltype:
        Skip cell types with fewer than this many cells.

    Returns
    -------
    anchor_indices:
        Global cell indices selected as anchors.
    """
    rng = np.random.default_rng(random_state)
    anchor_indices: List[int] = []

    # Convert label-keyed dict to numeric-keyed dict if needed.
    anchor_counts: Optional[Dict[int, int]] = None

    if isinstance(n_anchors_per_celltype, dict):
        anchor_counts = {}
        for key, value in n_anchors_per_celltype.items():
            if isinstance(key, (int, np.integer)):
                anchor_counts[int(key)] = int(value)
            else:
                if ctype_map is None or key not in ctype_map:
                    raise ValueError(
                        "n_anchors_per_celltype contains cell-type labels, but "
                        "ctype_map is missing or does not contain one of the labels."
                    )
                anchor_counts[int(ctype_map[key])] = int(value)

    for ctype in np.unique(y_ctypes):
        ctype = int(ctype)
        idx = np.where(y_ctypes == ctype)[0]

        if len(idx) < min_cells_per_celltype:
            continue

        if anchor_counts is None:
            n_anchors = int(n_anchors_per_celltype)
        else:
            n_anchors = int(anchor_counts.get(ctype, 0))

        if n_anchors <= 0:
            continue

        n_anchors = min(n_anchors, len(idx))
        X = latent_matrix[idx]

        # First anchor chosen randomly within cell type.
        first_local = int(rng.integers(0, len(idx)))
        selected_local = [first_local]

        # Squared distance to nearest selected anchor.
        min_dist = np.sum((X - X[first_local]) ** 2, axis=1)

        for _ in range(1, n_anchors):
            next_local = int(np.argmax(min_dist))
            selected_local.append(next_local)

            new_dist = np.sum((X - X[next_local]) ** 2, axis=1)
            min_dist = np.minimum(min_dist, new_dist)

        anchor_indices.extend(idx[selected_local].tolist())

    return np.array(anchor_indices, dtype=int)


# =============================================================================
# Partial-list construction
# =============================================================================

def _build_partial_directional_ranks(
    diff: np.ndarray,
    direction: str,
    top_k: int,
    n_genes_universe: int,
) -> np.ndarray:
    """
    Builds one partial RRA ranked list.

    Missing / omitted genes are assigned rank 1.0.

    Parameters
    ----------
    diff:
        Disease-control expression difference for one matched biological pair.

    direction:
        "up" or "down".
        up ranks genes by largest positive disease-control difference.
        down ranks genes by most negative disease-control difference.

    top_k:
        Number of genes retained in the partial list.

    n_genes_universe:
        Total number of genes in the tested universe after filtering.

    Returns
    -------
    ranks:
        Normalized ranks in [0, 1], shape (n_genes,).
        Smaller values indicate stronger rank.
    """
    n_genes = len(diff)
    ranks = np.ones(n_genes, dtype=float)

    if direction == "up":
        candidate_genes = np.where(diff > 0)[0]
        if len(candidate_genes) == 0:
            return ranks
        ordered_genes = candidate_genes[np.argsort(-diff[candidate_genes])]

    elif direction == "down":
        candidate_genes = np.where(diff < 0)[0]
        if len(candidate_genes) == 0:
            return ranks
        ordered_genes = candidate_genes[np.argsort(diff[candidate_genes])]

    else:
        raise ValueError("direction must be 'up' or 'down'.")

    selected = ordered_genes[: min(top_k, len(ordered_genes))]

    normalized_positions = (
        np.arange(1, len(selected) + 1, dtype=float) /
        float(n_genes_universe)
    )

    ranks[selected] = normalized_positions

    return ranks


# =============================================================================
# Per-focal-cell / per-anchor method
# =============================================================================

def process_focal_cell_rra(
    focal_idx: int,
    focal_ctype_idx: int,
    X_matrix: Any,
    y_ctypes: np.ndarray,
    y_batches: np.ndarray,
    latent_matrix: np.ndarray,
    k: int,
    batch_conditions_map: Dict[int, int],
    batch_sample_map: Union[Dict[int, Any], None] = None,
    celltype_calipers: Union[Dict[int, float], None] = None,
    min_cells_per_batch_neighborhood: int = 3,
    min_matched_pairs: int = 3,
    top_genes_per_list: Optional[int] = 100,
    fdr_threshold: float = 0.05,
) -> Union[Dict[int, Dict[str, Any]], None]:
    """
    Processes one focal anchor cell using local matched partial-list RRA.

    Main inferential structure
    --------------------------
    For each focal cell / anchor:
        1. Build batch-specific same-cell-type local neighborhoods.
        2. Match control and disease biological-replicate neighborhoods.
        3. For each matched biological pair, compute disease-control mean
           expression difference.
        4. Build partial ranked lists per biological pair:
               up-list   = top K positive differences
               down-list = top K negative differences
        5. Run Bonferroni-corrected RRA separately for up and down.
        6. Compute gene-level FDR within direction.
        7. Gene score:
               avg_expr_diff * -log10(direction-specific FDR)

    Notes
    -----
    top_genes_per_list controls the number of retained genes per ranked list.
    If top_genes_per_list is None, all genes are allowed into each directional
    list, subject to sign.
    """
    focal_latent = latent_matrix[focal_idx]
    unique_batches = np.unique(y_batches)

    batch_neighborhoods: Dict[int, np.ndarray] = {}
    batch_centroids: Dict[int, np.ndarray] = {}

    # ---------------------------------------------------------------------
    # SECTION 1: BATCH-FIRST ANCHOR NEIGHBORHOODS
    # ---------------------------------------------------------------------
    for b in unique_batches:
        b = int(b)
        mask = (y_batches == b) & (y_ctypes == focal_ctype_idx)

        # Exclude the focal cell if it is in this batch/cell-type subset.
        mask[focal_idx] = False

        b_indices = np.where(mask)[0]

        if len(b_indices) < min_cells_per_batch_neighborhood:
            continue

        b_latent = latent_matrix[b_indices]
        dists_sq = np.sum((b_latent - focal_latent) ** 2, axis=1)

        top_k_cells = min(k, len(b_indices))
        top_k_idx = np.argsort(dists_sq)[:top_k_cells]
        selected_cells = b_indices[top_k_idx]

        if len(selected_cells) < min_cells_per_batch_neighborhood:
            continue

        batch_neighborhoods[b] = selected_cells
        batch_centroids[b] = np.mean(latent_matrix[selected_cells], axis=0)

    control_batches = [
        b for b in batch_neighborhoods
        if batch_conditions_map[b] == 0
    ]

    disease_batches = [
        b for b in batch_neighborhoods
        if batch_conditions_map[b] == 1
    ]

    if not control_batches or not disease_batches:
        return None

    # ---------------------------------------------------------------------
    # SECTION 2: BATCH PAIRING + CALIPER FILTERING
    # ---------------------------------------------------------------------
    paired_c_batches: List[int] = []
    paired_d_batches: List[int] = []
    paired_distances: List[float] = []

    caliper = None
    if celltype_calipers is not None:
        caliper = celltype_calipers.get(int(focal_ctype_idx), None)

    if batch_sample_map is not None:
        # Strict paired design is validated globally in process_cell_type.
        # Here we keep only pairs where both batches have valid focal-cell
        # neighborhoods and pass the latent caliper.
        control_by_sample = {
            batch_sample_map[b]: b
            for b in control_batches
        }

        disease_by_sample = {
            batch_sample_map[b]: b
            for b in disease_batches
        }

        shared_samples = sorted(
            set(control_by_sample.keys()) & set(disease_by_sample.keys())
        )

        for sample_id in shared_samples:
            c_b = int(control_by_sample[sample_id])
            d_b = int(disease_by_sample[sample_id])

            pair_dist = float(np.linalg.norm(
                batch_centroids[c_b] - batch_centroids[d_b]
            ))

            if caliper is not None and pair_dist > caliper:
                continue

            paired_c_batches.append(c_b)
            paired_d_batches.append(d_b)
            paired_distances.append(pair_dist)

    else:
        # Distance-based one-to-one optimal matching between available
        # control and disease biological replicates.
        dist_matrix = np.zeros((len(control_batches), len(disease_batches)))

        for i, c_b in enumerate(control_batches):
            for j, d_b in enumerate(disease_batches):
                dist_matrix[i, j] = np.linalg.norm(
                    batch_centroids[c_b] - batch_centroids[d_b]
                )

        row_ind, col_ind = linear_sum_assignment(dist_matrix)

        for i, j in zip(row_ind, col_ind):
            c_b = int(control_batches[i])
            d_b = int(disease_batches[j])
            pair_dist = float(dist_matrix[i, j])

            if caliper is not None and pair_dist > caliper:
                continue

            paired_c_batches.append(c_b)
            paired_d_batches.append(d_b)
            paired_distances.append(pair_dist)

    if len(paired_c_batches) < min_matched_pairs:
        return None

    # ---------------------------------------------------------------------
    # SECTION 3: PAIRWISE EXPRESSION DIFFERENCES + PARTIAL RANKED LISTS
    # ---------------------------------------------------------------------
    n_genes = X_matrix.shape[1]
    n_pairs = len(paired_c_batches)

    if top_genes_per_list is None:
        top_k_genes = n_genes
    else:
        top_k_genes = int(top_genes_per_list)

    top_k_genes = max(1, min(top_k_genes, n_genes))

    ranks_up = np.ones((n_pairs, n_genes), dtype=float)
    ranks_down = np.ones((n_pairs, n_genes), dtype=float)
    mean_expr_diffs = np.zeros((n_pairs, n_genes), dtype=float)

    for p_idx, (c_b, d_b) in enumerate(zip(paired_c_batches, paired_d_batches)):
        c_cells = batch_neighborhoods[c_b]
        d_cells = batch_neighborhoods[d_b]

        Xc = X_matrix[c_cells]
        Xd = X_matrix[d_cells]

        if issparse(Xd):
            mean_d = Xd.mean(axis=0).A1
        else:
            mean_d = np.mean(Xd, axis=0)

        if issparse(Xc):
            mean_c = Xc.mean(axis=0).A1
        else:
            mean_c = np.mean(Xc, axis=0)

        diff = mean_d - mean_c
        mean_expr_diffs[p_idx, :] = diff

        ranks_up[p_idx, :] = _build_partial_directional_ranks(
            diff=diff,
            direction="up",
            top_k=top_k_genes,
            n_genes_universe=n_genes,
        )

        ranks_down[p_idx, :] = _build_partial_directional_ranks(
            diff=diff,
            direction="down",
            top_k=top_k_genes,
            n_genes_universe=n_genes,
        )

    # ---------------------------------------------------------------------
    # SECTION 4: DIRECTIONAL RRA P-VALUES AND FDR
    # ---------------------------------------------------------------------
    p_up = rra_bonferroni_pvalues(ranks_up)
    p_down = rra_bonferroni_pvalues(ranks_down)

    fdr_up = apply_bh_fdr(p_up)
    fdr_down = apply_bh_fdr(p_down)

    sig_up = fdr_up < fdr_threshold
    sig_down = fdr_down < fdr_threshold
    sig_mask = sig_up | sig_down

    # Optional two-sided-style minimum directional p-value, for convenience.
    p_combined = np.minimum(2.0 * np.minimum(p_up, p_down), 1.0)
    fdr_combined = apply_bh_fdr(p_combined)

    # ---------------------------------------------------------------------
    # SECTION 5: EFFECT SIZE, DIRECTION CONSISTENCY, AND GENE SCORE
    # ---------------------------------------------------------------------
    avg_diff = np.mean(mean_expr_diffs, axis=0)

    avg_sign = np.sign(avg_diff)
    pair_signs = np.sign(mean_expr_diffs)

    direction_consistency = np.mean(
        pair_signs == avg_sign[None, :],
        axis=0,
    )
    direction_consistency[avg_sign == 0] = 0.0

    directional_fdr = np.where(avg_diff >= 0, fdr_up, fdr_down)
    directional_p = np.where(avg_diff >= 0, p_up, p_down)

    gene_score = avg_diff * -np.log10(directional_fdr + 1e-300)

    cell_divergence_score = np.sum(np.abs(gene_score))
    n_sig_genes = int(np.sum(sig_mask))
    n_sig_up = int(np.sum(sig_up))
    n_sig_down = int(np.sum(sig_down))

    return {
        int(focal_idx): {
            "p_up": p_up,
            "p_down": p_down,
            "fdr_up": fdr_up,
            "fdr_down": fdr_down,
            "p_combined": p_combined,
            "fdr_combined": fdr_combined,
            "directional_p": directional_p,
            "directional_fdr": directional_fdr,
            "avg_expr_diff": avg_diff,
            "direction_consistency": direction_consistency,
            "gene_score": gene_score,
            "cell_divergence_score": cell_divergence_score,
            "n_sig_genes": n_sig_genes,
            "n_sig_up": n_sig_up,
            "n_sig_down": n_sig_down,
            "n_pairs": int(n_pairs),
            "paired_control_batches": np.array(paired_c_batches, dtype=int),
            "paired_disease_batches": np.array(paired_d_batches, dtype=int),
            "paired_distances": np.array(paired_distances, dtype=float),
            "top_genes_per_list_used": int(top_k_genes),
        }
    }



# =============================================================================
# Main driver
# =============================================================================

def process_cell_type(
    adata: sc.AnnData,
    conditions: Tuple[str, str],
    condition_col: str = "condition",
    batch_col: str = "batch",
    sample_col: Union[str, None] = None,
    n_neighborhoods: int = 20,
    cell_type_col: str = "Cell_Type",
    n_jobs: int = 4,
    latent_key: str = "X_pca_harmony",
    caliper_percentile: float = 75.0,
    min_cells_per_batch_neighborhood: int = 3,
    min_matched_pairs: int = 3,
    top_genes_per_list: Optional[int] = 100,
    fdr_threshold: float = 0.05,
    n_anchors_per_celltype: Optional[Union[int, Dict[Union[int, str], int]]] = None,
    anchor_indices: Optional[np.ndarray] = None,
    anchor_random_state: int = 0,
) -> Tuple[Dict[str, Dict[str, Any]], np.ndarray, Dict[str, Any]]:
    """
    Main processing function.

    Implements latent-neighborhood matched biological-replicate partial-list RRA
    for local condition-associated expression shifts.

    Run modes
    ---------
    1. All-cell mode:
       Set n_anchors_per_celltype=None and anchor_indices=None.
       The model runs over range(adata.n_obs).

    2. Anchor mode by automatic anchor selection:
       Set n_anchors_per_celltype to an int or dict.
       The model selects evenly spaced anchors within each cell type and runs
       RRA only on those anchors.

    3. Anchor mode by manual anchors:
       Provide anchor_indices.
       The model runs RRA only on those anchors.

    Anchor-mode results are only for the tested anchors.

    Parameters
    ----------
    adata:
        AnnData object.

    conditions:
        Tuple of condition labels:
            (disease_label, control_label)

    sample_col:
        If provided, strict paired design is enforced:
            each sample must contain exactly one control batch and one disease batch.

    top_genes_per_list:
        Number of top genes retained in each pair-level partial ranked list.
        If None, all genes are used subject to directional sign.

    n_anchors_per_celltype:
        None:
            run over every cell.
        int:
            select that many evenly spaced anchors per cell type.
        dict:
            mapping numeric cell-type id or original cell-type label -> anchor count.

    anchor_indices:
        Optional manual list/array of global cell indices to use as anchors.
        If provided, this triggers anchor mode regardless of n_anchors_per_celltype.

    Returns
    -------
    cell_importances:
        Dict mapping tested focal cell / anchor IDs to result dictionaries.

    gene_names:
        Gene names after filtering.

    run_info:
        Dict containing mode information, selected anchor indices, successful
        focal indices, mappings, and calipers.
    """

    # ---------------------------------------------------------------------
    # STEP 0: BASIC VALIDATION
    # ---------------------------------------------------------------------
    if latent_key not in adata.obsm:
        raise ValueError(f"Latent key '{latent_key}' not found in adata.obsm.")

    if batch_col not in adata.obs.columns:
        raise ValueError(f"Batch column '{batch_col}' not found in adata.obs.")

    if condition_col not in adata.obs.columns:
        raise ValueError(f"Condition column '{condition_col}' not found in adata.obs.")

    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"Cell-type column '{cell_type_col}' not found in adata.obs.")

    if top_genes_per_list is not None and int(top_genes_per_list) < 1:
        raise ValueError("top_genes_per_list must be None or a positive integer.")

    # ---------------------------------------------------------------------
    # STEP 1: GENE FILTERING
    # ---------------------------------------------------------------------
    min_cells = int(adata.n_obs * 0.005)
    min_cells = max(min_cells, 3)

    # Modifies adata in place. Use adata.copy() before calling if needed.
    sc.pp.filter_genes(adata, min_cells=min_cells)
    gene_names = adata.var_names.to_numpy()

    print(f"Processing {adata.n_obs} cells and {adata.n_vars} genes.")
    print(f"Using latent key: {latent_key}")
    print(f"Using top_genes_per_list: {top_genes_per_list}")
    print(f"Using {caliper_percentile}th percentile cell-type-specific caliper.")

    # ---------------------------------------------------------------------
    # STEP 2: PREPARE MATRICES
    # ---------------------------------------------------------------------
    latent_matrix = adata.obsm[latent_key]

    if latent_matrix.shape[0] != adata.n_obs:
        raise ValueError(
            f"Latent matrix '{latent_key}' has {latent_matrix.shape[0]} rows, "
            f"but adata has {adata.n_obs} cells."
        )

    if issparse(adata.X):
        input_X = adata.X.tocsr()
    else:
        input_X = adata.X.copy()

    # ---------------------------------------------------------------------
    # STEP 3: ENCODE CELL TYPES
    # ---------------------------------------------------------------------
    unique_ctypes = adata.obs[cell_type_col].unique()
    ctype_map = {ct: i for i, ct in enumerate(unique_ctypes)}
    reverse_ctype_map = {v: k for k, v in ctype_map.items()}
    y_ctypes = adata.obs[cell_type_col].map(ctype_map).to_numpy()

    # ---------------------------------------------------------------------
    # STEP 4: ENCODE CONDITIONS
    # ---------------------------------------------------------------------
    disease_label, control_label = conditions

    y_map = {
        control_label: 0,
        disease_label: 1,
    }

    observed_conditions = set(adata.obs[condition_col].unique())
    expected_conditions = set(y_map.keys())

    if not observed_conditions.issubset(expected_conditions):
        raise ValueError(
            f"Observed condition labels {observed_conditions} are not fully "
            f"contained in expected labels {expected_conditions}."
        )

    y_conditions = adata.obs[condition_col].map(y_map).to_numpy()

    # ---------------------------------------------------------------------
    # STEP 5: ENCODE BATCHES
    # ---------------------------------------------------------------------
    unique_batches = adata.obs[batch_col].unique()
    batch_map = {b: i for i, b in enumerate(unique_batches)}
    reverse_batch_map = {v: k for k, v in batch_map.items()}
    y_batches = adata.obs[batch_col].map(batch_map).to_numpy()

    batch_conditions_map: Dict[int, int] = {}

    for b_id in np.unique(y_batches):
        b_id = int(b_id)
        conds_in_batch = np.unique(y_conditions[y_batches == b_id])

        if len(conds_in_batch) != 1:
            original_batch = unique_batches[b_id]
            raise ValueError(
                f"Batch '{original_batch}' contains multiple conditions: "
                f"{conds_in_batch}. Each batch must belong to one condition only."
            )

        batch_conditions_map[b_id] = int(conds_in_batch[0])

    # ---------------------------------------------------------------------
    # STEP 6: OPTIONAL STRICT SAMPLE-PAIRED DESIGN
    # ---------------------------------------------------------------------
    batch_sample_map = None

    if sample_col is not None:
        if sample_col not in adata.obs.columns:
            raise ValueError(f"Sample column '{sample_col}' not found in adata.obs.")

        y_samples = adata.obs[sample_col].to_numpy()
        batch_sample_map = {}

        for b_id in np.unique(y_batches):
            b_id = int(b_id)
            samples_in_batch = np.unique(y_samples[y_batches == b_id])

            if len(samples_in_batch) != 1:
                original_batch = unique_batches[b_id]
                raise ValueError(
                    f"Batch '{original_batch}' spans multiple samples. "
                    f"Ensure batches are nested properly inside samples."
                )

            batch_sample_map[b_id] = samples_in_batch[0]

        sample_to_control_batches: Dict[Any, List[int]] = {}
        sample_to_disease_batches: Dict[Any, List[int]] = {}

        for b_id, sample_id in batch_sample_map.items():
            cond = batch_conditions_map[b_id]

            if cond == 0:
                sample_to_control_batches.setdefault(sample_id, []).append(b_id)
            elif cond == 1:
                sample_to_disease_batches.setdefault(sample_id, []).append(b_id)
            else:
                raise ValueError(
                    f"Unexpected condition code {cond} for batch id {b_id}."
                )

        all_samples = (
            set(sample_to_control_batches.keys()) |
            set(sample_to_disease_batches.keys())
        )

        invalid_samples = []

        for sample_id in all_samples:
            n_control = len(sample_to_control_batches.get(sample_id, []))
            n_disease = len(sample_to_disease_batches.get(sample_id, []))

            if n_control != 1 or n_disease != 1:
                invalid_samples.append((sample_id, n_control, n_disease))

        if invalid_samples:
            msg = "\n".join(
                [
                    f"  sample={sample_id}: "
                    f"{n_control} control batch(es), "
                    f"{n_disease} disease batch(es)"
                    for sample_id, n_control, n_disease in invalid_samples
                ]
            )

            raise ValueError(
                "Strict paired design violated. Each sample must contain "
                "exactly one control batch and exactly one disease batch.\n"
                f"Invalid samples:\n{msg}"
            )

        print(
            "Using strict one-to-one sample-paired batch matching "
            "with latent-distance caliper."
        )

    else:
        print(
            "Using distance-based bipartite batch matching "
            "with latent-distance caliper."
        )

    # ---------------------------------------------------------------------
    # STEP 7: COMPUTE CELL-TYPE-SPECIFIC CALIPERS
    # ---------------------------------------------------------------------
    celltype_calipers = compute_celltype_local_calipers(
        latent_matrix=latent_matrix,
        y_ctypes=y_ctypes,
        k_ref=n_neighborhoods,
        percentile=caliper_percentile,
        min_cells=max(30, n_neighborhoods + 1),
    )

    if len(celltype_calipers) == 0:
        raise ValueError(
            "No cell-type calipers could be computed. "
            "Check whether cell types have enough cells."
        )

    print("Computed cell-type-specific calipers:")
    for ctype_id, caliper in celltype_calipers.items():
        print(f"  {reverse_ctype_map[ctype_id]}: {caliper:.4f}")

    # ---------------------------------------------------------------------
    # STEP 8: CHOOSE ALL-CELL MODE OR ANCHOR MODE
    # ---------------------------------------------------------------------
    anchor_mode = anchor_indices is not None or n_anchors_per_celltype is not None

    if anchor_mode:
        if anchor_indices is None:
            anchor_indices = select_evenly_spaced_anchors(
                latent_matrix=latent_matrix,
                y_ctypes=y_ctypes,
                n_anchors_per_celltype=n_anchors_per_celltype,  # type: ignore[arg-type]
                ctype_map=ctype_map,
                random_state=anchor_random_state,
                min_cells_per_celltype=max(30, n_neighborhoods + 1),
            )
        else:
            anchor_indices = np.asarray(anchor_indices, dtype=int)

        if len(anchor_indices) == 0:
            raise ValueError(
                "No anchors were selected. Increase n_anchors_per_celltype or "
                "check cell-type sizes."
            )

        if np.any(anchor_indices < 0) or np.any(anchor_indices >= adata.n_obs):
            raise ValueError("anchor_indices contains values outside [0, adata.n_obs).")

        focal_indices_to_run = np.asarray(anchor_indices, dtype=int)

        print(
            f"Running focal-cell RRA on {len(focal_indices_to_run)} selected "
            f"anchors instead of all {adata.n_obs} cells."
        )

    else:
        focal_indices_to_run = np.arange(adata.n_obs, dtype=int)
        anchor_indices = None

        print(f"Running focal-cell RRA over range(adata.n_obs) = {adata.n_obs} cells.")

    # ---------------------------------------------------------------------
    # STEP 9: RUN FOCAL-CELL / ANCHOR RRA
    # ---------------------------------------------------------------------
    print("Starting parallel focal-cell partial-list RRA...")

    results = Parallel(n_jobs=n_jobs, verbose=50)(
        delayed(process_focal_cell_rra)(
            focal_idx=int(i),
            focal_ctype_idx=int(y_ctypes[int(i)]),
            X_matrix=input_X,
            y_ctypes=y_ctypes,
            y_batches=y_batches,
            latent_matrix=latent_matrix,
            k=n_neighborhoods,
            batch_conditions_map=batch_conditions_map,
            batch_sample_map=batch_sample_map,
            celltype_calipers=celltype_calipers,
            min_cells_per_batch_neighborhood=min_cells_per_batch_neighborhood,
            min_matched_pairs=min_matched_pairs,
            top_genes_per_list=top_genes_per_list,
            fdr_threshold=fdr_threshold,
        )
        for i in focal_indices_to_run
    )

    # ---------------------------------------------------------------------
    # STEP 10: COLLECT RESULTS
    # ---------------------------------------------------------------------
    cell_importances: Dict[str, Dict[str, Any]] = {}
    focal_ids = adata.obs.index

    n_success = 0
    n_failed = 0
    successful_focal_indices: List[int] = []

    for res in results:
        if res is not None:
            idx = int(list(res.keys())[0])
            cell_id = str(focal_ids[idx])
            cell_importances[cell_id] = res[idx]
            cell_importances[cell_id]["focal_idx"] = idx
            cell_importances[cell_id]["focal_id"] = cell_id
            cell_importances[cell_id]["is_anchor_mode"] = bool(anchor_mode)
            successful_focal_indices.append(idx)
            n_success += 1
        else:
            n_failed += 1

    print("Finished focal-cell partial-list RRA.")
    print(f"Successful focal cells/anchors: {n_success}")
    print(f"Skipped focal cells/anchors: {n_failed}")

    run_info: Dict[str, Any] = {
        "anchor_mode": bool(anchor_mode),
        "n_anchors_per_celltype": n_anchors_per_celltype,
        "anchor_indices": None if anchor_indices is None else np.asarray(anchor_indices, dtype=int),
        "successful_focal_indices": np.asarray(successful_focal_indices, dtype=int),
        "ctype_map": ctype_map,
        "reverse_ctype_map": reverse_ctype_map,
        "batch_map": batch_map,
        "reverse_batch_map": reverse_batch_map,
        "batch_conditions_map": batch_conditions_map,
        "celltype_calipers": celltype_calipers,
        "gene_names": gene_names,
    }

    return cell_importances, gene_names, run_info


# =============================================================================
# AnnData integration helpers
# =============================================================================

def add_rra_results_to_adata(
    adata: sc.AnnData,
    cell_importances: Dict[str, Dict[str, Any]],
    prefix: str = "rra",
    gene_level_keys: Optional[List[str]] = None,
    obs_level_keys: Optional[List[str]] = None,
    overwrite: bool = True,
) -> sc.AnnData:
    """
    Adds RRA results to an AnnData object by matching result keys to
    adata.obs_names.

    Use this for all-cell mode. It also works for anchor mode, but then
    non-anchor cells will receive NaN values.

    Gene-level arrays are stored in adata.layers[prefix + '_' + key].
    Cell-level scalar arrays are stored in adata.obs[prefix + '_' + key].
    """
    if gene_level_keys is None:
        gene_level_keys = [
            "avg_expr_diff",
            "directional_p",
            "directional_fdr",
            "p_up",
            "p_down",
            "fdr_up",
            "fdr_down",
            "p_combined",
            "fdr_combined",
            "gene_score",
            "direction_consistency",
        ]

    if obs_level_keys is None:
        obs_level_keys = [
            "cell_divergence_score",
            "n_sig_genes",
            "n_sig_up",
            "n_sig_down",
            "n_pairs",
            "top_genes_per_list_used",
            "focal_idx",
            "is_anchor_mode",
        ]

    n_cells, n_genes = adata.shape
    cell_to_idx = {str(name): i for i, name in enumerate(adata.obs_names)}

    def out_name(key: str) -> str:
        return f"{prefix}_{key}" if prefix else key

    gene_matrices = {
        key: np.full((n_cells, n_genes), np.nan, dtype=np.float32)
        for key in gene_level_keys
    }

    obs_arrays = {
        key: np.full(n_cells, np.nan, dtype=np.float32)
        for key in obs_level_keys
    }

    # Object dtype arrays for string metadata.
    obs_object_arrays: Dict[str, np.ndarray] = {
        "focal_id": np.full(n_cells, None, dtype=object),
    }

    n_mapped = 0

    for cell_id, metrics in cell_importances.items():
        cell_id = str(cell_id)
        if cell_id not in cell_to_idx:
            continue

        idx = cell_to_idx[cell_id]
        n_mapped += 1

        for key in gene_level_keys:
            if key not in metrics:
                continue
            value = np.asarray(metrics[key], dtype=float)
            if value.ndim != 1 or value.shape[0] != n_genes:
                raise ValueError(
                    f"Result key '{key}' for '{cell_id}' has shape {value.shape}; "
                    f"expected ({n_genes},)."
                )
            gene_matrices[key][idx, :] = value.astype(np.float32)

        for key in obs_level_keys:
            if key not in metrics:
                continue
            obs_arrays[key][idx] = float(metrics[key])

        if "focal_id" in metrics:
            obs_object_arrays["focal_id"][idx] = str(metrics["focal_id"])

    for key, matrix in gene_matrices.items():
        name = out_name(key)
        if not overwrite and name in adata.layers:
            raise ValueError(f"adata.layers['{name}'] already exists.")
        adata.layers[name] = matrix

    for key, values in obs_arrays.items():
        name = out_name(key)
        if not overwrite and name in adata.obs.columns:
            raise ValueError(f"adata.obs['{name}'] already exists.")
        adata.obs[name] = values

    for key, values in obs_object_arrays.items():
        name = out_name(key)
        if not overwrite and name in adata.obs.columns:
            raise ValueError(f"adata.obs['{name}'] already exists.")
        adata.obs[name] = values

    print(f"Mapped {n_mapped} RRA result entries to AnnData.")
    return adata


def make_results_adata(
    adata: sc.AnnData,
    cell_importances: Dict[str, Dict[str, Any]],
    run_info: Optional[Dict[str, Any]] = None,
    prefix: str = "rra",
    subset_to_results: bool = True,
    copy: bool = True,
    gene_level_keys: Optional[List[str]] = None,
    obs_level_keys: Optional[List[str]] = None,
) -> sc.AnnData:
    """
    Creates an AnnData object containing RRA results.

    For anchor mode, the recommended use is:
        anchor_adata = make_results_adata(
            adata,
            cell_importances,
            run_info,
            subset_to_results=True,
        )

    This returns an anchor-only AnnData object with RRA gene-level outputs in
    .layers and anchor-level summaries in .obs.

    For all-cell mode, subset_to_results=False keeps the original adata and
    adds the results to all cells.
    """
    if subset_to_results:
        if run_info is not None and "successful_focal_indices" in run_info:
            idx = np.asarray(run_info["successful_focal_indices"], dtype=int)
        else:
            name_to_idx = {str(name): i for i, name in enumerate(adata.obs_names)}
            idx = np.array(
                [name_to_idx[str(k)] for k in cell_importances if str(k) in name_to_idx],
                dtype=int,
            )

        result_adata = adata[idx, :].copy() if copy else adata[idx, :]
        result_adata.obs[f"{prefix}_original_cell_idx"] = idx
        result_adata.obs[f"{prefix}_is_result_cell"] = True
    else:
        result_adata = adata.copy() if copy else adata

    result_adata = add_rra_results_to_adata(
        adata=result_adata,
        cell_importances=cell_importances,
        prefix=prefix,
        gene_level_keys=gene_level_keys,
        obs_level_keys=obs_level_keys,
        overwrite=True,
    )

    if run_info is not None:
        result_adata.uns[f"{prefix}_run_info"] = {
            "anchor_mode": bool(run_info.get("anchor_mode", False)),
            "n_anchors_per_celltype": run_info.get("n_anchors_per_celltype", None),
            "n_successful_focal_indices": int(len(run_info.get("successful_focal_indices", []))),
            "n_requested_anchor_indices": (
                None if run_info.get("anchor_indices", None) is None
                else int(len(run_info.get("anchor_indices")))
            ),
            "gene_names": np.asarray(run_info.get("gene_names", result_adata.var_names.to_numpy())).astype(str),
        }

    return result_adata




# =============================================================================
# Five-dataset, one-core runtime benchmark driver
# =============================================================================
DATA_DIR = Path(os.environ.get("RUNTIME_DATA_DIR", "/home/anilprakash/labs/Mei/projects/anil/srda/notebooks/data/scrna_seq/simulation"))
MANIFEST = Path(os.environ.get("RUNTIME_MANIFEST", DATA_DIR / "simulation_manifest_runtime_5_cell_counts.csv"))
RESULTS_DIR = Path(os.environ.get("SCCST_RUNTIME_RESULTS", DATA_DIR / "sccst_runtime_5_cell_counts"))
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
N_HVGS = 2000
N_PCS = 50
NHOOD_SIZE = int(os.environ.get("SCCST_NHOOD_SIZE", "20"))
DEFAULT_HARMONY_THETA = 2.0
DEFAULT_RESPONDER_PERCENT = 40
DEFAULT_DROPOUT_RATE = 0.4
sc.settings.n_jobs = 1


def preprocess_for_sccst(adata: sc.AnnData, theta: float = DEFAULT_HARMONY_THETA) -> sc.AnnData:
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    run = sc.AnnData(X=counts.copy(), obs=adata.obs.copy(), var=adata.var.copy())
    if sp.issparse(run.X):
        run.X = run.X.tocsr()
    sc.pp.normalize_total(run, target_sum=1e4)
    sc.pp.log1p(run)
    sc.pp.highly_variable_genes(run, n_top_genes=min(N_HVGS, run.n_vars), flavor="seurat")
    if "highly_variable" in run.var and int(run.var["highly_variable"].sum()) > 0:
        run = run[:, run.var["highly_variable"]].copy()
    n_comps = min(N_PCS, run.n_vars - 1, run.n_obs - 1)
    sc.tl.pca(run, n_comps=n_comps, svd_solver="arpack")
    sc.external.pp.harmony_integrate(run, key="sim_batch", theta=theta, max_iter_harmony=50)
    return run


def resolve_file(path_value: str) -> Path:
    p = Path(path_value)
    if p.exists(): return p
    q = DATA_DIR / p.name
    if q.exists(): return q
    raise FileNotFoundError(path_value)


def main() -> None:
    manifest = pd.read_csv(MANIFEST)
    if "status" in manifest:
        manifest = manifest.loc[manifest["status"].astype(str).str.lower().eq("success")]
    manifest = manifest.sort_values("base_n_cells").head(5)
    if len(manifest) != 5:
        raise ValueError(f"Expected exactly five successful runtime datasets, found {len(manifest)}")
    rows = []
    for _, row in manifest.iterrows():
        input_file = resolve_file(str(row["file"]))
        loaded = sc.read_h5ad(input_file)
        t0 = time.perf_counter()
        run = preprocess_for_sccst(loaded)
        preprocessing_seconds = time.perf_counter() - t0
        t1 = time.perf_counter()
        cell_importances, gene_names, run_info = process_cell_type(
            adata=run, conditions=("Disease", "Control"), condition_col="condition",
            batch_col="sim_batch", sample_col=None, n_neighborhoods=NHOOD_SIZE,
            cell_type_col="Cell_Type", n_jobs=1, latent_key="X_pca_harmony",
            caliper_percentile=75.0, min_cells_per_batch_neighborhood=3,
            min_matched_pairs=3, top_genes_per_list=100, fdr_threshold=0.05,
            n_anchors_per_celltype=None, anchor_indices=None, anchor_random_state=0)
        analysis_seconds = time.perf_counter() - t1
        rows.append({"method":"scCST","input_file":str(input_file),
            "base_n_cells":int(row["base_n_cells"]),"n_cells":int(run.n_obs),
            "n_genes_after_hvg":int(run.n_vars),"nhood_size":NHOOD_SIZE,"n_cores":1,
            "responder_percent":DEFAULT_RESPONDER_PERCENT,"dropout_rate":DEFAULT_DROPOUT_RATE,
            "harmony_theta":DEFAULT_HARMONY_THETA,"preprocessing_seconds":preprocessing_seconds,
            "analysis_seconds":analysis_seconds,"analysis_component":"process_cell_type",
            "benchmark_total_seconds":preprocessing_seconds+analysis_seconds,
            "n_successful_focal_cells":len(cell_importances)})
        pd.DataFrame(rows).to_csv(RESULTS_DIR/"sccst_runtime_summary.csv",index=False)
        del loaded, run, cell_importances
    print(pd.DataFrame(rows))

if __name__ == "__main__":
    main()
