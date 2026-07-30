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
# Configured for the datasets created by the parameter-sweep simulation notebook.
# 
# Expected simulation directory:
# 
# ```python
# /path/to/data/scrna_seq/simulation
# ```
# 
# It reads `simulation_manifest_de_parameter_sweeps.csv` when available, runs scCST on each generated `.h5ad`, and writes one result object plus per-cell summaries for each dataset and neighborhood-size setting. The output filenames include `sweep_name` and `replicate` so the repeated default-parameter datasets from different sweeps do not overwrite each other.
# 

# In[3]:


#sccst performance on simulated data.

import scanpy as sc
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


# In[ ]:


# ============================================================
# Run scCST on parameter-sweep replicated simulation datasets
# ============================================================
# This cell is configured for outputs from:
# create_sim_parameter_sweeps_replicates.ipynb
#
# Expected inputs:
#   simulation_manifest_de_parameter_sweeps.csv
#   simulation_adata_de_<sweep_name>_replicate_<replicate>_responder_percent_<...>_dropout_<...>_harmony_theta_<...>.h5ad
#
# Expected total datasets from the requested design:
#   3 responder-percent values x 3 replicates = 9
#   3 dropout-rate values x 3 replicates = 9
#   3 Harmony-theta values x 3 replicates = 9
#   total = 27

from pathlib import Path
import re
import numpy as np
import pandas as pd
import scanpy as sc


# ------------------------------------------------------------
# Input/output directories
# ------------------------------------------------------------
data_dir = Path(
    "/path/to/data/"
    "scrna_seq/simulation"
)
results_dir = data_dir / "sccst_results"
results_dir.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# scCST run settings
# ------------------------------------------------------------
conditions = ("Disease", "Control")

condition_col = "condition"
batch_col = "sim_batch"
cell_type_col = "Cell_Type"
latent_key = "X_pca_harmony"

# You can test multiple neighborhood sizes if desired.
# For a first full run on all 27 datasets, start with [10].
nhood_sizes = [10, 20, 30]

caliper_percentile = 30.0
top_genes_per_list = 100
fdr_threshold = 0.05

# For 2,400-cell simulations, all-cell mode is usually okay.
# If runtime is too slow, set n_anchors_per_celltype to an integer, e.g. 100.
n_anchors_per_celltype = None

# Inner parallelism: focal-cell jobs inside process_cell_type.
# Override with SCCST_INNER_JOBS or --inner-jobs.
inner_jobs = int(os.environ.get("SCCST_INNER_JOBS", "16"))
n_jobs = inner_jobs

# Outer parallelism: number of dataset x neighborhood-size runs in parallel.
# Total CPU pressure is roughly outer_workers x inner_jobs.
# Start conservatively, e.g. outer=2 and inner=8 on a 16-core node.
outer_workers = int(os.environ.get("SCCST_OUTER_WORKERS", "1"))
outer_workers = max(1, outer_workers)

# Set to True to skip result files that already exist.
skip_existing = True

# Safer resume mode: when skip_existing=True, verify existing .h5ad can be opened.
# If verification fails, the file is treated as incomplete and rerun.
verify_existing = True

# Optional quick test.
test_mode = False
max_runs_for_test = 2


def _parse_bool(x):
    return str(x).lower() in {"true", "t", "1", "yes", "y"}


def parse_cli_args():
    """Small CLI parser for HPC batch scripts without adding dependencies."""
    global data_dir, results_dir, manifest_file, nhood_sizes
    global outer_workers, inner_jobs, n_jobs, skip_existing, verify_existing
    global test_mode, max_runs_for_test, n_anchors_per_celltype

    args = os.sys.argv[1:]

    def get_arg(prefix):
        hits = [a for a in args if a.startswith(prefix)]
        if not hits:
            return None
        return hits[-1][len(prefix):]

    x = get_arg("--data-dir=")
    if x:
        data_dir = Path(x)
        manifest_file = data_dir / "simulation_manifest_de_parameter_sweeps.csv"

    x = get_arg("--results-dir=")
    if x:
        results_dir = Path(x)

    x = get_arg("--nhood-sizes=")
    if x:
        nhood_sizes = [int(v.strip()) for v in x.split(",") if v.strip()]

    x = get_arg("--outer-workers=")
    if x:
        outer_workers = max(1, int(x))

    x = get_arg("--inner-jobs=")
    if x:
        inner_jobs = int(x)
        n_jobs = inner_jobs

    x = get_arg("--skip-existing=")
    if x is not None:
        skip_existing = _parse_bool(x)

    x = get_arg("--verify-existing=")
    if x is not None:
        verify_existing = _parse_bool(x)

    x = get_arg("--test-mode=")
    if x is not None:
        test_mode = _parse_bool(x)

    x = get_arg("--max-runs-for-test=")
    if x:
        max_runs_for_test = int(x)

    x = get_arg("--n-anchors-per-celltype=")
    if x:
        if x.lower() in {"none", "null", "na"}:
            n_anchors_per_celltype = None
        else:
            n_anchors_per_celltype = int(x)


parse_cli_args()
results_dir.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Locate parameter-sweep simulation files
# ------------------------------------------------------------
manifest_file = data_dir / "simulation_manifest_de_parameter_sweeps.csv"

if manifest_file.exists():
    manifest_df = pd.read_csv(manifest_file)
    manifest_df["file"] = manifest_df["file"].astype(str)

    # Prefer paths recorded in the manifest, but repair stale absolute paths
    # by resolving them against the current data_dir. This allows the SCCST
    # notebook to run after changing the simulation notebook's output_dir.
    input_files = []
    for recorded_file in manifest_df["file"].tolist():
        recorded_path = Path(recorded_file)
        if recorded_path.exists():
            input_files.append(recorded_path)
        else:
            relocated_path = data_dir / recorded_path.name
            input_files.append(relocated_path)
else:
    input_files = sorted(
        data_dir.glob(
            "simulation_adata_de_*_replicate_*_responder_percent_*_dropout_*_harmony_theta_*.h5ad"
        )
    )
    manifest_df = pd.DataFrame({"file": [str(x) for x in input_files]})

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def _float_tag_to_value(x):
    if x is None:
        return None
    return float(str(x).replace("p", "."))


def _value_to_tag(x):
    return str(x).replace(".", "p")


def parse_simulation_params_from_filename(path):
    """
    Parse filenames like:
    simulation_adata_de_responder_percent_sweep_replicate_1_responder_percent_40_dropout_0p2_harmony_theta_2.h5ad
    """
    name = path.name

    pattern = (
        r"simulation_adata_(?P<simulation_type>.+?)"
        r"_(?P<sweep_name>.+?)"
        r"_replicate_(?P<replicate>\d+)"
        r"_responder_percent_(?P<responder_percent>\d+)"
        r"_dropout_(?P<dropout>[\dp]+)"
        r"_harmony_theta_(?P<harmony_theta>[\dp]+)"
        r"\.h5ad"
    )

    m = re.match(pattern, name)
    if m is None:
        return {
            "simulation_type": None,
            "sweep_name": None,
            "swept_parameter": None,
            "swept_value": None,
            "replicate": None,
            "responder_percent": None,
            "dropout_rate": None,
            "harmony_theta": None,
        }

    return {
        "simulation_type": m.group("simulation_type"),
        "sweep_name": m.group("sweep_name"),
        "swept_parameter": None,
        "swept_value": None,
        "replicate": int(m.group("replicate")),
        "responder_percent": int(m.group("responder_percent")),
        "dropout_rate": _float_tag_to_value(m.group("dropout")),
        "harmony_theta": _float_tag_to_value(m.group("harmony_theta")),
    }


def get_manifest_params(input_file, manifest_df):
    """
    Prefer the manifest metadata because it includes swept_parameter and swept_value.
    Fall back to filename parsing if the manifest is unavailable or incomplete.
    """
    parsed = parse_simulation_params_from_filename(input_file)

    if manifest_df is None or manifest_df.empty or "file" not in manifest_df.columns:
        parsed["input_file"] = str(input_file)
        return parsed

    matches = manifest_df.loc[manifest_df["file"].astype(str) == str(input_file)]

    # Also try basename matching in case the manifest was moved with the data directory.
    if matches.empty:
        matches = manifest_df.loc[
            manifest_df["file"].astype(str).map(lambda x: Path(x).name) == input_file.name
        ]

    if matches.empty:
        parsed["input_file"] = str(input_file)
        return parsed

    row = matches.iloc[0].to_dict()

    params = parsed.copy()
    for key in [
        "simulation_type",
        "sweep_name",
        "swept_parameter",
        "swept_value",
        "replicate",
        "responder_percent",
        "dropout_rate",
        "harmony_theta",
        "seed",
        "n_cells",
        "n_genes",
        "n_responders",
        "n_leiden_clusters",
    ]:
        if key in row and pd.notna(row[key]):
            params[key] = row[key]

    # Normalize numeric types where possible.
    for key in ["replicate", "responder_percent", "seed", "n_cells", "n_genes", "n_responders", "n_leiden_clusters"]:
        if params.get(key) is not None and pd.notna(params.get(key)):
            params[key] = int(params[key])

    for key in ["dropout_rate", "harmony_theta"]:
        if params.get(key) is not None and pd.notna(params.get(key)):
            params[key] = float(params[key])

    params["input_file"] = str(input_file)
    return params


def sanitize_adata_for_h5ad(adata):
    """
    Avoid common AnnData write errors from object columns or mixed dtypes.
    """
    for col in adata.obs.columns:
        if adata.obs[col].dtype == "object":
            adata.obs[col] = adata.obs[col].astype(str)

    for col in adata.var.columns:
        if adata.var[col].dtype == "object":
            adata.var[col] = adata.var[col].astype(str)

    return adata


def summarize_sccst_result(result_adata, params, nhood_size, out_file):
    """
    Create one-row summary for each scCST result object.
    """
    row = {
        "input_file": params.get("input_file"),
        "result_file": str(out_file),
        "simulation_type": params.get("simulation_type"),
        "sweep_name": params.get("sweep_name"),
        "swept_parameter": params.get("swept_parameter"),
        "swept_value": params.get("swept_value"),
        "replicate": params.get("replicate"),
        "responder_percent": params.get("responder_percent"),
        "dropout_rate": params.get("dropout_rate"),
        "harmony_theta": params.get("harmony_theta"),
        "seed": params.get("seed"),
        "nhood_size": nhood_size,
        "n_result_cells": int(result_adata.n_obs),
        "n_genes": int(result_adata.n_vars),
    }

    if "rra_cell_divergence_score" in result_adata.obs.columns:
        div = pd.to_numeric(
            result_adata.obs["rra_cell_divergence_score"],
            errors="coerce",
        )

        row.update(
            {
                "mean_cell_divergence_score": float(np.nanmean(div)),
                "median_cell_divergence_score": float(np.nanmedian(div)),
                "max_cell_divergence_score": float(np.nanmax(div)),
            }
        )

    if "rra_n_sig_genes" in result_adata.obs.columns:
        n_sig = pd.to_numeric(
            result_adata.obs["rra_n_sig_genes"],
            errors="coerce",
        )

        row.update(
            {
                "mean_n_sig_genes": float(np.nanmean(n_sig)),
                "median_n_sig_genes": float(np.nanmedian(n_sig)),
                "max_n_sig_genes": float(np.nanmax(n_sig)),
            }
        )

    if cell_type_col in result_adata.obs.columns:
        row["n_cell_types"] = int(result_adata.obs[cell_type_col].nunique())

    return row


def build_output_file(results_dir, params, nhood_size):
    sweep_name = params.get("sweep_name") or "unknown_sweep"
    replicate = params.get("replicate") if params.get("replicate") is not None else "NA"
    responder_tag = params.get("responder_percent") if params.get("responder_percent") is not None else "NA"
    dropout_tag = _value_to_tag(params.get("dropout_rate") if params.get("dropout_rate") is not None else "NA")
    theta_tag = _value_to_tag(params.get("harmony_theta") if params.get("harmony_theta") is not None else "NA")
    swept_parameter = params.get("swept_parameter") or "unknown_parameter"
    swept_value = _value_to_tag(params.get("swept_value") if params.get("swept_value") is not None else "NA")

    return (
        results_dir
        / f"sccst_result"
          f"_{sweep_name}"
          f"_{swept_parameter}_{swept_value}"
          f"_replicate_{replicate}"
          f"_nhood_{nhood_size}"
          f"_responder_percent_{responder_tag}"
          f"_dropout_{dropout_tag}"
          f"_harmony_theta_{theta_tag}"
          f".h5ad"
    )


# ------------------------------------------------------------
# Main loop helpers with outer parallelism
# ------------------------------------------------------------
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback


def existing_h5ad_is_valid(path: Path) -> bool:
    """Return True if an existing .h5ad appears readable."""
    if not path.exists():
        return False
    if not verify_existing:
        return True
    try:
        tmp = sc.read_h5ad(path, backed="r")
        tmp.file.close()
        return True
    except Exception as e:
        print(f"Existing output appears invalid and will be rerun: {path}\n  {repr(e)}")
        return False


def build_run_tasks(input_files, manifest_df):
    tasks = []
    missing_rows = []

    for input_file in input_files:
        input_file = Path(input_file)

        if not input_file.exists():
            print(f"Skipping missing file: {input_file}")
            missing_rows.append({
                "input_file": str(input_file),
                "nhood_size": None,
                "error": "Input file does not exist",
            })
            continue

        params = get_manifest_params(input_file, manifest_df)

        for nhood_size in nhood_sizes:
            out_file = build_output_file(results_dir, params, nhood_size)
            tasks.append({
                "input_file": str(input_file),
                "params": params,
                "nhood_size": int(nhood_size),
                "out_file": str(out_file),
            })

    if test_mode:
        tasks = tasks[:max_runs_for_test]

    return tasks, missing_rows


def run_one_sccst_task(task):
    """Run one dataset x neighborhood-size combination."""
    input_file = Path(task["input_file"])
    params = task["params"]
    nhood_size = int(task["nhood_size"])
    out_file = Path(task["out_file"])

    try:
        if skip_existing and existing_h5ad_is_valid(out_file):
            return {
                "status": "skipped",
                "skipped": {
                    "input_file": str(input_file),
                    "result_file": str(out_file),
                    "nhood_size": nhood_size,
                    "reason": "valid result already exists",
                },
            }

        print("\n============================================================")
        print(f"Running scCST on: {input_file.name}")
        print(f"nhood_size = {nhood_size}")
        print(f"outer_workers = {outer_workers}; inner_jobs = {inner_jobs}")
        print(f"params = {params}")
        print("============================================================")

        # Read one simulated dataset
        adata = sc.read_h5ad(input_file)

        # Basic validation
        required_obs = [condition_col, batch_col, cell_type_col]
        for col in required_obs:
            if col not in adata.obs.columns:
                raise ValueError(f"Missing required adata.obs column: {col}")

        if latent_key not in adata.obsm:
            raise ValueError(f"Missing required adata.obsm key: {latent_key}")

        print(adata)
        print("Condition counts:")
        print(adata.obs[condition_col].value_counts())

        print("Batch counts:")
        print(adata.obs[batch_col].value_counts())

        print("Leiden-derived Cell_Type counts:")
        print(adata.obs[cell_type_col].value_counts())

        # ------------------------------------------------
        # Run scCST
        # Important: process_cell_type filters genes in place.
        # Use adata.copy() so the original loaded object is preserved.
        # ------------------------------------------------
        run_adata = adata.copy()

        cell_importances, gene_names, run_info = process_cell_type(
            adata=run_adata,
            conditions=conditions,
            condition_col=condition_col,
            batch_col=batch_col,
            sample_col=None,
            n_neighborhoods=nhood_size,
            cell_type_col=cell_type_col,
            n_jobs=inner_jobs,
            latent_key=latent_key,
            caliper_percentile=caliper_percentile,
            min_cells_per_batch_neighborhood=3,
            min_matched_pairs=3,
            top_genes_per_list=top_genes_per_list,
            fdr_threshold=fdr_threshold,
            n_anchors_per_celltype=n_anchors_per_celltype,
            anchor_indices=None,
            anchor_random_state=0,
        )

        # ------------------------------------------------
        # Make result AnnData
        # subset_to_results=True keeps only successful focal cells/anchors.
        # ------------------------------------------------
        result_adata = make_results_adata(
            adata=run_adata,
            cell_importances=cell_importances,
            run_info=run_info,
            prefix="rra",
            subset_to_results=True,
            copy=True,
            gene_level_keys=[
                "gene_score",
                "avg_expr_diff",
                "directional_fdr",
            ],
            obs_level_keys=[
                "cell_divergence_score",
                "n_sig_genes",
                "n_pairs",
            ],
        )

        # Store run metadata safely
        result_adata.uns["sccst_run_params"] = {
            "input_file": str(input_file),
            "conditions_disease": str(conditions[0]),
            "conditions_control": str(conditions[1]),
            "condition_col": str(condition_col),
            "batch_col": str(batch_col),
            "cell_type_col": str(cell_type_col),
            "latent_key": str(latent_key),
            "nhood_size": int(nhood_size),
            "caliper_percentile": float(caliper_percentile),
            "top_genes_per_list": int(top_genes_per_list),
            "fdr_threshold": float(fdr_threshold),
            "n_anchors_per_celltype": (
                "None" if n_anchors_per_celltype is None
                else str(n_anchors_per_celltype)
            ),
            "inner_jobs": int(inner_jobs),
            "outer_workers": int(outer_workers),
            "simulation_type": str(params.get("simulation_type")),
            "sweep_name": str(params.get("sweep_name")),
            "swept_parameter": str(params.get("swept_parameter")),
            "swept_value": str(params.get("swept_value")),
            "replicate": -1 if params.get("replicate") is None else int(params.get("replicate")),
            "responder_percent": -1 if params.get("responder_percent") is None else int(params.get("responder_percent")),
            "dropout_rate": -1.0 if params.get("dropout_rate") is None else float(params.get("dropout_rate")),
            "harmony_theta": -1.0 if params.get("harmony_theta") is None else float(params.get("harmony_theta")),
            "seed": -1 if params.get("seed") is None else int(params.get("seed")),
        }

        result_adata = sanitize_adata_for_h5ad(result_adata)

        # ------------------------------------------------
        # Save result atomically: write temp file, then rename.
        # This reduces the chance that cancelled runs leave a corrupt file
        # that skip_existing would later skip.
        # ------------------------------------------------
        tmp_out_file = out_file.with_suffix(out_file.suffix + ".tmp")
        if tmp_out_file.exists():
            tmp_out_file.unlink()
        result_adata.write_h5ad(tmp_out_file, compression="gzip")
        tmp_out_file.replace(out_file)

        # ------------------------------------------------
        # Save per-cell summary
        # ------------------------------------------------
        per_cell_summary_file = out_file.with_suffix(".per_cell_summary.csv")

        per_cell_cols = [
            condition_col,
            batch_col,
            cell_type_col,
            "sim_group_truth",
            "is_responder",
            "rra_cell_divergence_score",
            "rra_n_sig_genes",
            "rra_n_pairs",
        ]

        available_cols = [
            col for col in per_cell_cols
            if col in result_adata.obs.columns
        ]

        per_cell_summary = result_adata.obs[available_cols].copy()
        for key in [
            "simulation_type",
            "sweep_name",
            "swept_parameter",
            "swept_value",
            "replicate",
            "responder_percent",
            "dropout_rate",
            "harmony_theta",
            "seed",
            "nhood_size",
        ]:
            per_cell_summary[key] = params.get(key, nhood_size if key == "nhood_size" else None)
        per_cell_summary["nhood_size"] = nhood_size
        per_cell_summary.to_csv(per_cell_summary_file)

        # ------------------------------------------------
        # Add one-row run summary
        # ------------------------------------------------
        summary_row = summarize_sccst_result(
            result_adata=result_adata,
            params=params,
            nhood_size=nhood_size,
            out_file=out_file,
        )

        summary_row["per_cell_summary_file"] = str(per_cell_summary_file)
        summary_row["inner_jobs"] = int(inner_jobs)
        summary_row["outer_workers"] = int(outer_workers)

        print(f"Saved scCST result: {out_file}")
        print(f"Saved per-cell summary: {per_cell_summary_file}")

        return {"status": "success", "summary": summary_row}

    except Exception as e:
        print(f"FAILED: {input_file.name}, nhood_size={nhood_size}")
        print(repr(e))

        failed_row = {
            "input_file": str(input_file),
            "nhood_size": nhood_size,
            "error": repr(e),
            "traceback": traceback.format_exc(),
        }
        failed_row.update(params)
        return {"status": "failed", "failed": failed_row}


def main():
    print(f"Simulation input directory: {data_dir}")
    print(f"scCST results directory: {results_dir}")
    print(f"Found {len(input_files)} simulation files.")
    print(f"Neighborhood sizes: {nhood_sizes}")
    print(f"Outer workers: {outer_workers}")
    print(f"Inner focal-cell jobs per run: {inner_jobs}")
    print(f"Approx active CPU pressure: {outer_workers * max(1, inner_jobs)}")
    print(f"skip_existing: {skip_existing}; verify_existing: {verify_existing}")
    if len(input_files) != 27:
        print(
            "Note: the requested parameter-sweep design should create 27 datasets. "
            "If this count differs, check that the simulation notebook completed successfully "
            "and that data_dir points to the same output directory."
        )

    tasks, missing_rows = build_run_tasks(input_files, manifest_df)
    print(f"Total dataset x nhood tasks to consider: {len(tasks)}")

    all_summary_rows = []
    failed_runs = list(missing_rows)
    skipped_runs = []

    if outer_workers == 1:
        results = [run_one_sccst_task(task) for task in tasks]
    else:
        results = []
        with ProcessPoolExecutor(max_workers=outer_workers) as ex:
            future_to_task = {ex.submit(run_one_sccst_task, task): task for task in tasks}
            for fut in as_completed(future_to_task):
                results.append(fut.result())

    for result in results:
        if result["status"] == "success":
            all_summary_rows.append(result["summary"])
        elif result["status"] == "skipped":
            skipped_runs.append(result["skipped"])
        elif result["status"] == "failed":
            failed_runs.append(result["failed"])

    # ------------------------------------------------------------
    # Save full run summary
    # ------------------------------------------------------------
    summary_df = pd.DataFrame(all_summary_rows)
    summary_file = results_dir / "sccst_parameter_sweeps_run_summary.csv"
    summary_df.to_csv(summary_file, index=False)

    failed_df = pd.DataFrame(failed_runs)
    failed_file = results_dir / "sccst_parameter_sweeps_failed_runs.csv"
    failed_df.to_csv(failed_file, index=False)

    skipped_df = pd.DataFrame(skipped_runs)
    skipped_file = results_dir / "sccst_parameter_sweeps_skipped_runs.csv"
    skipped_df.to_csv(skipped_file, index=False)

    # Optional aggregate summaries by sweep/value/replicate.
    if not summary_df.empty:
        aggregate_cols = ["sweep_name", "swept_parameter", "swept_value", "nhood_size"]
        metric_cols = [
            col for col in [
                "mean_cell_divergence_score",
                "median_cell_divergence_score",
                "max_cell_divergence_score",
                "mean_n_sig_genes",
                "median_n_sig_genes",
                "max_n_sig_genes",
                "n_result_cells",
                "n_genes",
            ]
            if col in summary_df.columns
        ]

        aggregate_summary = (
            summary_df
            .groupby(aggregate_cols, dropna=False)[metric_cols]
            .agg(["mean", "std", "min", "max"])
        )
        aggregate_summary.columns = ["_".join(col).strip("_") for col in aggregate_summary.columns]
        aggregate_summary = aggregate_summary.reset_index()

        aggregate_file = results_dir / "sccst_parameter_sweeps_aggregate_summary.csv"
        aggregate_summary.to_csv(aggregate_file, index=False)
    else:
        aggregate_summary = pd.DataFrame()
        aggregate_file = results_dir / "sccst_parameter_sweeps_aggregate_summary.csv"
        aggregate_summary.to_csv(aggregate_file, index=False)

    print("\nDone.")
    print(f"Successful runs: {len(summary_df)}")
    print(f"Failed runs: {len(failed_df)}")
    print(f"Skipped existing runs: {len(skipped_df)}")
    print(f"Summary saved to: {summary_file}")
    print(f"Failed-run log saved to: {failed_file}")
    print(f"Skipped-run log saved to: {skipped_file}")
    print(f"Aggregate summary saved to: {aggregate_file}")

    return summary_df


if __name__ == "__main__":
    summary_df = main()
