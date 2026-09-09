"""Behaviour tests for the core parameters of process_cell_type."""

import numpy as np
import pytest

pytest.importorskip("scanpy")

import anndata as ad
import scanpy as sc

from sccst import make_results_adata, process_cell_type
from sccst.core import check_log_normalized


def _toy_adata(
    seed: int = 0,
    n_per_group: int = 30,
    n_genes: int = 60,
    n_de_genes: int = 8,
    de_shift: float = 0.25,
):
    """Two cell types, two conditions, two replicate batches per condition.

    The disease effect is deliberately confined to a few genes and kept small
    next to the cell-type signal. A condition effect spread across every gene
    would push case and control cells apart in the embedding, so no case cell
    would ever share a neighborhood with a control cell and zero matched pairs
    would form -- a silent empty result rather than a test failure.
    """
    rng = np.random.default_rng(seed)
    blocks, obs = [], []

    for ctype in ("A", "B"):
        for cond, batches in (("disease", ("d1", "d2")), ("control", ("c1", "c2"))):
            for batch in batches:
                centre = np.full(n_genes, 1.0 if ctype == "A" else 2.0)
                if cond == "disease":
                    centre[:n_de_genes] += de_shift
                blocks.append(np.clip(rng.normal(centre, 0.3, (n_per_group, n_genes)), 0, None))
                obs.extend([(ctype, cond, batch, batch[-1])] * n_per_group)

    adata = ad.AnnData(np.vstack(blocks).astype(np.float64))
    adata.obs["Cell_Type"] = [o[0] for o in obs]
    adata.obs["condition"] = [o[1] for o in obs]
    adata.obs["batch"] = [o[2] for o in obs]
    adata.obs["sample"] = [o[3] for o in obs]
    adata.obs_names = [f"cell{i}" for i in range(adata.n_obs)]
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    sc.pp.pca(adata, n_comps=10)
    return adata


RUN_KWARGS = dict(
    conditions=("disease", "control"),
    condition_col="condition",
    batch_col="batch",
    cell_type_col="Cell_Type",
    latent_key="X_pca",
    n_jobs=1,
)


# --- the log-normalisation heuristic -------------------------------------

def test_lognorm_check_passes_on_lognormalised_values():
    rng = np.random.default_rng(0)
    assert check_log_normalized(rng.random((50, 20)) * 5.0) is None


def test_lognorm_check_flags_raw_counts():
    rng = np.random.default_rng(0)
    message = check_log_normalized(rng.integers(0, 500, (50, 20)).astype(float))
    assert message is not None and "RAW COUNTS" in message


def test_lognorm_check_flags_normalised_but_unlogged():
    rng = np.random.default_rng(0)
    message = check_log_normalized(rng.random((50, 20)) * 1000.0)
    assert message is not None and "log-transformed" in message


def test_process_cell_type_warns_on_raw_counts():
    adata = _toy_adata()
    adata.X = np.round(adata.X * 200.0)
    with pytest.warns(UserWarning, match="RAW COUNTS"):
        process_cell_type(adata, neighborhood_size=10, **RUN_KWARGS)


# --- min_matched_pairs now defaults to 2 ---------------------------------

def test_two_replicates_per_condition_produce_results():
    """Exactly two matched pairs are available, so the default must be <= 2."""
    adata = _toy_adata()
    importances, _, _ = process_cell_type(adata, neighborhood_size=10, **RUN_KWARGS)
    assert len(importances) > 0
    assert all(m["n_pairs"] >= 2 for m in importances.values())


def test_requiring_three_pairs_yields_nothing_on_two_replicates():
    adata = _toy_adata()
    importances, _, _ = process_cell_type(
        adata, neighborhood_size=10, min_matched_pairs=3, **RUN_KWARGS
    )
    assert len(importances) == 0


# --- the neighborhood_size rename ----------------------------------------

def test_deprecated_n_neighborhoods_still_works_and_warns():
    adata = _toy_adata()
    with pytest.warns(DeprecationWarning, match="n_neighborhoods"):
        importances, _, _ = process_cell_type(adata, n_neighborhoods=10, **RUN_KWARGS)
    assert len(importances) > 0


def test_passing_both_names_raises():
    adata = _toy_adata()
    with pytest.raises(ValueError, match="only one of"):
        process_cell_type(
            adata, neighborhood_size=10, n_neighborhoods=10, **RUN_KWARGS
        )


def test_neighborhood_size_must_be_positive():
    adata = _toy_adata()
    with pytest.raises(ValueError, match="positive integer"):
        process_cell_type(adata, neighborhood_size=0, **RUN_KWARGS)


# --- float64 results ------------------------------------------------------

def test_results_are_float64():
    adata = _toy_adata()
    importances, _, info = process_cell_type(adata, neighborhood_size=10, **RUN_KWARGS)
    result = make_results_adata(adata, importances, info, prefix="rra", subset_to_results=True)

    for layer in ("rra_gene_score", "rra_avg_expr_diff", "rra_directional_fdr"):
        assert result.layers[layer].dtype == np.float64, layer
    for col in ("rra_cell_divergence_score", "rra_n_sig_genes", "rra_n_pairs"):
        assert result.obs[col].dtype == np.float64, col
    assert np.isfinite(result.layers["rra_gene_score"]).any()
