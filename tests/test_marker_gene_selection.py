"""Tests for up-regulated marker gene selection used by GO enrichment."""

import numpy as np
import pytest

pytest.importorskip("scanpy")

from sccst.downstream import select_up_marker_genes


def test_keeps_genes_up_here_and_enriched_over_rest():
    genes = select_up_marker_genes(
        gene_names=["good"],
        mean_target=[1.5],
        mean_rest=[0.2],
        n_genes=10,
    )
    assert genes == ["good"]


def test_rejects_gene_that_is_only_less_down_than_the_rest():
    """The original bug: a big positive difference from two negative means."""
    genes = select_up_marker_genes(
        gene_names=["less_down"],
        mean_target=[-0.1],   # not up in this cell type at all
        mean_rest=[-2.0],     # difference is +1.9, which used to rank it first
        n_genes=10,
    )
    assert genes == []


def test_rejects_gene_up_here_but_lower_than_the_rest():
    """Up-regulated, but more so everywhere else, so not a marker for this type."""
    genes = select_up_marker_genes(
        gene_names=["ubiquitous"],
        mean_target=[0.5],
        mean_rest=[2.0],      # difference is negative
        n_genes=10,
    )
    assert genes == []


def test_both_failure_modes_are_dropped_from_a_mixed_set():
    genes = select_up_marker_genes(
        gene_names=["good_a", "less_down", "ubiquitous", "good_b"],
        mean_target=[1.5, -0.1, 0.5, 3.0],
        mean_rest=[0.2, -2.0, 2.0, 0.1],
        n_genes=10,
    )
    # good_b has the larger difference (2.9 against 1.3), so it ranks first.
    assert genes == ["good_b", "good_a"]


def test_ranks_by_difference_and_respects_n_genes():
    genes = select_up_marker_genes(
        gene_names=["a", "b", "c"],
        mean_target=[1.0, 1.0, 1.0],
        mean_rest=[0.9, 0.1, 0.5],   # differences: 0.1, 0.9, 0.5
        n_genes=2,
    )
    assert genes == ["b", "c"]


def test_returns_empty_when_nothing_qualifies():
    assert select_up_marker_genes(
        gene_names=["x", "y"],
        mean_target=[-1.0, -2.0],
        mean_rest=[-3.0, -4.0],
        n_genes=10,
    ) == []


def test_accepts_numpy_inputs():
    genes = select_up_marker_genes(
        gene_names=np.array(["a", "b"]),
        mean_target=np.array([2.0, -1.0]),
        mean_rest=np.array([0.0, -5.0]),
        n_genes=10,
    )
    assert genes == ["a"]
