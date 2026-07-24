import pytest

pytest.importorskip("scanpy")

import numpy as np

from sccst.core import apply_bh_fdr, rra_bonferroni_pvalues


def test_apply_bh_fdr_shape_and_bounds():
    p = np.array([0.01, 0.04, 0.2, np.nan])
    q = apply_bh_fdr(p)
    assert q.shape == p.shape
    assert np.isnan(q[-1])
    assert np.all((q[:-1] >= 0) & (q[:-1] <= 1))


def test_rra_bonferroni_pvalues_shape_and_bounds():
    ranks = np.array([[0.01, 0.2, 1.0], [0.05, 0.4, 1.0]])
    p = rra_bonferroni_pvalues(ranks)
    assert p.shape == (3,)
    assert np.all((p >= 0) & (p <= 1))
