"""scCST: single-cell Cell State Transition.

Reusable Python implementation extracted from the development notebooks.
"""

try:
    from .core import (
        apply_bh_fdr,
        rra_bonferroni_pvalues,
        compute_celltype_local_calipers,
        select_evenly_spaced_anchors,
        process_focal_cell_rra,
        process_cell_type,
        add_rra_results_to_adata,
        make_results_adata,
    )
except Exception:  # pragma: no cover
    # Allows metadata inspection even when optional scientific dependencies are absent.
    pass

__version__ = "0.1.0"
