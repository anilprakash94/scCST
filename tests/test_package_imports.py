import pytest

pytest.importorskip("scanpy")


def test_import_core_symbols():
    import sccst

    assert hasattr(sccst, "process_cell_type")
    assert hasattr(sccst, "make_results_adata")
