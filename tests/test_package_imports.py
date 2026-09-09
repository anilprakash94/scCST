import pytest

pytest.importorskip("scanpy")


def test_import_core_symbols():
    import sccst

    assert hasattr(sccst, "process_cell_type")
    assert hasattr(sccst, "make_results_adata")


def test_public_plotting_names():
    from sccst import downstream, plotting

    # visualize_gene_stir was renamed to visualize_gene_map; both stay importable.
    assert downstream.visualize_gene_map.__name__ == "visualize_gene_map"
    assert plotting.visualize_gene_map is downstream.visualize_gene_map
    assert plotting.visualize_gene_stir is downstream.visualize_gene_stir
    assert {"visualize_gene_map", "visualize_gene_stir"} <= set(plotting.__all__)


def test_deprecated_visualize_gene_stir_warns_and_forwards(monkeypatch):
    from sccst import downstream

    seen = {}

    def fake(*args, **kwargs):
        seen["args"] = args
        seen["kwargs"] = kwargs
        return "called"

    monkeypatch.setattr(downstream, "visualize_gene_map", fake)

    with pytest.warns(DeprecationWarning, match="visualize_gene_stir"):
        result = downstream.visualize_gene_stir("adata", "GENE", layer="rra_gene_score")

    assert result == "called"
    assert seen["args"] == ("adata", "GENE")
    assert seen["kwargs"] == {"layer": "rra_gene_score"}
