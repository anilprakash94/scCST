# scCST — single-cell Cell State Transition

Local, replicate-aware case–control differential analysis for scRNA-seq. For every cell
(or a set of representative *anchor* cells) and gene, scCST scores how strongly expression
shifts between a condition (e.g. disease) and control **within matched biological-replicate
neighbourhoods**, and summarises each cell with a *divergence score*. It is robust to batch
effects and to cell-type composition differences between conditions.

## Installation

```bash
pip install git+https://github.com/<your-org>/scCST.git         # from GitHub
pip install ./scCST                                             # from a local clone
# optional extras for the downstream helpers (leiden clustering + GO enrichment):
pip install "sccst[downstream] @ git+https://github.com/<your-org>/scCST.git"
```

Requires Python ≥ 3.9 (numpy, scipy, scikit-learn, scanpy, anndata, joblib, matplotlib, seaborn).

## Input

scCST operates on an `AnnData` object with:

- **log-normalised** expression in `adata.X`
- a **latent embedding** in `adata.obsm[latent_key]` (e.g. `X_pca`, or a batch-integrated
  `X_pca_harmony`)
- three `adata.obs` columns: **condition**, **biological batch/replicate** (donor/patient),
  and **cell type**

```python
import scanpy as sc
adata = sc.read_h5ad("your_data.h5ad")
adata                       # inspect obs columns and obsm keys
```

## Pipeline

### Step 1 — run scCST

`process_cell_type()` builds matched case/control neighbourhoods per replicate and runs the
local rank-aggregation (RRA) test for every cell (or anchor).

```python
from sccst import process_cell_type, make_results_adata

cell_importances, gene_names, run_info = process_cell_type(
    adata,
    conditions=("stim", "ctrl"),   # (case, control) — CASE FIRST
    condition_col="label",         # obs column with the two states
    batch_col="batch",             # biological replicate / donor
    cell_type_col="cell_type",     # cell-type labels
    latent_key="X_pca",            # obsm embedding key
    n_neighborhoods=50,
    n_jobs=8,
    # n_anchors_per_celltype=100,  # anchor mode for large datasets (see Notes)
)
```

### Step 2 — assemble results into AnnData

`make_results_adata()` writes the per-cell/gene outputs back into an `AnnData`
(prefix `rra`). With `subset_to_results=True` it returns only the analysed cells/anchors.

```python
result = make_results_adata(
    adata, cell_importances, run_info,
    prefix="rra", subset_to_results=True,
)

result.layers["rra_gene_score"]          # signed gene score (cells x genes)
result.layers["rra_avg_expr_diff"]       # case - control mean expression difference
result.layers["rra_directional_fdr"]     # direction-specific FDR
result.obs["rra_cell_divergence_score"]  # per-cell divergence  (sum of |gene score|)
result.obs["rra_n_sig_genes"]            # per-cell number of significant genes
result.obs["rra_n_pairs"]                # matched replicate pairs used

result.write_h5ad("sccst_result.h5ad")
```

## Downstream analysis

Each analysis is a single function call on the `result` object. Every plotting function
takes `save_path=` (writes a file) and `figsize=`. (Clustering and GO enrichment need the
`sccst[downstream]` extra.)

```python
from sccst.downstream import (
    analyze_global_drivers, plot_ordered_violin, cluster_on_gene_scores,
    plot_cell_type_specific_gene_score_markers, plot_single_celltype_lollipop,
    visualize_gene_stir, run_go_celltype, plot_go_enrichment,
)

# 1. Genes shifting most across the whole dataset (bar plot + table)
drivers = analyze_global_drivers(result, top_n=20, save_path="global_drivers.pdf")

# 2. Divergence score per cell type (sorted, min-max normalized violins)
plot_ordered_violin(result, "rra_cell_divergence_score", "cell_type",
                    normalize=True, save_path="divergence_by_celltype.pdf")

# 3. Cluster cells by their gene-score profile (adds leiden_gene_score_* to obs)
result = cluster_on_gene_scores(result, n_top_genes=2000, resolutions=[0.5, 1.0])

# 4. Cell-type-specific marker genes (mean-gene-score matrix plot)
plot_cell_type_specific_gene_score_markers(result, cell_type_col="cell_type",
                    top_k=10, save_path="celltype_markers.pdf")

# 5. Top up/down driver genes for one cell type (lollipop)
plot_single_celltype_lollipop(result, "CD14+ Monocytes", cell_type_col="cell_type",
                    top_k=20, save_path="lollipop_CD14_Monocytes.pdf")

# 6. Map one gene's disease-relevance across cells
visualize_gene_stir(result, "ISG15", condition_col="label",
                    layer="rra_gene_score", save_path="ISG15_map.pdf")

# 7. GO enrichment per cell type (needs sccst[downstream] + internet)
go = run_go_celltype(result, cell_type_col="cell_type")
plot_go_enrichment(go, save_path="go_enrichment.pdf")
```

## Same two steps from the command line

```bash
# Step 1
python scripts/run_sccst.py --input your_data.h5ad --output sccst_result.h5ad \
    --condition-col label --case stim --control ctrl \
    --batch-col batch --cell-type-col cell_type --latent-key X_pca --n-jobs 8
    # add: --n-anchors-per-celltype 100   for large datasets

# Step 2
python scripts/run_downstream.py --input sccst_result.h5ad --output sccst_annotated.h5ad
```

## Outputs

| where | key | meaning |
|---|---|---|
| `layers` | `rra_gene_score` | signed per-cell, per-gene score (effect size × significance) |
| `layers` | `rra_avg_expr_diff` | case − control mean expression difference |
| `layers` | `rra_directional_fdr` | direction-specific FDR |
| `obs` | `rra_cell_divergence_score` | per-cell divergence = Σ\|gene score\| |
| `obs` | `rra_n_sig_genes` | number of significant genes for the cell |
| `obs` | `rra_n_pairs` | matched control/disease replicate pairs used |

## Notes

- **Large datasets:** set `n_anchors_per_celltype=100` in Step 1 to score ~100
  representative anchors per cell type instead of every cell (a 1.26M-cell dataset runs in
  ~35 min this way). Results are then per-anchor rather than per-cell.
- **Case order matters:** `conditions[0]` is the case/disease state; the gene score is
  case − control.
- **`batch_col`** must be the biological replicate (donor/patient/sample) that scCST matches
  across conditions — not a purely technical label.
- **No latent embedding / raw counts?** Compute `sc.pp.normalize_total` + `sc.pp.log1p`, then
  `sc.pp.pca` (and, if there are batch effects, a batch-integrated embedding such as Harmony
  into `X_pca_harmony`), and point `latent_key` at it.

## Examples

Worked, end-to-end notebooks on the Kang IFN-β PBMC data are in [`examples/`](examples/):
`kang_ifnb_run_sccst.ipynb` (Steps 1–2) and `kang_ifnb_downstream_analysis.ipynb`
(downstream).

## Citation

If you use scCST, please cite _(citation to be added)_.

## License

MIT — see [LICENSE](LICENSE).
