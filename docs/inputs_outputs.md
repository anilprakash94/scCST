# Inputs and outputs

## Inputs

scCST expects an `AnnData` object with normalized/log-transformed expression in `adata.X`, a latent representation in `adata.obsm[latent_key]`, cell-type labels, condition labels, and biological batch labels. Optional sample-pair labels can be used for strict paired matching.

## Outputs

The main gene-level output is a signed gene score for each analyzed focal cell or anchor and gene. The main cell-level output is the cell divergence score, calculated as the sum of absolute signed gene scores for a cell or anchor.

Common AnnData outputs include:

```python
adata.layers["rra_gene_score"]
adata.layers["rra_avg_expr_diff"]
adata.layers["rra_directional_fdr"]
adata.obs["rra_cell_divergence_score"]
adata.obs["rra_n_sig_genes"]
adata.obs["rra_n_pairs"]
```
