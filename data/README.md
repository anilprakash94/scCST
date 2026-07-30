# Data

scCST was evaluated on four published case–control or perturbation scRNA-seq datasets.
None of the raw data are redistributed here; download them from the sources below and
convert each to an `AnnData` (`.h5ad`) object with log-normalised `X`, a latent embedding
in `obsm` (e.g. `X_pca`), and `obs` columns for condition, biological replicate/donor, and
cell type (see the top-level [README](../README.md) for the expected input format).

| Dataset | Design | Condition (case vs control) | Source |
|---|---|---|---|
| **Kang et al.** — IFN-β-stimulated PBMCs | Multiplexed droplet scRNA-seq of cell-type-specific and inter-individual responses to IFN-β | `stim` vs `ctrl` | [figshare 19397624](https://doi.org/10.6084/m9.figshare.19397624) |
| **Ji et al.** — cutaneous squamous cell carcinoma (cSCC) | Single-cell, spatial, and imaging assays of cSCC tumors and matched normal skin | tumor vs normal skin | [GEO GSE144236](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144236) |
| **Habermann et al.** — pulmonary fibrosis | 114,396 cells from 20 pulmonary fibrosis lungs and 10 non-fibrotic controls | fibrotic vs control lung | [GEO GSE135893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893) |
| **Perez et al.** — systemic lupus erythematosus (SLE) PBMCs | >1.2 million PBMCs from 162 SLE cases and 99 healthy controls | SLE vs healthy | [CELLxGENE collection](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2) |
