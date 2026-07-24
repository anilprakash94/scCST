# Method overview

scCST quantifies local phenotype-associated expression shifts by comparing matched same-cell-type neighborhoods across biological replicates. For each focal cell or selected anchor, cells from each batch are searched in a latent representation to construct batch-specific local neighborhoods. Control and case neighborhoods are matched by sample identity or by latent centroid distance and filtered using a cell-type-specific caliper. Disease-control mean expression differences are computed for each matched replicate pair, converted into directional partial ranked gene lists, and aggregated using an RRA-style order-statistic procedure. The method reports signed gene scores and cell divergence scores that summarize local disease-associated transcriptional shifts.

See `scCST_full_context_README.md` for the complete method notes used to prepare this repository.
