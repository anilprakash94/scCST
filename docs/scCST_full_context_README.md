# scCST README for New Chats

## Method name

**scCST**  
**Full name:** single-cell Cell State Transition  
**Manuscript title:** *Quantifying phenotype-associated cell state differences from case-control single-cell data with scCST*  
**Full descriptive method phrase:** latent-neighborhood matched biological-replicate rank aggregation for local condition-associated expression shifts.

---

## 1. Purpose of this README

This README is intended to be uploaded into new chats or shared with manuscript-writing agents so they understand the scCST method, the statistical design, the expected manuscript framing, the benchmarking strategy, and the real-data results already generated.

The key point is that scCST is not a conventional whole-cell-type differential-expression method and not a cell-level pseudoreplicated test. It is a local, replicate-aware method for finding genes and cell states that show consistent phenotype-associated expression shifts in matched case-control neighborhoods.

---

## 2. Short summaries

**One sentence:** scCST identifies local phenotype-associated cell-state shifts by comparing matched same-cell-type neighborhoods across biological replicate batches in latent space and aggregating replicate-level directional expression differences into signed gene scores and cell divergence scores.

**One paragraph:** scCST is a method for replicated case-control single-cell RNA-seq studies. For each focal cell or selected anchor cell within a cell type, scCST builds same-cell-type neighborhoods separately within each biological batch, matches control and disease neighborhoods in a latent representation, computes disease-control mean expression differences for each matched biological-replicate pair, converts these differences into up- and down-regulated ranked gene lists, and aggregates the ranked lists using an RRA-style procedure. The output is a signed gene score for each anchor-gene pair and a cell divergence score for each anchor. Positive gene scores indicate genes locally increased in disease/case neighborhoods, negative scores indicate genes locally decreased in disease/case neighborhoods, and cell divergence scores summarize the magnitude of local phenotype-associated transcriptional deviation.

---

## 3. Recommended terminology

| Concept | Preferred term |
|---|---|
| Method | scCST |
| Full name | single-cell Cell State Transition |
| Long technical description | latent-neighborhood matched biological-replicate rank aggregation |
| Input design | replicated case-control scRNA-seq |
| Local analysis unit | focal cell or focal anchor |
| Statistical replicate unit | biological batch, sample, or matched sample pair |
| Gene-level output | signed gene score |
| Cell-level output | cell divergence score |
| Main biological interpretation | local phenotype-associated cell-state shift |
| Main statistical interpretation | local replicate-aware rank-aggregation evidence |

Avoid saying that scCST performs ordinary cell-level differential expression. It does not treat each cell as an independent replicate.

---

## 4. Scientific motivation

Single-cell RNA-seq can resolve disease-associated transcriptional changes across cell types, subtypes, activation states, and continuous cellular manifolds. However, standard differential-expression approaches often either compare individual cells directly or aggregate all cells within a broad cell type.

Direct cell-level testing can inflate significance if cells from the same donor or batch are treated as independent biological replicates. Pseudobulk and mixed-model methods improve replicate-aware inference by aggregating cells by sample and cell type or modeling sample-level effects, but whole-cell-type aggregation can dilute local expression shifts that occur only in a subset of a cell type or along a continuous cell-state trajectory.

scCST addresses this gap by performing local comparisons around focal anchors while preserving the biological replicate as the unit of evidence. It tests whether genes are consistently shifted in disease/case neighborhoods relative to matched control neighborhoods across biological replicate pairs.

---

## 5. Primary scientific question

For a given focal cell or anchor within a cell type, scCST asks:

> Which genes show consistent local disease-control expression shifts across matched biological-replicate neighborhoods?

scCST does not primarily ask which genes are globally differentially expressed across an entire cell type. Instead, it asks whether case cells near a local cell state are transcriptionally shifted relative to matched control cells, and which genes drive that local shift.

---

## 6. Inputs required by scCST

scCST expects an annotated single-cell object, usually an `AnnData` object, with:

1. **Expression matrix:** `adata.X`, usually normalized/log-transformed for analysis. scCST computes mean expression differences between local neighborhoods.
2. **Latent representation:** `adata.obsm[latent_key]`, for example `X_pca_harmony`, PCA, Harmony-corrected PCA, scVI, or another cell-state embedding.
3. **Cell-type labels:** for example `adata.obs["Cell_Type"]`. scCST only compares cells within the same annotated cell type.
4. **Condition labels:** for example `adata.obs["condition"]`. Two conditions are expected, such as disease and control. In the implementation, control is encoded as 0 and disease/case as 1.
5. **Batch labels:** for example `adata.obs["batch"]` or `adata.obs["sim_batch"]`. Each batch must belong to only one condition.
6. **Optional sample pairing labels:** for example `adata.obs["sample"]`. If provided, scCST enforces strict paired matching, where each sample must have exactly one control batch and one disease/case batch.

---

## 7. Main assumptions

scCST is most appropriate when the dataset has multiple biological replicates per condition, each batch belongs to only one condition, matched case-control local neighborhoods are biologically comparable in latent space, cell-type annotations are sufficiently accurate, and the latent representation captures meaningful cell-state structure.

The ranked lists used for aggregation should be approximately independent because they come from different biological-replicate pairs. There must also be enough cells per batch-specific neighborhood and enough matched pairs per anchor.

Important caveat: because the latent space is derived from expression, scCST estimates disease-associated shifts conditional on local similarity in that latent representation. This is useful for local matching but should not be described as an unconditional global DE test.

---

## 8. Core algorithm

### Step 1. Gene filtering

Genes detected in too few cells are removed before testing. In the provided implementation, the minimum-cell threshold is 0.5% of cells, with a lower bound of three cells.

### Step 2. Encode metadata

scCST encodes cell types, conditions, batches, and optional sample-pairing metadata. The method validates that each batch contains cells from only one condition.

### Step 3. Compute cell-type-specific local latent calipers

For each cell type, scCST estimates a latent-distance caliper from same-cell-type nearest-neighbor distances. The caliper is used to reject case-control neighborhood matches that are too far apart in latent space.

### Step 4. Choose focal cells or anchors

scCST supports all-cell mode and anchor mode. In all-cell mode, every cell is used as a focal cell. In anchor mode, a subset of representative anchors is used. Anchors can be manually supplied or automatically selected by farthest-point sampling within each cell type. Anchor mode reduces redundant testing of highly overlapping neighborhoods.

### Step 5. Build batch-specific neighborhoods

For each focal anchor and each biological batch, scCST selects the nearest same-cell-type cells to the focal anchor in latent space. This creates one local neighborhood per batch. The focal cell is excluded from its own batch-specific neighborhood when applicable.

### Step 6. Match control and disease/case neighborhoods

In strict paired mode, if `sample_col` is supplied, scCST pairs control and disease/case batches by sample identity. Each sample must contain exactly one control batch and one disease/case batch.

In distance-based unpaired mode, if no `sample_col` is supplied, scCST performs one-to-one optimal bipartite matching between control and disease/case batch neighborhoods using latent centroid distances.

### Step 7. Apply caliper filtering

Matched neighborhood pairs are retained only if their centroid distance is below the cell-type-specific latent caliper.

### Step 8. Compute replicate-level expression differences

For each matched biological-replicate neighborhood pair and gene:

```text
Delta[p, g] = mean(case neighborhood p, gene g) - mean(control neighborhood p, gene g)
```

where `p` is the matched pair and `g` is the gene.

### Step 9. Build directional partial ranked lists

For each matched pair, scCST creates an up list and a down list. The up list contains genes with positive case-control differences, ranked from largest positive difference downward. The down list contains genes with negative case-control differences, ranked from most negative difference upward.

Only the top `K` genes per direction are retained by default. Omitted genes receive normalized rank 1.0.

### Step 10. Perform RRA-style aggregation

scCST applies directional RRA-style rank aggregation separately for up and down lists. The implementation computes beta order-statistic rho scores and converts them into conservative Bonferroni-style RRA p-values.

Important wording: say “RRA-style” or “conservative RRA rho p-values.” Do not claim exact `RobustRankAggreg` package p-values unless that exact implementation is used.

### Step 11. Multiple testing

Benjamini-Hochberg FDR correction is applied within each anchor and direction.

Important caveat: the current implementation controls local anchor-level directional FDR. If reporting significant genes across all anchors, additional global or hierarchical correction may be needed.

### Step 12. Compute gene scores and cell divergence scores

For each anchor-gene pair, scCST computes:

```text
S[g] = average_delta[g] * -log10(FDR[g] + epsilon)
```

where `average_delta[g]` is the average disease-control expression difference across matched pairs.

The cell divergence score is:

```text
D[cell] = sum(abs(S[cell, gene])) across genes
```

---

## 9. Main outputs

### Gene-level outputs

The gene-level outputs are `avg_expr_diff`, `directional_fdr`, and `gene_score`.

### Cell-level or anchor-level outputs

The cell-level (or anchor-level) outputs are `cell_divergence_score`, `n_sig_genes`, and `n_pairs`.

### AnnData integration

Gene-level matrices are usually stored in layers, for example:

```python
adata.layers["rra_gene_score"]
adata.layers["rra_avg_expr_diff"]
adata.layers["rra_directional_fdr"]
```

Cell-level metrics are usually stored in obs, for example:

```python
adata.obs["rra_cell_divergence_score"]
adata.obs["rra_n_sig_genes"]
adata.obs["rra_n_pairs"]
```

---

## 10. Downstream analyses of scCST outputs

### 10.1 Clustering on gene scores

Cells can be clustered using the scCST gene-score matrix rather than raw expression. This identifies groups of cells with similar phenotype-associated effect profiles.

Workflow:

1. Extract `adata.layers["rra_gene_score"]`.
2. Replace NaNs with zero.
3. Select top genes by variance of gene scores across cells.
4. Run PCA on the selected gene-score matrix.
5. Build a neighborhood graph.
6. Run Leiden clustering at one or more resolutions.

This creates clusters based on disease-associated gene-score patterns, not baseline expression similarity.

### 10.2 Visualizing cell divergence score

Cell divergence scores can be min-max scaled and plotted on UMAP to map the magnitude of disease-associated local transcriptional deviation across the cell-state manifold.

### 10.3 Global gene-score drivers

Global disease-associated drivers can be identified by averaging signed gene scores across cells. Genes with large positive mean scores are globally increased in disease/case local neighborhoods. Genes with large negative mean scores are globally decreased.

Penetrance can be computed as the fraction of cells with large absolute gene score, for example `abs(gene_score) > 1.96`.

### 10.4 Cell-type-specific gene-score markers

Cell-type-specific gene-score markers can be identified by comparing each cell type against the strongest competing other cell type.

```python
ranking_score = target_mean - max_other_mean
```

where `target_mean` is the mean gene score in the target cell type and `max_other_mean` is the maximum mean gene score across all other eligible cell types.

### 10.5 One-cell-type marker analysis

For a selected target cell type, the same formula can be used. If two marker functions give different genes, check that the same layer, same cell-type column, same `top_k`, same `min_cells`, and same target label are used. Compare `marker_dict[target_cell_type]`, not the full dictionary or a matrixplot containing markers for multiple cell types.

### 10.6 Gene ontology enrichment

GO enrichment can be performed per cell type using top positively scored genes. In the downstream analysis, genes were selected based on directional FDR significance in a minimum number of target cells, mean target score greater than mean rest score, and top positive score. GO enrichment was run using Enrichr through `gseapy`, for example with `GO_Biological_Process_2023`.

---

## 11. Simulated benchmark datasets

### 11.1 Simulation 1: heterogeneous responder and batch-artifact benchmark

Purpose: test whether scCST recovers local cell-type-specific disease marker genes in responder disease cells despite dropout, responder heterogeneity, and batch-specific artifacts.

Main design:

- 1,200 baseline cells.
- 2,000 genes.
- Counts generated from gene-specific negative-binomial distributions.
- Dropout rate: 0.4.
- Cells divided into groups of 120 cells.
- Each group became a simulated cell type.
- For each cell type, a control subset was copied to create a disease subset.
- Identity genes were increased in both control and disease cells to preserve cell-type structure.
- Disease marker genes were injected only into disease responder cells.
- Responders were sampled within disease batches so response was distributed across disease batches.
- Batch-specific artifact genes were injected independently into each batch.
- Final dataset had 2,400 cells and 2,000 genes.

Important note: the code comment says 60% responders, but the code used 80% responders per disease batch.

Preprocessing:

- Normalize total counts.
- Log-transform.
- Select highly variable genes.
- PCA.
- Harmony integration by simulated batch.
- Neighbors and UMAP using `X_pca_harmony`.

Benchmarking:

- scCST was run using the Harmony-corrected latent space.
- Marker-set boxplots showed increased gene scores for the correct marker genes in the corresponding cell type.
- Heatmaps showed block/diagonal structure across cell types.
- scCST gene scores were compared against CACOA and miloDE z-scores.
- scCST achieved higher cell-macro and gene-macro AUPRC in the simulated disease-cell benchmark.

AUPRC metrics:

- Cell-macro AUPRC: for each responder disease cell, test whether true marker genes are ranked above other genes.
- Gene-macro AUPRC: for each true marker gene, test whether the correct responder disease cells receive higher scores than other disease cells.

### 11.2 Simulation 2: graded cell-divergence benchmark

Purpose: test whether the cell divergence score increases with the magnitude of the simulated disease perturbation.

Main design:

- 2,400 baseline cells.
- 2,000 genes.
- Negative-binomial counts with mean 2.0 and overdispersion 2.0.
- Cells divided into groups of 240 cells.
- Each group became a simulated cell type.
- Control subset was copied to create disease subset.
- Four control batches and four disease batches.
- Low-level shared noise added.
- Each cell type had 100 identity genes increased in both conditions.
- Disease marker genes increased progressively by cell type: CellType_0 had 5 disease marker genes, CellType_1 had 10, CellType_2 had 15, and so on.
- Disease marker counts were increased in disease cells.
- Control cells were set to zero for disease marker genes.

Results:

- Boxplots showed higher scCST gene scores for marker genes in the corresponding simulated cell type.
- Heatmaps showed cell-type-specific block structure.
- Cell divergence score strongly tracked the number of disease marker genes per cell type.
- Spearman correlation between normalized divergence score and number of disease marker genes was approximately 0.995.

Interpretation: scCST gene scores recovered the correct local disease markers, and scCST cell divergence score quantitatively reflected the magnitude of disease-associated transcriptional burden.

---

## 12. Real datasets analyzed

### 12.1 IFN-beta-stimulated PBMC dataset

Dataset: IFN-beta-stimulated PBMC vs control from Kang et al. This is often used as a benchmark for perturbation and integration analyses.

Original finding: IFN-beta induced broad transcriptomic shifts across immune cells, including shared antiviral-response genes and more cell-type-specific programs in myeloid and lymphoid populations. Myeloid cells showed strong IFN-beta response signatures, including antiviral response, chemokine signaling, and systemic lupus erythematosus pathways.

scCST results:

- Cell divergence scores were highest in myeloid cells.
- CD14+ monocytes had the highest divergence.
- FCGR3A+ monocytes and dendritic cells also showed high divergence.
- Lower divergence was seen in CD4 T cells, CD8 T cells, NK cells, and B cells.

Top global upregulated genes included ISG15, IFI6, ISG20, IFIT3, IFIT1, LY6E, MX1, CXCL10, IFIT2, TNFSF10, RSAD2, OAS1, IFITM3, and IRF7.

GO enrichment showed defense response to virus, defense response to symbiont, negative regulation of viral process, negative regulation of viral genome replication, response to interferon-beta, cellular response to cytokine stimulus, regulation of cytokine production, and translation-related processes in some lymphoid populations.

Interpretation: scCST recovered canonical IFN-beta response programs and showed that myeloid populations had the strongest treatment-associated transcriptomic shifts.

### 12.2 Cutaneous squamous cell carcinoma dataset

Dataset: cutaneous squamous cell carcinoma from Ji et al., Cell 2020, with tumor and matched normal skin.

Original finding: normal epidermal keratinocytes separate into basal, cycling, and differentiating states. Tumor keratinocytes largely recapitulate these programs but also contain a cancer-specific tumor-specific keratinocyte population, TSK. Key TSK markers include MMP10, VIM, ITGA5, ITGB1, PLAU, MMP9, TNC, SERPINE2, CALM1, and TGFB1. Tumor keratinocytes show altered metabolic and stress programs including glycolysis and hypoxia genes.

scCST results:

- Cell divergence scores were high in fibroblasts and keratinocyte populations.
- High-scoring keratinocyte populations included KC_Diff, KC_Cyc, KC_Basal, and TSK.
- Global upregulated genes included KRT6B, KRT6A, S100A9, FABP5, GJB2, IFITM3, S100A8, S100A2, SERPINB3, TYMP, TMSB10, SERPINB4, CTSB, TXNDC17, GJB6, KRT6C, S100A7, CLCA2, and CTSC.
- Global downregulated genes included KRT1, KRT10, DUSP1, EGR1, FOS, FOSB, DMKN, CCL27, KLF4, JUNB, and ATF3.
- Gene-score Leiden cluster 18 recapitulated the TSK signature.
- GO enrichment showed epithelium development, epidermis development, intermediate filament organization, skin development, peptide cross-linking, extracellular matrix organization, collagen fibril organization, and glycolytic process.

Interpretation: scCST recovered tumor-associated keratinocyte programs, stromal remodeling, and a TSK-like state from local tumor-normal expression shifts.

### 12.3 Pulmonary fibrosis / IPF dataset

Dataset: pulmonary fibrosis / idiopathic pulmonary fibrosis from Habermann et al., Science Advances 2020. It contains 114,396 cells from 20 pulmonary fibrosis lungs and 10 non-fibrotic control lungs.

Original finding: major IPF/PF markers examined at single-cell resolution included MUC5B, MMP7, FN1, COL1A1, CDKN2A, and SMAD3.

Related finding: Morse et al. reported that highly proliferative SPP1-high macrophages contribute to activation of IPF myofibroblasts in lung fibrosis.

scCST results:

- Global upregulated genes included FN1, CRIP1, S100A10, CCL18, SPP1, SEPW1, SCGB3A1, GSN, MT2A, PPDPF, TIMP1, LGALS1, TYMP, EMP3, SERPINA1, HCST, TMSB4X, and CALM1.
- Global downregulated genes included C1QC, CD163, C1QA, C1QB, VSIG4, SLPI, APOD, GPX3, TXNIP, and JUN.
- High divergence was seen in fibroblasts, MUC5B+ epithelial cells, transitional AT2 cells, macrophages, AT2 cells, and endothelial cells.
- scCST identified a macrophage subpopulation with increased shift and high SPP1, interpreted as an SPP1+ profibrotic macrophage state.
- In the macrophage-associated gene-score driver plot, SPP1 was the strongest positive driver. Other positive drivers included CCL18, CSTB, FN1, C15orf48, PPDPF, HCST, LILRB4, FABP5, and CCL2.
- GO enrichment showed antigen processing and presentation via MHC class II, neutrophil migration, granulocyte chemotaxis, neutrophil chemotaxis, positive regulation of cell adhesion, regulation of epithelial cell proliferation, regulation of angiogenesis, and blood vessel morphogenesis.

Interpretation: scCST recovered profibrotic markers, identified shifted epithelial and stromal populations, and highlighted SPP1-high profibrotic macrophages consistent with known IPF biology.

### 12.4 Systemic lupus erythematosus 1.2 million PBMC dataset

Dataset: systemic lupus erythematosus PBMCs from Perez et al., Science 2022, containing more than 1.2 million PBMCs from 162 SLE cases and 99 controls.

Special mode: because the dataset is very large, scCST used evenly spaced anchor cells instead of all cells.

Original finding: classical monocytes had the strongest expression of pan-cell-type and myeloid-specific interferon-stimulated genes. Representative/global SLE markers included ISG15, IFI27, IFI6, IFI44L, LY6E, HLA-C, MT2A, IFITM2, PSME2, and PSMB9. The myeloid-specific module was upregulated mainly in myeloid cells, especially classical monocytes, and included IFITM1, IFITM3, APOBEC3A, RNASE2, and IFIT2.

scCST divergence results:

- Classical monocytes showed the highest cell divergence scores.
- Non-classical monocytes and conventional dendritic cells also showed high divergence.
- Natural killer cells, plasmacytoid dendritic cells, B cells, and CD8-positive alpha-beta T cells had intermediate divergence.
- Plasmablasts and progenitor cells had low divergence.

Global upregulated genes included ISG15, IFITM3, IFI6, LY6E, HLA-A, IFI44L, HLA-C, JUN, EPSTI1, MT2A, HLA-B, PSME2, NFKBIA, BST2, TYMP, HLA-DRB5, CXCR4, CD69, KLF6, and PSMB9.

GO enrichment showed defense response to virus, defense response to symbiont, negative regulation of viral process, negative regulation of viral genome replication, positive regulation of cytokine production, positive regulation of inflammatory response, cellular defense response, natural killer cell mediated immunity, and translation-related terms in some immune populations.

Gene-score cluster result:

- scCST recapitulated the myeloid-specific interferon module.
- Gene-score Leiden cluster 17 showed strong positive scores for interferon-stimulated and myeloid-associated genes.
- Genes included IFITM3, IFI6, LY6E, MT2A, ISG15, NFKBIA, TNFSF13B, IFIT2, and VAMP5.

Candidate novel marker result:

- scCST identified PATL2 as a potential disease-associated marker in a terminal effector CD8-positive alpha-beta T-cell state.
- CD8-positive alpha-beta T-cell subcluster 10 showed PATL2 signal.
- Dot plots showed PATL2 enriched in SLE compared with normal PBMCs.
- Sex-stratified plots showed strongest PATL2 signal in SLE male samples.
- This should be described as a candidate or potential novel marker, not as a validated marker unless follow-up validation is performed.

Interpretation: scCST recovered known SLE interferon biology, identified high divergence in classical monocytes, recapitulated a myeloid-specific ISG module, and highlighted a candidate PATL2-positive CD8 T-cell substate.

---

## 13. Manuscript structure suggestions

### Introduction: five-paragraph structure

1. Single-cell RNA-seq in disease studies and the need for replicated case-control design.
2. Existing tools for identifying important genes between case and control cells, including Wilcoxon/Seurat tests, MAST, edgeR, DESeq2, limma-voom, muscat, and NEBULA.
3. Shortcomings of existing methods, including cell-level pseudoreplication, dilution of local effects in pseudobulk, dependence on cluster annotation granularity, and difficulty detecting localized cell-state shifts.
4. Neighborhood-based methods such as Milo, CACOA, and miloDE, including their advantages and remaining challenges.
5. scCST, emphasizing local latent-neighborhood matching, biological-replicate comparisons, directional RRA-style aggregation, gene scores, cell divergence scores, and anchor mode.

### Methods section structure

Suggested sections:

1. Overview of scCST.
2. Input data and preprocessing.
3. Cell-type-stratified analysis.
4. Focal-cell and anchor selection.
5. Cell-type-specific latent calipers.
6. Batch-specific local neighborhood construction.
7. Case-control neighborhood matching.
8. Replicate-level expression differences.
9. Directional partial ranked lists.
10. RRA-style rank aggregation.
11. Multiple-testing correction.
12. Gene scores and cell divergence scores.
13. AnnData output structure.
14. Downstream gene-score analyses.
15. Simulated datasets.
16. Real datasets.
17. Benchmarking against CACOA and miloDE.
18. Statistical interpretation and limitations.

### Results section structure

Suggested sections:

1. scCST method overview.
2. Simulation 1: scCST recovers local responder disease markers and outperforms CACOA/miloDE AUPRC.
3. Simulation 2: cell divergence score tracks number of injected disease markers.
4. IFN-beta PBMC: scCST detects myeloid-dominant interferon response.
5. cSCC: scCST identifies fibroblast/keratinocyte shifts and TSK-like programs.
6. IPF: scCST identifies profibrotic programs and SPP1-high macrophages.
7. SLE: anchor-mode scCST recovers classical monocyte interferon module and PATL2-positive CD8 T-cell candidate state.

---

## 14. Recommended Discussion points

Advantages:

1. Replicate-aware: uses biological batches or matched samples as evidence units.
2. Local: detects cell-state-specific shifts that broad pseudobulk can dilute.
3. Directional: separates up and down disease-associated shifts.
4. Interpretable: produces gene scores and cell divergence scores.
5. Scalable: anchor mode enables use on very large datasets.
6. Complementary: can be used alongside pseudobulk, mixed models, CACOA, Milo, and miloDE.
7. Robust to noisy replicate-specific rankings through RRA-style aggregation.

Limitations:

1. Results depend on the chosen embedding/integration method.
2. Current FDR is within anchor and direction, not global across all anchors.
3. RRA p-values are conservative RRA-style rho p-values, not exact RobustRankAggreg package p-values.
4. Requires enough biological replicates and matched neighborhoods.
5. Rare cell types or poorly overlapping case-control states may be skipped.
6. Anchor mode provides representative local measurements, not exhaustive per-cell inference.
7. Candidate novel markers require external validation.

Conclusion:

scCST provides a local, replicate-aware framework for quantifying phenotype-associated cell-state transitions in case-control scRNA-seq. Simulated and real-data analyses show that scCST recovers known disease biology, identifies local disease-associated gene programs, and provides a cell divergence score that summarizes the magnitude of phenotype-associated transcriptional deviation.

---

## 15. Suggested citations to use in manuscript

Verify all citations before submission.

### General scRNA-seq and preprocessing

- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018.
- Virshup I et al. The scverse project provides a computational ecosystem for single-cell omics data analysis. *Nature Biotechnology*. 2023.
- Butler A et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species. *Nature Biotechnology*. 2018.
- Stuart T et al. Comprehensive integration of single-cell data. *Cell*. 2019.

### Differential expression and replicate-aware inference

- Soneson C, Robinson MD. Bias, robustness and scalability in single-cell differential expression analysis. *Nature Methods*. 2018.
- Crowell HL et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nature Communications*. 2020.
- Squair JW et al. Confronting false discoveries in single-cell differential expression. *Nature Communications*. 2021.
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014.
- Robinson MD, McCarthy DJ, Smyth GK. edgeR. *Bioinformatics*. 2010.
- Ritchie ME et al. limma powers differential expression analyses for RNA-seq and microarray studies. *Nucleic Acids Research*. 2015.
- Finak G et al. MAST. *Genome Biology*. 2015.
- He L et al. NEBULA. *Communications Biology*. 2021.

### Rank aggregation

- Kolde R, Laur S, Adler P, Vilo J. Robust rank aggregation for gene list integration and meta-analysis. *Bioinformatics*. 2012.

### Neighborhood-based single-cell methods

- Dann E et al. Differential abundance testing on single-cell data using k-nearest neighbor graphs. *Nature Biotechnology*. 2022.
- CACOA preprint: case-control analysis of single-cell RNA-seq cohorts. bioRxiv 2022.
- miloDE: Leveraging neighborhood representations of single-cell data to achieve sensitive DE testing. 2024.

### Real datasets

- Kang HM et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nature Biotechnology*. 2018.
- Ji AL et al. Multimodal analysis of composition and spatial architecture in human squamous cell carcinoma. *Cell*. 2020.
- Habermann AC et al. Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis. *Science Advances*. 2020.
- Morse C et al. Proliferating SPP1/MERTK-expressing macrophages in idiopathic pulmonary fibrosis. *European Respiratory Journal*. 2019.
- Perez RK et al. Single-cell RNA-seq reveals cell type-specific molecular and genetic associations to lupus. *Science*. 2022.

---

## 16. Safe claim strength

Strong claims:

- scCST identifies local, replicate-consistent phenotype-associated expression shifts.
- scCST maps disease-associated transcriptional divergence across cell-state space.
- scCST recovers known interferon, tumor, fibrotic, and lupus-associated programs in real datasets.
- scCST cell divergence score tracks simulated perturbation magnitude.

Moderate claims:

- scCST identifies candidate disease-associated cell states.
- scCST highlights candidate novel markers such as PATL2 in SLE CD8 T-cell substates.

Avoid unless externally validated:

- scCST proves causal cell-state transitions.
- scCST discovers definitive disease mechanisms without validation.
- scCST globally controls FDR across all cells, anchors, genes, and directions.
- scCST eliminates all confounding from batch, donor, or latent structure.

---

## 17. Minimal pseudo-code

```python
cell_importances, gene_names, run_info = process_cell_type(
    adata=adata,
    conditions=("Disease", "Control"),
    condition_col="condition",
    batch_col="batch",
    sample_col=None,
    n_neighborhoods=20,
    cell_type_col="Cell_Type",
    n_jobs=4,
    latent_key="X_pca_harmony",
    caliper_percentile=75.0,
    min_cells_per_batch_neighborhood=3,
    min_matched_pairs=3,
    top_genes_per_list=100,
    fdr_threshold=0.05,
    n_anchors_per_celltype=None
)

result_adata = make_results_adata(
    adata,
    cell_importances,
    run_info,
    prefix="rra",
    subset_to_results=True
)
```

For large datasets:

```python
cell_importances, gene_names, run_info = process_cell_type(
    adata=adata,
    conditions=("Disease", "Control"),
    condition_col="condition",
    batch_col="batch",
    cell_type_col="Cell_Type",
    latent_key="X_pca_harmony",
    n_anchors_per_celltype=500,
    anchor_random_state=0
)
```

---

## 18. Best short description for a new chat

Use this paragraph when starting a new chat:

> I am developing scCST, single-cell Cell State Transition, a method for replicated case-control scRNA-seq. For each focal cell or selected anchor within a cell type, scCST builds same-cell-type local neighborhoods separately within each biological batch, matches case and control neighborhoods in latent space, computes replicate-level disease-control expression differences, converts these differences into directional ranked gene lists, and aggregates them using an RRA-style procedure. It outputs signed anchor-gene scores and cell divergence scores. The goal is to identify local, replicate-consistent phenotype-associated expression shifts and map disease-associated cell-state transitions.

---

## 19. Best brief methods wording

> scCST quantifies local phenotype-associated expression shifts by comparing matched same-cell-type neighborhoods across biological replicates. For each focal anchor, cells from each batch are searched in a latent representation to construct batch-specific local neighborhoods. Control and case neighborhoods are then matched by sample identity or by latent centroid distance and filtered using a cell-type-specific caliper. Disease-control mean expression differences are computed for each matched replicate pair, converted into directional partial ranked gene lists, and aggregated using an RRA-style order-statistic procedure. The method reports signed gene scores and cell divergence scores that summarize local disease-associated transcriptional shifts.

---

## 20. Best brief Results conclusion

> Across simulated and real datasets, scCST recovered known phenotype-associated gene programs and localized them to biologically relevant cell states. In simulations, scCST recovered the injected cell-type-specific disease markers and its cell divergence score tracked the number of simulated disease genes. In IFN-beta PBMCs and SLE PBMCs, scCST recovered interferon-stimulated gene programs with strongest shifts in myeloid populations. In cSCC, scCST identified keratinocyte and fibroblast shifts and recapitulated a tumor-specific keratinocyte-like program. In IPF, scCST identified profibrotic gene programs and SPP1-high macrophage states.
