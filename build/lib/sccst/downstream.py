"""Downstream analysis and plotting utilities extracted from the scCST analysis notebook.

These functions are optional helpers and are not required to run the core scCST method.
"""


# %% Extracted from sc_stir_analysis_kang notebook cell 4

def cluster_on_gene_scores(adata, n_top_genes=2000, resolutions=[0.5, 1.0, 1.5]):
    """
    Clusters cells using the 'gene_score' layer, applying custom 
    feature selection to isolate the most dynamic effect sizes.
    """
    print("Preparing data for clustering based on 'gene_score'...")
    
    # 1. Extract the gene_score layer and handle NaNs
    gene_scores = adata.layers['rra_gene_score'].copy()
    gene_scores = np.nan_to_num(gene_scores, nan=0.0)
    
    # 2. Create the temporary AnnData object
    temp_adata = sc.AnnData(X=gene_scores, obs=adata.obs.copy())
    temp_adata.var_names = adata.var_names # Ensure gene names match
    
    # ==========================================
    # NEW: Custom Feature Selection
    # ==========================================
    print(f"Selecting the top {n_top_genes} genes by gene_score variance...")
    
    # Calculate variance for each gene across all cells
    gene_variances = np.var(temp_adata.X, axis=0)
    
    # Get the indices of the top N most variable genes
    top_indices = np.argsort(gene_variances)[::-1][:n_top_genes]
    
    # Subset the temporary AnnData to only these highly dynamic genes
    temp_adata = temp_adata[:, top_indices].copy()
    print(f"Filtered down to {temp_adata.n_vars} highly variable effect genes.")
    # ==========================================
    
    # 3. Compute PCA on the FILTERED gene scores
    print("Computing PCA...")
    sc.tl.pca(temp_adata, svd_solver='arpack')
    
    # 4. Compute neighborhood graph
    print("Computing neighborhood graph...")
    sc.pp.neighbors(temp_adata)
    
    # 5. Run Leiden clustering
    for res in resolutions:
        print(f"Running Leiden clustering at resolution {res}...")
        col_name = f'leiden_gene_score_{res}'
        sc.tl.leiden(temp_adata, resolution=res, key_added=col_name)
        adata.obs[col_name] = temp_adata.obs[col_name]
        
    print("Clustering complete.")
    return adata


# %% Extracted from sc_stir_analysis_kang notebook cell 10

import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

def plot_ordered_violin(
    adata,
    score_key,
    celltype_col,
    method="mean",
    save_path=None,
    figsize=(16, 8),
    normalize=False,
):
    """
    Sorts cell types by the mean or median of a score and generates a violin plot
    using the same colors as sc.pl.umap(..., color=celltype_col).
    """

    # Optionally min-max normalize the score to [0, 1] before plotting.
    if normalize:
        _v = adata.obs[score_key].astype(float)
        _rng = (_v.max() - _v.min()) or 1.0
        score_key = f"{score_key}_normalized"
        adata.obs[score_key] = (_v - _v.min()) / _rng

    # 1. Calculate statistics and get the sorted order
    stats = (
        adata.obs.groupby(celltype_col, observed=False)[score_key]
        .agg(method)
        .sort_values(ascending=False)
    )
    new_order = stats.index.tolist()

    # 2. Ensure UMAP colors exist for this categorical column
    color_key = f"{celltype_col}_colors"

    if color_key not in adata.uns:
        # This forces Scanpy to create/store the palette in adata.uns[color_key]
        sc.pl.umap(
            adata,
            color=celltype_col,
            show=False
        )
        plt.close()

    # 3. Get the original category order and colors
    original_categories = list(adata.obs[celltype_col].cat.categories)
    original_colors = list(adata.uns[color_key])

    color_map = dict(zip(original_categories, original_colors))

    # 4. Re-cast the column as categorical with the NEW order
    adata.obs[celltype_col] = pd.Categorical(
        adata.obs[celltype_col],
        categories=new_order,
        ordered=True
    )

    # 5. Reorder colors to match the new category order
    adata.uns[color_key] = [
        color_map[celltype] for celltype in new_order
    ]

    # 6. Plotting
    fig, ax = plt.subplots(figsize=figsize)

    sc.pl.violin(
        adata,
        keys=score_key,
        groupby=celltype_col,
        stripplot=False,
        order=new_order,
        ax=ax,
        show=False
    )

    # 7. Refine aesthetics
    plt.xticks(rotation=45, ha="right", fontsize=8)

    title_str = f"{method.capitalize()} {score_key.replace('_', ' ')}"
    plt.title(title_str, pad=20)

    # 8. Save and cleanup
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.close()


# %% Extracted from sc_stir_analysis_kang notebook cell 12

def analyze_global_drivers(adata, layer='rra_gene_score', top_n=20, save_path=None, figsize=(12, 6)):
    """
    Identifies genes with the highest/lowest average t_stat-scores across the entire dataset.
    """
    print("--- Computing Global Drivers ---")
    
    # 1. Get the matrix (ignoring NaNs from skipped cells)
    t_stat_matrix = adata.layers[layer]
    
    # 2. Compute stats ignoring NaNs
    global_means = np.nanmean(t_stat_matrix, axis=0)
    
    # Penetrance: % of total cells with significant score
    is_sig = np.abs(t_stat_matrix) > 1.96
    penetrance = np.nanmean(is_sig, axis=0)
    
    # 3. Create DataFrame
    df_global = pd.DataFrame({
        'gene': adata.var_names,
        'mean_score': global_means,
        'penetrance': penetrance
    }).set_index('gene')
    
    # 4. Get Top Positive and Negative
    top_pos = df_global.nlargest(top_n, 'mean_score')
    top_neg = df_global.nsmallest(top_n, 'mean_score').iloc[::-1] 
    
    # 5. Plotting (FIXED: Added hue and legend=False)
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # Positive Plot
    sns.barplot(
        data=top_pos, 
        x='mean_score', 
        y=top_pos.index, 
        hue=top_pos.index,  # Fix: Assign y to hue
        legend=False,       # Fix: Hide redundant legend
        ax=axes[0], 
        palette='Reds_r'
    )
    axes[0].set_title(f"Top {top_n} Upregulated in Disease (Global)")
    axes[0].set_xlabel("Mean Signed Score")
    
    # Negative Plot
    sns.barplot(
        data=top_neg, 
        x='mean_score', 
        y=top_neg.index, 
        hue=top_neg.index,  # Fix: Assign y to hue
        legend=False,       # Fix: Hide redundant legend
        ax=axes[1], 
        palette='Blues_r'
    )
    axes[1].set_title(f"Top {top_n} Downregulated in Disease (Global)")
    axes[1].set_xlabel("Mean Signed Score")
    
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return df_global


# %% Extracted from sc_stir_analysis_kang notebook cell 14

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import warnings

from scipy.sparse import issparse


def plot_cell_type_specific_gene_score_markers(
    adata,
    cell_type_col="Cell_Type",
    layer="rra_gene_score",
    top_k=10,
    min_cells=10,
    clip=None,
    vmin=-5,
    vmax=5,
    cmap="RdBu_r",
    figsize=None,
    save_path=None,
    is_dendrogram=True,
    plot_cell_types=None,
):
    """
    Plot cell-type-specific gene-score markers.

    Marker genes are identified by comparing each cell type against the
    strongest competing other cell type.

    Ranking score:

        target_mean - max_other_mean

    No filtering is applied. The top_k genes per cell type are selected only
    by this ranking score.

    If plot_cell_types is provided, only those cell types are shown in the plot,
    but marker selection is still performed using all eligible cell types.
    """

    print(f"--- Finding and plotting cell-type-specific markers from `{layer}` ---")

    if layer not in adata.layers:
        raise ValueError(f"Layer `{layer}` not found in adata.layers.")

    if cell_type_col not in adata.obs:
        raise ValueError(f"Column `{cell_type_col}` not found in adata.obs.")

    valid_mask = adata.obs[cell_type_col].notna().values
    cell_types = adata.obs.loc[valid_mask, cell_type_col].astype(str)

    counts = cell_types.value_counts()
    keep_cell_types = counts[counts >= min_cells].index.astype(str).tolist()

    if len(keep_cell_types) < 2:
        raise ValueError("Need at least two cell types with enough cells.")

    valid_mask = (
        valid_mask
        & adata.obs[cell_type_col].astype(str).isin(keep_cell_types).values
    )

    X = adata.layers[layer][valid_mask, :]

    if issparse(X):
        X = X.toarray()

    X = np.asarray(X, dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    ct_values = adata.obs.loc[valid_mask, cell_type_col].astype(str).values
    categories = pd.Index(sorted(pd.unique(ct_values)))

    codes = pd.Categorical(ct_values, categories=categories).codes

    n_cells = X.shape[0]
    n_types = len(categories)

    design = np.zeros((n_cells, n_types), dtype=np.float32)
    design[np.arange(n_cells), codes] = 1.0

    group_counts = design.sum(axis=0)[:, None]

    mean_scores = (design.T @ X) / group_counts

    max_other_mean = np.empty_like(mean_scores)

    for i in range(n_types):
        others = np.arange(n_types) != i
        max_other_mean[i, :] = mean_scores[others, :].max(axis=0)

    ranking_scores = mean_scores - max_other_mean

    marker_score_df = pd.DataFrame(
        ranking_scores,
        index=categories,
        columns=adata.var_names,
    )

    marker_dict = {}

    for ct in categories:
        scores = marker_score_df.loc[ct].replace([np.inf, -np.inf], np.nan).dropna()
        marker_dict[ct] = scores.nlargest(top_k).index.tolist()

    if plot_cell_types is not None:
        if isinstance(plot_cell_types, str):
            plot_cell_types = [plot_cell_types]

        plot_cell_types = [str(ct) for ct in plot_cell_types]

        missing = sorted(set(plot_cell_types) - set(categories.astype(str)))
        if len(missing) > 0:
            raise ValueError(
                "These plot_cell_types were not found among eligible cell types: "
                f"{missing}"
            )

        marker_dict_plot = {
            ct: marker_dict[ct]
            for ct in plot_cell_types
            if ct in marker_dict
        }

    else:
        marker_dict_plot = marker_dict.copy()

    all_markers = []
    for genes in marker_dict_plot.values():
        all_markers.extend(genes)

    all_markers = list(dict.fromkeys(all_markers))

    if len(all_markers) == 0:
        print("No markers found for selected plotted cell types.")
        return marker_dict, marker_score_df

    ad_plot = sc.AnnData(
        X=adata.layers[layer].copy(),
        obs=adata.obs[[cell_type_col]].copy(),
        var=adata.var.copy(),
    )

    ad_plot.obs[cell_type_col] = ad_plot.obs[cell_type_col].astype(str)

    if plot_cell_types is not None:
        plot_mask = ad_plot.obs[cell_type_col].isin(plot_cell_types).values
        ad_plot = ad_plot[plot_mask, :].copy()

        if ad_plot.n_obs == 0:
            raise ValueError("No cells found for selected plot_cell_types.")

    ad_plot = ad_plot[:, all_markers].copy()

    if issparse(ad_plot.X):
        ad_plot.X = ad_plot.X.toarray()

    ad_plot.X = np.asarray(ad_plot.X, dtype=np.float32)
    ad_plot.X = np.nan_to_num(ad_plot.X, nan=0.0, posinf=0.0, neginf=0.0)

    if clip is not None:
        ad_plot.X = np.clip(ad_plot.X, clip[0], clip[1])

    # Precompute the dendrogram with a dimension-safe PCA. scanpy's default dendrogram
    # PCA uses n_comps=50, which errors when there are few marker genes or cell types.
    use_dendro = is_dendrogram
    if is_dendrogram:
        try:
            npc = max(2, min(50, ad_plot.n_vars - 1, ad_plot.n_obs - 1))
            sc.pp.pca(ad_plot, n_comps=npc)
            sc.tl.dendrogram(ad_plot, groupby=cell_type_col, use_rep="X_pca")
        except Exception:
            use_dendro = False

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        if figsize is not None:
            plt.figure(figsize=figsize)

        sc.pl.matrixplot(
            ad_plot,
            var_names=marker_dict_plot,
            groupby=cell_type_col,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            dendrogram=use_dendro,
            colorbar_title="Mean gene score",
            title="",
            show=False,
        )

        plt.tight_layout()

        if save_path is not None:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")

        plt.close()

    return marker_dict, marker_score_df


# %% Extracted from sc_stir_analysis_kang notebook cell 16

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.sparse import issparse


def plot_one_cell_type_specific_marker_gene_scores(
    adata,
    target_cell_type,
    cell_type_col="Cell_Type",
    layer="rra_gene_score",
    top_k=20,
    min_cells=10,
    normalize_plot_score=None,
    figsize=None,
    save_path=None,
    return_markers=True,
):
    """
    Identify marker genes specific to one target cell type and plot their
    target-cell-type aggregate gene scores as a horizontal bar plot.

    Marker selection is based on:

        ranking_score = target_mean - max_other_mean

    No filtering is applied.

    The plot always shows target_mean.

    Optionally, target_mean can be min-max normalized across the selected
    top marker genes for plotting.

    Parameters
    ----------
    adata:
        AnnData object.

    target_cell_type:
        Cell type for which markers should be identified.

    cell_type_col:
        Column in adata.obs containing cell type labels.

    layer:
        adata.layers key containing gene scores.

    top_k:
        Number of top marker genes to plot.

    min_cells:
        Minimum number of cells required per cell type.

    normalize_plot_score:
        How to normalize plotted target_mean values.

        Options:
            None     : plot raw target-cell-type mean gene score
            "minmax" : min-max scale target_mean across selected genes to [0, 1]

        Default:
            None

    figsize:
        Optional figure size.

    save_path:
        Optional path to save the plot.

    return_markers:
        If True, return marker genes and full marker statistics.

    Returns
    -------
    If return_markers=True:
        marker_genes, marker_stats_df

    If return_markers=False:
        None
    """

    print(
        f"--- Finding `{target_cell_type}`-specific markers from `{layer}` "
        "and plotting target mean gene scores ---"
    )

    if layer not in adata.layers:
        raise ValueError(f"Layer `{layer}` not found in adata.layers.")

    if cell_type_col not in adata.obs:
        raise ValueError(f"Column `{cell_type_col}` not found in adata.obs.")

    if normalize_plot_score not in {None, "minmax"}:
        raise ValueError("normalize_plot_score must be None or 'minmax'.")

    target_cell_type = str(target_cell_type)

    valid_mask = adata.obs[cell_type_col].notna().values
    cell_types = adata.obs.loc[valid_mask, cell_type_col].astype(str)

    counts = cell_types.value_counts()
    keep_cell_types = counts[counts >= min_cells].index.astype(str).tolist()

    if target_cell_type not in keep_cell_types:
        raise ValueError(
            f"Target cell type `{target_cell_type}` not found or has fewer than "
            f"{min_cells} cells."
        )

    if len(keep_cell_types) < 2:
        raise ValueError("Need at least two eligible cell types for specificity.")

    valid_mask = (
        valid_mask
        & adata.obs[cell_type_col].astype(str).isin(keep_cell_types).values
    )

    X = adata.layers[layer][valid_mask, :]

    if issparse(X):
        X = X.toarray()

    X = np.asarray(X, dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    ct_values = adata.obs.loc[valid_mask, cell_type_col].astype(str).values
    categories = pd.Index(sorted(pd.unique(ct_values)))

    codes = pd.Categorical(ct_values, categories=categories).codes

    n_cells = X.shape[0]
    n_types = len(categories)

    design = np.zeros((n_cells, n_types), dtype=np.float32)
    design[np.arange(n_cells), codes] = 1.0

    group_counts = design.sum(axis=0)[:, None]

    mean_scores = (design.T @ X) / group_counts

    target_idx = categories.get_loc(target_cell_type)
    other_idx = np.arange(n_types) != target_idx

    target_mean = mean_scores[target_idx, :]
    max_other_mean = mean_scores[other_idx, :].max(axis=0)

    strongest_competitor_idx = mean_scores[other_idx, :].argmax(axis=0)
    other_categories = categories[other_idx]

    strongest_competitor = other_categories[
        strongest_competitor_idx
    ].astype(str)

    ranking_scores = target_mean - max_other_mean

    marker_stats_df = pd.DataFrame(
        {
            "gene": adata.var_names.astype(str),
            "ranking_score": ranking_scores,
            "target_mean": target_mean,
            "max_other_mean": max_other_mean,
            "strongest_competitor": strongest_competitor,
        }
    )

    marker_stats_df = marker_stats_df.replace([np.inf, -np.inf], np.nan)

    selected_df = (
        marker_stats_df
        .dropna(subset=["ranking_score"])
        .sort_values("ranking_score", ascending=False)
        .head(top_k)
        .copy()
    )

    marker_genes = selected_df["gene"].tolist()

    if len(marker_genes) == 0:
        print(f"No markers found for `{target_cell_type}`.")

        if return_markers:
            return marker_genes, marker_stats_df

        return None

    plot_df = selected_df.copy()

    if normalize_plot_score is None:
        plot_df["plot_value"] = plot_df["target_mean"]
        xlabel = "Aggregated gene score"

    else:
        values = plot_df["target_mean"].astype(float).values
        value_min = values.min()
        value_max = values.max()

        if value_max == value_min:
            plot_df["plot_value"] = 0.0
        else:
            plot_df["plot_value"] = (values - value_min) / (value_max - value_min)

        xlabel = "Normalized aggregate gene score"

    plot_df = plot_df.sort_values("plot_value", ascending=True)

    if figsize is None:
        figsize = (6, max(3, 0.35 * len(plot_df)))

    plt.figure(figsize=figsize)

    plt.barh(
        plot_df["gene"],
        plot_df["plot_value"],
    )

    plt.xlabel(xlabel)
    plt.ylabel("")
    plt.title(f"{target_cell_type}-specific marker genes")

    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    plt.close()

    print(f"Selected {len(marker_genes)} markers for `{target_cell_type}`:")
    print(marker_genes)

    if return_markers:
        return marker_genes, marker_stats_df

    return None


# %% Extracted from sc_stir_analysis_kang notebook cell 18

import matplotlib.pyplot as plt
import seaborn as sns

def plot_single_celltype_lollipop(adata, cell_type, cell_type_col='Cell_Type', layer='rra_gene_score', top_k=20, save_path=None, figsize=(6, 8)):
    # 1. Extract and average data for the specific cell type
    mask = adata.obs[cell_type_col] == cell_type
    gene_score = np.nanmean(adata.layers[layer][mask, :], axis=0)
    
    df = pd.DataFrame({'gene': adata.var_names, 'gene_score': gene_score})
    
    # 2. Get top positive and negative drivers
    top_df = pd.concat([
        df.nlargest(top_k // 2, 'gene_score'),
        df.nsmallest(top_k // 2, 'gene_score')
    ]).sort_values('gene_score', ascending=True)

    # 3. Plotting
    plt.figure(figsize=figsize)
    colors = ['#d73027' if x > 0 else '#4575b4' for x in top_df['gene_score']]
    
    plt.hlines(y=top_df['gene'], xmin=0, xmax=top_df['gene_score'], color=colors, alpha=0.6)
    plt.scatter(top_df['gene_score'], top_df['gene'], color=colors, s=80, edgecolors='white')
    
    # Add significance lines
    plt.axvline(0, color='black', linewidth=0.8)
    
    plt.title(f"Top Gene score Drivers: {cell_type}")
    plt.xlabel("Mean Score")
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


# %% Extracted from sc_stir_analysis_kang notebook cell 20

import pandas as pd
import numpy as np
import warnings
import os
import time  # Added for retry delay

def run_go_celltype(
    adata, 
    cell_type_col='Cell_Type', 
    layer='rra_avg_expr_diff', 
    layer_p='rra_directional_fdr',
    n_genes=50,                  
    min_sig_cells=5,               
    organism='human',     
    db='GO_Biological_Process_2023',
    output_csv='go_enrichment_results.csv',
    max_retries=3,       # Added max retries
    retry_delay=5        # Added delay between retries (in seconds)
):
    """
    Computes GO enrichment for up-regulated marker genes of each cell type.
    Includes a retry mechanism for API failures.
    """
    import gseapy as gp  # optional dependency; imported lazily so `import sccst.downstream` works without it
    genes = adata.var_names.to_numpy()
    unique_types = adata.obs[cell_type_col].dropna().unique()
    all_enrichment_results = []
    
    print(f"--- Processing GO Enrichment: {cell_type_col} ---")
    
    for ct in unique_types:
        is_target = (adata.obs[cell_type_col] == ct).values
        is_rest = (adata.obs[cell_type_col] != ct) & (adata.obs[cell_type_col].notna())
        
        if np.sum(is_target) < 3: 
            continue
        
        # 1. Extract Data
        score_target = adata.layers[layer][is_target, :]
        score_rest = adata.layers[layer][is_rest, :]
        p_target = adata.layers[layer_p][is_target, :]

        # 2. Filter: Significant in at least X cells
        sig_counts = np.sum(p_target < 0.05, axis=0)
        valid_gene_mask = sig_counts >= min_sig_cells
        
        if np.sum(valid_gene_mask) < 10:
            continue

        # 3. Ranking: Focus on Enrichment (Target > Rest)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_target = np.nan_to_num(np.nanmean(score_target[:, valid_gene_mask], axis=0), nan=0.0)
            mean_rest = np.nan_to_num(np.nanmean(score_rest[:, valid_gene_mask], axis=0), nan=0.0)
        
        diff_score = mean_target - mean_rest
        valid_genes_names = genes[valid_gene_mask]
        
        df_diff = pd.DataFrame({'gene': valid_genes_names, 'score': diff_score})
        
        # Select only top N genes with the highest positive scores
        up_genes = df_diff.sort_values('score', ascending=False).head(n_genes)['gene'].tolist()

        # 4. Enrichment Helper with Retry Mechanism
        if len(up_genes) >= 10:
            for attempt in range(max_retries):
                try:
                    enr = gp.enrichr(gene_list=up_genes, gene_sets=db, organism=organism, cutoff=0.05)
                    res = enr.results
                    
                    if res is not None and not res.empty: 
                        res = res[res['Adjusted P-value'] < 0.05].copy()
                        
                        if not res.empty:
                            res['Cell_Type'] = ct
                            res['LogP'] = -np.log10(res['Adjusted P-value'])
                            res['Term'] = res['Term'].str.split(r' \(GO').str[0]
                            
                            all_enrichment_results.append(res)
                            
                    # Break the retry loop if successful
                    break 
                    
                except Exception as e:
                    if attempt < max_retries - 1:
                        print(f"Error processing {ct}: {e}. Retrying in {retry_delay} seconds... (Attempt {attempt + 1} of {max_retries})")
                        time.sleep(retry_delay)
                    else:
                        print(f"Failed to process {ct} after {max_retries} attempts. Final error: {e}")

    # 5. Finalize and Save
    if not all_enrichment_results:
        print("No significant enrichment found.")
        return None

    final_df = pd.concat(all_enrichment_results, ignore_index=True)
    final_df = final_df.sort_values(['Cell_Type', 'Adjusted P-value'])

    final_df.to_csv(output_csv, index=False)
    print(f"Enrichment results saved to: {output_csv}")

    return final_df


# %% Extracted from sc_stir_analysis_kang notebook cell 23

import textwrap
import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram

def plot_go_enrichment(df, cell_types=None, top_k=5, title_suffix="", data_dir=None, figsize=None, save_path=None):
    """
    Plots GO enrichment with Gene Hits.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain columns: 'Cell_Type', 'Term', 'Overlap', and 'LogP'.
    cell_types : list[str] or None
        Cell types to plot, e.g. ['Neuron', 'Astrocyte'].
        If None, plots all cell types.
    top_k : int
        Number of top GO terms to show per cell type.
    title_suffix : str
        Optional suffix for plot title and output filename.
    data_dir : pathlib.Path or str or None
        Directory where the PDF should be saved.
        If None, the plot is shown but not saved.
    """

    if df is None or df.empty:
        print("No data to plot.")
        return

    required_cols = {'Cell_Type', 'Term', 'Overlap', 'LogP'}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    plot_df = df.copy()

    # Filter by selected cell types
    if cell_types is not None:
        available_types = plot_df['Cell_Type'].unique()
        valid_types = [ct for ct in cell_types if ct in available_types]

        if not valid_types:
            print(f"None of the provided cell types {cell_types} found in data.")
            return

        plot_df = plot_df[plot_df['Cell_Type'].isin(valid_types)]

    # Parse 'Overlap' to get Gene Hits, e.g. "11/189" -> 11
    plot_df['Gene_Hits'] = (
        plot_df['Overlap']
        .astype(str)
        .str.split('/')
        .str[0]
        .astype(int)
    )

    # Select Top K per Cell Type
    plot_df = (
        plot_df.sort_values(['Cell_Type', 'LogP'], ascending=[True, False])
        .groupby('Cell_Type', group_keys=False)
        .head(top_k)
    )

    if plot_df.empty:
        print("No data left after filtering.")
        return

    # Cluster cell types
    pivot_df = (
        plot_df
        .pivot_table(index='Term', columns='Cell_Type', values='LogP')
        .fillna(0)
    )

    if pivot_df.shape[1] > 1:
        Z = linkage(pivot_df.T, method='ward')
        cluster_order = [
            pivot_df.columns[i]
            for i in dendrogram(Z, no_plot=True)['leaves']
        ]
        plot_df['Cell_Type'] = pd.Categorical(
            plot_df['Cell_Type'],
            categories=cluster_order,
            ordered=True
        )

    # Wrap long GO Term names
    wrap_width = 35
    plot_df['Term'] = plot_df['Term'].apply(
        lambda x: '\n'.join(textwrap.wrap(str(x), width=wrap_width))
    )

    # Final sort for visualization
    plot_df = plot_df.sort_values(['Cell_Type', 'LogP'], ascending=[True, False])

    # Dynamic height that accounts for wrapped, multi-line y-axis labels
    unique_terms = plot_df['Term'].drop_duplicates()
    total_term_lines = sum(term.count('\n') + 1 for term in unique_terms)

    fig_height = max(10, total_term_lines * 0.35)

    plt.figure(figsize=figsize if figsize is not None else (12, fig_height))

    scatter = plt.scatter(
        x=plot_df['Cell_Type'],
        y=plot_df['Term'],
        s=plot_df['Gene_Hits'] * 20,
        c=plot_df['LogP'],
        cmap='Reds',
        alpha=0.9,
        edgecolors='black',
        linewidth=0.5
    )

    # Reduce y-axis font size to avoid overlap
    plt.yticks(fontsize=7)

    # Colorbar
    cbar = plt.colorbar(scatter, fraction=0.03, pad=0.1)
    cbar.set_label('-log10(Adjusted P-value)', rotation=270, labelpad=15)

    # Legend for Gene Hits
    handles, labels = [], []
    sizes = [5, 10, 20, 50]

    for size in sizes:
        handles.append(
            plt.scatter(
                [],
                [],
                s=size * 20,
                color='gray',
                alpha=0.5,
                edgecolors='black'
            )
        )
        labels.append(str(size))

    plt.legend(
        handles,
        labels,
        title="Gene Hits",
        bbox_to_anchor=(1.25, 0.5),
        loc='center left',
        labelspacing=2.5,
        borderpad=1.5,
        handletextpad=1.5,
        frameon=True
    )

    plt.title(
        f"GO Enrichment: {title_suffix}" if title_suffix else "GO Enrichment",
        fontsize=14,
        pad=20
    )

    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Gene Ontology Terms")
    plt.xlabel("Cell Clusters")
    plt.grid(True, axis='both', linestyle=':', alpha=0.3)

    plt.tight_layout()

    # Save to save_path (a file, like the other plotting functions); fall back to data_dir.
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved plot to: {save_path}")
    elif data_dir is not None:
        from pathlib import Path

        data_dir = Path(data_dir)
        data_dir.mkdir(parents=True, exist_ok=True)

        safe_suffix = title_suffix.replace(" ", "_") if title_suffix else ""
        output_file = data_dir / f'go_enrichment_plot{("_" + safe_suffix) if safe_suffix else ""}.pdf'

        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved plot to: {output_file}")

    plt.close()


# %% Extracted from sc_stir_analysis_kang notebook cell 25

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import issparse

def visualize_gene_stir(
    adata, 
    gene_name: str, 
    condition_col: str = 'condition', 
    layer: str = None,  # None for Expression (.X), 'stir_z_scores_signed' for Z-scores
    save_path: str = None,
    figsize=None,
):
    """
    Intelligent plotter:
    - If layer is None (Expression): Plots 2 side-by-side UMAPs (Condition A vs Condition B).
    - If layer is 'gene score': Plots 1 single UMAP (Disease Relevance Map).
    """
    
    # 1. Validation
    if gene_name not in adata.var_names:
        print(f"Error: {gene_name} not found.")
        return

    # 2. Extract Data for Scaling (Global vmin/vmax)
    if layer:
        data_vec = adata[:, gene_name].layers[layer]
    else:
        data_vec = adata[:, gene_name].X
        
    if issparse(data_vec): data_vec = data_vec.toarray()
    
    # ---------------------------------------------------------
    # MODE A: pi-scores (Single Plot, Divergent Color)
    # ---------------------------------------------------------
    if layer != None:
        print(f"Plotting Gene Score Map for '{gene_name}'...")
        
        # Centered scaling for Z-scores (e.g., -3 to +3)
        max_val = np.nanmax(np.abs(data_vec))
        # Cap visual noise at 5.0 usually
        limit = min(max_val, 5.0) 
        
        plt.figure(figsize=figsize if figsize is not None else (7, 6))
        
        sc.pl.umap(
            adata,
            color=gene_name,
            layer=layer,
            cmap='RdBu_r',       # Red = High Disease, Blue = High Control
            vmin=-limit,
            vmax=limit,
            title=f"{gene_name}\n(Gene Importance Score)",
            frameon=False,
            show=False
        )

    # ---------------------------------------------------------
    # MODE B: GENE EXPRESSION (Two Plots, Sequential Color)
    # ---------------------------------------------------------
    else:
        conditions = adata.obs[condition_col].unique()
        if len(conditions) != 2:
            print("Error: Need exactly 2 conditions for split plotting.")
            return
            
        print(f"Plotting Expression Split for '{gene_name}' ({conditions[0]} vs {conditions[1]})...")
        
        # 0 to Max scaling for expression
        vmax = np.nanmax(data_vec)
        vmin = np.nanmin(data_vec)
        
        fig, axes = plt.subplots(1, 2, figsize=figsize if figsize is not None else (12, 5.5))
        
        for i, cond in enumerate(conditions):
            # Subset
            subset = adata[adata.obs[condition_col] == cond]
            
            sc.pl.umap(
                subset, 
                color=gene_name,
                layer=layer,
                ax=axes[i],
                vmin=vmin, 
                vmax=vmax,
                cmap='viridis',  # Standard expression colormap
                title=f"{gene_name}\n({cond})",
                frameon=False,
                show=False,
                use_raw=False,
                colorbar_loc=None if i==0 else 'right'
            )
        
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
