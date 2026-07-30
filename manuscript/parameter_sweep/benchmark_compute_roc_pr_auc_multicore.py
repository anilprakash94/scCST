#!/usr/bin/env python
# ============================================================
# Standalone benchmark: ROC-AUC + PR-AUC (AUPRC) for scCST / CACOA / miloDE
#
# Computes, for every simulation dataset x neighborhood size x method:
#   - cell_macro   (mean over responder disease cells)
#   - gene_macro   (mean over true marker genes)
#   - global       (all disease-cell/gene pairs pooled)
# for BOTH roc_auc and pr_auc in a single pass (each score file loaded once).
#
# Also exports the responder vs non-responder DE-gene-score table used by the
# boxplot notebook: per (responder disease cell x its group's marker gene) score
# and per (non-responder disease cell x its group's marker gene) score.
#
# This is the heavy-compute step. Run it once via sbatch; the plotting notebooks
# then read the saved tables only. Logic (file matching, orientation, truth
# construction, alignment) mirrors the AUPRC/ROC benchmark notebooks exactly.
# ============================================================

from pathlib import Path
import os
import re
import sys
import time
import warnings

# Avoid oversubscription when running many benchmark workers.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import roc_auc_score, average_precision_score

from joblib import Parallel, delayed

warnings.filterwarnings("ignore")


# ============================================================
# 1. Configuration
# ============================================================
simulation_roots = [
    Path("/path/to/data/scrna_seq/simulation")
]
simulation_root_configs = []
default_base_dir = Path("/path/to/data/scrna_seq/simulation")

# Sweep/dataset selection (empty list -> use all values found in manifest).
requested_sweeps = {
    "responder_percent": [20, 40, 60],
    "dropout_rate": [0.2, 0.4, 0.6],
    "harmony_theta": [1, 2, 3],
}
requested_replicates = [1, 2, 3]
benchmark_nhood_sizes = [10, 20, 30]

method_order = ["scCST", "CACOA", "miloDE"]

# Output directory for benchmark tables and figures.
benchmark_output_dir = "/path/to/data/scrna_seq/benchmark_roc_pr_auc_combined"

# Parallelism. Workers taken from env (SLURM), else this default.
n_benchmark_workers = int(os.environ.get("WORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "16")))
parallel_backend = "loky"
parallel_verbose = 10
worker_batch_size = 1

# DE-gene-score export: cap sampled pairs per (dataset x nhood x method x status)
# group so the boxplot table stays a manageable size while preserving shape.
MAX_DE_PAIRS_PER_GROUP = 3000

metric_families = {
    "roc_auc": roc_auc_score,
    "pr_auc": average_precision_score,
}


# ============================================================
# 2. Resolve roots and load manifests
# ============================================================
def _as_path(x):
    return Path(str(x)).expanduser().resolve()


def _default_milode_dir(sim_dir):
    new_dir = sim_dir / "milode_parameter_sweeps_results_from_hvg_rds"
    old_dir = sim_dir / "milode_parameter_sweeps_results"
    if new_dir.exists() or not old_dir.exists():
        return new_dir
    return old_dir


def normalize_root_configs(simulation_roots, simulation_root_configs, default_base_dir):
    configs = []
    for root in simulation_roots:
        root = _as_path(root)
        configs.append({
            "name": root.name,
            "sim_dir": root,
            "sccst_dir": root / "sccst_results",
            "cacoa_dir": root / "cacoa_multi_results_parameter_sweeps",
            "milode_dir": _default_milode_dir(root),
        })
    for cfg in simulation_root_configs:
        sim_dir = _as_path(cfg["sim_dir"])
        configs.append({
            "name": str(cfg.get("name", sim_dir.name)),
            "sim_dir": sim_dir,
            "sccst_dir": _as_path(cfg.get("sccst_dir", sim_dir / "sccst_results")),
            "cacoa_dir": _as_path(cfg.get("cacoa_dir", sim_dir / "cacoa_multi_results_parameter_sweeps")),
            "milode_dir": _as_path(cfg.get("milode_dir", _default_milode_dir(sim_dir))),
        })
    if len(configs) == 0:
        root = _as_path(default_base_dir)
        configs.append({
            "name": root.name,
            "sim_dir": root,
            "sccst_dir": root / "sccst_results",
            "cacoa_dir": root / "cacoa_multi_results_parameter_sweeps",
            "milode_dir": _default_milode_dir(root),
        })
    seen = {}
    for cfg in configs:
        name = cfg["name"]
        seen[name] = seen.get(name, 0) + 1
        if seen[name] > 1:
            cfg["name"] = f"{name}_{seen[name]}"
    return configs


root_configs = normalize_root_configs(simulation_roots, simulation_root_configs, default_base_dir)

if benchmark_output_dir is None:
    benchmark_dir = root_configs[0]["sim_dir"] / "benchmark_roc_pr_auc_combined"
else:
    benchmark_dir = _as_path(benchmark_output_dir)

fig_dir = benchmark_dir / "figures"
table_dir = benchmark_dir / "tables"
for d in [benchmark_dir, fig_dir, table_dir]:
    d.mkdir(parents=True, exist_ok=True)

print("Benchmark output:", benchmark_dir, flush=True)
for cfg in root_configs:
    print("-", cfg["name"], flush=True)
    print("  sim_dir  :", cfg["sim_dir"], "exists=", cfg["sim_dir"].exists(), flush=True)
    print("  scCST    :", cfg["sccst_dir"], "exists=", cfg["sccst_dir"].exists(), flush=True)
    print("  CACOA    :", cfg["cacoa_dir"], "exists=", cfg["cacoa_dir"].exists(), flush=True)
    print("  miloDE   :", cfg["milode_dir"], "exists=", cfg["milode_dir"].exists(), flush=True)


required_cols = [
    "simulation_type", "sweep_name", "swept_parameter", "swept_value",
    "replicate", "responder_percent", "dropout_rate", "harmony_theta",
]


def _find_h5ad_column(df):
    for col in ["file", "input_h5ad", "h5ad_file", "adata_file"]:
        if col in df.columns:
            return col
    return None


def resolve_sim_file(path_string, sim_dir):
    p = Path(str(path_string))
    if p.exists():
        return p
    alt = Path(sim_dir) / p.name
    if alt.exists():
        return alt
    return p


manifest_parts = []
manifest_errors = []
for cfg in root_configs:
    manifest_file = cfg["sim_dir"] / "simulation_manifest_de_parameter_sweeps.csv"
    if not manifest_file.exists():
        manifest_errors.append({"source_name": cfg["name"], "manifest_file": str(manifest_file), "reason": "manifest not found"})
        continue
    df = pd.read_csv(manifest_file)
    h5ad_col = _find_h5ad_column(df)
    missing_cols = [c for c in required_cols if c not in df.columns]
    if h5ad_col is None:
        missing_cols.append("file/input_h5ad")
    if missing_cols:
        manifest_errors.append({"source_name": cfg["name"], "manifest_file": str(manifest_file), "reason": f"missing columns: {missing_cols}"})
        continue
    df = df.copy()
    df["file"] = df[h5ad_col].astype(str)
    df["source_name"] = cfg["name"]
    df["sim_dir"] = str(cfg["sim_dir"])
    df["sccst_dir"] = str(cfg["sccst_dir"])
    df["cacoa_dir"] = str(cfg["cacoa_dir"])
    df["milode_dir"] = str(cfg["milode_dir"])
    df["manifest_file"] = str(manifest_file)
    df["resolved_file"] = df["file"].map(lambda x: str(resolve_sim_file(x, cfg["sim_dir"])))
    df = df[df["resolved_file"].map(lambda x: Path(x).exists())].copy()
    manifest_parts.append(df)

if manifest_errors:
    pd.DataFrame(manifest_errors).to_csv(table_dir / "manifest_loading_errors.csv", index=False)

if not manifest_parts:
    raise FileNotFoundError("No valid simulation manifests were found.")

manifest_df = pd.concat(manifest_parts, ignore_index=True)
manifest_df["swept_value_num"] = pd.to_numeric(manifest_df["swept_value"], errors="coerce")
manifest_df["replicate"] = pd.to_numeric(manifest_df["replicate"], errors="coerce").astype("Int64")
manifest_df = manifest_df[manifest_df["swept_parameter"].astype(str).isin(requested_sweeps.keys())].copy()


def _keep_requested_sweep_value(row):
    param = str(row["swept_parameter"])
    wanted = requested_sweeps.get(param, [])
    if wanted is None or len(wanted) == 0:
        return True
    return float(row["swept_value_num"]) in [float(x) for x in wanted]


manifest_df = manifest_df[manifest_df.apply(_keep_requested_sweep_value, axis=1)].copy()
manifest_df = manifest_df[manifest_df["replicate"].isin(requested_replicates)].copy()
manifest_df = manifest_df.sort_values(["source_name", "swept_parameter", "swept_value_num", "replicate"]).reset_index(drop=True)

print(f"Using {len(manifest_df)} simulation files from {len(root_configs)} root(s).", flush=True)
manifest_df.to_csv(table_dir / "benchmark_input_manifest.csv", index=False)


# ============================================================
# 3. Filename/tag utilities and method-output finders
# ============================================================
def safe_value_tag(x):
    text = str(x)
    if text.endswith(".0"):
        text = text[:-2]
    text = text.replace(".", "p")
    text = re.sub(r"[^A-Za-z0-9_-]+", "_", text)
    return text


def sccst_value_tag(x):
    text = str(x).replace(".", "p")
    text = re.sub(r"[^A-Za-z0-9_-]+", "_", text)
    return text


def make_dataset_tag(row):
    return (
        f"{row['simulation_type']}"
        f"_{row['sweep_name']}"
        f"_{row['swept_parameter']}_{safe_value_tag(row['swept_value'])}"
        f"_replicate_{int(row['replicate'])}"
        f"_responder_percent_{int(row['responder_percent'])}"
        f"_dropout_{safe_value_tag(row['dropout_rate'])}"
        f"_harmony_theta_{safe_value_tag(row['harmony_theta'])}"
    )


def make_sccst_tag(row, nhood_size, tag_style="sccst"):
    tagger = sccst_value_tag if tag_style == "sccst" else safe_value_tag
    return (
        f"{row['sweep_name']}"
        f"_{row['swept_parameter']}_{tagger(row['swept_value'])}"
        f"_replicate_{int(row['replicate'])}"
        f"_nhood_{int(nhood_size)}"
        f"_responder_percent_{int(row['responder_percent'])}"
        f"_dropout_{tagger(row['dropout_rate'])}"
        f"_harmony_theta_{tagger(row['harmony_theta'])}"
    )


def dataset_label_from_row(row):
    return (
        f"{row['source_name']}|"
        f"{row['sweep_name']}|"
        f"{row['swept_parameter']}={float(row['swept_value_num']):g}|"
        f"rep={int(row['replicate'])}"
    )


def find_sccst_file(row, nhood_size):
    sccst_dir = Path(row["sccst_dir"])
    candidates = [
        sccst_dir / f"sccst_result_{make_sccst_tag(row, nhood_size, tag_style='sccst')}.h5ad",
        sccst_dir / f"sccst_result_{make_sccst_tag(row, nhood_size, tag_style='safe')}.h5ad",
    ]
    for exact in candidates:
        if exact.exists():
            return exact
    sweep_tags = sorted(set([sccst_value_tag(row["swept_value"]), safe_value_tag(row["swept_value"])]))
    dropout_tags = sorted(set([sccst_value_tag(row["dropout_rate"]), safe_value_tag(row["dropout_rate"])]))
    theta_tags = sorted(set([sccst_value_tag(row["harmony_theta"]), safe_value_tag(row["harmony_theta"])]))
    matches = []
    for sv in sweep_tags:
        for dr in dropout_tags:
            for th in theta_tags:
                pattern = (
                    f"sccst_result_{row['sweep_name']}_{row['swept_parameter']}_{sv}"
                    f"_replicate_{int(row['replicate'])}_nhood_{int(nhood_size)}"
                    f"_responder_percent_{int(row['responder_percent'])}"
                    f"_dropout_{dr}_harmony_theta_{th}.h5ad"
                )
                matches.extend(sccst_dir.glob(pattern))
    matches = sorted(set(matches))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise RuntimeError("Ambiguous scCST files:\n" + "\n".join(str(x) for x in matches))
    return None


def find_cacoa_score_file(row, nhood_size):
    cacoa_dir = Path(row["cacoa_dir"])
    tag = f"nhood_{int(nhood_size)}_{make_dataset_tag(row)}"
    exact = cacoa_dir / f"cacoa_cluster_free_de_z_{tag}.csv"
    if exact.exists():
        return exact
    sweep_tags = sorted(set([safe_value_tag(row["swept_value"]), sccst_value_tag(row["swept_value"])]))
    dropout_tags = sorted(set([safe_value_tag(row["dropout_rate"]), sccst_value_tag(row["dropout_rate"])]))
    theta_tags = sorted(set([safe_value_tag(row["harmony_theta"]), sccst_value_tag(row["harmony_theta"])]))
    matches = []
    for sv in sweep_tags:
        for dr in dropout_tags:
            for th in theta_tags:
                pattern = (
                    f"cacoa_cluster_free_de_z_nhood_{int(nhood_size)}_{row['simulation_type']}"
                    f"_{row['sweep_name']}_{row['swept_parameter']}_{sv}"
                    f"_replicate_{int(row['replicate'])}"
                    f"_responder_percent_{int(row['responder_percent'])}"
                    f"_dropout_{dr}_harmony_theta_{th}.csv"
                )
                matches.extend(cacoa_dir.glob(pattern))
    matches = sorted(set(matches))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise RuntimeError("Ambiguous CACOA score files:\n" + "\n".join(str(x) for x in matches))
    return None


def find_milode_score_file(row, nhood_size):
    milode_dir = Path(row["milode_dir"])
    tag = f"nhood_{int(nhood_size)}_{make_dataset_tag(row)}"
    exact = milode_dir / f"milode_cell_z_scores_{tag}.csv"
    if exact.exists():
        return exact
    search_dirs = [milode_dir]
    old_dir = Path(row["sim_dir"]) / "milode_parameter_sweeps_results"
    if old_dir not in search_dirs:
        search_dirs.append(old_dir)
    sweep_tags = sorted(set([safe_value_tag(row["swept_value"]), sccst_value_tag(row["swept_value"])]))
    dropout_tags = sorted(set([safe_value_tag(row["dropout_rate"]), sccst_value_tag(row["dropout_rate"])]))
    theta_tags = sorted(set([safe_value_tag(row["harmony_theta"]), sccst_value_tag(row["harmony_theta"])]))
    matches = []
    for d in search_dirs:
        if not d.exists():
            continue
        for sv in sweep_tags:
            for dr in dropout_tags:
                for th in theta_tags:
                    pattern = (
                        f"milode_cell_z_scores_nhood_{int(nhood_size)}_{row['simulation_type']}"
                        f"_{row['sweep_name']}_{row['swept_parameter']}_{sv}"
                        f"_replicate_{int(row['replicate'])}"
                        f"_responder_percent_{int(row['responder_percent'])}"
                        f"_dropout_{dr}_harmony_theta_{th}.csv"
                    )
                    matches.extend(d.glob(pattern))
    matches = sorted(set(matches))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise RuntimeError("Ambiguous miloDE score files:\n" + "\n".join(str(x) for x in matches))
    return None


# ============================================================
# 4. Score and truth loading helpers
# ============================================================
def read_csv_matrix(path):
    path = Path(path)
    df = pd.read_csv(path, index_col=0)
    empty_cols = [c for c in df.columns if str(c).startswith("Unnamed")]
    if empty_cols:
        df = df.drop(columns=empty_cols)
    df_num = df.apply(pd.to_numeric, errors="coerce")
    non_numeric = df_num.columns[df_num.isna().all()].tolist()
    if non_numeric:
        df_num = df_num.drop(columns=non_numeric)
    return df_num


def orient_score_matrix(score_df, obs_names, var_names):
    obs_set = set(map(str, obs_names))
    var_set = set(map(str, var_names))
    row_names = list(map(str, score_df.index))
    col_names = list(map(str, score_df.columns))
    n_rows_obs = sum(x in obs_set for x in row_names)
    n_cols_obs = sum(x in obs_set for x in col_names)
    n_rows_var = sum(x in var_set for x in row_names)
    n_cols_var = sum(x in var_set for x in col_names)
    if n_rows_obs >= max(5, n_cols_obs) and n_cols_var >= max(5, n_rows_var):
        out = score_df.copy()
    elif n_rows_var >= max(5, n_cols_var) and n_cols_obs >= max(5, n_rows_obs):
        out = score_df.T.copy()
    elif n_cols_var >= n_rows_var:
        out = score_df.copy()
    else:
        out = score_df.T.copy()
    out.index = out.index.astype(str)
    out.columns = out.columns.astype(str)
    return out


def clean_gene_names(names):
    return [str(x).replace(".", "-") for x in names]


def load_truth_from_sim(sim_file):
    # Backed read: obs/var/uns load into memory but the (unused) expression matrix
    # .X stays on disk. This is a large memory saving for the 20k x 20k sims.
    adata = sc.read_h5ad(sim_file, backed="r")
    adata.var_names = clean_gene_names(adata.var_names)
    marker_dict = adata.uns.get("hidden_group_disease_markers", None)
    if marker_dict is None:
        marker_dict = adata.uns.get("cell_type_markers", None)
    if marker_dict is None:
        raise ValueError("No hidden_group_disease_markers or cell_type_markers in adata.uns")
    marker_dict = {str(k): set(clean_gene_names(v)) for k, v in marker_dict.items()}
    obs = adata.obs.copy()
    obs.index = obs.index.astype(str)
    obs["cell_barcode"] = obs.index
    return adata, obs, marker_dict


def make_group_marker_matrix(obs, genes, marker_dict, responders_only):
    """Mark (cell, gene)=1 when gene is a marker of the cell's hidden group.

    responders_only=True  -> only responder cells marked (this is the benchmark truth).
    responders_only=False -> every cell marked (used to also locate non-responder
                             marker-gene scores for the DE-score boxplots).
    """
    genes = list(map(str, genes))
    mat = np.zeros((obs.shape[0], len(genes)), dtype=np.int8)
    if "sim_group_truth" not in obs.columns:
        raise ValueError("obs must contain sim_group_truth")
    if "is_responder" not in obs.columns:
        raise ValueError("obs must contain is_responder")
    is_resp = obs["is_responder"].astype(str).str.lower().isin(["true", "1", "yes"]).to_numpy()
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    for i, (_, row) in enumerate(obs.iterrows()):
        if responders_only and not is_resp[i]:
            continue
        markers = marker_dict.get(str(row["sim_group_truth"]), set())
        for g in markers:
            j = gene_to_idx.get(g)
            if j is not None:
                mat[i, j] = 1
    return mat, is_resp


def align_scores_to_truth(score_df, truth_obs, marker_dict):
    score_df = score_df.copy()
    score_df.index = score_df.index.astype(str)
    score_df.columns = [str(x).replace(".", "-") for x in score_df.columns]
    truth_obs = truth_obs.copy()
    truth_obs.index = truth_obs.index.astype(str)
    common_cells = [c for c in score_df.index if c in set(truth_obs.index)]
    common_genes = [
        g for g in score_df.columns
        if any(g in markers for markers in marker_dict.values()) or g.startswith("Gene")
    ]
    if len(common_genes) < 10:
        common_genes = list(score_df.columns)
    if len(common_cells) == 0:
        raise ValueError("No overlapping cells between score matrix and truth metadata.")
    obs_aligned = truth_obs.loc[common_cells].copy()
    score_aligned = score_df.loc[common_cells, common_genes].copy()
    if "condition" in obs_aligned.columns:
        disease_mask = obs_aligned["condition"].astype(str).eq("Disease")
        obs_aligned = obs_aligned.loc[disease_mask].copy()
        score_aligned = score_aligned.loc[obs_aligned.index].copy()
    genes = list(score_aligned.columns)
    truth, is_resp = make_group_marker_matrix(obs_aligned, genes, marker_dict, responders_only=True)
    scores = score_aligned.to_numpy(dtype=np.float32)
    scores = np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)
    return obs_aligned, genes, scores, truth, is_resp


def load_sccst_scores(sccst_file, sim_adata=None):
    # Backed read avoids loading the scCST .X matrix; only the rra_gene_score
    # layer (the scores we need) is materialized, as float32 to halve memory.
    a = sc.read_h5ad(sccst_file, backed="r")
    a.var_names = clean_gene_names(a.var_names)
    if "rra_gene_score" not in a.layers:
        raise ValueError(f"rra_gene_score not found in {sccst_file}")
    X = a.layers["rra_gene_score"]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X).astype(np.float32, copy=False)
    score_index = a.obs_names.astype(str)
    if sim_adata is not None and "rra_original_cell_idx" in a.obs.columns:
        original_idx = pd.to_numeric(a.obs["rra_original_cell_idx"], errors="coerce")
        valid = original_idx.notna()
        if int(valid.sum()) == a.n_obs:
            original_idx = original_idx.astype(int).to_numpy()
            if np.all(original_idx >= 0) and np.all(original_idx < sim_adata.n_obs):
                score_index = sim_adata.obs_names[original_idx].astype(str)
    return pd.DataFrame(X, index=score_index, columns=a.var_names.astype(str))


def load_csv_scores(score_file, sim_adata):
    raw_df = read_csv_matrix(score_file)
    raw_df.index = raw_df.index.astype(str)
    raw_df.columns = raw_df.columns.astype(str)
    return orient_score_matrix(
        raw_df,
        obs_names=sim_adata.obs_names.astype(str),
        var_names=sim_adata.var_names.astype(str),
    )


# ============================================================
# 5. Metric functions (both roc_auc and pr_auc)
# ============================================================
def safe_metric(y_true, y_score, kind):
    """kind in {'roc_auc','pr_auc'} with guards for degenerate label sets."""
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    mask = np.isfinite(y_score)
    y_true = y_true[mask]
    y_score = y_score[mask]
    if y_true.size == 0:
        return np.nan
    n_pos = int(np.sum(y_true == 1))
    n_neg = int(np.sum(y_true == 0))
    if n_pos == 0:
        return np.nan
    if n_neg == 0:
        # ROC-AUC undefined with no negatives; AUPRC is 1.0 by convention.
        return np.nan if kind == "roc_auc" else 1.0
    return float(metric_families[kind](y_true, y_score))


def compute_distribution_metrics(scores, truth):
    """Return summary dict + per-cell/per-gene rows carrying both metrics."""
    scores = np.asarray(scores, dtype=float)
    truth = np.asarray(truth, dtype=int)

    cell_rows = []
    for i in range(scores.shape[0]):
        if truth[i, :].sum() > 0:
            cell_rows.append({
                "unit": f"cell_{i}",
                "metric": "cell_macro",
                "roc_auc": safe_metric(truth[i, :], scores[i, :], "roc_auc"),
                "pr_auc": safe_metric(truth[i, :], scores[i, :], "pr_auc"),
            })

    gene_rows = []
    for j in range(scores.shape[1]):
        if truth[:, j].sum() > 0:
            gene_rows.append({
                "unit": f"gene_{j}",
                "metric": "gene_macro",
                "roc_auc": safe_metric(truth[:, j], scores[:, j], "roc_auc"),
                "pr_auc": safe_metric(truth[:, j], scores[:, j], "pr_auc"),
            })

    flat_true = truth.ravel()
    flat_score = scores.ravel()
    summary = {}
    for kind in metric_families:
        cvals = [r[kind] for r in cell_rows if np.isfinite(r[kind])]
        gvals = [r[kind] for r in gene_rows if np.isfinite(r[kind])]
        summary[f"cell_macro_{kind}"] = float(np.nanmean(cvals)) if len(cvals) else np.nan
        summary[f"gene_macro_{kind}"] = float(np.nanmean(gvals)) if len(gvals) else np.nan
        summary[f"global_{kind}"] = safe_metric(flat_true, flat_score, kind)
    summary["n_cell_units"] = int(sum(np.isfinite(r["roc_auc"]) or np.isfinite(r["pr_auc"]) for r in cell_rows))
    summary["n_gene_units"] = int(sum(np.isfinite(r["roc_auc"]) or np.isfinite(r["pr_auc"]) for r in gene_rows))
    return summary, cell_rows, gene_rows


def extract_de_gene_scores(obs_aligned, genes, scores, marker_dict, rng):
    """Responder vs non-responder DE-gene scores: score at (cell, its group's
    marker gene) split by the cell's responder status. Subsampled per group."""
    group_mat, is_resp = make_group_marker_matrix(obs_aligned, genes, marker_dict, responders_only=False)
    out = {}
    for status, sel in [("responder", is_resp), ("non_responder", ~is_resp)]:
        rows = np.where(sel)[0]
        if rows.size == 0:
            out[status] = np.array([], dtype=float)
            continue
        sub_mask = group_mat[rows, :] == 1
        vals = scores[rows, :][sub_mask]
        vals = vals[np.isfinite(vals)]
        if vals.size > MAX_DE_PAIRS_PER_GROUP:
            idx = rng.choice(vals.size, size=MAX_DE_PAIRS_PER_GROUP, replace=False)
            vals = vals[idx]
        out[status] = vals
    return out


# ============================================================
# 6. Per-task benchmark
# ============================================================
def _record_for_worker(row, nhood_size):
    d = row.to_dict()
    out = {}
    for k, v in d.items():
        out[k] = str(v) if isinstance(v, Path) else v
    out["nhood_size"] = int(nhood_size)
    return out


PARAM_KEEP = [
    "source_name", "simulation_type", "sweep_name", "swept_parameter", "swept_value",
    "swept_value_num", "replicate", "responder_percent", "dropout_rate", "harmony_theta",
]


def benchmark_one_dataset_nhood(row_record):
    row = pd.Series(row_record)
    params = {k: row_record.get(k) for k in PARAM_KEEP}
    nhood_size = int(row_record["nhood_size"])
    sim_file = Path(row_record.get("resolved_file", row_record["file"]))
    dataset_label = dataset_label_from_row(row)
    # Deterministic per-task RNG seed for reproducible DE-score subsampling.
    seed = (abs(hash((dataset_label, nhood_size))) % (2**32))
    rng = np.random.default_rng(seed)

    summary_rows, distribution_rows, de_rows, missing_rows, log_lines = [], [], [], [], []

    def log(msg):
        log_lines.append(str(msg))

    log(f"[{params['source_name']}] {sim_file.name} | nhood={nhood_size} | {dataset_label}")

    try:
        sim_adata, truth_obs, marker_dict = load_truth_from_sim(sim_file)
        sim_adata.var_names = clean_gene_names(sim_adata.var_names)
    except Exception as e:
        missing_rows.append({**params, "nhood_size": nhood_size, "method": "truth", "reason": repr(e)})
        return {"summary_rows": summary_rows, "distribution_rows": distribution_rows,
                "de_rows": de_rows, "missing_rows": missing_rows, "log": "\n".join(log_lines)}

    finders = {"scCST": find_sccst_file, "CACOA": find_cacoa_score_file, "miloDE": find_milode_score_file}
    method_files = {}
    for method, finder in finders.items():
        try:
            method_files[method] = finder(row, nhood_size)
        except Exception as e:
            method_files[method] = None
            missing_rows.append({**params, "nhood_size": nhood_size, "method": method, "reason": repr(e)})

    for method, score_file in method_files.items():
        if score_file is None:
            missing_rows.append({**params, "nhood_size": nhood_size, "method": method, "reason": "missing score output"})
            continue
        try:
            if method == "scCST":
                score_df = load_sccst_scores(score_file, sim_adata=sim_adata)
            else:
                score_df = load_csv_scores(score_file, sim_adata=sim_adata)

            obs_aligned, genes_aligned, scores, truth, is_resp = align_scores_to_truth(
                score_df=score_df, truth_obs=truth_obs, marker_dict=marker_dict,
            )
            if truth.sum() == 0:
                raise ValueError("Aligned truth has no positive labels (marker genes absent from output).")

            summary, cell_rows, gene_rows = compute_distribution_metrics(scores, truth)
            summary_rows.append({
                **params, "dataset_label": dataset_label, "method": method, "nhood_size": nhood_size,
                "score_file": str(score_file),
                "n_cells_used": int(scores.shape[0]), "n_genes_used": int(scores.shape[1]),
                "n_positive_pairs": int(truth.sum()), **summary,
            })
            for r in cell_rows + gene_rows:
                distribution_rows.append({
                    **params, "dataset_label": dataset_label, "method": method, "nhood_size": nhood_size,
                    "metric": r["metric"], "unit": r["unit"],
                    "roc_auc": r["roc_auc"], "pr_auc": r["pr_auc"],
                })

            de_scores = extract_de_gene_scores(obs_aligned, genes_aligned, scores, marker_dict, rng)
            for status, vals in de_scores.items():
                for v in vals:
                    de_rows.append({
                        **params, "dataset_label": dataset_label, "method": method,
                        "nhood_size": nhood_size, "responder_status": status, "score": float(v),
                    })

            log(f"  {method}: cell_roc={summary['cell_macro_roc_auc']:.3f} gene_roc={summary['gene_macro_roc_auc']:.3f} "
                f"cell_pr={summary['cell_macro_pr_auc']:.3f} gene_pr={summary['gene_macro_pr_auc']:.3f} "
                f"cells={scores.shape[0]} genes={scores.shape[1]}")
        except Exception as e:
            missing_rows.append({**params, "nhood_size": nhood_size, "method": method,
                                 "score_file": str(score_file), "reason": repr(e)})
            log(f"  {method}: failed -> {repr(e)}")

    del sim_adata, truth_obs, marker_dict
    return {"summary_rows": summary_rows, "distribution_rows": distribution_rows,
            "de_rows": de_rows, "missing_rows": missing_rows, "log": "\n".join(log_lines)}


# ============================================================
# 7. Run
# ============================================================
def main():
    t0 = time.time()
    benchmark_records = [
        _record_for_worker(row, nhood_size)
        for _, row in manifest_df.iterrows()
        for nhood_size in benchmark_nhood_sizes
    ]
    n_tasks = len(benchmark_records)
    n_jobs = max(1, min(int(n_benchmark_workers), n_tasks))
    print(f"Running {n_tasks} dataset x nhood tasks with n_jobs={n_jobs}.", flush=True)

    if n_jobs <= 1:
        results = [benchmark_one_dataset_nhood(rec) for rec in benchmark_records]
    else:
        results = Parallel(n_jobs=n_jobs, backend=parallel_backend,
                           verbose=parallel_verbose, batch_size=worker_batch_size)(
            delayed(benchmark_one_dataset_nhood)(rec) for rec in benchmark_records
        )

    summary_rows, distribution_rows, de_rows, missing_rows = [], [], [], []
    for res in results:
        if res.get("log"):
            print(res["log"], flush=True)
        summary_rows.extend(res.get("summary_rows", []))
        distribution_rows.extend(res.get("distribution_rows", []))
        de_rows.extend(res.get("de_rows", []))
        missing_rows.extend(res.get("missing_rows", []))

    summary_df = pd.DataFrame(summary_rows)
    distribution_df = pd.DataFrame(distribution_rows)
    de_df = pd.DataFrame(de_rows)
    missing_df = pd.DataFrame(missing_rows)

    summary_out = table_dir / "benchmark_summary_roc_pr_auc.csv"
    distribution_out = table_dir / "benchmark_distribution_long_roc_pr_auc.csv"
    de_out = table_dir / "de_gene_scores_responder_status_long.csv"
    missing_out = table_dir / "benchmark_missing_or_failed.csv"

    summary_df.to_csv(summary_out, index=False)
    distribution_df.to_csv(distribution_out, index=False)
    de_df.to_csv(de_out, index=False)
    missing_df.to_csv(missing_out, index=False)

    print("\n=== DONE ===", flush=True)
    print("Summary rows      :", len(summary_df), "->", summary_out, flush=True)
    print("Distribution rows :", len(distribution_df), "->", distribution_out, flush=True)
    print("DE-score rows     :", len(de_df), "->", de_out, flush=True)
    print("Missing/failed    :", len(missing_df), "->", missing_out, flush=True)
    print(f"Elapsed: {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
