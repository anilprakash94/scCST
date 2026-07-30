#!/usr/bin/env python
# coding: utf-8
"""Parallel fast export of simulation AnnData files for CACOA and miloDE.

This is a multicore version of convert_parameter_sweeps_adata_to_csv_fast.py.

It exports only the files consumed by the current CACOA/miloDE R scripts:
  - cells x genes expression TSV
  - metadata CSV
  - gene annotation CSV
  - export manifest

It does not create latent_csv, genes-by-cells TSV, or other optional files because
those are not used by the current CACOA/miloDE scripts and add substantial I/O.

Run:
  python convert_parameter_sweeps_adata_to_csv_fast_multicore.py

Tune worker count:
  EXPORT_N_WORKERS=4 python convert_parameter_sweeps_adata_to_csv_fast_multicore.py

Memory note:
  Each worker may materialize one dense cells x genes matrix before writing TSV.
  For 10k-20k cells x ~2k genes, start with 2-4 workers, then increase if RAM and
  disk bandwidth are not saturated.
"""

from __future__ import annotations

# Keep BLAS/HDF5-backed libraries from spawning many threads inside each worker.
# This avoids N_WORKERS x BLAS_THREADS oversubscription.
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import re
import time
import gc
import traceback

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# ------------------------------------------------------------
# Input/output directories
# ------------------------------------------------------------
data_dir = Path(
    "/path/to/data/"
    "scrna_seq/simulation"
)

manifest_file = data_dir / "simulation_manifest_de_parameter_sweeps.csv"

export_dir = data_dir / "cacoa_milode_exports_parameter_sweeps"
expr_dir = export_dir / "expression_tsv"
obs_dir = export_dir / "metadata_csv"
var_dir = export_dir / "gene_annotations_csv"

# Fast default: skip rewriting outputs that already exist.
OVERWRITE = False

# Default to a conservative worker count because each worker can hold a dense matrix.
# Override with EXPORT_N_WORKERS=<int>.
DEFAULT_WORKERS = min(4, os.cpu_count() or 1)
N_WORKERS = int(os.environ.get("EXPORT_N_WORKERS", DEFAULT_WORKERS))
N_WORKERS = max(1, N_WORKERS)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def parse_float_tag(value):
    return float(str(value).replace("p", "."))


def float_to_tag(value):
    text = str(value)
    if text.endswith(".0"):
        text = text[:-2]
    return text.replace(".", "p")


def safe_token(x):
    try:
        if pd.isna(x):
            return "NA"
    except Exception:
        pass
    return re.sub(r"[^A-Za-z0-9_-]+", "_", float_to_tag(x))


def parse_parameter_sweep_filename(path):
    name = Path(path).name
    pattern = (
        r"simulation_adata_(?P<simulation_type>.+?)_"
        r"(?P<sweep_name>(?:responder_percent|dropout_rate|harmony_theta)_sweep)_"
        r"replicate_(?P<replicate>\d+)_"
        r"responder_percent_(?P<responder_percent>\d+)_"
        r"dropout_(?P<dropout>[\dp]+)_"
        r"harmony_theta_(?P<harmony_theta>[\dp]+)"
        r"\.h5ad$"
    )
    m = re.match(pattern, name)
    if m is None:
        # Runtime benchmark filenames or other simple files.
        base_cells_match = re.search(r"(?:base_cells|base_n_cells|n_cells)_(\d+)", name)
        base_cells = int(base_cells_match.group(1)) if base_cells_match else np.nan
        return {
            "simulation_type": "runtime" if "runtime" in name else "unknown",
            "sweep_name": "runtime_cell_count" if not pd.isna(base_cells) else "unknown",
            "swept_parameter": "base_n_cells" if not pd.isna(base_cells) else "unknown",
            "swept_value": base_cells,
            "replicate": 1,
            "responder_percent": np.nan,
            "dropout_rate": np.nan,
            "harmony_theta": np.nan,
        }

    sweep_name = m.group("sweep_name")
    swept_parameter = sweep_name.replace("_sweep", "")
    responder_percent = int(m.group("responder_percent"))
    dropout_rate = parse_float_tag(m.group("dropout"))
    harmony_theta = parse_float_tag(m.group("harmony_theta"))
    swept_value = {
        "responder_percent": responder_percent,
        "dropout_rate": dropout_rate,
        "harmony_theta": harmony_theta,
    }[swept_parameter]

    return {
        "simulation_type": m.group("simulation_type"),
        "sweep_name": sweep_name,
        "swept_parameter": swept_parameter,
        "swept_value": swept_value,
        "replicate": int(m.group("replicate")),
        "responder_percent": responder_percent,
        "dropout_rate": dropout_rate,
        "harmony_theta": harmony_theta,
    }


def params_from_manifest_row_or_filename(h5ad_file: str, manifest_row: dict | None = None):
    parsed = parse_parameter_sweep_filename(h5ad_file)
    if manifest_row is None:
        return parsed

    out = parsed.copy()
    for key in [
        "simulation_type", "sweep_name", "swept_parameter", "swept_value", "replicate",
        "responder_percent", "dropout_rate", "harmony_theta", "seed", "cluster_composition_file",
        "base_n_cells", "n_cells", "base_cells",
    ]:
        if key in manifest_row and not pd.isna(manifest_row[key]):
            out[key] = manifest_row[key]

    # Fill runtime fields when present.
    for base_key in ["base_n_cells", "base_cells", "n_cells"]:
        if base_key in out and not pd.isna(out[base_key]) and out.get("swept_parameter", "unknown") == "unknown":
            out["simulation_type"] = "runtime"
            out["sweep_name"] = "runtime_cell_count"
            out["swept_parameter"] = "base_n_cells"
            out["swept_value"] = int(out[base_key])
            out["replicate"] = 1
            break

    for key in ["replicate", "responder_percent", "seed"]:
        if key in out and not pd.isna(out[key]):
            out[key] = int(out[key])
    for key in ["swept_value", "dropout_rate", "harmony_theta"]:
        if key in out and not pd.isna(out[key]):
            out[key] = float(out[key])
    return out


def tag_from_params(params):
    sim = safe_token(params.get("simulation_type", "unknown"))
    sweep = safe_token(params.get("sweep_name", "unknown"))
    param = safe_token(params.get("swept_parameter", "unknown"))
    value = safe_token(params.get("swept_value", "NA"))
    rep = int(params.get("replicate", 1) if not pd.isna(params.get("replicate", 1)) else 1)
    return f"{sim}_{sweep}_{param}_{value}_replicate_{rep}"


def make_unique_names(names):
    seen = {}
    out = []
    for name in names:
        name = str(name).replace(".", "-")
        if name not in seen:
            seen[name] = 0
            out.append(name)
        else:
            seen[name] += 1
            out.append(f"{name}-{seen[name]}")
    return out


def stringify_categories(df):
    out = df.copy()
    for col in out.columns:
        if isinstance(out[col].dtype, pd.CategoricalDtype) or out[col].dtype == "object":
            out[col] = out[col].astype(str)
    return out


def get_counts_matrix(adata):
    if "counts" in adata.layers:
        return adata.layers["counts"], "layers['counts']"
    if adata.raw is not None:
        return adata.raw.X, "raw.X"
    return adata.X, "X"


def add_alias_columns(obs_df):
    out = obs_df.copy()
    if "condition" in out.columns:
        out["group"] = out["condition"].astype(str)
        out["condition_label"] = out["condition"].astype(str)
    if "sim_batch" in out.columns:
        out["batch"] = out["sim_batch"].astype(str)
        out["sample"] = out["sim_batch"].astype(str)
        out["sample_id"] = out["sim_batch"].astype(str)
    if "Cell_Type" in out.columns:
        out["cell_type"] = out["Cell_Type"].astype(str)
        out["cluster"] = out["Cell_Type"].astype(str)
    return out


def write_expression_tsv(X, cell_names, gene_names, path: Path, overwrite: bool):
    """Write cells x genes TSV. Current R scripts consume this format."""
    if path.exists() and not overwrite:
        return "skipped_existing"

    if sp.issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)

    # Counts are usually integer. Integer text is smaller/faster than float text.
    if np.issubdtype(X.dtype, np.floating) and np.all(np.isfinite(X)) and np.allclose(X, np.rint(X)):
        X = X.astype(np.int32, copy=False)

    df = pd.DataFrame(X, index=pd.Index(cell_names, name="cell_barcode"), columns=gene_names)
    df.to_csv(path, sep="\t", index=True)
    del df
    return "written"



def export_one(task: dict):
    """Worker function. Returns one manifest row or one failure row."""
    h5ad_file = Path(task["h5ad_file"])
    manifest_row = task.get("manifest_row")
    overwrite = bool(task.get("overwrite", False))

    t0 = time.perf_counter()
    try:
        if not h5ad_file.exists():
            raise FileNotFoundError(f"Missing input h5ad: {h5ad_file}")

        params = params_from_manifest_row_or_filename(str(h5ad_file), manifest_row)
        file_tag = tag_from_params(params)

        expr_file = expr_dir / f"sim_adata_expr_{file_tag}.tsv"
        obs_file = obs_dir / f"sim_adata_obs_{file_tag}.csv"
        var_file = var_dir / f"sim_adata_var_{file_tag}.csv"

        required_exist = expr_file.exists() and obs_file.exists() and var_file.exists()
        if required_exist and not overwrite:
            elapsed = time.perf_counter() - t0
            return {
                "ok": True,
                "row": {
                    "input_h5ad": str(h5ad_file),
                    "simulation_type": params.get("simulation_type", "unknown"),
                    "sweep_name": params.get("sweep_name", "unknown"),
                    "swept_parameter": params.get("swept_parameter", "unknown"),
                    "swept_value": params.get("swept_value", np.nan),
                    "replicate": params.get("replicate", 1),
                    "responder_percent": params.get("responder_percent", np.nan),
                    "dropout_rate": params.get("dropout_rate", np.nan),
                    "harmony_theta": params.get("harmony_theta", np.nan),
                    "seed": params.get("seed", np.nan),
                    "base_n_cells": params.get("base_n_cells", params.get("base_cells", np.nan)),
                    "counts_source": "not_reloaded_existing_output",
                    "expression_cells_by_genes_tsv": str(expr_file),
                    "metadata_csv": str(obs_file),
                    "gene_annotations_csv": str(var_file),
                    "n_cells": np.nan,
                    "n_genes": np.nan,
                    "export_seconds": elapsed,
                    "worker_pid": os.getpid(),
                    "status": "skipped_existing",
                },
            }

        adata = sc.read_h5ad(h5ad_file)
        gene_names = make_unique_names(adata.var_names)
        cell_names = adata.obs_names.astype(str)
        X_counts, counts_source = get_counts_matrix(adata)

        expression_status = write_expression_tsv(X_counts, cell_names, gene_names, expr_file, overwrite=overwrite)

        obs_df = stringify_categories(adata.obs)
        obs_df.index = obs_df.index.astype(str)
        obs_df.index.name = "cell_barcode"
        obs_df = add_alias_columns(obs_df)
        for key in [
            "simulation_type", "sweep_name", "swept_parameter", "swept_value", "replicate",
            "responder_percent", "dropout_rate", "harmony_theta", "seed", "base_n_cells",
        ]:
            if key in params:
                obs_df[key] = params[key]
        if (not obs_file.exists()) or overwrite:
            obs_df.to_csv(obs_file)
            metadata_status = "written"
        else:
            metadata_status = "skipped_existing"

        var_df = stringify_categories(adata.var)
        var_df.index = pd.Index(gene_names, name="gene")
        var_df["gene"] = gene_names
        if (not var_file.exists()) or overwrite:
            var_df.to_csv(var_file)
            var_status = "written"
        else:
            var_status = "skipped_existing"


        n_cells = int(adata.n_obs)
        n_genes = int(adata.n_vars)
        del adata, obs_df, var_df
        gc.collect()

        elapsed = time.perf_counter() - t0
        return {
            "ok": True,
            "row": {
                "input_h5ad": str(h5ad_file),
                "simulation_type": params.get("simulation_type", "unknown"),
                "sweep_name": params.get("sweep_name", "unknown"),
                "swept_parameter": params.get("swept_parameter", "unknown"),
                "swept_value": params.get("swept_value", np.nan),
                "replicate": params.get("replicate", 1),
                "responder_percent": params.get("responder_percent", np.nan),
                "dropout_rate": params.get("dropout_rate", np.nan),
                "harmony_theta": params.get("harmony_theta", np.nan),
                "seed": params.get("seed", np.nan),
                "base_n_cells": params.get("base_n_cells", params.get("base_cells", np.nan)),
                "counts_source": counts_source,
                "expression_cells_by_genes_tsv": str(expr_file),
                "metadata_csv": str(obs_file),
                "gene_annotations_csv": str(var_file),
                "n_cells": n_cells,
                "n_genes": n_genes,
                "export_seconds": elapsed,
                "worker_pid": os.getpid(),
                "status": "completed",
                "expression_status": expression_status,
                "metadata_status": metadata_status,
                "var_status": var_status,
            },
        }
    except Exception as e:
        elapsed = time.perf_counter() - t0
        return {
            "ok": False,
            "row": {
                "input_h5ad": str(h5ad_file),
                "reason": repr(e),
                "traceback": traceback.format_exc(),
                "export_seconds": elapsed,
                "worker_pid": os.getpid(),
            },
        }


def find_h5ad_tasks():
    if manifest_file.exists():
        manifest_df = pd.read_csv(manifest_file)
        if "file" not in manifest_df.columns:
            raise ValueError(f"Manifest exists but does not contain a 'file' column: {manifest_file}")

        tasks = []
        for _, row in manifest_df.iterrows():
            p = Path(row["file"])
            if not p.exists():
                repaired = data_dir / p.name
                if repaired.exists():
                    p = repaired
            tasks.append({
                "h5ad_file": str(p),
                "manifest_row": row.to_dict(),
                "overwrite": OVERWRITE,
            })
        return tasks

    print("Manifest not found. Falling back to filename glob.")
    return [
        {"h5ad_file": str(p), "manifest_row": None, "overwrite": OVERWRITE}
        for p in sorted(data_dir.glob("simulation_adata*.h5ad"))
    ]


def main():
    for d in [export_dir, expr_dir, obs_dir, var_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print(f"Data directory:   {data_dir}")
    print(f"Manifest file:    {manifest_file}")
    print(f"Export directory: {export_dir}")
    print(f"N_WORKERS={N_WORKERS}")
    print("Exporting only files required by CACOA/miloDE: expression TSV, metadata CSV, gene annotation CSV, manifest")
    print(f"OVERWRITE={OVERWRITE}")

    tasks = find_h5ad_tasks()
    print(f"Found {len(tasks)} simulated h5ad files.")
    if len(tasks) == 0:
        raise FileNotFoundError(f"No .h5ad files found in {data_dir}")

    start_all = time.perf_counter()
    export_rows = []
    failed_rows = []

    if N_WORKERS == 1:
        for i, task in enumerate(tasks, start=1):
            result = export_one(task)
            if result["ok"]:
                export_rows.append(result["row"])
                print(f"[{i}/{len(tasks)}] OK {Path(task['h5ad_file']).name} in {result['row']['export_seconds']:.2f}s")
            else:
                failed_rows.append(result["row"])
                print(f"[{i}/{len(tasks)}] FAILED {Path(task['h5ad_file']).name}: {result['row']['reason']}")
    else:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
            future_to_task = {ex.submit(export_one, task): task for task in tasks}
            for i, fut in enumerate(as_completed(future_to_task), start=1):
                task = future_to_task[fut]
                result = fut.result()
                if result["ok"]:
                    export_rows.append(result["row"])
                    print(f"[{i}/{len(tasks)}] OK {Path(task['h5ad_file']).name} in {result['row']['export_seconds']:.2f}s")
                else:
                    failed_rows.append(result["row"])
                    print(f"[{i}/{len(tasks)}] FAILED {Path(task['h5ad_file']).name}: {result['row']['reason']}")

    export_manifest = pd.DataFrame(export_rows)
    if not export_manifest.empty:
        sort_cols = [c for c in ["sweep_name", "swept_parameter", "swept_value", "replicate", "input_h5ad"] if c in export_manifest.columns]
        export_manifest = export_manifest.sort_values(sort_cols).reset_index(drop=True)

    export_manifest_file = export_dir / "cacoa_milode_export_manifest_parameter_sweeps.csv"
    export_manifest.to_csv(export_manifest_file, index=False)

    failed_df = pd.DataFrame(failed_rows)
    failed_file = export_dir / "cacoa_milode_export_failed_parameter_sweeps.csv"
    failed_df.to_csv(failed_file, index=False)

    total_elapsed = time.perf_counter() - start_all
    print("\nAll exports complete.")
    print(f"Successful exports: {len(export_manifest)}")
    print(f"Failed/skipped files: {len(failed_df)}")
    print(f"Total elapsed seconds: {total_elapsed:.2f}")
    print(f"Export manifest saved to: {export_manifest_file}")
    print(f"Failed log saved to: {failed_file}")

    if not export_manifest.empty:
        group_cols = [c for c in ["sweep_name", "swept_parameter", "swept_value"] if c in export_manifest.columns]
        if group_cols:
            print("\nExports per dataset group:")
            print(export_manifest.groupby(group_cols).size().reset_index(name="n_exports").to_string(index=False))


if __name__ == "__main__":
    main()