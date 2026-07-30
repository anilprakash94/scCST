# Generate exactly five single-core runtime benchmark datasets with increasing cell counts
# Runtime sweep over base cell count only
# Fixed defaults: responder_percent=40, dropout_rate=0.4, harmony_theta=2
# One dataset per cell count; no replicates
#
# Important design:
# - Cell types are NOT assigned during simulation.
# - Hidden groups are used only internally to inject structure and truth labels.
# - Observed Cell_Type is derived from Leiden clusters after Harmony integration.
#
# Dropout design:
# - Dropout is applied after adding:
#   1. hidden-group identity genes
#   2. disease-response genes
#   3. batch-specific artifact genes
# - Dropout is NOT applied at the beginning.
#
# Boost design:
# - Disease-response, hidden-group identity, and batch-specific artifact
#   boosts are all sampled from disease_boost_range.

import os
# Avoid N_WORKERS x BLAS/Scanpy thread oversubscription when running simulations in parallel.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import random
import gc
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse
import anndata as ad
import scanpy as sc
from joblib import Parallel, delayed

# Keep per-simulation Scanpy work single-threaded. Parallelism happens across datasets.
sc.settings.n_jobs = 1

# -----------------------------
# Output directory
# -----------------------------
output_dir = Path(
    "/path/to/data/"
    "scrna_seq/simulation"
)
output_dir.mkdir(parents=True, exist_ok=True)


# -----------------------------
# Simulation grid
# -----------------------------
simulation_type = "de"

# Default parameters
default_responder_percent = 40
default_dropout_rate = 0.4
default_harmony_theta = 2

# One-parameter-at-a-time sweeps requested for benchmarking
runtime_base_cell_counts = [2000, 5000, 10000, 15000, 20000]

n_replicates = 1
base_seed = 123

# -----------------------------
# Multicore simulation settings
# -----------------------------
# Each worker generates and writes one full AnnData object. These simulations can be
# memory-heavy, so the default is conservative. Increase SIM_N_JOBS after checking RAM.
# Example: SIM_N_JOBS=4 when launching the notebook/kernel environment.
n_simulation_jobs = 1


# -----------------------------
# Core simulation settings
# -----------------------------
n_base_cells = 10000
n_genes = 20000
group_size = 1000

n_control_batches = 5
n_disease_batches = 5

identity_genes_per_group = 200
disease_genes_per_group = 100
batch_artifact_genes_per_batch = 100

disease_boost_range = (4.0, 10.0)
# Draw identity and batch-artifact boosts from the same range as disease boosts.
identity_boost_range = disease_boost_range
batch_artifact_boost_range = disease_boost_range

k_neighbors = 40
leiden_resolution = 0.8


# -----------------------------
# Helper functions
# -----------------------------
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)


def generate_noisy_nb_counts(
    n_cells,
    n_genes,
    mean_range=(1, 5),
    overdispersion_range=(2, 10),
):
    """
    Generate noisy negative-binomial count matrix.
    """
    means = np.random.uniform(mean_range[0], mean_range[1], size=n_genes)
    dispersions = np.random.uniform(
        overdispersion_range[0],
        overdispersion_range[1],
        size=n_genes,
    )

    counts = np.zeros((n_cells, n_genes), dtype=np.float32)

    for i in range(n_genes):
        m = means[i]
        v = m * dispersions[i]

        # Negative-binomial parameterization
        p = m / v
        n = (m * p) / (1 - p)

        counts[:, i] = np.random.negative_binomial(
            n=n,
            p=p,
            size=n_cells,
        )

    return counts


def apply_dropout(counts, dropout_rate=0.4):
    """
    Randomly zero out entries to mimic dropout.
    """
    mask = np.random.binomial(
        1,
        1 - dropout_rate,
        size=counts.shape,
    )
    return counts * mask


def allocate_gene_block(shuffled_genes, gene_cursor, n):
    """
    Allocate n genes from shuffled gene list.
    """
    genes = shuffled_genes[gene_cursor : gene_cursor + n]
    gene_cursor += n

    if len(genes) < n:
        raise ValueError(
            "Ran out of genes while allocating marker/artifact genes. "
            "Increase n_genes or reduce marker sizes."
        )

    return genes, gene_cursor


def make_safe_gene_dict(gene_dict):
    """
    Convert dict[str, list] into h5ad-safe dict[str, list[str]].
    """
    return {
        str(k): [str(x) for x in v]
        for k, v in gene_dict.items()
    }


def make_one_simulation(
    responder_percent,
    dropout_rate,
    harmony_theta,
    seed,
    sweep_name,
    replicate,
):
    set_seed(seed)

    print(
        f"\n--- Generating simulation: "
        f"sweep={sweep_name}, "
        f"replicate={replicate}, "
        f"responder={responder_percent}, "
        f"dropout={dropout_rate}, "
        f"harmony_theta={harmony_theta}, "
        f"seed={seed} ---"
    )

    # -----------------------------
    # Step 1: Generate base counts
    # -----------------------------
    nb_matrix = generate_noisy_nb_counts(
        n_cells=n_base_cells,
        n_genes=n_genes,
    )

    # IMPORTANT:
    # Do NOT apply dropout here.
    # Dropout is applied later after identity, disease,
    # and batch-specific signals have all been injected.
    control_adata = ad.AnnData(nb_matrix.astype(np.float32))
    control_adata.var_names = [f"Gene{i}" for i in range(n_genes)]

    n_cells = control_adata.shape[0]
    num_full_groups = n_cells // group_size

    group_indices_list = [
        np.arange(i * group_size, (i + 1) * group_size)
        for i in range(num_full_groups)
    ]

    # -----------------------------
    # Step 2: Gene allocation
    # -----------------------------
    shuffled_genes = random.sample(
        list(control_adata.var_names),
        len(control_adata.var_names),
    )

    gene_cursor = 0

    control_batches = [f"Ctrl_B{i + 1}" for i in range(n_control_batches)]
    disease_batches = [f"Dis_B{i + 1}" for i in range(n_disease_batches)]
    all_batches = control_batches + disease_batches

    batch_specific_genes = {}

    for batch in all_batches:
        genes, gene_cursor = allocate_gene_block(
            shuffled_genes,
            gene_cursor,
            batch_artifact_genes_per_batch,
        )
        batch_specific_genes[batch] = genes

    hidden_group_identity_genes = {}
    hidden_group_disease_markers = {}

    adatas_to_concat = []

    # -----------------------------
    # Step 3: Generate hidden groups
    # -----------------------------
    # IMPORTANT:
    # We do not assign Cell_Type here.
    # Hidden groups are only used to inject structured signals and truth labels.
    # Cell_Type will be assigned later from Leiden clusters.
    # -----------------------------
    for group_idx, g_idx in enumerate(group_indices_list):
        hidden_group = f"HiddenGroup_{group_idx}"

        sub_control = control_adata[g_idx].copy()
        sub_control.obs_names = [
            f"{hidden_group}_control_{i}"
            for i in range(sub_control.n_obs)
        ]

        sub_control.obs["condition"] = "Control"
        sub_control.obs["is_responder"] = False
        sub_control.obs["sim_group_truth"] = hidden_group

        sub_disease = sub_control.copy()
        sub_disease.obs_names = [
            f"{hidden_group}_disease_{i}"
            for i in range(sub_disease.n_obs)
        ]

        sub_disease.obs["condition"] = "Disease"
        sub_disease.obs["sim_group_truth"] = hidden_group

        # Assign batches
        sub_control.obs["sim_batch"] = np.random.choice(
            control_batches,
            size=sub_control.n_obs,
        )

        sub_disease.obs["sim_batch"] = np.random.choice(
            disease_batches,
            size=sub_disease.n_obs,
        )

        # Add hidden-group identity signal to both control and disease
        identity_genes, gene_cursor = allocate_gene_block(
            shuffled_genes,
            gene_cursor,
            identity_genes_per_group,
        )

        identity_idx = [
            sub_control.var_names.get_loc(g)
            for g in identity_genes
        ]

        control_identity_boosts = np.random.uniform(
            low=identity_boost_range[0],
            high=identity_boost_range[1],
            size=(sub_control.n_obs, len(identity_idx)),
        )

        disease_identity_boosts = np.random.uniform(
            low=identity_boost_range[0],
            high=identity_boost_range[1],
            size=(sub_disease.n_obs, len(identity_idx)),
        )

        sub_control.X[:, identity_idx] += (
            np.round(control_identity_boosts).astype(np.float32)
        )
        sub_disease.X[:, identity_idx] += (
            np.round(disease_identity_boosts).astype(np.float32)
        )

        hidden_group_identity_genes[hidden_group] = identity_genes

        # Add disease marker signal only to responder disease cells
        disease_genes, gene_cursor = allocate_gene_block(
            shuffled_genes,
            gene_cursor,
            disease_genes_per_group,
        )

        disease_idx = [
            sub_disease.var_names.get_loc(g)
            for g in disease_genes
        ]

        num_disease_cells = sub_disease.n_obs
        is_responder_mask = np.zeros(num_disease_cells, dtype=bool)
        responder_indices_list = []

        # Uniform responder allocation within each disease batch
        for dbatch in disease_batches:
            batch_indices = np.where(
                sub_disease.obs["sim_batch"].values == dbatch
            )[0]

            n_batch_responders = int(
                len(batch_indices) * responder_percent / 100
            )

            if n_batch_responders > 0:
                batch_responder_idx = np.random.choice(
                    batch_indices,
                    size=n_batch_responders,
                    replace=False,
                )
                responder_indices_list.extend(batch_responder_idx)

        responder_indices = np.array(responder_indices_list, dtype=int)

        if len(responder_indices) > 0:
            is_responder_mask[responder_indices] = True

            random_boosts = np.random.uniform(
                low=disease_boost_range[0],
                high=disease_boost_range[1],
                size=(len(responder_indices), len(disease_idx)),
            )

            sub_disease.X[np.ix_(responder_indices, disease_idx)] += (
                np.round(random_boosts).astype(np.float32)
            )

        sub_disease.obs["is_responder"] = is_responder_mask

        # Disease marker genes are absent in controls
        sub_control.X[:, disease_idx] = 0

        hidden_group_disease_markers[hidden_group] = disease_genes

        adatas_to_concat.extend([sub_control, sub_disease])

    # -----------------------------
    # Step 4: Concatenate
    # -----------------------------
    sim_adata = ad.concat(
        adatas_to_concat,
        axis=0,
        join="inner",
        merge="same",
    )

    # Make dense float32 for simple modification and preprocessing
    if scipy.sparse.issparse(sim_adata.X):
        sim_adata.X = sim_adata.X.toarray()

    sim_adata.X = sim_adata.X.astype(np.float32)

    # -----------------------------
    # Step 5: Inject batch-specific artifacts
    # -----------------------------
    for batch_name, b_genes in batch_specific_genes.items():
        b_idx = [
            sim_adata.var_names.get_loc(g)
            for g in b_genes
        ]

        batch_rows = np.where(
            sim_adata.obs["sim_batch"].values == batch_name
        )[0]

        if len(batch_rows) > 0:
            artifact_boosts = np.random.uniform(
                low=batch_artifact_boost_range[0],
                high=batch_artifact_boost_range[1],
                size=(len(batch_rows), len(b_idx)),
            )

            sim_adata.X[np.ix_(batch_rows, b_idx)] += (
                np.round(artifact_boosts).astype(np.float32)
            )

    # -----------------------------
    # Step 6: Apply dropout after all simulated signals
    # -----------------------------
    sim_adata.X = apply_dropout(
        sim_adata.X,
        dropout_rate=dropout_rate,
    ).astype(np.float32)

    # Store raw counts after signal injection and dropout,
    # before normalization/log transformation.
    sim_adata.layers["counts"] = sim_adata.X.copy()

    # -----------------------------
    # Step 7: Store h5ad-safe simulation truth
    # -----------------------------
    sim_adata.uns["simulation_params"] = {
        "simulation_type": str(simulation_type),
        "sweep_name": str(sweep_name),
        "replicate": int(replicate),
        "responder_percent": int(responder_percent),
        "dropout_rate": float(dropout_rate),
        "harmony_theta": float(harmony_theta),
        "seed": int(seed),
        "n_base_cells": int(n_base_cells),
        "n_genes": int(n_genes),
        "group_size": int(group_size),
        "n_control_batches": int(n_control_batches),
        "n_disease_batches": int(n_disease_batches),
        "identity_genes_per_group": int(identity_genes_per_group),
        "disease_genes_per_group": int(disease_genes_per_group),
        "batch_artifact_genes_per_batch": int(batch_artifact_genes_per_batch),
        "disease_boost_min": float(disease_boost_range[0]),
        "disease_boost_max": float(disease_boost_range[1]),
        "identity_boost_min": float(identity_boost_range[0]),
        "identity_boost_max": float(identity_boost_range[1]),
        "batch_artifact_boost_min": float(batch_artifact_boost_range[0]),
        "batch_artifact_boost_max": float(batch_artifact_boost_range[1]),
        "k_neighbors": int(k_neighbors),
        "leiden_resolution": float(leiden_resolution),
    }

    sim_adata.uns["hidden_group_identity_genes"] = make_safe_gene_dict(
        hidden_group_identity_genes
    )

    sim_adata.uns["hidden_group_disease_markers"] = make_safe_gene_dict(
        hidden_group_disease_markers
    )

    sim_adata.uns["batch_genes"] = make_safe_gene_dict(
        batch_specific_genes
    )

    print(
        f"Simulation count matrix complete: "
        f"{sim_adata.n_obs} cells x {sim_adata.n_vars} genes"
    )

    print(
        "Total disease responders:",
        int(sim_adata.obs["is_responder"].sum()),
    )

    # -----------------------------
    # Step 8: Preprocessing
    # -----------------------------
    sc.pp.normalize_total(sim_adata)
    sc.pp.log1p(sim_adata)

    sim_adata.raw = sim_adata

    sc.pp.highly_variable_genes(
        sim_adata,
        n_top_genes=n_genes,
        flavor="seurat",
    )

    sc.tl.pca(
        sim_adata,
        svd_solver="arpack",
        use_highly_variable=True,
    )

    # Harmony integration
    sc.external.pp.harmony_integrate(
        sim_adata,
        key="sim_batch",
        theta=harmony_theta,
        max_iter_harmony=50,
    )

    # -----------------------------
    # Step 9: Leiden-derived Cell_Type
    # -----------------------------
    sc.pp.neighbors(
        sim_adata,
        n_neighbors=k_neighbors,
        use_rep="X_pca_harmony",
    )

    sc.tl.umap(sim_adata)

    sc.tl.leiden(
        sim_adata,
        resolution=leiden_resolution,
        key_added="leiden",
        random_state=seed,
    )

    # Observed cell-type label used by scCST
    sim_adata.obs["Cell_Type"] = (
        "Leiden_" + sim_adata.obs["leiden"].astype(str)
    )

    # Save cluster composition externally to avoid nested uns write issues
    cluster_composition = (
        sim_adata.obs
        .groupby(["Cell_Type", "sim_group_truth"], observed=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["Cell_Type", "n_cells"], ascending=[True, False])
    )

    dropout_tag = str(dropout_rate).replace(".", "p")
    theta_tag = str(harmony_theta).replace(".", "p")
    sweep_tag = str(sweep_name).replace(" ", "_")

    composition_file = (
        output_dir
        / f"cluster_composition_runtime_base_cells_{n_base_cells}"
          f"_{sweep_tag}"
          f"_replicate_{replicate}"
          f"_responder_percent_{responder_percent}"
          f"_dropout_{dropout_tag}"
          f"_harmony_theta_{theta_tag}"
          f".csv"
    )

    cluster_composition.to_csv(composition_file, index=False)

    sim_adata.uns["cluster_composition_file"] = str(composition_file)

    print("Leiden-derived Cell_Type counts:")
    print(sim_adata.obs["Cell_Type"].value_counts().sort_index())

    return sim_adata, composition_file


# -----------------------------
# Build exactly five runtime datasets
# -----------------------------
simulation_plan = [
    {
        "sweep_name": "runtime_cell_count",
        "swept_parameter": "base_n_cells",
        "swept_value": int(n_cells),
        "replicate": 1,
        "responder_percent": default_responder_percent,
        "dropout_rate": default_dropout_rate,
        "harmony_theta": default_harmony_theta,
        "base_n_cells": int(n_cells),
    }
    for n_cells in runtime_base_cell_counts
]

print(f"Prepared {len(simulation_plan)} runtime simulations.")
print(pd.DataFrame(simulation_plan))


# -----------------------------
# Run requested simulations in parallel
# -----------------------------
def _safe_tag(value):
    return str(value).replace(".", "p")


def run_one_simulation_from_plan(dataset_counter, params):
    """Generate one simulation, write files inside the worker, and return one manifest row."""
    seed = base_seed + dataset_counter

    try:
        global n_base_cells
        n_base_cells = int(params["base_n_cells"])
        sim_adata, composition_file = make_one_simulation(
            responder_percent=params["responder_percent"],
            dropout_rate=params["dropout_rate"],
            harmony_theta=params["harmony_theta"],
            seed=seed,
            sweep_name=params["sweep_name"],
            replicate=params["replicate"],
        )

        dropout_tag = _safe_tag(params["dropout_rate"])
        theta_tag = _safe_tag(params["harmony_theta"])
        sweep_tag = str(params["sweep_name"]).replace(" ", "_")

        out_file = (
            output_dir
            / f"simulation_adata_runtime_base_cells_{params['base_n_cells']}"
              f"_responder_percent_{params['responder_percent']}"
              f"_dropout_{dropout_tag}"
              f"_harmony_theta_{theta_tag}"
              f".h5ad"
        )

        sim_adata.write_h5ad(out_file, compression="gzip")

        summary_row = {
            "file": str(out_file),
            "cluster_composition_file": str(composition_file),
            "simulation_type": simulation_type,
            "sweep_name": str(params["sweep_name"]),
            "swept_parameter": str(params["swept_parameter"]),
            "swept_value": params["swept_value"],
            "replicate": int(params["replicate"]),
            "responder_percent": int(params["responder_percent"]),
            "dropout_rate": float(params["dropout_rate"]),
            "harmony_theta": float(params["harmony_theta"]),
            "seed": int(seed),
            "base_n_cells": int(params["base_n_cells"]),
            "n_cells": int(sim_adata.n_obs),
            "n_genes": int(sim_adata.n_vars),
            "n_responders": int(sim_adata.obs["is_responder"].sum()),
            "n_leiden_clusters": int(sim_adata.obs["Cell_Type"].nunique()),
            "status": "success",
            "error": "",
        }

        print(f"Saved h5ad: {out_file}")
        print(f"Saved cluster composition: {composition_file}")

        del sim_adata
        gc.collect()
        return summary_row

    except Exception as e:
        error_row = {
            "file": "",
            "cluster_composition_file": "",
            "simulation_type": simulation_type,
            "sweep_name": str(params.get("sweep_name", "")),
            "swept_parameter": str(params.get("swept_parameter", "")),
            "swept_value": params.get("swept_value", np.nan),
            "replicate": int(params.get("replicate", -1)),
            "responder_percent": int(params.get("responder_percent", -1)),
            "dropout_rate": float(params.get("dropout_rate", np.nan)),
            "harmony_theta": float(params.get("harmony_theta", np.nan)),
            "seed": int(seed),
            "n_cells": np.nan,
            "n_genes": np.nan,
            "n_responders": np.nan,
            "n_leiden_clusters": np.nan,
            "status": "failed",
            "error": repr(e),
            "traceback": traceback.format_exc(),
        }
        print("FAILED simulation:", error_row)
        gc.collect()
        return error_row


print(f"Running {len(simulation_plan)} simulations serially with one core")
print("Tip: set SIM_N_JOBS=1,2,4,... before launching the kernel to tune parallelism.")

if True:
    summary_rows = [
        run_one_simulation_from_plan(i, params)
        for i, params in enumerate(simulation_plan, start=1)
    ]
else:
    # Use loky/cloudpickle so notebook-defined functions are serializable.
    summary_rows = Parallel(
        n_jobs=n_simulation_jobs,
        backend="loky",
        verbose=10,
    )(
        delayed(run_one_simulation_from_plan)(i, params)
        for i, params in enumerate(simulation_plan, start=1)
    )


# -----------------------------
# Save simulation manifest
# -----------------------------
summary_df = pd.DataFrame(summary_rows)

manifest_file = output_dir / "simulation_manifest_runtime_5_cell_counts.csv"
summary_df.to_csv(manifest_file, index=False)

print("\nAll simulations complete.")
print(f"Manifest saved to: {manifest_file}")
print(summary_df)