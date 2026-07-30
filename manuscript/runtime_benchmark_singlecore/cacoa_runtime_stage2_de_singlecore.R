#!/usr/bin/env Rscript
# ============================================================
# CACOA Stage 2: load Stage-1 HVG+Harmony Seurat RDS files and
# run FindNeighbors + CACOA cluster-free DE for each nhood size.
#
# Parallelism:
#   - outer workers: multiple dataset/nhood-size runs in parallel
#   - n_cores: cores passed to cao$estimateClusterFreeDE()
#
# Keep outer_workers * n_cores within available CPU/RAM.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(cacoa)
  library(dplyr)
  library(glue)
  library(future)
  library(parallel)
})

# -----------------------------
# Defaults
# -----------------------------
base_dir <- "/path/to/data/scrna_seq/simulation"
preprocess_dir <- file.path(base_dir, "cacoa_runtime_hvg_harmony_rds")
preprocess_manifest <- file.path(preprocess_dir, "cacoa_runtime_preprocess_manifest.csv")

# Keep this default compatible with the existing benchmark code.
results_dir <- file.path(base_dir, "cacoa_runtime_5_cell_counts_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

nhood_sizes <- c(20)
harmony_dims <- 1:30
ref.level <- "Control"
target.level <- "Disease"
condition_col <- "condition"
sample_id_col <- "sim_batch"
cell_type_col <- "Cell_Type"
reduction_name <- "harmony"
min_expr_frac <- 0.0001
adjust_pvalues <- TRUE
smooth_de <- TRUE

outer_workers <- 1L
outer_backend <- "serial"
n_cores <- 1L

# Keep future sequential by default to avoid exporting multi-GB Seurat/CACOA objects.
future_plan_mode <- "sequential" # sequential or multisession for tiny tests only
future_max_size_gb <- 100

overwrite <- FALSE
test_mode <- FALSE
max_runs_for_test <- 2
save_cacoa_rds <- FALSE
save_hvg_file <- TRUE
save_metadata <- FALSE

# Optional filters
run_only_sweep_names <- NULL
run_only_swept_parameters <- NULL
run_only_swept_values <- NULL
run_only_replicates <- NULL

# -----------------------------
# CLI parsing
# -----------------------------
parse_cli_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(prefix) {
    hit <- args[startsWith(args, prefix)]
    if (length(hit) == 0) return(NULL)
    sub(prefix, "", hit[[1]], fixed = TRUE)
  }
  split_arg <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  parse_bool <- function(x) tolower(x) %in% c("true", "t", "1", "yes", "y")

  x <- get_arg("--base-dir="); if (!is.null(x)) base_dir <<- x
  x <- get_arg("--preprocess-dir="); if (!is.null(x)) preprocess_dir <<- x
  x <- get_arg("--preprocess-manifest="); if (!is.null(x)) preprocess_manifest <<- x
  x <- get_arg("--results-dir="); if (!is.null(x)) results_dir <<- x
  x <- get_arg("--nhood-sizes="); if (!is.null(x) && nzchar(x)) nhood_sizes <<- as.integer(split_arg(x))
  x <- get_arg("--outer-workers="); if (!is.null(x)) outer_workers <<- as.integer(x)
  x <- get_arg("--outer-backend="); if (!is.null(x)) outer_backend <<- x
  x <- get_arg("--n-cores="); if (!is.null(x)) n_cores <<- as.integer(x)
  x <- get_arg("--future-plan="); if (!is.null(x)) future_plan_mode <<- x
  x <- get_arg("--future-max-size-gb="); if (!is.null(x)) future_max_size_gb <<- as.numeric(x)
  x <- get_arg("--overwrite="); if (!is.null(x)) overwrite <<- parse_bool(x)
  x <- get_arg("--test-mode="); if (!is.null(x)) test_mode <<- parse_bool(x)
  x <- get_arg("--save-cacoa-rds="); if (!is.null(x)) save_cacoa_rds <<- parse_bool(x)
  x <- get_arg("--save-hvg-file="); if (!is.null(x)) save_hvg_file <<- parse_bool(x)
  x <- get_arg("--save-metadata="); if (!is.null(x)) save_metadata <<- parse_bool(x)
  invisible(NULL)
}

safe_value_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.0$", "", x)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x
}

make_run_tag <- function(row, nhood_size) {
  dataset_id <- if ("dataset_id" %in% names(row) && !is.na(row[["dataset_id"]]) && nzchar(as.character(row[["dataset_id"]]))) {
    safe_value_tag(row[["dataset_id"]])
  } else {
    paste0(
      as.character(row[["simulation_type"]]),
      "_", as.character(row[["sweep_name"]]),
      "_", as.character(row[["swept_parameter"]]), "_", safe_value_tag(row[["swept_value"]]),
      "_replicate_", row[["replicate"]],
      "_responder_percent_", row[["responder_percent"]],
      "_dropout_", safe_value_tag(row[["dropout_rate"]]),
      "_harmony_theta_", safe_value_tag(row[["harmony_theta"]])
    )
  }
  paste0("nhood_", nhood_size, "_", dataset_id)
}

choose_nn_graph <- function(so) {
  graph_names <- names(so@graphs)
  if ("RNA_nn" %in% graph_names) return("RNA_nn")
  nn_graphs <- graph_names[grepl("_nn$", graph_names)]
  if (length(nn_graphs) > 0) return(nn_graphs[1])
  if (length(graph_names) > 0) {
    warning("No *_nn graph found. Using first graph: ", graph_names[1])
    return(graph_names[1])
  }
  stop("No Seurat neighbor graph found in so@graphs.")
}

timed <- function(label, expr) {
  cat("[TIMER START]", label, "\n")
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  seconds <- round(proc.time()[["elapsed"]] - t0, 3)
  cat("[TIMER DONE]", label, ":", seconds, "seconds\n")
  list(value = value, seconds = seconds)
}

setup_future <- function() {
  if (identical(future_plan_mode, "multisession")) {
    future::plan("multisession", workers = n_cores)
  } else {
    future_plan_mode <<- "sequential"
    future::plan("sequential")
  }
  options(future.globals.maxSize = future_max_size_gb * 1024^3)
}

load_preprocess_manifest <- function() {
  if (!file.exists(preprocess_manifest)) stop("Preprocess manifest not found: ", preprocess_manifest)
  df <- read.csv(preprocess_manifest, stringsAsFactors = FALSE, check.names = FALSE)

  if (!("hvg_seurat_rds" %in% colnames(df))) stop("Preprocess manifest is missing hvg_seurat_rds.")
  status_col <- if ("preprocess_status" %in% colnames(df)) "preprocess_status" else if ("status.1" %in% colnames(df)) "status.1" else if ("status" %in% colnames(df)) "status" else NA_character_
  if (is.na(status_col)) stop("Preprocess manifest has no status column.")
  cat("Using preprocessing status column:", status_col, "\n")

  df <- df[df[[status_col]] %in% c("success", "skipped_existing"), , drop = FALSE]
  df <- df[file.exists(df$hvg_seurat_rds), , drop = FALSE]

  if (!is.null(run_only_sweep_names)) df <- df[df$sweep_name %in% run_only_sweep_names, , drop = FALSE]
  if (!is.null(run_only_swept_parameters)) df <- df[df$swept_parameter %in% run_only_swept_parameters, , drop = FALSE]
  if (!is.null(run_only_swept_values)) df <- df[df$swept_value %in% run_only_swept_values, , drop = FALSE]
  if (!is.null(run_only_replicates)) df <- df[df$replicate %in% run_only_replicates, , drop = FALSE]

  if (nrow(df) == 0) stop("No successful preprocessed CACOA datasets found in: ", preprocess_manifest)
  df
}

run_one <- function(task) {
  row <- as.list(task$row)
  nhood_size <- as.integer(task$nhood_size)
  run_tag <- make_run_tag(row, nhood_size)
  hvg_seurat_rds <- as.character(row$hvg_seurat_rds)

  z_file <- file.path(results_dir, glue("cacoa_cluster_free_de_z_{run_tag}.csv"))
  rds_file <- file.path(results_dir, glue("cacoa_results_{run_tag}.rds"))
  hvg_file <- file.path(results_dir, glue("cacoa_hvg_genes_{run_tag}.csv"))
  meta_file <- file.path(results_dir, glue("cacoa_cell_metadata_{run_tag}.csv"))

  if (file.exists(z_file) && !overwrite) {
    return(data.frame(as.data.frame(row), nhood_size = nhood_size, run_tag = run_tag, run_status = "skipped_existing", error = NA_character_, z_file = z_file, rds_file = if (file.exists(rds_file)) rds_file else NA_character_, stringsAsFactors = FALSE))
  }

  t_all <- proc.time()[["elapsed"]]
  tryCatch({
    cat("\n============================================================\n")
    cat("Running CACOA from HVG/Harmony RDS:", run_tag, "\n")
    cat("============================================================\n")

    if (!file.exists(hvg_seurat_rds)) stop("Missing hvg_seurat_rds: ", hvg_seurat_rds)
    x <- timed("read_hvg_harmony_seurat_rds", readRDS(hvg_seurat_rds)); so <- x$value; read_rds_seconds <- x$seconds
    cat("HVG Seurat object:", nrow(so), "genes x", ncol(so), "cells\n")

    if (!(reduction_name %in% names(so@reductions))) stop("Reduction not found in Seurat object: ", reduction_name)
    dims_available <- seq_len(ncol(Embeddings(so, reduction_name)))
    dims_use <- harmony_dims[harmony_dims %in% dims_available]
    if (length(dims_use) == 0) stop("No valid Harmony dimensions available for FindNeighbors.")

    x <- timed("find_neighbors", FindNeighbors(
      so,
      reduction = reduction_name,
      dims = dims_use,
      k.param = nhood_size,
      verbose = FALSE
    )); so <- x$value; find_neighbors_seconds <- x$seconds

    graph_name <- choose_nn_graph(so)
    cat("Using graph:", graph_name, "\n")

    cell.groups <- setNames(as.character(so@meta.data[[cell_type_col]]), rownames(so@meta.data))
    sample.per.cell <- setNames(as.character(so@meta.data[[sample_id_col]]), rownames(so@meta.data))
    unique_samples <- unique(sample.per.cell)
    sample_conditions <- sapply(unique_samples, function(s) {
      unique(so@meta.data[[condition_col]][so@meta.data[[sample_id_col]] == s])[1]
    }, USE.NAMES = FALSE)
    sample.groups <- setNames(sample_conditions, unique_samples)
    sample.groups <- factor(sample.groups, levels = c(ref.level, target.level))
    if (any(is.na(sample.groups))) stop("sample.groups contains NA after factor conversion. Check condition labels.")

    x <- timed("create_cacoa_object", Cacoa$new(
      so,
      sample.groups = sample.groups,
      cell.groups = cell.groups,
      sample.per.cell = sample.per.cell,
      ref.level = ref.level,
      target.level = target.level,
      graph.name = graph_name,
      data.slot = "data"
    )); cao <- x$value; create_cacoa_seconds <- x$seconds

    hvg_genes <- rownames(so)
    x <- timed("estimate_cluster_free_de", cao$estimateClusterFreeDE(
      genes = hvg_genes,
      min.expr.frac = min_expr_frac,
      adjust.pvalues = adjust_pvalues,
      smooth = smooth_de,
      verbose = TRUE,
      n.cores = n_cores
    )); de_seconds <- x$seconds

    x <- timed("write_outputs", {
      z <- cao$test.results$cluster.free.de$z
      write.csv(z, z_file)
      if (save_cacoa_rds) saveRDS(cao, rds_file)
      if (save_hvg_file) write.csv(data.frame(gene = hvg_genes), hvg_file, row.names = FALSE)
      if (save_metadata) write.csv(so@meta.data, meta_file)
      TRUE
    }); write_seconds <- x$seconds

    total_seconds <- round(proc.time()[["elapsed"]] - t_all, 3)
    data.frame(
      as.data.frame(row),
      nhood_size = nhood_size,
      run_tag = run_tag,
      run_status = "success",
      error = NA_character_,
      run_n_cells = as.integer(ncol(so)),
      run_n_genes = as.integer(nrow(so)),
      n_hvg = as.integer(length(hvg_genes)),
      n_samples = as.integer(length(unique_samples)),
      n_cell_groups = as.integer(length(unique(cell.groups))),
      graph_name = graph_name,
      read_rds_seconds = read_rds_seconds,
      find_neighbors_seconds = find_neighbors_seconds,
      create_cacoa_seconds = create_cacoa_seconds,
      de_seconds = de_seconds,
      preprocessing_seconds = if ("preprocess_seconds" %in% names(row)) as.numeric(row$preprocess_seconds) else NA_real_,
      benchmark_total_seconds = if ("preprocess_seconds" %in% names(row)) as.numeric(row$preprocess_seconds) + de_seconds else NA_real_,
      benchmark_definition = "preprocessing + estimateClusterFreeDE",
      write_seconds = write_seconds,
      total_seconds = total_seconds,
      z_file = z_file,
      rds_file = if (save_cacoa_rds) rds_file else NA_character_,
      hvg_file = if (save_hvg_file) hvg_file else NA_character_,
      metadata_file = if (save_metadata) meta_file else NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(as.data.frame(row), nhood_size = nhood_size, run_tag = run_tag, run_status = "failed", error = conditionMessage(e), z_file = z_file, stringsAsFactors = FALSE)
  })
}

# -----------------------------
# Main
# -----------------------------
parse_cli_args()
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
setup_future()

cat("Base directory:", base_dir, "\n")
cat("Preprocess manifest:", preprocess_manifest, "\n")
cat("Results directory:", results_dir, "\n")
cat("nhood_sizes:", paste(nhood_sizes, collapse = ","), "\n")
cat("outer_workers:", outer_workers, "outer_backend:", outer_backend, "\n")
cat("CACOA DE n_cores:", n_cores, "\n")
cat("future plan:", future_plan_mode, "\n")
cat("Estimated active cores if all nested parallelism is active:", outer_workers * n_cores, "\n")

prep_df <- load_preprocess_manifest()
cat("Preprocessed datasets:", nrow(prep_df), "\n")

tasks <- list()
for (i in seq_len(nrow(prep_df))) {
  for (k in nhood_sizes) {
    tasks[[length(tasks) + 1]] <- list(row = prep_df[i, , drop = FALSE], nhood_size = k)
  }
}
if (test_mode) tasks <- head(tasks, max_runs_for_test)
cat("CACOA Stage 2 runs:", length(tasks), "\n")

if (outer_workers <= 1 || identical(outer_backend, "serial")) {
  rows <- lapply(tasks, run_one)
} else {
  rows <- parallel::mclapply(tasks, run_one, mc.cores = outer_workers, mc.preschedule = FALSE)
}

summary_df <- bind_rows(rows)
summary_file <- file.path(results_dir, "cacoa_parameter_sweeps_run_summary_from_hvg_rds.csv")
failed_file <- file.path(results_dir, "cacoa_parameter_sweeps_failed_runs_from_hvg_rds.csv")
aggregate_file <- file.path(results_dir, "cacoa_parameter_sweeps_aggregate_summary_from_hvg_rds.csv")

write.csv(summary_df, summary_file, row.names = FALSE)
write.csv(summary_df[summary_df$run_status == "failed", , drop = FALSE], failed_file, row.names = FALSE)

if (nrow(summary_df) > 0 && "run_status" %in% colnames(summary_df)) {
  aggregate_df <- summary_df %>%
    filter(run_status %in% c("success", "skipped_existing")) %>%
    group_by(sweep_name, swept_parameter, swept_value, nhood_size) %>%
    summarize(
      n_successful_replicates = n_distinct(replicate),
      n_runs = n(),
      mean_n_cells = mean(run_n_cells, na.rm = TRUE),
      mean_n_genes = mean(run_n_genes, na.rm = TRUE),
      mean_de_seconds = mean(de_seconds, na.rm = TRUE),
      mean_total_seconds = mean(total_seconds, na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(sweep_name, swept_value, nhood_size)
} else {
  aggregate_df <- data.frame()
}
write.csv(aggregate_df, aggregate_file, row.names = FALSE)

cat("\nDone CACOA Stage 2.\n")
cat("Successful/skipped:", sum(summary_df$run_status %in% c("success", "skipped_existing")), "\n")
cat("Failed:", sum(summary_df$run_status == "failed"), "\n")
cat("Summary file:", summary_file, "\n")
cat("Failed file:", failed_file, "\n")
cat("Aggregate file:", aggregate_file, "\n")
