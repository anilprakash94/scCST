#!/usr/bin/env Rscript
# ============================================================
# miloDE Stage 2: load Stage-1 HVG+Harmony SCE RDS files and
# run neighbourhood assignment, AUC filtering, and miloDE DE.
#
# Parallelism:
#   - outer workers: multiple dataset/nhood-size runs in parallel
#   - inner cores: BiocParallel workers inside AUC and DE
#
# Keep outer_workers * inner_cores within your available CPU/RAM.
# ============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(miloDE)
  library(miloR)
  library(dplyr)
  library(BiocParallel)
  library(glue)
  library(Matrix)
})

# -----------------------------
# Defaults
# -----------------------------
base_dir <- "/path/to/data/scrna_seq/simulation"
preprocess_dir <- file.path(base_dir, "milode_runtime_hvg_harmony_rds")
preprocess_manifest <- file.path(preprocess_dir, "milode_runtime_preprocess_manifest.csv")
results_dir <- file.path(base_dir, "milode_runtime_5_cell_counts_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

nhood_sizes <- c(20)
order_value <- 2
filtering_value <- TRUE
reduced_dim_name <- "HARMONY"
sample_id_col <- "sim_batch"
condition_col <- "condition"
control_label <- "Control"
disease_label <- "Disease"

outer_workers <- 1L
inner_cores <- 1L
outer_backend <- "serial"
inner_backend <- "serial"

de_verbose <- FALSE
save_cell_z <- TRUE
save_nhood_z <- FALSE
save_auc <- FALSE
save_de_stat_rds <- FALSE
save_sce_milo_rds <- FALSE
overwrite <- FALSE
test_mode <- FALSE
max_runs_for_test <- 2

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
  x <- get_arg("--inner-cores="); if (!is.null(x)) inner_cores <<- as.integer(x)
  x <- get_arg("--outer-backend="); if (!is.null(x)) outer_backend <<- x
  x <- get_arg("--inner-backend="); if (!is.null(x)) inner_backend <<- x
  x <- get_arg("--de-verbose="); if (!is.null(x)) de_verbose <<- parse_bool(x)
  x <- get_arg("--save-cell-z="); if (!is.null(x)) save_cell_z <<- parse_bool(x)
  x <- get_arg("--save-nhood-z="); if (!is.null(x)) save_nhood_z <<- parse_bool(x)
  x <- get_arg("--save-auc="); if (!is.null(x)) save_auc <<- parse_bool(x)
  x <- get_arg("--save-de-stat-rds="); if (!is.null(x)) save_de_stat_rds <<- parse_bool(x)
  x <- get_arg("--save-sce-milo-rds="); if (!is.null(x)) save_sce_milo_rds <<- parse_bool(x)
  x <- get_arg("--overwrite="); if (!is.null(x)) overwrite <<- parse_bool(x)
  x <- get_arg("--test-mode="); if (!is.null(x)) test_mode <<- parse_bool(x)
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

make_bpparam <- function(workers, backend) {
  workers <- max(1, as.integer(workers))
  if (workers == 1 || identical(backend, "serial")) return(SerialParam())
  if (identical(backend, "snow")) return(SnowParam(workers = workers, type = "SOCK"))
  MulticoreParam(workers = workers)
}

timed <- function(label, expr) {
  cat("[TIMER START]", label, "\n")
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  seconds <- round(proc.time()[["elapsed"]] - t0, 3)
  cat("[TIMER DONE]", label, ":", seconds, "seconds\n")
  list(value = value, seconds = seconds)
}

save_matrix_csv <- function(mat, path) {
  df <- as.data.frame(as.matrix(mat))
  write.csv(df, path, quote = TRUE)
}

project_z_to_cells <- function(sce_milo, de_stat, tested_nhoods) {
  p_mat <- assay(de_stat, "pval_corrected_across_nhoods")
  p_mat[is.na(p_mat)] <- 1
  logfc_mat <- assay(de_stat, "logFC")
  logfc_mat[is.na(logfc_mat)] <- 0

  z_mat <- sign(logfc_mat) * qnorm(p_mat / 2, lower.tail = FALSE)
  z_mat[is.na(z_mat)] <- 0
  z_mat[is.infinite(z_mat)] <- sign(z_mat[is.infinite(z_mat)]) * 38

  nhood_membership <- nhoods(sce_milo)[, tested_nhoods, drop = FALSE]
  z_mat_t <- t(z_mat) # neighborhoods x genes
  common_nhoods <- intersect(colnames(nhood_membership), rownames(z_mat_t))

  if (length(common_nhoods) == 0) {
    if (ncol(nhood_membership) == nrow(z_mat_t)) {
      common_nhoods <- colnames(nhood_membership)
      rownames(z_mat_t) <- common_nhoods
    } else {
      stop("Could not align neighbourhood membership columns with DE z-score neighbourhoods.")
    }
  }

  nhood_membership <- nhood_membership[, common_nhoods, drop = FALSE]
  z_mat_t <- z_mat_t[common_nhoods, , drop = FALSE]

  cell_z <- nhood_membership %*% z_mat_t
  n_nhoods_per_cell <- rowSums(nhood_membership)
  cell_z <- sweep(cell_z, 1, n_nhoods_per_cell, "/")
  cell_z[n_nhoods_per_cell == 0, ] <- 0

  list(cell_z = cell_z, z_mat = z_mat)
}

run_one <- function(task) {
  row <- task$row
  nhood_size <- as.integer(task$nhood_size)
  run_tag <- make_run_tag(row, nhood_size)
  hvg_sce_rds <- as.character(row$hvg_sce_rds)

  cell_z_file <- file.path(results_dir, glue("milode_cell_z_scores_{run_tag}.csv"))
  nhood_z_file <- file.path(results_dir, glue("milode_nhood_z_scores_{run_tag}.csv"))
  auc_file <- file.path(results_dir, glue("milode_neighbourhood_auc_{run_tag}.csv"))
  sce_milo_file <- file.path(results_dir, glue("milode_sce_milo_{run_tag}.rds"))
  de_stat_file <- file.path(results_dir, glue("milode_de_stat_{run_tag}.rds"))

  if (save_cell_z && file.exists(cell_z_file) && !overwrite) {
    return(data.frame(row, nhood_size = nhood_size, run_tag = run_tag, status = "skipped_existing", error = NA_character_, cell_z_file = cell_z_file, stringsAsFactors = FALSE))
  }

  t_all <- proc.time()[["elapsed"]]
  tryCatch({
    cat("\n------------------------------------------------------------\n")
    cat("Running miloDE from HVG RDS:", run_tag, "\n")
    cat("------------------------------------------------------------\n")

    if (!file.exists(hvg_sce_rds)) stop("Missing hvg_sce_rds: ", hvg_sce_rds)
    x <- timed("read_hvg_harmony_sce_rds", readRDS(hvg_sce_rds)); sce <- x$value; read_rds_seconds <- x$seconds
    cat("HVG SCE:", nrow(sce), "genes x", ncol(sce), "cells\n")

    if (!(reduced_dim_name %in% reducedDimNames(sce))) stop("ReducedDim not found in SCE: ", reduced_dim_name)
    colData(sce)[[condition_col]] <- factor(colData(sce)[[condition_col]], levels = c(control_label, disease_label))

    inner_param <- make_bpparam(inner_cores, inner_backend)

    set.seed(42)
    x <- timed("assign_neighbourhoods", assign_neighbourhoods(
      sce,
      k = nhood_size,
      order = order_value,
      filtering = filtering_value,
      reducedDim_name = reduced_dim_name,
      verbose = TRUE
    )); sce_milo <- x$value; assign_seconds <- x$seconds

    x <- timed("auc_neighbourhood_filtering", suppressWarnings(calc_AUC_per_neighbourhood(
      sce_milo,
      sample_id = sample_id_col,
      condition_id = condition_col,
      min_n_cells_per_sample = 1,
      BPPARAM = inner_param
    ))); stat_auc <- x$value; auc_seconds <- x$seconds

    tested_nhoods <- stat_auc$Nhood[!is.na(stat_auc$auc)]
    if (length(tested_nhoods) == 0) stop("No neighborhoods passed AUC filtering.")

    x <- timed("milode_de_test_neighbourhoods", de_test_neighbourhoods(
      sce_milo,
      sample_id = sample_id_col,
      design = as.formula(glue("~{condition_col}")),
      covariates = c(condition_col),
      subset_nhoods = tested_nhoods,
      output_type = "SCE",
      plot_summary_stat = FALSE,
      BPPARAM = inner_param,
      verbose = de_verbose
    )); de_stat <- x$value; de_seconds <- x$seconds
    cat("miloDE DE testing finished in", de_seconds, "seconds.\n")

    x <- timed("project_nhood_z_to_cell_z", project_z_to_cells(sce_milo, de_stat, tested_nhoods)); projected <- x$value; project_seconds <- x$seconds
    cell_z <- projected$cell_z; z_mat <- projected$z_mat

    x <- timed("write_outputs", {
      if (save_cell_z) {
        cell_z_df <- as.data.frame(as.matrix(cell_z))
        cell_z_df$cell_barcode <- rownames(cell_z_df)
        cell_z_df <- cell_z_df[, c("cell_barcode", setdiff(colnames(cell_z_df), "cell_barcode")), drop = FALSE]
        write.csv(cell_z_df, cell_z_file, row.names = FALSE)
      } else {
        cell_z_file <- NA_character_
      }
      if (save_nhood_z) save_matrix_csv(z_mat, nhood_z_file) else nhood_z_file <- NA_character_
      if (save_auc) write.csv(as.data.frame(stat_auc), auc_file, row.names = FALSE) else auc_file <- NA_character_
      if (save_sce_milo_rds) saveRDS(sce_milo, sce_milo_file, compress = FALSE) else sce_milo_file <- NA_character_
      if (save_de_stat_rds) saveRDS(de_stat, de_stat_file, compress = FALSE) else de_stat_file <- NA_character_
      TRUE
    }); write_seconds <- x$seconds

    total_seconds <- round(proc.time()[["elapsed"]] - t_all, 3)
    data.frame(
      row,
      nhood_size = nhood_size,
      run_tag = run_tag,
      status = "success",
      error = NA_character_,
      n_cells = as.integer(ncol(sce_milo)),
      n_hvgs = as.integer(nrow(sce_milo)),
      n_neighbourhoods_total = as.integer(ncol(nhoods(sce_milo))),
      n_neighbourhoods_tested = as.integer(length(tested_nhoods)),
      cell_z_file = cell_z_file,
      nhood_z_file = nhood_z_file,
      auc_file = auc_file,
      sce_milo_rds = sce_milo_file,
      de_stat_rds = de_stat_file,
      read_rds_seconds = read_rds_seconds,
      assign_seconds = assign_seconds,
      auc_seconds = auc_seconds,
      de_seconds = de_seconds,
      preprocessing_seconds = if ("preprocess_seconds" %in% names(row)) as.numeric(row$preprocess_seconds) else NA_real_,
      benchmark_total_seconds = if ("preprocess_seconds" %in% names(row)) as.numeric(row$preprocess_seconds) + de_seconds else NA_real_,
      benchmark_definition = "preprocessing + de_test_neighbourhoods",
      project_seconds = project_seconds,
      write_seconds = write_seconds,
      total_seconds = total_seconds,
      outer_workers = outer_workers,
      inner_cores = inner_cores,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(row, nhood_size = nhood_size, run_tag = run_tag, status = "failed", error = conditionMessage(e), cell_z_file = cell_z_file, total_seconds = round(proc.time()[["elapsed"]] - t_all, 3), stringsAsFactors = FALSE)
  })
}

# -----------------------------
# Main
# -----------------------------
parse_cli_args()
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(preprocess_manifest)) stop("Preprocess manifest not found: ", preprocess_manifest, "\nRun milode_stage1_preprocess_hvg_harmony_parallel.R first.")

prep_df <- read.csv(preprocess_manifest, stringsAsFactors = FALSE, check.names = FALSE)

# Stage-1 output may contain two status columns:
#   - status: original conversion/export status, e.g. "completed"
#   - status.1: Stage-1 preprocessing status, e.g. "success"
# This happens when the conversion manifest already had a column named status and
# Stage 1 appended its own status column. Prefer the Stage-1 status column when present.
find_preprocess_status_col <- function(df) {
  candidate_cols <- grep("^status(\\.[0-9]+)?$", colnames(df), value = TRUE)
  if (length(candidate_cols) == 0) {
    return(NA_character_)
  }

  # Prefer the column containing Stage-1 values.
  for (col in candidate_cols) {
    vals <- unique(as.character(df[[col]]))
    if (any(vals %in% c("success", "skipped_existing", "failed"), na.rm = TRUE)) {
      return(col)
    }
  }

  # Fallback: use the last duplicate status-like column.
  tail(candidate_cols, 1)
}

preprocess_status_col <- find_preprocess_status_col(prep_df)
if (is.na(preprocess_status_col)) {
  stop("Could not find a preprocessing status column in: ", preprocess_manifest)
}

cat("Using preprocessing status column:", preprocess_status_col, "\n")
prep_df$preprocess_status <- as.character(prep_df[[preprocess_status_col]])

prep_df <- prep_df[prep_df$preprocess_status %in% c("success", "skipped_existing"), , drop = FALSE]

# Extra safety: keep only rows with an existing Stage-1 RDS file.
if (!("hvg_sce_rds" %in% colnames(prep_df))) {
  stop("Preprocess manifest is missing hvg_sce_rds column: ", preprocess_manifest)
}
prep_df <- prep_df[!is.na(prep_df$hvg_sce_rds) & file.exists(prep_df$hvg_sce_rds), , drop = FALSE]

if (nrow(prep_df) == 0) {
  stop(
    "No successful preprocessed datasets with existing hvg_sce_rds files found in: ",
    preprocess_manifest,
    "\nDetected preprocessing status column: ", preprocess_status_col,
    "\nCheck that hvg_sce_rds paths exist and that Stage 1 completed successfully."
  )
}

if (!is.null(run_only_sweep_names)) prep_df <- prep_df[prep_df$sweep_name %in% run_only_sweep_names, , drop = FALSE]
if (!is.null(run_only_swept_parameters)) prep_df <- prep_df[prep_df$swept_parameter %in% run_only_swept_parameters, , drop = FALSE]
if (!is.null(run_only_swept_values)) prep_df <- prep_df[prep_df$swept_value %in% run_only_swept_values, , drop = FALSE]
if (!is.null(run_only_replicates)) prep_df <- prep_df[prep_df$replicate %in% run_only_replicates, , drop = FALSE]

run_grid <- expand.grid(row_i = seq_len(nrow(prep_df)), nhood_size = nhood_sizes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
if (test_mode) run_grid <- head(run_grid, max_runs_for_test)

tasks <- lapply(seq_len(nrow(run_grid)), function(i) {
  list(row = prep_df[run_grid$row_i[i], , drop = FALSE], nhood_size = run_grid$nhood_size[i])
})

cat("Preprocess manifest:", preprocess_manifest, "\n")
cat("Results dir:", results_dir, "\n")
cat("Runs:", length(tasks), "\n")
cat("Outer workers:", outer_workers, outer_backend, "\n")
cat("Inner cores per run:", inner_cores, inner_backend, "\n")
cat("Approx max active cores:", outer_workers * inner_cores, "\n")

outer_param <- make_bpparam(outer_workers, outer_backend)
summary_list <- bplapply(tasks, run_one, BPPARAM = outer_param)
summary_df <- bind_rows(summary_list)

summary_file <- file.path(results_dir, "milode_hvg_rds_run_summary.csv")
failed_file <- file.path(results_dir, "milode_hvg_rds_failed_runs.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
write.csv(summary_df[summary_df$status == "failed", , drop = FALSE], failed_file, row.names = FALSE)

aggregate_file <- file.path(results_dir, "milode_hvg_rds_aggregate_summary.csv")
if (nrow(summary_df) > 0 && any(summary_df$status == "success")) {
  aggregate_df <- summary_df %>%
    filter(status == "success") %>%
    group_by(sweep_name, swept_parameter, swept_value, nhood_size) %>%
    summarize(
      n_runs = n(),
      mean_n_cells = mean(n_cells),
      mean_n_hvgs = mean(n_hvgs),
      mean_n_neighbourhoods_total = mean(n_neighbourhoods_total),
      mean_n_neighbourhoods_tested = mean(n_neighbourhoods_tested),
      mean_assign_seconds = mean(assign_seconds),
      mean_auc_seconds = mean(auc_seconds),
      mean_de_seconds = mean(de_seconds),
      mean_total_seconds = mean(total_seconds),
      .groups = "drop"
    ) %>% arrange(sweep_name, swept_value, nhood_size)
} else {
  aggregate_df <- data.frame()
}
write.csv(aggregate_df, aggregate_file, row.names = FALSE)

cat("\nDone Stage 2.\n")
cat("Successful runs:", sum(summary_df$status == "success"), "\n")
cat("Skipped existing:", sum(summary_df$status == "skipped_existing"), "\n")
cat("Failed runs:", sum(summary_df$status == "failed"), "\n")
cat("Summary:", summary_file, "\n")
cat("Failed:", failed_file, "\n")
cat("Aggregate:", aggregate_file, "\n")
