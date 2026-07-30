#!/usr/bin/env Rscript
# ============================================================
# miloDE Stage 1: read TSV/CSV once, normalize, select HVGs,
# subset to HVGs, run PCA + Harmony, and save compact RDS files.
#
# This avoids repeating expensive I/O/HVG/Harmony work for every
# neighbourhood size in the miloDE benchmark.
# ============================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(BiocParallel)
  library(dplyr)
  library(glue)
  library(harmony)
  library(Matrix)
})

has_data_table <- requireNamespace("data.table", quietly = TRUE)
if (has_data_table) suppressPackageStartupMessages(library(data.table))

# -----------------------------
# Defaults
# -----------------------------
base_dir <- "/home/anilprakash/labs/Mei/projects/anil/srda/notebooks/data/scrna_seq/simulation"
export_dir <- file.path(base_dir, "cacoa_milode_exports_runtime_5_cell_counts")
manifest_path <- file.path(export_dir, "cacoa_milode_export_manifest_runtime_5_cell_counts.csv")

preprocess_dir <- file.path(base_dir, "milode_runtime_hvg_harmony_rds")
dir.create(preprocess_dir, recursive = TRUE, showWarnings = FALSE)

input_mode <- "manifest"        # manifest, selected, direct
selected_manifest_path <- NULL
selected_files <- NULL
selected_files_file <- NULL
expression_files <- NULL
metadata_files <- NULL
gene_annotation_files <- NULL
input_h5ad_files <- NULL
selected_dataset_ids <- NULL

n_hvgs <- 2000
n_pcs <- 50
harmony_batch_col <- "sim_batch"
reduced_dim_name <- "HARMONY"
default_harmony_theta <- 2
condition_col <- "condition"
sample_id_col <- "sim_batch"
control_label <- "Control"
disease_label <- "Disease"

preprocess_workers <- 1L
preprocess_backend <- "serial"
overwrite <- FALSE
test_mode <- FALSE
max_datasets_for_test <- 2

# Optional manifest filters
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

  x <- get_arg("--base-dir="); if (!is.null(x)) base_dir <<- x
  x <- get_arg("--export-dir="); if (!is.null(x)) export_dir <<- x
  x <- get_arg("--manifest-path="); if (!is.null(x)) manifest_path <<- x
  x <- get_arg("--preprocess-dir="); if (!is.null(x)) preprocess_dir <<- x
  x <- get_arg("--input-mode="); if (!is.null(x)) input_mode <<- x
  x <- get_arg("--selected-manifest="); if (!is.null(x)) { selected_manifest_path <<- x; input_mode <<- "selected" }
  x <- get_arg("--selected-files="); if (!is.null(x) && nzchar(x)) { selected_files <<- split_arg(x); input_mode <<- "selected" }
  x <- get_arg("--selected-files-file="); if (!is.null(x)) { selected_files_file <<- x; input_mode <<- "selected" }
  x <- get_arg("--expression-files="); if (!is.null(x) && nzchar(x)) { expression_files <<- split_arg(x); input_mode <<- "direct" }
  x <- get_arg("--metadata-files="); if (!is.null(x) && nzchar(x)) { metadata_files <<- split_arg(x); input_mode <<- "direct" }
  x <- get_arg("--gene-annotation-files="); if (!is.null(x) && nzchar(x)) { gene_annotation_files <<- split_arg(x); input_mode <<- "direct" }
  x <- get_arg("--input-h5ad-files="); if (!is.null(x) && nzchar(x)) input_h5ad_files <<- split_arg(x)
  x <- get_arg("--dataset-ids="); if (!is.null(x) && nzchar(x)) selected_dataset_ids <<- split_arg(x)
  x <- get_arg("--n-hvgs="); if (!is.null(x)) n_hvgs <<- as.integer(x)
  x <- get_arg("--n-pcs="); if (!is.null(x)) n_pcs <<- as.integer(x)
  x <- get_arg("--default-harmony-theta="); if (!is.null(x)) default_harmony_theta <<- as.numeric(x)
  x <- get_arg("--workers="); if (!is.null(x)) preprocess_workers <<- as.integer(x)
  x <- get_arg("--backend="); if (!is.null(x)) preprocess_backend <<- x
  x <- get_arg("--overwrite="); if (!is.null(x)) overwrite <<- tolower(x) %in% c("true","t","1","yes","y")
  x <- get_arg("--test-mode="); if (!is.null(x)) test_mode <<- tolower(x) %in% c("true","t","1","yes","y")
  invisible(NULL)
}

# -----------------------------
# Input helpers
# -----------------------------
safe_value_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.0$", "", x)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x
}

make_dataset_id <- function(row) {
  if ("dataset_id" %in% colnames(row) && !is.na(row[["dataset_id"]]) && nzchar(as.character(row[["dataset_id"]]))) {
    return(safe_value_tag(row[["dataset_id"]]))
  }
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

read_selected_files <- function() {
  vals <- character(0)
  if (!is.null(selected_files)) vals <- c(vals, selected_files)
  if (!is.null(selected_files_file)) {
    if (!file.exists(selected_files_file)) stop("selected_files_file does not exist: ", selected_files_file)
    file_vals <- trimws(readLines(selected_files_file, warn = FALSE))
    file_vals <- file_vals[nzchar(file_vals)]
    file_vals <- file_vals[!(tolower(file_vals) %in% c("file", "files", "input_h5ad", "expression_cells_by_genes_tsv"))]
    vals <- c(vals, file_vals)
  }
  unique(trimws(vals[nzchar(vals)]))
}

filter_manifest_by_selected_files <- function(df, selected) {
  if (length(selected) == 0) stop("input_mode='selected' but no selected files were provided.")
  selected_basenames <- basename(selected)
  candidate_cols <- intersect(c("input_h5ad", "expression_cells_by_genes_tsv", "metadata_csv", "gene_annotations_csv"), colnames(df))
  keep <- rep(FALSE, nrow(df))
  for (col in candidate_cols) {
    values <- as.character(df[[col]])
    keep <- keep | values %in% selected | basename(values) %in% selected_basenames
  }
  out <- df[keep, , drop = FALSE]
  if (nrow(out) == 0) stop("Selected files did not match any manifest rows.")
  out
}

resolve_file_path <- function(path, preferred_dir, label) {
  path <- as.character(path)
  if (is.na(path) || path == "") return(path)
  if (file.exists(path)) return(normalizePath(path, mustWork = TRUE))
  candidates <- c(file.path(preferred_dir, basename(path)), file.path(export_dir, basename(path)), file.path(base_dir, basename(path)))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))
  stop("Could not find ", label, " file: ", path, "\nTried as given and by basename under: ", preferred_dir)
}

make_direct_manifest <- function() {
  if (is.null(expression_files) || is.null(metadata_files) || is.null(gene_annotation_files)) {
    stop("Direct mode requires --expression-files, --metadata-files, and --gene-annotation-files.")
  }
  n <- length(expression_files)
  if (length(metadata_files) != n || length(gene_annotation_files) != n) stop("Direct-mode file vectors must have the same length.")
  input_h5ad <- if (is.null(input_h5ad_files)) rep(NA_character_, n) else input_h5ad_files
  if (length(input_h5ad) != n) stop("input_h5ad_files must be absent or have the same length as expression_files.")
  dataset_id <- if (is.null(selected_dataset_ids)) tools::file_path_sans_ext(basename(expression_files)) else selected_dataset_ids
  if (length(dataset_id) != n) stop("dataset_ids must be absent or have the same length as expression_files.")
  data.frame(
    input_h5ad = input_h5ad,
    simulation_type = "selected_files",
    sweep_name = "selected_files",
    swept_parameter = "selected_file",
    swept_value = seq_len(n),
    replicate = seq_len(n),
    responder_percent = NA_integer_,
    dropout_rate = NA_real_,
    harmony_theta = rep(default_harmony_theta, n),
    dataset_id = dataset_id,
    expression_cells_by_genes_tsv = expression_files,
    metadata_csv = metadata_files,
    gene_annotations_csv = gene_annotation_files,
    stringsAsFactors = FALSE
  )
}

load_manifest <- function() {
  active_manifest_path <- if (!is.null(selected_manifest_path)) selected_manifest_path else manifest_path
  if (identical(input_mode, "direct")) {
    df <- make_direct_manifest()
  } else {
    if (!file.exists(active_manifest_path)) stop("Manifest not found: ", active_manifest_path)
    df <- read.csv(active_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  required <- c("expression_cells_by_genes_tsv", "metadata_csv", "gene_annotations_csv", "simulation_type", "sweep_name", "swept_parameter", "swept_value", "replicate", "responder_percent", "dropout_rate", "harmony_theta")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("Manifest is missing required columns: ", paste(missing, collapse = ", "))

  if (!identical(input_mode, "direct") && identical(input_mode, "selected") && is.null(selected_manifest_path)) {
    df <- filter_manifest_by_selected_files(df, read_selected_files())
  }

  if ("input_h5ad" %in% colnames(df)) {
    df$input_h5ad <- vapply(df$input_h5ad, function(x) {
      if (is.na(x) || x == "") return(as.character(x))
      if (file.exists(x)) return(normalizePath(x, mustWork = TRUE))
      hit <- file.path(base_dir, basename(x))
      if (file.exists(hit)) normalizePath(hit, mustWork = TRUE) else as.character(x)
    }, character(1))
  } else {
    df$input_h5ad <- NA_character_
  }
  df$expression_cells_by_genes_tsv <- vapply(df$expression_cells_by_genes_tsv, resolve_file_path, character(1), preferred_dir = file.path(export_dir, "expression_tsv"), label = "expression TSV")
  df$metadata_csv <- vapply(df$metadata_csv, resolve_file_path, character(1), preferred_dir = file.path(export_dir, "metadata_csv"), label = "metadata CSV")
  df$gene_annotations_csv <- vapply(df$gene_annotations_csv, resolve_file_path, character(1), preferred_dir = file.path(export_dir, "gene_annotations_csv"), label = "gene annotation CSV")

  df$replicate <- suppressWarnings(as.integer(df$replicate))
  df$responder_percent <- suppressWarnings(as.integer(df$responder_percent))
  df$dropout_rate <- suppressWarnings(as.numeric(df$dropout_rate))
  df$harmony_theta <- suppressWarnings(as.numeric(df$harmony_theta))
  df$swept_value <- suppressWarnings(as.numeric(df$swept_value))

  if (!is.null(run_only_sweep_names)) df <- df[df$sweep_name %in% run_only_sweep_names, , drop = FALSE]
  if (!is.null(run_only_swept_parameters)) df <- df[df$swept_parameter %in% run_only_swept_parameters, , drop = FALSE]
  if (!is.null(run_only_swept_values)) df <- df[df$swept_value %in% run_only_swept_values, , drop = FALSE]
  if (!is.null(run_only_replicates)) df <- df[df$replicate %in% run_only_replicates, , drop = FALSE]

  df <- df %>% arrange(sweep_name, swept_parameter, swept_value, replicate)
  if (test_mode) df <- head(df, max_datasets_for_test)
  if (nrow(df) == 0) stop("Zero datasets after filtering.")
  df$dataset_id <- vapply(seq_len(nrow(df)), function(i) make_dataset_id(df[i, , drop = FALSE]), character(1))
  df
}

# -----------------------------
# Reading and SCE helpers
# -----------------------------
safe_read_expression <- function(expr_path) {
  if (has_data_table) {
    expr_df <- as.data.frame(data.table::fread(expr_path, sep = "\t", data.table = FALSE, check.names = FALSE))
    cell_ids <- expr_df[[1]]
    expr_df[[1]] <- NULL
    count_matrix <- as.matrix(expr_df)
    rownames(count_matrix) <- cell_ids
  } else {
    expr_df <- read.delim(expr_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    count_matrix <- as.matrix(expr_df)
  }
  storage.mode(count_matrix) <- "numeric"
  t(count_matrix) # genes x cells
}

safe_read_metadata <- function(obs_path) {
  if (has_data_table) {
    obs <- as.data.frame(data.table::fread(obs_path, data.table = FALSE, check.names = FALSE))
    rownames(obs) <- obs[[1]]; obs[[1]] <- NULL; obs
  } else {
    read.csv(obs_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

safe_read_var <- function(var_path) {
  if (has_data_table) {
    var <- as.data.frame(data.table::fread(var_path, data.table = FALSE, check.names = FALSE))
    rownames(var) <- var[[1]]; var[[1]] <- NULL; var
  } else {
    read.csv(var_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

validate_metadata <- function(obs_df) {
  needed <- c(condition_col, sample_id_col, "Cell_Type")
  missing <- setdiff(needed, colnames(obs_df))
  if (length(missing) > 0) stop("Metadata is missing required columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

align_inputs <- function(count_matrix, obs_df, var_df) {
  common_cells <- intersect(colnames(count_matrix), rownames(obs_df))
  common_genes <- intersect(rownames(count_matrix), rownames(var_df))
  if (length(common_cells) == 0) stop("No overlapping cell barcodes between expression and metadata.")
  if (length(common_genes) == 0) stop("No overlapping genes between expression and gene annotation.")
  count_matrix <- count_matrix[common_genes, common_cells, drop = FALSE]
  obs_df <- obs_df[common_cells, , drop = FALSE]
  var_df <- var_df[common_genes, , drop = FALSE]
  list(count_matrix = count_matrix, obs_df = obs_df, var_df = var_df)
}

make_sce <- function(count_matrix, obs_df, var_df) {
  validate_metadata(obs_df)
  obs_df[[condition_col]] <- factor(obs_df[[condition_col]], levels = c(control_label, disease_label))
  SingleCellExperiment(assays = list(counts = count_matrix), colData = obs_df, rowData = var_df)
}

timed <- function(label, expr) {
  cat("[TIMER START]", label, "\n")
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  seconds <- round(proc.time()[["elapsed"]] - t0, 3)
  cat("[TIMER DONE]", label, ":", seconds, "seconds\n")
  list(value = value, seconds = seconds)
}

preprocess_one <- function(row) {
  row <- as.data.frame(row, stringsAsFactors = FALSE)
  dataset_id <- as.character(row$dataset_id)
  out_rds <- file.path(preprocess_dir, paste0("milode_hvg_harmony_sce_", dataset_id, ".rds"))
  hvg_file <- file.path(preprocess_dir, paste0("milode_hvg_genes_", dataset_id, ".csv"))

  if (file.exists(out_rds) && file.exists(hvg_file) && !overwrite) {
    return(data.frame(row, hvg_sce_rds = out_rds, hvg_file = hvg_file, status = "skipped_existing", error = NA_character_, preprocess_seconds = 0, stringsAsFactors = FALSE))
  }

  t_all <- proc.time()[["elapsed"]]
  tryCatch({
    cat("\n------------------------------------------------------------\n")
    cat("Preprocessing miloDE dataset:", dataset_id, "\n")
    cat("------------------------------------------------------------\n")

    x <- timed("read_expression_tsv", safe_read_expression(row$expression_cells_by_genes_tsv)); count_matrix <- x$value; read_expression_seconds <- x$seconds
    x <- timed("read_metadata_csv", safe_read_metadata(row$metadata_csv)); obs_df <- x$value; read_metadata_seconds <- x$seconds
    x <- timed("read_gene_annotation_csv", safe_read_var(row$gene_annotations_csv)); var_df <- x$value; read_var_seconds <- x$seconds
    x <- timed("align_expression_metadata_genes", align_inputs(count_matrix, obs_df, var_df)); aligned <- x$value; align_seconds <- x$seconds
    rm(count_matrix, obs_df, var_df); gc()

    count_matrix <- aligned$count_matrix; obs_df <- aligned$obs_df; var_df <- aligned$var_df
    n_genes_original <- nrow(count_matrix); n_cells <- ncol(count_matrix)
    cat("Counts matrix:", n_genes_original, "genes x", n_cells, "cells\n")

    x <- timed("make_sce", make_sce(count_matrix, obs_df, var_df)); sce <- x$value; make_sce_seconds <- x$seconds
    rm(count_matrix, obs_df, var_df, aligned); gc()

    x <- timed("log_normalization", logNormCounts(sce)); sce <- x$value; lognorm_seconds <- x$seconds
    x <- timed("hvg_selection_modelGeneVar", modelGeneVar(sce)); dec.sce <- x$value; hvg_seconds <- x$seconds
    hvg.genes <- getTopHVGs(dec.sce, n = min(n_hvgs, nrow(sce)))
    if (length(hvg.genes) < 3) stop("Too few HVGs selected.")

    rowData(sce)$is_hvg <- rownames(sce) %in% hvg.genes
    x <- timed("subset_sce_to_hvgs", sce[hvg.genes, , drop = FALSE]); sce <- x$value; subset_seconds <- x$seconds
    cat("Subsetted SCE to HVGs:", nrow(sce), "genes x", ncol(sce), "cells\n")

    x <- timed("pca_on_hvgs", runPCA(sce, ncomponents = min(n_pcs, length(hvg.genes) - 1), subset_row = rownames(sce))); sce <- x$value; pca_seconds <- x$seconds

    harmony_theta_value <- as.numeric(row$harmony_theta)
    if (is.na(harmony_theta_value)) harmony_theta_value <- default_harmony_theta
    cat("Running Harmony with theta =", harmony_theta_value, "\n")
    x <- timed("harmony", {
      emb <- RunHarmony(
        data_mat = reducedDim(sce, "PCA"),
        meta_data = as.data.frame(colData(sce)),
        vars_use = harmony_batch_col,
        theta = harmony_theta_value,
        verbose = FALSE
      )
      reducedDim(sce, reduced_dim_name) <- emb
      sce
    }); sce <- x$value; harmony_seconds <- x$seconds

    metadata(sce)$milode_preprocess <- list(
      dataset_id = dataset_id,
      n_genes_original = n_genes_original,
      n_hvgs = length(hvg.genes),
      n_cells = n_cells,
      harmony_theta = harmony_theta_value,
      reduced_dim_name = reduced_dim_name
    )

    dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
    x <- timed("save_hvg_harmony_rds", saveRDS(sce, out_rds, compress = FALSE)); save_seconds <- x$seconds
    write.csv(data.frame(gene = hvg.genes), hvg_file, row.names = FALSE)

    total_seconds <- round(proc.time()[["elapsed"]] - t_all, 3)
    data.frame(
      row,
      hvg_sce_rds = out_rds,
      hvg_file = hvg_file,
      status = "success",
      error = NA_character_,
      n_genes_original = as.integer(n_genes_original),
      n_hvgs = as.integer(length(hvg.genes)),
      n_cells = as.integer(n_cells),
      read_expression_seconds = read_expression_seconds,
      read_metadata_seconds = read_metadata_seconds,
      read_var_seconds = read_var_seconds,
      align_seconds = align_seconds,
      make_sce_seconds = make_sce_seconds,
      lognorm_seconds = lognorm_seconds,
      hvg_seconds = hvg_seconds,
      subset_seconds = subset_seconds,
      pca_seconds = pca_seconds,
      harmony_seconds = harmony_seconds,
      save_seconds = save_seconds,
      preprocess_seconds = total_seconds,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(row, hvg_sce_rds = out_rds, hvg_file = hvg_file, status = "failed", error = conditionMessage(e), preprocess_seconds = round(proc.time()[["elapsed"]] - t_all, 3), stringsAsFactors = FALSE)
  })
}

make_bpparam <- function(workers, backend) {
  workers <- max(1, as.integer(workers))
  if (workers == 1 || identical(backend, "serial")) return(SerialParam())
  if (identical(backend, "snow")) return(SnowParam(workers = workers, type = "SOCK"))
  MulticoreParam(workers = workers)
}

# -----------------------------
# Main
# -----------------------------
parse_cli_args()
dir.create(preprocess_dir, recursive = TRUE, showWarnings = FALSE)
manifest_df <- load_manifest()

cat("Input mode:", input_mode, "\n")
cat("Datasets to preprocess:", nrow(manifest_df), "\n")
cat("Preprocess dir:", preprocess_dir, "\n")
cat("Workers:", preprocess_workers, "backend:", preprocess_backend, "\n")

rows <- split(manifest_df, seq_len(nrow(manifest_df)))
param <- make_bpparam(preprocess_workers, preprocess_backend)
summary_list <- bplapply(rows, preprocess_one, BPPARAM = param)
summary_df <- bind_rows(summary_list)

summary_file <- file.path(preprocess_dir, "milode_runtime_preprocess_manifest.csv")
failed_file <- file.path(preprocess_dir, "milode_hvg_harmony_preprocess_failed.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
write.csv(summary_df[summary_df$status == "failed", , drop = FALSE], failed_file, row.names = FALSE)

cat("\nDone preprocessing.\n")
cat("Successful:", sum(summary_df$status %in% c("success", "skipped_existing")), "\n")
cat("Failed:", sum(summary_df$status == "failed"), "\n")
cat("Preprocessed manifest:", summary_file, "\n")
cat("Failed file:", failed_file, "\n")
