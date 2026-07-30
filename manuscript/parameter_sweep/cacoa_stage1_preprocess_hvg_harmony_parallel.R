#!/usr/bin/env Rscript
# ============================================================
# CACOA Stage 1: read TSV/CSV once, create Seurat object,
# normalize, select HVGs, physically subset to HVGs, run PCA +
# Harmony, and save compact Seurat RDS files.
#
# This avoids repeating expensive TSV I/O/HVG/Harmony work for
# every neighbourhood size in the CACOA benchmark.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(glue)
  library(harmony)
  library(parallel)
})

has_data_table <- requireNamespace("data.table", quietly = TRUE)
if (has_data_table) suppressPackageStartupMessages(library(data.table))

# -----------------------------
# Defaults
# -----------------------------
base_dir <- "/path/to/data/scrna_seq/simulation"
export_dir <- file.path(base_dir, "cacoa_milode_exports_parameter_sweeps")
manifest_file <- file.path(export_dir, "cacoa_milode_export_manifest_parameter_sweeps.csv")

preprocess_dir <- file.path(base_dir, "cacoa_hvg_harmony_rds")
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

n_hvg <- 2000
n_pcs <- 50
harmony_dims <- 1:30
default_harmony_theta <- 2
harmony_batch_col <- "sim_batch"
condition_col <- "condition"
sample_id_col <- "sim_batch"
cell_type_col <- "Cell_Type"
control_label <- "Control"
disease_label <- "Disease"

preprocess_workers <- as.integer(Sys.getenv("CACOA_PREPROCESS_WORKERS", "2"))
preprocess_backend <- Sys.getenv("CACOA_PREPROCESS_BACKEND", "multicore") # multicore or serial
overwrite <- FALSE
test_mode <- FALSE
max_datasets_for_test <- 2

# Optional manifest filters
run_only_sweeps <- NULL
run_only_parameters <- NULL
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
  x <- get_arg("--export-dir="); if (!is.null(x)) export_dir <<- x
  x <- get_arg("--manifest-path="); if (!is.null(x)) manifest_file <<- x
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
  x <- get_arg("--n-hvg="); if (!is.null(x)) n_hvg <<- as.integer(x)
  x <- get_arg("--n-pcs="); if (!is.null(x)) n_pcs <<- as.integer(x)
  x <- get_arg("--default-harmony-theta="); if (!is.null(x)) default_harmony_theta <<- as.numeric(x)
  x <- get_arg("--workers="); if (!is.null(x)) preprocess_workers <<- as.integer(x)
  x <- get_arg("--backend="); if (!is.null(x)) preprocess_backend <<- x
  x <- get_arg("--overwrite="); if (!is.null(x)) overwrite <<- parse_bool(x)
  x <- get_arg("--test-mode="); if (!is.null(x)) test_mode <<- parse_bool(x)
  invisible(NULL)
}

# -----------------------------
# Helpers
# -----------------------------
safe_value_tag <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.0$", "", x)
  x <- gsub("\\.", "p", x)
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x
}

make_dataset_id <- function(row) {
  if ("dataset_id" %in% names(row) && !is.na(row[["dataset_id"]]) && nzchar(as.character(row[["dataset_id"]]))) {
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
  active_manifest <- if (!is.null(selected_manifest_path)) selected_manifest_path else manifest_file
  if (identical(input_mode, "direct")) {
    df <- make_direct_manifest()
  } else {
    if (!file.exists(active_manifest)) stop("Manifest not found: ", active_manifest)
    df <- read.csv(active_manifest, stringsAsFactors = FALSE, check.names = FALSE)
  }

  required <- c("input_h5ad", "simulation_type", "sweep_name", "swept_parameter", "swept_value", "replicate", "responder_percent", "dropout_rate", "harmony_theta", "expression_cells_by_genes_tsv", "metadata_csv", "gene_annotations_csv")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("Manifest is missing required columns: ", paste(missing, collapse = ", "))

  if (!identical(input_mode, "direct") && identical(input_mode, "selected") && is.null(selected_manifest_path)) {
    df <- filter_manifest_by_selected_files(df, read_selected_files())
  }

  df$input_h5ad <- vapply(df$input_h5ad, function(x) {
    if (is.na(x) || x == "") return(as.character(x))
    if (file.exists(x)) return(normalizePath(x, mustWork = TRUE))
    hit <- file.path(base_dir, basename(x))
    if (file.exists(hit)) normalizePath(hit, mustWork = TRUE) else as.character(x)
  }, character(1))
  df$expression_cells_by_genes_tsv <- vapply(df$expression_cells_by_genes_tsv, resolve_file_path, character(1), preferred_dir = file.path(export_dir, "expression_tsv"), label = "expression TSV")
  df$metadata_csv <- vapply(df$metadata_csv, resolve_file_path, character(1), preferred_dir = file.path(export_dir, "metadata_csv"), label = "metadata CSV")
  df$gene_annotations_csv <- vapply(df$gene_annotations_csv, resolve_file_path, character(1), preferred_dir = file.path(export_dir, "gene_annotations_csv"), label = "gene annotation CSV")

  if (!is.null(run_only_sweeps)) df <- df[df$sweep_name %in% run_only_sweeps, , drop = FALSE]
  if (!is.null(run_only_parameters)) df <- df[df$swept_parameter %in% run_only_parameters, , drop = FALSE]
  if (!is.null(run_only_replicates)) df <- df[df$replicate %in% run_only_replicates, , drop = FALSE]

  df$replicate <- as.integer(df$replicate)
  df$responder_percent <- suppressWarnings(as.integer(df$responder_percent))
  df$dropout_rate <- suppressWarnings(as.numeric(df$dropout_rate))
  df$harmony_theta <- suppressWarnings(as.numeric(df$harmony_theta))
  df$swept_value <- suppressWarnings(as.numeric(df$swept_value))

  df <- df %>% arrange(sweep_name, swept_parameter, swept_value, replicate)
  if (test_mode) df <- head(df, max_datasets_for_test)
  if (nrow(df) == 0) stop("No datasets to preprocess after filtering.")
  df
}

read_table_with_cell_barcode <- function(path, sep = "\t") {
  if (!file.exists(path)) stop("File not found: ", path)
  if (has_data_table) {
    df <- as.data.frame(data.table::fread(path, sep = sep, data.table = FALSE, check.names = FALSE))
  } else {
    df <- read.table(path, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!("cell_barcode" %in% colnames(df))) stop("Expected a 'cell_barcode' column in: ", path)
  rownames(df) <- df$cell_barcode
  df$cell_barcode <- NULL
  df
}

read_csv_with_rowname <- function(path, rowname_col) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (has_data_table) {
    df <- as.data.frame(data.table::fread(path, data.table = FALSE, check.names = FALSE))
  } else {
    df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (!(rowname_col %in% colnames(df))) stop("Expected column '", rowname_col, "' in: ", path)
  rownames(df) <- df[[rowname_col]]
  df[[rowname_col]] <- NULL
  df
}

validate_cacoa_metadata <- function(obs) {
  needed <- c(condition_col, sample_id_col, cell_type_col)
  missing <- setdiff(needed, colnames(obs))
  if (length(missing) > 0) stop("Metadata is missing required columns: ", paste(missing, collapse = ", "))
  if (length(unique(obs[[condition_col]])) != 2) stop("Expected exactly two conditions in metadata.")
  batch_condition_counts <- obs %>% group_by(.data[[sample_id_col]]) %>% summarize(n_conditions = n_distinct(.data[[condition_col]]), .groups = "drop")
  bad_batches <- batch_condition_counts[[sample_id_col]][batch_condition_counts$n_conditions != 1]
  if (length(bad_batches) > 0) stop("Each ", sample_id_col, " must belong to exactly one condition. Bad batches: ", paste(bad_batches, collapse = ", "))
  invisible(TRUE)
}

timed <- function(label, expr) {
  cat("[TIMER START]", label, "\n")
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  seconds <- round(proc.time()[["elapsed"]] - t0, 3)
  cat("[TIMER DONE]", label, ":", seconds, "seconds\n")
  list(value = value, seconds = seconds)
}

process_one <- function(row_df) {
  row <- as.list(row_df[1, , drop = FALSE])
  dataset_id <- make_dataset_id(row)
  seurat_rds <- file.path(preprocess_dir, glue("cacoa_hvg_harmony_seurat_{dataset_id}.rds"))
  hvg_file <- file.path(preprocess_dir, glue("cacoa_hvg_genes_{dataset_id}.csv"))

  if (file.exists(seurat_rds) && file.exists(hvg_file) && !overwrite) {
    return(data.frame(row_df, dataset_id = dataset_id, hvg_seurat_rds = seurat_rds, hvg_file = hvg_file, preprocess_status = "skipped_existing", error = NA_character_, stringsAsFactors = FALSE))
  }

  t_all <- proc.time()[["elapsed"]]
  tryCatch({
    cat("\n------------------------------------------------------------\n")
    cat("CACOA Stage 1 preprocessing:", dataset_id, "\n")
    cat("------------------------------------------------------------\n")

    x <- timed("read_expression_tsv", read_table_with_cell_barcode(as.character(row$expression_cells_by_genes_tsv), sep = "\t")); expr <- x$value; read_expression_seconds <- x$seconds
    x <- timed("read_metadata_csv", read_csv_with_rowname(as.character(row$metadata_csv), rowname_col = "cell_barcode")); obs <- x$value; read_metadata_seconds <- x$seconds
    x <- timed("read_gene_annotation_csv", read_csv_with_rowname(as.character(row$gene_annotations_csv), rowname_col = "gene")); var <- x$value; read_var_seconds <- x$seconds

    x <- timed("align_expression_metadata_genes", {
      non_numeric <- which(!vapply(expr, is.numeric, logical(1)))
      if (length(non_numeric) > 0) stop("Expression contains non-numeric columns: ", paste(colnames(expr)[non_numeric], collapse = ", "))
      common_cells <- intersect(rownames(expr), rownames(obs))
      if (length(common_cells) == 0) stop("No overlapping cell barcodes between expression and metadata.")
      expr <- expr[common_cells, , drop = FALSE]
      obs <- obs[common_cells, , drop = FALSE]
      validate_cacoa_metadata(obs)
      common_genes <- intersect(colnames(expr), rownames(var))
      if (length(common_genes) == 0) stop("No overlapping genes between expression and gene annotation.")
      expr <- expr[, common_genes, drop = FALSE]
      var <- var[common_genes, , drop = FALSE]
      list(expr = expr, obs = obs, var = var)
    }); aligned <- x$value; align_seconds <- x$seconds
    expr <- aligned$expr; obs <- aligned$obs; var <- aligned$var

    n_cells <- nrow(expr)
    n_genes_original <- ncol(expr)

    x <- timed("create_seurat_object", {
      expr_sparse <- Matrix(as.matrix(expr), sparse = TRUE)
      expr_sparse <- t(expr_sparse)
      so <- CreateSeuratObject(counts = expr_sparse, meta.data = obs)
      if (nrow(var) == nrow(so[["RNA"]])) {
        so[["RNA"]]@meta.features <- var
      }
      so
    }); so <- x$value; make_seurat_seconds <- x$seconds
    rm(expr, aligned); gc()

    x <- timed("normalize_data", NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)); so <- x$value; normalize_seconds <- x$seconds
    x <- timed("find_variable_features", FindVariableFeatures(so, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)); so <- x$value; hvg_seconds <- x$seconds

    hvg_genes <- VariableFeatures(so)
    if (length(hvg_genes) < 3) stop("Too few HVGs selected.")

    x <- timed("subset_to_hvgs", subset(so, features = hvg_genes)); so <- x$value; subset_seconds <- x$seconds
    cat("Subsetted Seurat object to HVGs:", nrow(so), "genes x", ncol(so), "cells\n")

    x <- timed("scale_data_hvgs", ScaleData(so, features = rownames(so), verbose = FALSE)); so <- x$value; scale_seconds <- x$seconds
    x <- timed("run_pca_hvgs", RunPCA(so, features = rownames(so), npcs = min(n_pcs, nrow(so) - 1), verbose = FALSE)); so <- x$value; pca_seconds <- x$seconds

    harmony_theta_value <- suppressWarnings(as.numeric(row$harmony_theta))
    if (is.na(harmony_theta_value)) harmony_theta_value <- default_harmony_theta
    pca_dims_available <- seq_len(ncol(Embeddings(so, "pca")))
    dims_use <- harmony_dims[harmony_dims %in% pca_dims_available]
    if (length(dims_use) == 0) stop("No valid Harmony dimensions available.")

    cat("Running Harmony with theta =", harmony_theta_value, "\n")
    x <- timed("run_harmony", RunHarmony(so, group.by.vars = harmony_batch_col, dims.use = dims_use, theta = harmony_theta_value, verbose = FALSE)); so <- x$value; harmony_seconds <- x$seconds

    x <- timed("save_hvg_harmony_rds", {
      write.csv(data.frame(gene = hvg_genes), hvg_file, row.names = FALSE)
      saveRDS(so, seurat_rds)
      TRUE
    }); save_seconds <- x$seconds

    elapsed <- round(proc.time()[["elapsed"]] - t_all, 3)
    data.frame(
      row_df,
      dataset_id = dataset_id,
      hvg_seurat_rds = seurat_rds,
      hvg_file = hvg_file,
      preprocess_status = "success",
      error = NA_character_,
      n_genes_original = as.integer(n_genes_original),
      n_hvgs = as.integer(length(hvg_genes)),
      n_cells = as.integer(n_cells),
      read_expression_seconds = read_expression_seconds,
      read_metadata_seconds = read_metadata_seconds,
      read_var_seconds = read_var_seconds,
      align_seconds = align_seconds,
      make_seurat_seconds = make_seurat_seconds,
      normalize_seconds = normalize_seconds,
      hvg_seconds = hvg_seconds,
      subset_seconds = subset_seconds,
      scale_seconds = scale_seconds,
      pca_seconds = pca_seconds,
      harmony_seconds = harmony_seconds,
      save_seconds = save_seconds,
      preprocess_seconds = elapsed,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(row_df, dataset_id = dataset_id, hvg_seurat_rds = seurat_rds, hvg_file = hvg_file, preprocess_status = "failed", error = conditionMessage(e), stringsAsFactors = FALSE)
  })
}

# -----------------------------
# Main
# -----------------------------
parse_cli_args()
dir.create(preprocess_dir, recursive = TRUE, showWarnings = FALSE)
cat("Base directory:", base_dir, "\n")
cat("Export directory:", export_dir, "\n")
cat("Manifest:", manifest_file, "\n")
cat("Preprocess directory:", preprocess_dir, "\n")
cat("Workers:", preprocess_workers, "backend:", preprocess_backend, "\n")

manifest <- load_manifest()
cat("Datasets to preprocess:", nrow(manifest), "\n")

tasks <- split(manifest, seq_len(nrow(manifest)))
if (preprocess_workers <= 1 || identical(preprocess_backend, "serial")) {
  rows <- lapply(tasks, process_one)
} else {
  rows <- parallel::mclapply(tasks, process_one, mc.cores = preprocess_workers, mc.preschedule = FALSE)
}

summary_df <- bind_rows(rows)
manifest_out <- file.path(preprocess_dir, "cacoa_hvg_harmony_preprocess_manifest.csv")
failed_out <- file.path(preprocess_dir, "cacoa_hvg_harmony_preprocess_failed.csv")
write.csv(summary_df, manifest_out, row.names = FALSE)
write.csv(summary_df[summary_df$preprocess_status == "failed", , drop = FALSE], failed_out, row.names = FALSE)

cat("\nDone CACOA Stage 1.\n")
cat("Successful/skipped:", sum(summary_df$preprocess_status %in% c("success", "skipped_existing")), "\n")
cat("Failed:", sum(summary_df$preprocess_status == "failed"), "\n")
cat("Preprocess manifest:", manifest_out, "\n")
cat("Failed file:", failed_out, "\n")
