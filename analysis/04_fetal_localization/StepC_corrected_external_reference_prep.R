#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(argv) {
  out <- list(
    base_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex",
    stepb_dir = NULL,
    outdir = NULL,
    preview_rows = 200L,
    preview_cols = 50L,
    nowakowski_use_stepb_rds = TRUE
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    val <- if (i < length(argv)) argv[[i + 1L]] else NA_character_
    if (key == "--base_dir") {
      out$base_dir <- val; i <- i + 2L
    } else if (key == "--stepb_dir") {
      out$stepb_dir <- val; i <- i + 2L
    } else if (key == "--outdir") {
      out$outdir <- val; i <- i + 2L
    } else if (key == "--preview_rows") {
      out$preview_rows <- as.integer(val); i <- i + 2L
    } else if (key == "--preview_cols") {
      out$preview_cols <- as.integer(val); i <- i + 2L
    } else if (key == "--nowakowski_use_stepb_rds") {
      out$nowakowski_use_stepb_rds <- tolower(val) %in% c("true","t","1","yes","y"); i <- i + 2L
    } else if (key %in% c("-h", "--help")) {
      cat(
"Usage:\n",
"  Rscript StepC_corrected_external_reference_prep.R [options]\n\n",
"Options:\n",
"  --base_dir <path>                   Base ASD_cortex directory\n",
"  --stepb_dir <path>                  StepB output directory (default: <base_dir>/stepB_external_reference_loading)\n",
"  --outdir <path>                     Output directory (default: <base_dir>/stepC_external_reference_prep)\n",
"  --preview_rows <int>                Gene rows to export in preview matrices (default: 200)\n",
"  --preview_cols <int>                Cell columns to export in preview matrices (default: 50)\n",
"  --nowakowski_use_stepb_rds <T/F>    Use StepB full RDS for Nowakowski if available (default: TRUE)\n",
sep = "")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
  }
  if (is.null(out$stepb_dir) || is.na(out$stepb_dir) || out$stepb_dir == "") {
    out$stepb_dir <- file.path(out$base_dir, "stepB_external_reference_loading")
  }
  if (is.null(out$outdir) || is.na(out$outdir) || out$outdir == "") {
    out$outdir <- file.path(out$base_dir, "stepC_external_reference_prep")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "StepC_corrected_external_reference_prep.log")
cat("", file = log_file)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_dir <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)

safe_fwrite <- function(x, path, sep = "\t") {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    fwrite(x, file = path, sep = sep, quote = FALSE, na = "NA", compress = "gzip")
  } else {
    fwrite(x, file = path, sep = sep, quote = FALSE, na = "NA")
  }
}

run_cmd <- function(cmd) {
  log_msg("CMD: ", cmd)
  out <- tryCatch(system(cmd, intern = TRUE), warning = function(w) system(cmd, intern = TRUE), error = function(e) e)
  if (inherits(out, "error")) stop(out$message)
  out
}

standardize_colnames <- function(nms) {
  out <- gsub("[^A-Za-z0-9]+", "_", nms)
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  out
}

read_gz_table <- function(path, header = TRUE) {
  fread(cmd = sprintf("zcat %s", shQuote(path)), sep = "\t", header = header,
        data.table = FALSE, check.names = FALSE, fill = TRUE)
}

split_tab_line <- function(x) {
  strsplit(x, "\t", fixed = TRUE)[[1]]
}

extract_matrix_header_cells <- function(matrix_path) {
  header_line <- run_cmd(sprintf("zcat %s | head -1", shQuote(matrix_path)))[1]
  second_line <- run_cmd(sprintf("zcat %s | head -2 | tail -1", shQuote(matrix_path)))[1]
  header_fields <- split_tab_line(header_line)
  second_fields <- split_tab_line(second_line)

  has_blank_gene_placeholder <- length(header_fields) == length(second_fields) && nzchar(second_fields[1]) && grepl("^(ENSG|[A-Za-z0-9_.-]+\\|)", second_fields[1]) && !nzchar(header_fields[1])

  if (has_blank_gene_placeholder) {
    cell_ids <- header_fields[-1]
    gene_col_name <- "gene_id"
  } else if (length(header_fields) == (length(second_fields) - 1L)) {
    cell_ids <- header_fields
    gene_col_name <- "gene_id"
  } else {
    # fallback: if first header token looks like a cell id, keep full header as cells
    cell_ids <- if (length(header_fields) > 1L && grepl("_", header_fields[1])) header_fields else header_fields[-1]
    gene_col_name <- "gene_id"
  }

  list(
    cell_ids = cell_ids,
    n_cells = length(cell_ids),
    has_blank_gene_placeholder = has_blank_gene_placeholder,
    first_header_token = if (length(header_fields)) header_fields[1] else NA_character_,
    first_gene_token = if (length(second_fields)) second_fields[1] else NA_character_
  )
}

count_fields_lines <- function(path) {
  nf <- run_cmd(sprintf("zcat %s | awk -F'\\t' 'NR==1{print NF}'", shQuote(path)))
  nr <- run_cmd(sprintf("zcat %s | awk 'END{print NR}'", shQuote(path)))
  list(header_fields = as.integer(nf[1]), total_lines = as.integer(nr[1]))
}

sample_value_type <- function(matrix_path, nlines = 201L, ncols = 200L) {
  dt <- fread(cmd = sprintf("zcat %s | head -%d", shQuote(matrix_path), nlines), sep = "\t", header = TRUE,
              data.table = FALSE, check.names = FALSE, fill = TRUE)
  if (ncol(dt) < 2L) return(list(value_type = "unknown", n_values = 0L, frac_noninteger = NA_real_, min = NA_real_, max = NA_real_))
  vals <- unlist(dt[seq_len(min(nrow(dt), nlines - 1L)), 2:min(ncol(dt), ncols + 1L), drop = FALSE], use.names = FALSE)
  suppressWarnings(nums <- as.numeric(vals))
  nums <- nums[is.finite(nums)]
  if (!length(nums)) return(list(value_type = "unknown", n_values = 0L, frac_noninteger = NA_real_, min = NA_real_, max = NA_real_))
  frac_noninteger <- mean(abs(nums - round(nums)) > 1e-8)
  value_type <- if (frac_noninteger > 0) "continuous_normalized_or_log" else "integer_counts_like"
  list(value_type = value_type, n_values = length(nums), frac_noninteger = frac_noninteger,
       min = min(nums), max = max(nums))
}

write_gene_ids <- function(matrix_path, out_path) {
  cmd <- sprintf("zcat %s | tail -n +2 | cut -f1", shQuote(matrix_path))
  gene_ids <- run_cmd(cmd)
  safe_fwrite(data.frame(gene_id = gene_ids, stringsAsFactors = FALSE), out_path)
}

write_matrix_preview <- function(matrix_path, out_path, nrows = 200L, ncols = 50L) {
  dt <- fread(cmd = sprintf("zcat %s | head -%d", shQuote(matrix_path), nrows + 1L), sep = "\t", header = TRUE,
              data.table = FALSE, check.names = FALSE, fill = TRUE)
  keep_cols <- seq_len(min(ncol(dt), ncols + 1L))
  dt <- dt[, keep_cols, drop = FALSE]
  colnames(dt)[1] <- "gene_id"
  safe_fwrite(dt, out_path)
}

extract_shared_ids <- function(a, b) {
  shared <- intersect(a, b)
  shared
}

standardize_age_weeks <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]+", "", x)))

prepare_nowakowski <- function(base_dir, stepb_dir, outdir, use_stepb_rds = TRUE) {
  ds_out <- file.path(outdir, "Nowakowski")
  safe_dir(ds_out)

  stepb_rds <- file.path(stepb_dir, "Nowakowski", "Nowakowski_reference_full.rds")
  stepb_meta <- file.path(stepb_dir, "Nowakowski", "Nowakowski_metadata.standardized.tsv.gz")
  stepb_manifest <- file.path(stepb_dir, "Nowakowski", "Nowakowski_manifest.tsv")
  stopifnot(file.exists(stepb_meta), file.exists(stepb_manifest))

  meta <- read_gz_table(stepb_meta, header = TRUE)
  manifest <- fread(stepb_manifest, sep = "\t", header = TRUE, data.table = FALSE)
  expr_loaded <- FALSE
  expr_obj <- NULL
  gene_table <- NULL

  if (isTRUE(use_stepb_rds) && file.exists(stepb_rds)) {
    log_msg("Loading StepB Nowakowski full RDS...")
    obj <- readRDS(stepb_rds)
    expr_obj <- obj$expr
    expr_loaded <- TRUE
    gene_table <- data.frame(
      gene_raw = rownames(expr_obj),
      gene_symbol = sub("\\|.*$", "", rownames(expr_obj)),
      stringsAsFactors = FALSE
    )
  } else {
    genes_path <- file.path(stepb_dir, "Nowakowski", "Nowakowski_genes.tsv.gz")
    if (file.exists(genes_path)) {
      gene_table <- read_gz_table(genes_path, header = TRUE)
    }
  }

  cell_alignment_ok <- FALSE
  if (expr_loaded) {
    cell_alignment_ok <- identical(meta$cell_id, colnames(expr_obj))
  }

  if (!is.null(gene_table)) safe_fwrite(gene_table, file.path(ds_out, "Nowakowski_gene_table.tsv.gz"))
  safe_fwrite(meta, file.path(ds_out, "Nowakowski_metadata.program_ready.tsv.gz"))

  out_obj <- list(
    dataset = "Nowakowski_UCSC",
    modality = "scRNAseq_browser_export",
    value_type = "continuous_normalized_or_log",
    matrix_storage = if (expr_loaded) "in_memory_full" else "stepB_external_file",
    expr = expr_obj,
    meta = meta,
    gene_table = gene_table,
    manifest = manifest,
    notes = "Primary fetal replication object; continuous browser-export expression."
  )
  saveRDS(out_obj, file.path(ds_out, "Nowakowski_program_ready.rds"), compress = "xz")

  summary <- data.frame(
    dataset = "Nowakowski_UCSC",
    n_cells = if (expr_loaded) ncol(expr_obj) else manifest$n_cells_from_matrix,
    n_genes = if (expr_loaded) nrow(expr_obj) else manifest$n_genes_from_matrix,
    n_meta_rows = nrow(meta),
    value_type = "continuous_normalized_or_log",
    matrix_header_n_cells = if (expr_loaded) ncol(expr_obj) else manifest$n_cells_from_matrix,
    matched_cells = if (expr_loaded) sum(meta$cell_id %in% colnames(expr_obj)) else NA_integer_,
    dropped_from_metadata = 0L,
    dropped_from_matrix = 0L,
    cell_alignment_ok = cell_alignment_ok,
    primary_role = "main_fetal_replication",
    notes = "Ready for immediate program scoring.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(summary, file.path(ds_out, "Nowakowski_summary.tsv"))
  summary
}

prepare_gse162170_rna <- function(base_dir, outdir, preview_rows = 200L, preview_cols = 50L) {
  ds_out <- file.path(outdir, "GSE162170_RNA")
  safe_dir(ds_out)

  matrix_path <- file.path(base_dir, "GSE162170", "GSE162170_rna_counts.tsv.gz")
  meta_path <- file.path(base_dir, "GSE162170", "GSE162170_rna_cell_metadata.txt.gz")
  cells_path <- file.path(base_dir, "GSE162170", "GSE162170_rna_cell_names.txt.gz")
  stopifnot(file.exists(matrix_path), file.exists(meta_path), file.exists(cells_path))

  meta <- read_gz_table(meta_path, header = TRUE)
  colnames(meta) <- standardize_colnames(colnames(meta))
  header_info <- extract_matrix_header_cells(matrix_path)
  matrix_cells <- header_info$cell_ids
  external_cells <- run_cmd(sprintf("zcat %s", shQuote(cells_path)))
  dims <- count_fields_lines(matrix_path)
  value_info <- sample_value_type(matrix_path)

  safe_fwrite(data.frame(cell_id = matrix_cells, stringsAsFactors = FALSE),
              file.path(ds_out, "GSE162170_RNA_matrix_header_cells.tsv.gz"))
  safe_fwrite(data.frame(cell_id = external_cells, stringsAsFactors = FALSE),
              file.path(ds_out, "GSE162170_RNA_external_cell_names.tsv.gz"))

  header_equals_external <- identical(matrix_cells, external_cells)
  shared_cells <- extract_shared_ids(matrix_cells, meta$Cell_ID)
  matched_cells <- matrix_cells[matrix_cells %in% shared_cells]
  matched_meta <- meta[match(matched_cells, meta$Cell_ID), , drop = FALSE]

  stopifnot(length(matched_cells) == nrow(matched_meta))
  stopifnot(identical(matched_cells, matched_meta$Cell_ID))

  meta_std <- data.frame(
    cell_id = matched_meta$Cell_ID,
    sample_id = if ("Sample_ID" %in% names(matched_meta)) matched_meta$Sample_ID else NA_character_,
    donor_id = if ("Tissue_ID" %in% names(matched_meta)) matched_meta$Tissue_ID else if ("Sample_ID" %in% names(matched_meta)) matched_meta$Sample_ID else NA_character_,
    dataset = "GSE162170_RNA",
    source_dataset = "GSE162170",
    modality = "scRNAseq",
    age_label = if ("Age" %in% names(matched_meta)) matched_meta$Age else NA_character_,
    age_weeks = standardize_age_weeks(if ("Age" %in% names(matched_meta)) matched_meta$Age else NA_character_),
    tissue = if ("Tissue_ID" %in% names(matched_meta)) matched_meta$Tissue_ID else NA_character_,
    region = if ("Tissue_ID" %in% names(matched_meta)) matched_meta$Tissue_ID else NA_character_,
    sample_type = if ("Sample_Type" %in% names(matched_meta)) matched_meta$Sample_Type else NA_character_,
    assay = if ("Assay" %in% names(matched_meta)) matched_meta$Assay else NA_character_,
    batch = if ("Batch" %in% names(matched_meta)) matched_meta$Batch else NA_character_,
    df_classification = if ("DF_classification" %in% names(matched_meta)) matched_meta$DF_classification else NA_character_,
    cell_type = if ("DF_classification" %in% names(matched_meta)) matched_meta$DF_classification else if ("seurat_clusters" %in% names(matched_meta)) as.character(matched_meta$seurat_clusters) else NA_character_,
    cluster_id = if ("seurat_clusters" %in% names(matched_meta)) as.character(matched_meta$seurat_clusters) else NA_character_,
    cell_barcode = if ("Cell_Barcode" %in% names(matched_meta)) matched_meta$Cell_Barcode else NA_character_,
    RNA_Counts = if ("RNA_Counts" %in% names(matched_meta)) matched_meta$RNA_Counts else NA_real_,
    RNA_Features = if ("RNA_Features" %in% names(matched_meta)) matched_meta$RNA_Features else NA_real_,
    stringsAsFactors = FALSE
  )

  safe_fwrite(meta_std, file.path(ds_out, "GSE162170_RNA_metadata.program_ready.tsv.gz"))
  write_gene_ids(matrix_path, file.path(ds_out, "GSE162170_RNA_gene_ids.tsv.gz"))
  write_matrix_preview(matrix_path, file.path(ds_out, sprintf("GSE162170_RNA_preview_first_%d_rows_%d_cells.tsv.gz", preview_rows, preview_cols)), preview_rows, preview_cols)

  alignment_qc <- data.frame(
    matrix_header_n_cells = length(matrix_cells),
    external_cell_names_n = length(external_cells),
    metadata_n_rows = nrow(meta),
    header_equals_external = header_equals_external,
    metadata_exact_match_to_header = identical(meta$Cell_ID, matrix_cells),
    shared_cells = length(shared_cells),
    dropped_from_metadata = nrow(meta) - length(shared_cells),
    dropped_from_matrix = length(matrix_cells) - length(shared_cells),
    value_type = value_info$value_type,
    frac_noninteger_preview = value_info$frac_noninteger,
    preview_min = value_info$min,
    preview_max = value_info$max,
    stringsAsFactors = FALSE
  )
  safe_fwrite(alignment_qc, file.path(ds_out, "GSE162170_RNA_alignment_qc.tsv"))

  manifest <- data.frame(
    dataset = "GSE162170_RNA",
    matrix_file = matrix_path,
    matrix_format = "dense_tsv_with_header_blank_gene_placeholder_gene_by_cell",
    gene_column = "gene_id",
    header_has_blank_gene_placeholder = header_info$has_blank_gene_placeholder,
    n_header_fields = dims$header_fields,
    n_total_lines = dims$total_lines,
    n_cells_from_matrix_header = length(matrix_cells),
    n_genes_from_matrix = dims$total_lines - 1L,
    n_external_cell_names = length(external_cells),
    n_meta_rows = nrow(meta),
    n_matched_cells = length(shared_cells),
    dropped_from_metadata = nrow(meta) - length(shared_cells),
    dropped_from_matrix = length(matrix_cells) - length(shared_cells),
    header_matches_external_cell_names = header_equals_external,
    cell_alignment_ok = identical(matched_cells, meta_std$cell_id),
    value_type = value_info$value_type,
    full_matrix_loaded = FALSE,
    notes = "Corrected StepC interpretation: matrix has a true header row of cell IDs; external cell-names file is concordant and retained as QC only.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(manifest, file.path(ds_out, "GSE162170_RNA_manifest.corrected.tsv"))

  light_obj <- list(
    dataset = "GSE162170_RNA",
    modality = "scRNAseq",
    value_type = value_info$value_type,
    matrix_path = matrix_path,
    matrix_has_header = TRUE,
    matrix_header_blank_gene_placeholder = header_info$has_blank_gene_placeholder,
    cell_ids = matched_cells,
    meta = meta_std,
    gene_ids_file = file.path(ds_out, "GSE162170_RNA_gene_ids.tsv.gz"),
    manifest = manifest,
    notes = "Program-ready light object. Use matrix_path + cell_ids for downstream streaming/extraction."
  )
  saveRDS(light_obj, file.path(ds_out, "GSE162170_RNA_program_ready_light.rds"), compress = "xz")

  summary <- data.frame(
    dataset = "GSE162170_RNA",
    n_cells = length(matched_cells),
    n_genes = dims$total_lines - 1L,
    n_meta_rows = nrow(meta),
    value_type = value_info$value_type,
    matrix_header_n_cells = length(matrix_cells),
    matched_cells = length(shared_cells),
    dropped_from_metadata = nrow(meta) - length(shared_cells),
    dropped_from_matrix = length(matrix_cells) - length(shared_cells),
    cell_alignment_ok = identical(matched_cells, meta_std$cell_id),
    primary_role = "main_geo_rna_replication",
    notes = if (header_equals_external) "Header-corrected and analysis-ready." else "Header-corrected; external cell-name file differs from matrix header.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(summary, file.path(ds_out, "GSE162170_RNA_summary.tsv"))
  summary
}

prepare_gse162170_multiome <- function(base_dir, outdir, preview_rows = 200L, preview_cols = 50L) {
  ds_out <- file.path(outdir, "GSE162170_multiome")
  safe_dir(ds_out)

  matrix_path <- file.path(base_dir, "GSE162170", "GSE162170_multiome_rna_logcounts.tsv.gz")
  meta_path <- file.path(base_dir, "GSE162170", "GSE162170_multiome_cell_metadata.txt.gz")
  cluster_path <- file.path(base_dir, "GSE162170", "GSE162170_multiome_cluster_names.txt.gz")
  stopifnot(file.exists(matrix_path), file.exists(meta_path), file.exists(cluster_path))

  meta <- read_gz_table(meta_path, header = TRUE)
  colnames(meta) <- standardize_colnames(colnames(meta))
  cluster_map <- read_gz_table(cluster_path, header = TRUE)
  colnames(cluster_map) <- standardize_colnames(colnames(cluster_map))

  cluster_lut <- NULL
  if (all(c("Cluster_ID", "Cluster_Name") %in% names(cluster_map))) {
    cluster_lut <- setNames(cluster_map$Cluster_Name, cluster_map$Cluster_ID)
  }

  header_info <- extract_matrix_header_cells(matrix_path)
  matrix_cells <- header_info$cell_ids
  dims <- count_fields_lines(matrix_path)
  value_info <- sample_value_type(matrix_path)

  shared_cells <- extract_shared_ids(matrix_cells, meta$Cell_ID)
  matched_cells <- matrix_cells[matrix_cells %in% shared_cells]
  matched_meta <- meta[match(matched_cells, meta$Cell_ID), , drop = FALSE]

  stopifnot(length(matched_cells) == nrow(matched_meta))
  stopifnot(identical(matched_cells, matched_meta$Cell_ID))

  cluster_name <- if (!is.null(cluster_lut) && "seurat_clusters" %in% names(matched_meta)) {
    unname(cluster_lut[as.character(matched_meta$seurat_clusters)])
  } else {
    rep(NA_character_, nrow(matched_meta))
  }

  meta_std <- data.frame(
    cell_id = matched_meta$Cell_ID,
    sample_id = if ("Sample_ID" %in% names(matched_meta)) matched_meta$Sample_ID else NA_character_,
    donor_id = if ("Dissociation_ID" %in% names(matched_meta)) matched_meta$Dissociation_ID else if ("Sample_ID" %in% names(matched_meta)) matched_meta$Sample_ID else NA_character_,
    dataset = "GSE162170_multiome",
    source_dataset = "GSE162170",
    modality = "multiome_RNA",
    age_label = if ("Sample_Age" %in% names(matched_meta)) matched_meta$Sample_Age else NA_character_,
    age_weeks = standardize_age_weeks(if ("Sample_Age" %in% names(matched_meta)) matched_meta$Sample_Age else NA_character_),
    sample_batch = if ("Sample_Batch" %in% names(matched_meta)) matched_meta$Sample_Batch else NA_character_,
    cell_barcode = if ("Cell_Barcode" %in% names(matched_meta)) matched_meta$Cell_Barcode else NA_character_,
    df_classification = if ("DF_classification" %in% names(matched_meta)) matched_meta$DF_classification else NA_character_,
    cluster_id = if ("seurat_clusters" %in% names(matched_meta)) as.character(matched_meta$seurat_clusters) else NA_character_,
    cluster_name = cluster_name,
    cell_type = ifelse(!is.na(cluster_name), cluster_name,
                       if ("DF_classification" %in% names(matched_meta)) matched_meta$DF_classification else NA_character_),
    RNA_Counts = if ("RNA_Counts" %in% names(matched_meta)) matched_meta$RNA_Counts else NA_real_,
    RNA_Features = if ("RNA_Features" %in% names(matched_meta)) matched_meta$RNA_Features else NA_real_,
    stringsAsFactors = FALSE
  )

  safe_fwrite(meta_std, file.path(ds_out, "GSE162170_multiome_metadata.program_ready.tsv.gz"))
  safe_fwrite(cluster_map, file.path(ds_out, "GSE162170_multiome_cluster_map.tsv.gz"))
  safe_fwrite(data.frame(cell_id = matrix_cells, stringsAsFactors = FALSE), file.path(ds_out, "GSE162170_multiome_matrix_header_cells.tsv.gz"))
  write_gene_ids(matrix_path, file.path(ds_out, "GSE162170_multiome_gene_ids.tsv.gz"))
  write_matrix_preview(matrix_path, file.path(ds_out, sprintf("GSE162170_multiome_preview_first_%d_rows_%d_cells.tsv.gz", preview_rows, preview_cols)), preview_rows, preview_cols)

  dropped_meta_cells <- setdiff(meta$Cell_ID, shared_cells)
  dropped_matrix_cells <- setdiff(matrix_cells, shared_cells)
  safe_fwrite(data.frame(cell_id = dropped_meta_cells, stringsAsFactors = FALSE), file.path(ds_out, "GSE162170_multiome_cells_dropped_from_metadata_side.tsv.gz"))
  safe_fwrite(data.frame(cell_id = dropped_matrix_cells, stringsAsFactors = FALSE), file.path(ds_out, "GSE162170_multiome_cells_dropped_from_matrix_side.tsv.gz"))

  alignment_qc <- data.frame(
    matrix_header_n_cells = length(matrix_cells),
    metadata_n_rows = nrow(meta),
    shared_cells = length(shared_cells),
    dropped_from_metadata = nrow(meta) - length(shared_cells),
    dropped_from_matrix = length(matrix_cells) - length(shared_cells),
    value_type = value_info$value_type,
    frac_noninteger_preview = value_info$frac_noninteger,
    preview_min = value_info$min,
    preview_max = value_info$max,
    stringsAsFactors = FALSE
  )
  safe_fwrite(alignment_qc, file.path(ds_out, "GSE162170_multiome_alignment_qc.tsv"))

  manifest <- data.frame(
    dataset = "GSE162170_multiome",
    matrix_file = matrix_path,
    matrix_format = "dense_tsv_with_header_blank_gene_placeholder_gene_by_cell",
    gene_column = "gene_id",
    header_has_blank_gene_placeholder = header_info$has_blank_gene_placeholder,
    n_header_fields = dims$header_fields,
    n_total_lines = dims$total_lines,
    n_cells_from_matrix_header = length(matrix_cells),
    n_genes_from_matrix = dims$total_lines - 1L,
    n_meta_rows = nrow(meta),
    n_matched_cells = length(shared_cells),
    dropped_from_metadata = nrow(meta) - length(shared_cells),
    dropped_from_matrix = length(matrix_cells) - length(shared_cells),
    cell_alignment_ok = identical(matched_cells, meta_std$cell_id),
    value_type = value_info$value_type,
    full_matrix_loaded = FALSE,
    notes = "Corrected StepC interpretation: matrix has a true header row of cell IDs. Metadata were intersected with matrix header before downstream use.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(manifest, file.path(ds_out, "GSE162170_multiome_manifest.corrected.tsv"))

  light_obj <- list(
    dataset = "GSE162170_multiome",
    modality = "multiome_RNA",
    value_type = value_info$value_type,
    matrix_path = matrix_path,
    matrix_has_header = TRUE,
    matrix_header_blank_gene_placeholder = header_info$has_blank_gene_placeholder,
    cell_ids = matched_cells,
    meta = meta_std,
    gene_ids_file = file.path(ds_out, "GSE162170_multiome_gene_ids.tsv.gz"),
    cluster_map = cluster_map,
    manifest = manifest,
    notes = "Program-ready light object after explicit matrix-header/metadata intersection."
  )
  saveRDS(light_obj, file.path(ds_out, "GSE162170_multiome_program_ready_light.rds"), compress = "xz")

  summary <- data.frame(
    dataset = "GSE162170_multiome",
    n_cells = length(matched_cells),
    n_genes = dims$total_lines - 1L,
    n_meta_rows = nrow(meta),
    value_type = value_info$value_type,
    matrix_header_n_cells = length(matrix_cells),
    matched_cells = length(shared_cells),
    dropped_from_metadata = nrow(meta) - length(shared_cells),
    dropped_from_matrix = length(matrix_cells) - length(shared_cells),
    cell_alignment_ok = identical(matched_cells, meta_std$cell_id),
    primary_role = "supportive_multiome_corroboration",
    notes = "Header-corrected; metadata intersected to remove 1-cell mismatch before analysis.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(summary, file.path(ds_out, "GSE162170_multiome_summary.tsv"))
  summary
}

log_msg("Starting StepC corrected external reference prep...")
log_msg("base_dir = ", args$base_dir)
log_msg("stepb_dir = ", args$stepb_dir)
log_msg("outdir   = ", args$outdir)

summaries <- list()

summaries[[1]] <- prepare_nowakowski(args$base_dir, args$stepb_dir, args$outdir, use_stepb_rds = args$nowakowski_use_stepb_rds)
summaries[[2]] <- prepare_gse162170_rna(args$base_dir, args$outdir, preview_rows = args$preview_rows, preview_cols = args$preview_cols)
summaries[[3]] <- prepare_gse162170_multiome(args$base_dir, args$outdir, preview_rows = args$preview_rows, preview_cols = args$preview_cols)

summary_dt <- rbindlist(summaries, fill = TRUE)
safe_fwrite(summary_dt, file.path(args$outdir, "StepC_dataset_correction_summary.tsv"))

sink(file.path(args$outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

log_msg("Finished StepC corrected external reference prep.")
