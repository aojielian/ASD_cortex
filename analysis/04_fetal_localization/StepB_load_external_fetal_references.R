#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Minimal argument parser
# -----------------------------
parse_args <- function(argv) {
  out <- list(
    base_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex",
    outdir = NULL,
    load_nowakowski_full = TRUE,
    preview_rows = 200L,
    gse132672_build_minimal_meta = TRUE
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    val <- if (i < length(argv)) argv[[i + 1L]] else NA_character_
    if (key == "--base_dir") {
      out$base_dir <- val; i <- i + 2L
    } else if (key == "--outdir") {
      out$outdir <- val; i <- i + 2L
    } else if (key == "--load_nowakowski_full") {
      out$load_nowakowski_full <- tolower(val) %in% c("true", "t", "1", "yes", "y"); i <- i + 2L
    } else if (key == "--preview_rows") {
      out$preview_rows <- as.integer(val); i <- i + 2L
    } else if (key == "--gse132672_build_minimal_meta") {
      out$gse132672_build_minimal_meta <- tolower(val) %in% c("true", "t", "1", "yes", "y"); i <- i + 2L
    } else if (key %in% c("-h", "--help")) {
      cat(
"Usage:\n",
"  Rscript StepB_load_external_fetal_references.R [options]\n\n",
"Options:\n",
"  --base_dir <path>                      Base ASD_cortex directory\n",
"  --outdir <path>                        Output directory (default: <base_dir>/stepB_external_reference_loading)\n",
"  --load_nowakowski_full <TRUE/FALSE>    Whether to load full Nowakowski matrix into RDS (default: TRUE)\n",
"  --preview_rows <int>                   Number of preview rows to export for large dense matrices (default: 200)\n",
"  --gse132672_build_minimal_meta <T/F>   Build minimal metadata from organoid cell IDs (default: TRUE)\n",
sep = "")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
  }
  if (is.null(out$outdir) || is.na(out$outdir) || out$outdir == "") {
    out$outdir <- file.path(out$base_dir, "stepB_external_reference_loading")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "StepB_load_external_fetal_references.log")
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

first_line_split <- function(path, gz = FALSE, n = 10L) {
  cmd <- if (gz) {
    sprintf("zcat %s | head -1", shQuote(path))
  } else {
    sprintf("head -1 %s", shQuote(path))
  }
  ln <- run_cmd(cmd)
  strsplit(ln, "\t", fixed = FALSE)[[1]][seq_len(min(length(strsplit(ln, "\t", fixed = FALSE)[[1]]), n))]
}

count_fields_lines <- function(path, gz = FALSE) {
  nf <- if (gz) run_cmd(sprintf("zcat %s | awk -F'\t' 'NR==1{print NF}'", shQuote(path))) else run_cmd(sprintf("awk -F'\t' 'NR==1{print NF}' %s", shQuote(path)))
  nr <- if (gz) run_cmd(sprintf("zcat %s | awk 'END{print NR}'", shQuote(path))) else run_cmd(sprintf("awk 'END{print NR}' %s", shQuote(path)))
  list(header_fields = as.integer(nf[1]), total_lines = as.integer(nr[1]))
}

read_text_lines <- function(path, gz = FALSE, n = 5L) {
  cmd <- if (gz) sprintf("zcat %s | head -%d", shQuote(path), n) else sprintf("head -%d %s", n, shQuote(path))
  run_cmd(cmd)
}

normalize_bool <- function(x) ifelse(is.na(x), FALSE, x)

guess_value_type <- function(x) {
  suppressWarnings(nums <- as.numeric(x))
  nums <- nums[is.finite(nums)]
  if (!length(nums)) return("unknown")
  nonint <- any(abs(nums - round(nums)) > 1e-8)
  if (nonint) return("continuous_normalized_or_log")
  if (all(nums >= 0)) return("integer_counts_like")
  "unknown"
}

standardize_colnames <- function(nms) {
  out <- gsub("[^A-Za-z0-9]+", "_", nms)
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  out
}

# -----------------------------
# Dataset-specific loaders
# -----------------------------
read_nowakowski <- function(base_dir, outdir, load_full = TRUE) {
  ds_out <- file.path(outdir, "Nowakowski")
  safe_dir(ds_out)

  expr_path <- file.path(base_dir, "Nowakowski", "exprMatrix.tsv.gz")
  meta_path <- file.path(base_dir, "Nowakowski", "meta.tsv")
  umap_path <- file.path(base_dir, "Nowakowski", "UMAP.coords.tsv.gz")

  stopifnot(file.exists(expr_path), file.exists(meta_path), file.exists(umap_path))
  log_msg("Loading Nowakowski metadata...")
  meta <- fread(meta_path, sep = "\t", header = TRUE, data.table = FALSE)
  colnames(meta) <- standardize_colnames(colnames(meta))

  log_msg("Loading Nowakowski UMAP with header=FALSE...")
  umap <- fread(cmd = sprintf("zcat %s", shQuote(umap_path)), sep = "\t", header = FALSE, data.table = FALSE)
  stopifnot(ncol(umap) >= 3)
  colnames(umap)[1:3] <- c("Cell", "UMAP1", "UMAP2")

  meta_std <- data.frame(
    cell_id = meta$Cell,
    sample_id = if ("Name" %in% names(meta)) meta$Name else NA_character_,
    donor_id = if ("Name" %in% names(meta)) meta$Name else NA_character_,
    dataset = "Nowakowski_UCSC",
    source_dataset = "Nowakowski",
    modality = "scRNAseq_browser_export",
    age_label = if ("Age_in_Weeks" %in% names(meta)) paste0("gw", meta$Age_in_Weeks) else NA_character_,
    age_weeks = if ("Age_in_Weeks" %in% names(meta)) suppressWarnings(as.numeric(meta$Age_in_Weeks)) else NA_real_,
    tissue = if ("RegionName" %in% names(meta)) meta$RegionName else NA_character_,
    region = if ("RegionName" %in% names(meta)) meta$RegionName else NA_character_,
    cell_type = if ("WGCNAcluster" %in% names(meta)) meta$WGCNAcluster else NA_character_,
    cluster_id = if ("WGCNAcluster" %in% names(meta)) meta$WGCNAcluster else NA_character_,
    laminae = if ("Laminae" %in% names(meta)) meta$Laminae else NA_character_,
    area = if ("Area" %in% names(meta)) meta$Area else NA_character_,
    stringsAsFactors = FALSE
  )

  meta_std <- merge(meta_std, umap[, c("Cell", "UMAP1", "UMAP2")], by.x = "cell_id", by.y = "Cell", all.x = TRUE, sort = FALSE)

  safe_fwrite(meta_std, file.path(ds_out, "Nowakowski_metadata.standardized.tsv.gz"))
  safe_fwrite(as.data.frame(umap), file.path(ds_out, "Nowakowski_umap.tsv.gz"))

  expr_dim <- count_fields_lines(expr_path, gz = TRUE)
  value_preview <- fread(cmd = sprintf("zcat %s | head -6", shQuote(expr_path)), sep = "\t", header = TRUE, data.table = FALSE)
  value_type <- guess_value_type(unlist(value_preview[1:min(5, nrow(value_preview)), 2:min(6, ncol(value_preview)), drop = FALSE]))

  manifest <- data.frame(
    dataset = "Nowakowski_UCSC",
    matrix_file = expr_path,
    matrix_format = "dense_tsv_with_header_gene_by_cell",
    gene_column = "gene",
    n_header_fields = expr_dim$header_fields,
    n_total_lines = expr_dim$total_lines,
    n_cells_from_matrix = expr_dim$header_fields - 1L,
    n_genes_from_matrix = expr_dim$total_lines - 1L,
    n_meta_rows = nrow(meta_std),
    cell_alignment_ok = identical(meta_std$cell_id, colnames(value_preview)[-1]),
    value_type = value_type,
    full_matrix_loaded = load_full,
    notes = "Browser-export matrix; treat as normalized continuous expression rather than raw counts.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(manifest, file.path(ds_out, "Nowakowski_manifest.tsv"))

  if (isTRUE(load_full)) {
    log_msg("Loading full Nowakowski expression matrix into memory...")
    expr_dt <- fread(cmd = sprintf("zcat %s", shQuote(expr_path)), sep = "\t", header = TRUE, data.table = FALSE)
    gene_ids <- expr_dt[[1]]
    expr_mat <- as.matrix(expr_dt[, -1, drop = FALSE])
    rownames(expr_mat) <- gene_ids
    rm(expr_dt)
    gc()

    saveRDS(
      list(
        expr = expr_mat,
        meta = meta_std,
        manifest = manifest
      ),
      file = file.path(ds_out, "Nowakowski_reference_full.rds"),
      compress = "xz"
    )

    gene_table <- data.frame(
      gene_raw = rownames(expr_mat),
      gene_symbol = sub("\\|.*$", "", rownames(expr_mat)),
      stringsAsFactors = FALSE
    )
    safe_fwrite(gene_table, file.path(ds_out, "Nowakowski_genes.tsv.gz"))
  }

  list(meta = meta_std, manifest = manifest)
}

read_gse162170_rna <- function(base_dir, outdir, preview_rows = 200L) {
  ds_out <- file.path(outdir, "GSE162170_RNA")
  safe_dir(ds_out)

  matrix_path <- file.path(base_dir, "GSE162170", "GSE162170_rna_counts.tsv.gz")
  meta_path <- file.path(base_dir, "GSE162170", "GSE162170_rna_cell_metadata.txt.gz")
  cells_path <- file.path(base_dir, "GSE162170", "GSE162170_rna_cell_names.txt.gz")
  stopifnot(file.exists(matrix_path), file.exists(meta_path), file.exists(cells_path))

  log_msg("Loading GSE162170 RNA metadata...")
  meta <- fread(cmd = sprintf("zcat %s", shQuote(meta_path)), sep = "\t", header = TRUE, data.table = FALSE)
  colnames(meta) <- standardize_colnames(colnames(meta))

  log_msg("Loading GSE162170 RNA cell names with header=FALSE...")
  cell_names <- fread(cmd = sprintf("zcat %s", shQuote(cells_path)), sep = "\t", header = FALSE, data.table = FALSE)
  stopifnot(ncol(cell_names) >= 1)
  cell_names_vec <- as.character(cell_names[[1]])

  meta_std <- data.frame(
    cell_id = meta$Cell_ID,
    sample_id = if ("Sample_ID" %in% names(meta)) meta$Sample_ID else NA_character_,
    donor_id = if ("Tissue_ID" %in% names(meta)) meta$Tissue_ID else if ("Sample_ID" %in% names(meta)) meta$Sample_ID else NA_character_,
    dataset = "GSE162170_RNA",
    source_dataset = "GSE162170",
    modality = "scRNAseq",
    age_label = if ("Age" %in% names(meta)) meta$Age else NA_character_,
    age_weeks = suppressWarnings(as.numeric(gsub("[^0-9.]", "", if ("Age" %in% names(meta)) meta$Age else NA_character_))),
    tissue = if ("Tissue_ID" %in% names(meta)) meta$Tissue_ID else NA_character_,
    region = if ("Tissue_ID" %in% names(meta)) meta$Tissue_ID else NA_character_,
    sample_type = if ("Sample_Type" %in% names(meta)) meta$Sample_Type else NA_character_,
    cell_type = if ("DF_classification" %in% names(meta)) meta$DF_classification else if ("seurat_clusters" %in% names(meta)) as.character(meta$seurat_clusters) else NA_character_,
    cluster_id = if ("seurat_clusters" %in% names(meta)) as.character(meta$seurat_clusters) else NA_character_,
    cell_barcode = if ("Cell_Barcode" %in% names(meta)) meta$Cell_Barcode else NA_character_,
    batch = if ("Batch" %in% names(meta)) meta$Batch else NA_character_,
    stringsAsFactors = FALSE
  )

  safe_fwrite(meta_std, file.path(ds_out, "GSE162170_RNA_metadata.standardized.tsv.gz"))
  safe_fwrite(data.frame(cell_id = cell_names_vec, stringsAsFactors = FALSE), file.path(ds_out, "GSE162170_RNA_cell_names.tsv.gz"))

  dims <- count_fields_lines(matrix_path, gz = TRUE)
  preview <- fread(cmd = sprintf("zcat %s | head -%d", shQuote(matrix_path), max(6L, preview_rows + 1L)), sep = "\t", header = FALSE, data.table = FALSE)
  colnames(preview) <- c("gene_id", paste0("cell_", seq_len(ncol(preview) - 1L)))
  gene_table <- data.frame(gene_id = preview[[1]], stringsAsFactors = FALSE)
  safe_fwrite(gene_table, file.path(ds_out, "GSE162170_RNA_gene_preview.tsv.gz"))

  preview_export <- preview[seq_len(min(preview_rows, nrow(preview))), , drop = FALSE]
  mapped_cell_names <- cell_names_vec[seq_len(min(length(cell_names_vec), ncol(preview_export) - 1L))]
  colnames(preview_export) <- c("gene_id", mapped_cell_names)
  safe_fwrite(preview_export, file.path(ds_out, sprintf("GSE162170_RNA_preview_first_%d_rows.tsv.gz", nrow(preview_export))))

  value_type <- guess_value_type(unlist(preview_export[1:min(5, nrow(preview_export)), 2:min(6, ncol(preview_export)), drop = FALSE]))
  alignment_ok <- identical(meta_std$cell_id, cell_names_vec)

  manifest <- data.frame(
    dataset = "GSE162170_RNA",
    matrix_file = matrix_path,
    matrix_format = "dense_tsv_no_header_gene_by_cell_with_external_cell_names",
    gene_column = "V1",
    external_cell_names_file = cells_path,
    n_header_fields = dims$header_fields,
    n_total_lines = dims$total_lines,
    n_cells_from_matrix = dims$header_fields - 1L,
    n_genes_from_matrix = dims$total_lines,
    n_meta_rows = nrow(meta_std),
    n_cell_names = length(cell_names_vec),
    cell_alignment_ok = alignment_ok,
    value_type = value_type,
    full_matrix_loaded = FALSE,
    notes = "Matrix appears to have no header. First column is ENSG gene IDs; columns 2..N map to external cell names in order.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(manifest, file.path(ds_out, "GSE162170_RNA_manifest.tsv"))

  list(meta = meta_std, manifest = manifest)
}

read_gse162170_multiome <- function(base_dir, outdir, preview_rows = 200L) {
  ds_out <- file.path(outdir, "GSE162170_multiome")
  safe_dir(ds_out)

  matrix_path <- file.path(base_dir, "GSE162170", "GSE162170_multiome_rna_logcounts.tsv.gz")
  meta_path <- file.path(base_dir, "GSE162170", "GSE162170_multiome_cell_metadata.txt.gz")
  cluster_path <- file.path(base_dir, "GSE162170", "GSE162170_multiome_cluster_names.txt.gz")
  stopifnot(file.exists(matrix_path), file.exists(meta_path), file.exists(cluster_path))

  log_msg("Loading GSE162170 multiome metadata...")
  meta <- fread(cmd = sprintf("zcat %s", shQuote(meta_path)), sep = "\t", header = TRUE, data.table = FALSE)
  colnames(meta) <- standardize_colnames(colnames(meta))

  log_msg("Loading GSE162170 multiome cluster map...")
  cluster_map <- fread(cmd = sprintf("zcat %s", shQuote(cluster_path)), sep = "\t", header = TRUE, data.table = FALSE)
  colnames(cluster_map) <- standardize_colnames(colnames(cluster_map))

  cluster_lut <- NULL
  if (all(c("Cluster_ID", "Cluster_Name") %in% names(cluster_map))) {
    cluster_lut <- setNames(cluster_map$Cluster_Name, cluster_map$Cluster_ID)
  }

  cluster_name <- if (!is.null(cluster_lut) && "seurat_clusters" %in% names(meta)) {
    unname(cluster_lut[as.character(meta$seurat_clusters)])
  } else {
    rep(NA_character_, nrow(meta))
  }

  meta_std <- data.frame(
    cell_id = meta$Cell_ID,
    sample_id = if ("Sample_ID" %in% names(meta)) meta$Sample_ID else NA_character_,
    donor_id = if ("Dissociation_ID" %in% names(meta)) meta$Dissociation_ID else if ("Sample_ID" %in% names(meta)) meta$Sample_ID else NA_character_,
    dataset = "GSE162170_multiome",
    source_dataset = "GSE162170",
    modality = "multiome_RNA",
    age_label = if ("Sample_Age" %in% names(meta)) meta$Sample_Age else NA_character_,
    age_weeks = suppressWarnings(as.numeric(gsub("[^0-9.]", "", if ("Sample_Age" %in% names(meta)) meta$Sample_Age else NA_character_))),
    tissue = if ("Sample_ID" %in% names(meta)) meta$Sample_ID else NA_character_,
    region = if ("Sample_ID" %in% names(meta)) meta$Sample_ID else NA_character_,
    cell_type = ifelse(!is.na(cluster_name), cluster_name,
                       if ("DF_classification" %in% names(meta)) meta$DF_classification else if ("seurat_clusters" %in% names(meta)) as.character(meta$seurat_clusters) else NA_character_),
    cluster_id = if ("seurat_clusters" %in% names(meta)) as.character(meta$seurat_clusters) else NA_character_,
    cluster_name = cluster_name,
    cell_barcode = if ("Cell_Barcode" %in% names(meta)) meta$Cell_Barcode else NA_character_,
    batch = if ("Sample_Batch" %in% names(meta)) meta$Sample_Batch else NA_character_,
    stringsAsFactors = FALSE
  )

  safe_fwrite(meta_std, file.path(ds_out, "GSE162170_multiome_metadata.standardized.tsv.gz"))
  safe_fwrite(cluster_map, file.path(ds_out, "GSE162170_multiome_cluster_map.tsv.gz"))

  dims <- count_fields_lines(matrix_path, gz = TRUE)
  preview <- fread(cmd = sprintf("zcat %s | head -%d", shQuote(matrix_path), max(6L, preview_rows + 1L)), sep = "\t", header = FALSE, data.table = FALSE)
  colnames(preview) <- c("gene_id", paste0("cell_", seq_len(ncol(preview) - 1L)))
  safe_fwrite(data.frame(gene_id = preview[[1]], stringsAsFactors = FALSE), file.path(ds_out, "GSE162170_multiome_gene_preview.tsv.gz"))

  # Unlike RNA counts, there is no separate cell-name file. Use metadata ordering as provisional matrix order.
  meta_cells <- meta_std$cell_id
  preview_export <- preview[seq_len(min(preview_rows, nrow(preview))), , drop = FALSE]
  mapped_cell_names <- meta_cells[seq_len(min(length(meta_cells), ncol(preview_export) - 1L))]
  colnames(preview_export) <- c("gene_id", mapped_cell_names)
  safe_fwrite(preview_export, file.path(ds_out, sprintf("GSE162170_multiome_preview_first_%d_rows.tsv.gz", nrow(preview_export))))

  value_type <- guess_value_type(unlist(preview_export[1:min(5, nrow(preview_export)), 2:min(6, ncol(preview_export)), drop = FALSE]))

  manifest <- data.frame(
    dataset = "GSE162170_multiome",
    matrix_file = matrix_path,
    matrix_format = "dense_tsv_no_header_gene_by_cell_provisional_metadata_order",
    gene_column = "V1",
    external_cell_names_file = NA_character_,
    n_header_fields = dims$header_fields,
    n_total_lines = dims$total_lines,
    n_cells_from_matrix = dims$header_fields - 1L,
    n_genes_from_matrix = dims$total_lines,
    n_meta_rows = nrow(meta_std),
    n_cell_names = nrow(meta_std),
    cell_alignment_ok = (nrow(meta_std) == (dims$header_fields - 1L)),
    value_type = value_type,
    full_matrix_loaded = FALSE,
    notes = "Matrix appears to have no header. Column order is provisionally assumed to match metadata row order; verify with downstream sanity checks before biological analysis.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(manifest, file.path(ds_out, "GSE162170_multiome_manifest.tsv"))

  list(meta = meta_std, manifest = manifest)
}

read_gse132672 <- function(base_dir, outdir, preview_rows = 200L, build_minimal_meta = TRUE) {
  ds_out <- file.path(outdir, "GSE132672")
  safe_dir(ds_out)

  matrix_path <- file.path(base_dir, "GSE132672", "GSE132672_allorganoids_withnew_matrix.txt.gz")
  stopifnot(file.exists(matrix_path))

  dims <- count_fields_lines(matrix_path, gz = TRUE)
  preview <- fread(cmd = sprintf("zcat %s | head -%d", shQuote(matrix_path), max(6L, preview_rows + 1L)), sep = "\t", header = TRUE, data.table = FALSE)
  colnames(preview)[1] <- "gene_id"
  safe_fwrite(data.frame(gene_id = preview$gene_id, stringsAsFactors = FALSE), file.path(ds_out, "GSE132672_gene_preview.tsv.gz"))
  safe_fwrite(preview[seq_len(min(preview_rows, nrow(preview))), , drop = FALSE], file.path(ds_out, sprintf("GSE132672_preview_first_%d_rows.tsv.gz", min(preview_rows, nrow(preview)))))

  cell_ids <- colnames(preview)[-1]
  minimal_meta <- NULL
  if (isTRUE(build_minimal_meta)) {
    sample_id <- sub("_[^_]+$", "", cell_ids)
    week_token <- ifelse(grepl("Week[0-9]+", sample_id, ignore.case = TRUE),
                         regmatches(sample_id, regexpr("Week[0-9]+", sample_id, ignore.case = TRUE)),
                         NA_character_)
    age_weeks <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", week_token)))

    minimal_meta <- data.frame(
      cell_id = cell_ids,
      sample_id = sample_id,
      donor_id = sample_id,
      dataset = "GSE132672_organoid",
      source_dataset = "GSE132672",
      modality = "organoid_scRNA_matrix_only",
      age_label = week_token,
      age_weeks = age_weeks,
      tissue = "organoid",
      region = NA_character_,
      cell_type = NA_character_,
      cluster_id = NA_character_,
      stringsAsFactors = FALSE
    )
    safe_fwrite(minimal_meta, file.path(ds_out, "GSE132672_minimal_metadata.standardized.tsv.gz"))
  }

  value_type <- guess_value_type(unlist(preview[1:min(5, nrow(preview)), 2:min(6, ncol(preview)), drop = FALSE]))
  manifest <- data.frame(
    dataset = "GSE132672_organoid",
    matrix_file = matrix_path,
    matrix_format = "dense_tsv_with_header_gene_by_cell",
    gene_column = "V1",
    external_cell_names_file = NA_character_,
    n_header_fields = dims$header_fields,
    n_total_lines = dims$total_lines,
    n_cells_from_matrix = dims$header_fields - 1L,
    n_genes_from_matrix = dims$total_lines - 1L,
    n_meta_rows = if (is.null(minimal_meta)) NA_integer_ else nrow(minimal_meta),
    n_cell_names = dims$header_fields - 1L,
    cell_alignment_ok = TRUE,
    value_type = value_type,
    full_matrix_loaded = FALSE,
    notes = "Expression matrix only; no curated metadata supplied in current download. Use as supportive/background context unless external metadata are added.",
    stringsAsFactors = FALSE
  )
  safe_fwrite(manifest, file.path(ds_out, "GSE132672_manifest.tsv"))

  list(meta = minimal_meta, manifest = manifest)
}

# -----------------------------
# Main
# -----------------------------
log_msg("Starting StepB external fetal reference loading...")
log_msg("base_dir = ", args$base_dir)
log_msg("outdir   = ", args$outdir)

all_manifests <- list()
all_summaries <- list()

nw <- read_nowakowski(args$base_dir, args$outdir, load_full = args$load_nowakowski_full)
all_manifests[[length(all_manifests) + 1L]] <- nw$manifest
all_summaries[[length(all_summaries) + 1L]] <- data.frame(
  dataset = nw$manifest$dataset,
  n_cells = nw$manifest$n_cells_from_matrix,
  n_genes = nw$manifest$n_genes_from_matrix,
  n_meta_rows = nw$manifest$n_meta_rows,
  value_type = nw$manifest$value_type,
  full_matrix_loaded = nw$manifest$full_matrix_loaded,
  primary_role = "main_fetal_replication",
  stringsAsFactors = FALSE
)

g1 <- read_gse162170_rna(args$base_dir, args$outdir, preview_rows = args$preview_rows)
all_manifests[[length(all_manifests) + 1L]] <- g1$manifest
all_summaries[[length(all_summaries) + 1L]] <- data.frame(
  dataset = g1$manifest$dataset,
  n_cells = g1$manifest$n_cells_from_matrix,
  n_genes = g1$manifest$n_genes_from_matrix,
  n_meta_rows = g1$manifest$n_meta_rows,
  value_type = g1$manifest$value_type,
  full_matrix_loaded = g1$manifest$full_matrix_loaded,
  primary_role = "main_geo_rna_replication_candidate",
  stringsAsFactors = FALSE
)

g2 <- read_gse162170_multiome(args$base_dir, args$outdir, preview_rows = args$preview_rows)
all_manifests[[length(all_manifests) + 1L]] <- g2$manifest
all_summaries[[length(all_summaries) + 1L]] <- data.frame(
  dataset = g2$manifest$dataset,
  n_cells = g2$manifest$n_cells_from_matrix,
  n_genes = g2$manifest$n_genes_from_matrix,
  n_meta_rows = g2$manifest$n_meta_rows,
  value_type = g2$manifest$value_type,
  full_matrix_loaded = g2$manifest$full_matrix_loaded,
  primary_role = "supportive_multiome_corroboration",
  stringsAsFactors = FALSE
)

g3 <- read_gse132672(args$base_dir, args$outdir, preview_rows = args$preview_rows,
                     build_minimal_meta = args$gse132672_build_minimal_meta)
all_manifests[[length(all_manifests) + 1L]] <- g3$manifest
all_summaries[[length(all_summaries) + 1L]] <- data.frame(
  dataset = g3$manifest$dataset,
  n_cells = g3$manifest$n_cells_from_matrix,
  n_genes = g3$manifest$n_genes_from_matrix,
  n_meta_rows = g3$manifest$n_meta_rows,
  value_type = g3$manifest$value_type,
  full_matrix_loaded = g3$manifest$full_matrix_loaded,
  primary_role = "supportive_background_organoid_only",
  stringsAsFactors = FALSE
)

summary_dt <- rbindlist(all_summaries, fill = TRUE)
manifests_dt <- rbindlist(all_manifests, fill = TRUE)

safe_fwrite(summary_dt, file.path(args$outdir, "StepB_dataset_loading_summary.tsv"))
safe_fwrite(manifests_dt, file.path(args$outdir, "StepB_all_manifests.tsv"))

sink(file.path(args$outdir, "sessionInfo.txt"))
print(sessionInfo())
sink()

log_msg("Finished StepB external fetal reference loading.")
