#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABLEDIR <- file.path(OUTDIR, "tables")
RDSDIR <- file.path(OUTDIR, "rds")
LOGDIR <- file.path(OUTDIR, "logs")

for (d in c(OUTDIR, METADIR, TABLEDIR, RDSDIR, LOGDIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(LOGDIR, "Step01A_inventory_and_standardize.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

# -----------------------------
# Input paths
# -----------------------------
sfari_file <- file.path(BASE, "sfari.genes.Sand1.txt")
gse102741_file <- file.path(BASE, "GSE102741_raw_counts_GRCh38.p13_NCBI.tsv.gz")
gse64018_file  <- file.path(BASE, "GSE64018_countlevel_12asd_12ctl.txt.gz")
brainspan_dir  <- file.path(BASE, "Gencode_v3c_summarized_to_genes")
brainspan_cols <- file.path(brainspan_dir, "columns_metadata.csv")
brainspan_rows <- file.path(brainspan_dir, "rows_metadata.csv")
brainspan_expr <- file.path(brainspan_dir, "expression_matrix.csv")

# -----------------------------
# 0. Check files exist
# -----------------------------
inputs <- data.table(
  label = c("sfari", "gse102741", "gse64018", "brainspan_cols", "brainspan_rows", "brainspan_expr"),
  path = c(sfari_file, gse102741_file, gse64018_file, brainspan_cols, brainspan_rows, brainspan_expr)
)
inputs[, exists := file.exists(path)]
safe_fwrite(inputs, file.path(METADIR, "00_input_file_check.tsv"))

if (any(!inputs$exists)) {
  stop("Missing required input files: ", paste(inputs[exists == FALSE, label], collapse = ", "))
}

log_msg("All required input files found.")

# -----------------------------
# 1. SFARI gene set (headerless, one symbol per line)
# -----------------------------
log_msg("Reading SFARI file: ", sfari_file)
sfari_raw <- fread(sfari_file, header = FALSE, sep = "\t", data.table = TRUE)
setnames(sfari_raw, "V1", "gene_symbol")
sfari_raw[, gene_symbol := toupper(trimws(gene_symbol))]
sfari_raw <- sfari_raw[gene_symbol != "" & !is.na(gene_symbol)]
sfari_std <- unique(sfari_raw)
setorder(sfari_std, gene_symbol)

log_msg("SFARI rows read: ", nrow(sfari_raw))
log_msg("Unique SFARI genes: ", nrow(sfari_std))

safe_fwrite(sfari_std, file.path(TABLEDIR, "01_SFARI_primary_Sand1_standardized.tsv"))
saveRDS(sfari_std, file.path(RDSDIR, "01_SFARI_primary_Sand1_standardized.rds"))

sfari_summary <- data.table(
  metric = c("n_rows_raw", "n_unique_genes"),
  value = c(nrow(sfari_raw), nrow(sfari_std))
)
safe_fwrite(sfari_summary, file.path(METADIR, "01_SFARI_summary.tsv"))

# -----------------------------
# 2. GSE102741 summary
# -----------------------------
log_msg("Inspecting GSE102741: ", gse102741_file)

gse102741_head <- fread(cmd = paste("zcat", shQuote(gse102741_file), "| head -5"), data.table = TRUE)
gse102741_dims <- fread(cmd = paste("zcat", shQuote(gse102741_file), "| awk 'END{print NR}'"), header = FALSE)

n_rows_total_102741 <- as.integer(gse102741_dims$V1) - 1L
n_cols_102741 <- ncol(gse102741_head)
n_samples_102741 <- n_cols_102741 - 1L

safe_fwrite(as.data.table(gse102741_head), file.path(METADIR, "02_GSE102741_head.tsv"))

# probe first 1000 rows for gene id properties
probe_102741 <- fread(cmd = paste("zcat", shQuote(gse102741_file), "| head -1001"), data.table = TRUE)
geneid_col_102741 <- colnames(probe_102741)[1]
probe_geneids_102741 <- probe_102741[[1]]
all_numeric_probe_102741 <- all(grepl("^[0-9]+$", probe_geneids_102741))

summary_102741 <- data.table(
  metric = c(
    "n_rows_total_excluding_header",
    "n_columns_total",
    "n_samples",
    "first_column_name",
    "first_column_numeric_probe_1000"
  ),
  value = c(
    n_rows_total_102741,
    n_cols_102741,
    n_samples_102741,
    geneid_col_102741,
    as.character(all_numeric_probe_102741)
  )
)
safe_fwrite(summary_102741, file.path(METADIR, "02_GSE102741_summary.tsv"))

# -----------------------------
# 3. GSE64018 summary
# -----------------------------
log_msg("Inspecting GSE64018: ", gse64018_file)

gse64018_head <- fread(cmd = paste("zcat", shQuote(gse64018_file), "| head -5"), data.table = TRUE)
gse64018_dims <- fread(cmd = paste("zcat", shQuote(gse64018_file), "| awk 'END{print NR}'"), header = FALSE)

n_rows_total_64018 <- as.integer(gse64018_dims$V1) - 1L
n_cols_64018 <- ncol(gse64018_head)
n_samples_64018 <- n_cols_64018 - 1L

safe_fwrite(as.data.table(gse64018_head), file.path(METADIR, "03_GSE64018_head.tsv"))

summary_64018 <- data.table(
  metric = c(
    "n_rows_total_excluding_header",
    "n_columns_total",
    "n_samples",
    "first_column_name"
  ),
  value = c(
    n_rows_total_64018,
    n_cols_64018,
    n_samples_64018,
    colnames(gse64018_head)[1]
  )
)
safe_fwrite(summary_64018, file.path(METADIR, "03_GSE64018_summary.tsv"))

# -----------------------------
# 4. BrainSpan inventory
# -----------------------------
log_msg("Reading BrainSpan files.")

brainspan_files <- data.table(
  file = c(brainspan_cols, brainspan_rows, brainspan_expr),
  basename = basename(c(brainspan_cols, brainspan_rows, brainspan_expr)),
  size_bytes = unname(file.info(c(brainspan_cols, brainspan_rows, brainspan_expr))$size)
)
safe_fwrite(brainspan_files, file.path(METADIR, "07_BrainSpan_inventory.tsv"))

cols_meta <- fread(brainspan_cols, data.table = TRUE)
rows_meta <- fread(brainspan_rows, data.table = TRUE)
expr_header <- fread(cmd = paste("head -2", shQuote(brainspan_expr)), data.table = TRUE)

# dictionary of BrainSpan columns/rows metadata
col_dict <- data.table(
  table = "columns_metadata",
  column_name = colnames(cols_meta)
)
row_dict <- data.table(
  table = "rows_metadata",
  column_name = colnames(rows_meta)
)
brainspan_dict <- rbindlist(list(col_dict, row_dict), use.names = TRUE)
safe_fwrite(brainspan_dict, file.path(METADIR, "08_BrainSpan_column_dictionary.tsv"))

# identify candidate mapping columns in BrainSpan row metadata
mapping_patterns <- c(
  "gene", "symbol", "ensembl", "entrez", "ncbi", "id"
)
map_candidates <- rows_meta[, lapply(.SD, function(x) class(x)[1]), .SDcols = colnames(rows_meta)]
map_candidates <- melt(map_candidates, measure.vars = colnames(rows_meta), variable.name = "column_name", value.name = "class")
map_candidates[, pattern_hit := vapply(column_name, function(z) any(grepl(paste(mapping_patterns, collapse = "|"), z, ignore.case = TRUE)), logical(1))]
map_candidates <- map_candidates[pattern_hit == TRUE]
safe_fwrite(map_candidates, file.path(METADIR, "09_BrainSpan_mapping_candidate_columns.tsv"))

# quick BrainSpan run summary
n_brainspan_samples <- nrow(cols_meta)
n_brainspan_genes <- nrow(rows_meta)
n_expr_columns <- ncol(expr_header)

run_summary <- rbindlist(list(
  data.table(section = "SFARI", metric = "n_unique_genes", value = as.character(nrow(sfari_std))),
  data.table(section = "GSE102741", metric = "n_samples", value = as.character(n_samples_102741)),
  data.table(section = "GSE102741", metric = "gene_id_style_numeric_probe_1000", value = as.character(all_numeric_probe_102741)),
  data.table(section = "GSE64018", metric = "n_samples", value = as.character(n_samples_64018)),
  data.table(section = "BrainSpan", metric = "n_column_metadata_rows", value = as.character(n_brainspan_samples)),
  data.table(section = "BrainSpan", metric = "n_row_metadata_rows", value = as.character(n_brainspan_genes)),
  data.table(section = "BrainSpan", metric = "n_expression_header_columns", value = as.character(n_expr_columns))
), use.names = TRUE)

safe_fwrite(run_summary, file.path(METADIR, "10_Step01A_run_summary.tsv"))

saveRDS(list(
  sfari = sfari_std,
  gse102741_head = gse102741_head,
  gse64018_head = gse64018_head,
  brainspan_columns_metadata = cols_meta[1:min(10, .N)],
  brainspan_rows_metadata = rows_meta[1:min(10, .N)]
), file.path(RDSDIR, "Step01A_preview_objects.rds"))

log_msg("Step01A completed successfully.")
