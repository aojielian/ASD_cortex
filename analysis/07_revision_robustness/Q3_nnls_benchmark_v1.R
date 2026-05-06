#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(nnls)
})

parse_args <- function(x) {
  out <- list(
    base = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    psych_dir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024",
    markers_per_class = 150L,
    min_cells_per_group = 20L,
    n_sims = 500L,
    seed = 20260331L,
    outdir = NA_character_
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    val <- if (i < length(x)) x[[i+1L]] else NA_character_
    if (key == "--base") out$base <- val
    if (key == "--psych_dir") out$psych_dir <- val
    if (key == "--markers_per_class") out$markers_per_class <- as.integer(val)
    if (key == "--min_cells_per_group") out$min_cells_per_group <- as.integer(val)
    if (key == "--n_sims") out$n_sims <- as.integer(val)
    if (key == "--seed") out$seed <- as.integer(val)
    if (key == "--outdir") out$outdir <- val
    i <- i + 2L
  }
  if (is.na(out$outdir)) out$outdir <- file.path(out$base, "reviewer_round2", "Q3_nnls_benchmark")
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))

dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "meta"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "Q3_nnls_benchmark.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
safe_fwrite <- function(x, file) fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
read_dt_no_header <- function(path) fread(path, header = FALSE)
clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}
collapse_sparse_by_symbol <- function(mtx, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mtx <- mtx[keep, , drop = FALSE]
  symbols <- symbols[keep]
  grp <- factor(symbols, levels = unique(symbols))
  P <- sparseMatrix(i = as.integer(grp), j = seq_along(grp), x = 1, dims = c(nlevels(grp), length(grp)))
  res <- P %*% mtx
  rownames(res) <- levels(grp)
  res
}
aggregate_sparse_by_group <- function(mtx, groups) {
  keep <- !is.na(groups) & groups != ""
  mtx <- mtx[, keep, drop = FALSE]
  groups <- groups[keep]
  grp <- factor(groups, levels = unique(groups))
  G <- sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1, dims = c(length(grp), nlevels(grp)))
  res <- mtx %*% G
  colnames(res) <- levels(grp)
  res
}
safe_log2p1 <- function(mat) log2(mat + 1)
dirichlet1 <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}
rmse <- function(x,y) sqrt(mean((x-y)^2, na.rm = TRUE))
mae <- function(x,y) mean(abs(x-y), na.rm = TRUE)

set.seed(args$seed)

meta_file <- file.path(args$psych_dir, "meta.tsv")
feature_file <- file.path(args$psych_dir, "counts_features.tsv.gz")
barcode_file <- file.path(args$psych_dir, "counts_barcodes.tsv.gz")
matrix_file <- file.path(args$psych_dir, "counts_matrix.mtx.gz")
for (f in c(meta_file, feature_file, barcode_file, matrix_file)) if (!file.exists(f)) stop("Missing required file: ", f)

log_msg("Reading PsychENCODE meta/features/barcodes/matrix")
meta <- fread(meta_file)
stopifnot("Cell_ID" %in% colnames(meta), "annotation" %in% colnames(meta))
broad_levels <- c("EXN", "INN", "AST", "OPC", "ODC", "MG", "END")
meta <- meta[annotation %in% broad_levels]

donor_candidates <- c("individual_ID","individual","donor_id","donor","sample_id")
donor_col <- donor_candidates[donor_candidates %in% colnames(meta)][1]
if (is.na(donor_col)) stop("Could not identify donor column in PsychENCODE meta")

feat <- read_dt_no_header(feature_file)
bc <- read_dt_no_header(barcode_file)
gene_symbol <- if (ncol(feat) >= 2L) feat[[2L]] else feat[[1L]]
gene_symbol <- clean_symbol(gene_symbol)

mtx <- readMM(matrix_file)
mtx <- as(mtx, "CsparseMatrix")

nrow_target <- min(nrow(mtx), length(gene_symbol))
if (nrow_target < nrow(mtx) || nrow_target < length(gene_symbol)) {
  mtx <- mtx[seq_len(nrow_target), , drop = FALSE]
  gene_symbol <- gene_symbol[seq_len(nrow_target)]
}
ncol_target <- min(ncol(mtx), nrow(meta))
if (ncol_target < ncol(mtx) || ncol_target < nrow(meta)) {
  mtx <- mtx[, seq_len(ncol_target), drop = FALSE]
  meta <- meta[seq_len(ncol_target)]
}
rownames(mtx) <- gene_symbol
colnames(mtx) <- meta$Cell_ID

log_msg("Collapsing duplicated genes and building class basis")
mtx_gene <- collapse_sparse_by_symbol(mtx, rownames(mtx))
meta_use <- copy(meta)
meta_use$annotation <- factor(meta_use$annotation, levels = broad_levels)

basis_counts <- aggregate_sparse_by_group(mtx_gene, as.character(meta_use$annotation))
basis_counts <- as.matrix(basis_counts)
lib_sizes <- colSums(basis_counts); lib_sizes[lib_sizes <= 0] <- 1
basis_cpm <- sweep(basis_counts, 2, lib_sizes, "/", check.margin = FALSE) * 1e6
basis_logcpm <- safe_log2p1(basis_cpm)

marker_list <- list()
for (cls in broad_levels) {
  this <- basis_logcpm[, cls]
  other_cols <- setdiff(colnames(basis_logcpm), cls)
  others <- rowMeans(basis_logcpm[, other_cols, drop = FALSE])
  dt <- data.table(
    gene_symbol = rownames(basis_logcpm),
    broad_class = cls,
    logFC = as.numeric(this - others),
    avg_expr = as.numeric(this),
    avg_others = as.numeric(others)
  )
  dt <- dt[is.finite(logFC) & !is.na(gene_symbol) & gene_symbol != ""]
  setorder(dt, -logFC, -avg_expr)
  dt <- dt[logFC > 0]
  if (nrow(dt) > args$markers_per_class) dt <- dt[seq_len(args$markers_per_class)]
  marker_list[[cls]] <- dt
}
marker_dt <- rbindlist(marker_list, use.names = TRUE, fill = TRUE)
if (nrow(marker_dt) == 0L) stop("No valid reference markers identified from PsychENCODE.")
basis <- basis_logcpm[unique(marker_dt$gene_symbol), broad_levels, drop = FALSE]
safe_fwrite(marker_dt, file.path(args$outdir, "tables", "Q3_reference_marker_summary.tsv"))

# donor-class pseudobulks
meta_use[, donor_id := as.character(get(donor_col))]
meta_use <- meta_use[!is.na(donor_id) & donor_id != ""]
meta_use[, donor_class_id := paste(donor_id, annotation, sep = "||")]
group_counts <- meta_use[, .N, by = donor_class_id]
valid_groups <- group_counts[N >= args$min_cells_per_group, donor_class_id]
meta_pb <- meta_use[donor_class_id %in% valid_groups, .(Cell_ID, donor_id, broad_class = annotation, donor_class_id)]
if (nrow(meta_pb) == 0) stop("No donor-class groups pass min_cells_per_group")

agg_counts <- aggregate_sparse_by_group(mtx_gene[, meta_pb$Cell_ID, drop = FALSE], meta_pb$donor_class_id)
pb_meta <- unique(meta_pb[, .(donor_class_id, donor_id, broad_class)])
pb_meta <- pb_meta[match(colnames(agg_counts), donor_class_id)]
lib_pb <- colSums(agg_counts); lib_pb[lib_pb <= 0] <- 1
pb_cpm <- sweep(as.matrix(agg_counts), 2, lib_pb, "/", check.margin = FALSE) * 1e6
pb_logcpm <- safe_log2p1(pb_cpm)
pb_logcpm <- pb_logcpm[rownames(pb_logcpm) %in% rownames(basis), , drop = FALSE]
basis <- basis[rownames(pb_logcpm), , drop = FALSE]

pool_dt <- pb_meta[, .N, by = broad_class][order(broad_class)]
safe_fwrite(pool_dt, file.path(args$outdir, "meta", "Q3_donorClass_pool_summary.tsv"))

# simulate mixtures
log_msg("Running synthetic mixture benchmark with n_sims=", args$n_sims)
mix_rows <- vector("list", args$n_sims)
truth_long <- vector("list", args$n_sims)
est_long <- vector("list", args$n_sims)
for (i in seq_len(args$n_sims)) {
  truth <- dirichlet1(rep(1, length(broad_levels)))
  names(truth) <- broad_levels
  picked <- pb_meta[, .SD[sample(.N, 1)], by = broad_class]
  picked_ids <- picked$donor_class_id
  target <- pb_logcpm[, picked_ids, drop = FALSE] %*% matrix(truth[match(picked$broad_class, broad_levels)], ncol = 1)
  fit <- nnls::nnls(as.matrix(basis), as.numeric(target))
  est <- coef(fit)
  est[!is.finite(est) | est < 0] <- 0
  if (sum(est) <= 0) est[] <- 1 / length(est) else est <- est / sum(est)
  names(est) <- colnames(basis)

  mix_rows[[i]] <- data.table(sim_id = i, donor_class_ids = paste(picked_ids, collapse = ";"), rss = fit$deviance)
  truth_long[[i]] <- data.table(sim_id = i, broad_class = names(truth), truth_fraction = as.numeric(truth))
  est_long[[i]] <- data.table(sim_id = i, broad_class = names(est), est_fraction = as.numeric(est))
}
mix_meta <- rbindlist(mix_rows)
truth_dt <- rbindlist(truth_long)
est_dt <- rbindlist(est_long)
bench <- merge(truth_dt, est_dt, by = c("sim_id","broad_class"))
bench <- merge(bench, mix_meta, by = "sim_id", all.x = TRUE)
bench[, abs_error := abs(est_fraction - truth_fraction)]
bench[, sq_error := (est_fraction - truth_fraction)^2]

by_ct <- bench[, .(
  n = .N,
  pearson_r = suppressWarnings(cor(truth_fraction, est_fraction, method = "pearson")),
  spearman_r = suppressWarnings(cor(truth_fraction, est_fraction, method = "spearman")),
  MAE = mae(est_fraction, truth_fraction),
  RMSE = rmse(est_fraction, truth_fraction),
  mean_truth = mean(truth_fraction),
  mean_est = mean(est_fraction)
), by = broad_class]

overall <- data.table(
  n_sims = args$n_sims,
  n_broad_classes = length(broad_levels),
  pearson_r_all = suppressWarnings(cor(bench$truth_fraction, bench$est_fraction, method = "pearson")),
  spearman_r_all = suppressWarnings(cor(bench$truth_fraction, bench$est_fraction, method = "spearman")),
  MAE_all = mae(bench$est_fraction, bench$truth_fraction),
  RMSE_all = rmse(bench$est_fraction, bench$truth_fraction)
)

safe_fwrite(bench, file.path(args$outdir, "tables", "Q3_nnls_benchmark_per_mixture.tsv.gz"))
safe_fwrite(by_ct, file.path(args$outdir, "tables", "Q3_nnls_benchmark_summary_by_celltype.tsv"))
safe_fwrite(overall, file.path(args$outdir, "tables", "Q3_nnls_benchmark_overall.tsv"))

run_summary <- data.table(
  item = c("n_cells_meta", "n_genes_basis", "n_marker_rows", "n_valid_donor_class_groups", "n_sims", "markers_per_class", "min_cells_per_group"),
  value = c(nrow(meta_use), nrow(basis), nrow(marker_dt), nrow(pb_meta), args$n_sims, args$markers_per_class, args$min_cells_per_group)
)
safe_fwrite(run_summary, file.path(args$outdir, "meta", "Q3_run_summary.tsv"))
log_msg("Completed Q3 benchmark.")
