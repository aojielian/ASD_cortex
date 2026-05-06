#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(nnls)
})

parse_args <- function(x){
  out <- list(
    base = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    psych_dir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024",
    gandal_rdata = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData",
    markers_per_class = 150L,
    min_cells_per_group = 20L,
    drop_fraction_class = "END",
    outdir = NA_character_
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    val <- if (i < length(x)) x[[i+1L]] else NA_character_
    if (key == "--base") out$base <- val
    if (key == "--psych_dir") out$psych_dir <- val
    if (key == "--gandal_rdata") out$gandal_rdata <- val
    if (key == "--markers_per_class") out$markers_per_class <- as.integer(val)
    if (key == "--min_cells_per_group") out$min_cells_per_group <- as.integer(val)
    if (key == "--drop_fraction_class") out$drop_fraction_class <- val
    if (key == "--outdir") out$outdir <- val
    i <- i + 2L
  }
  if (is.na(out$outdir)) out$outdir <- file.path(out$base, "reviewer_round3", "MA_round1B_prepost_CI")
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))

dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "meta"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "MA_round1B_prepost_CI.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
safe_fwrite <- function(x, file) fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x <- sub("\\..*$", "", x)
  x[nchar(x) == 0L] <- NA_character_
  toupper(x)
}
normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"+|"+$', '', x)
  x <- gsub("^'+|'+$", "", x)
  x <- sub("^(ASD|CTL)_", "", x)
  gsub("\\s+", "", x)
}
extract_attr <- function(x, key) {
  m <- regexec(paste0(key, ' "([^"]+)"'), x)
  regmatches(x, m) |> lapply(function(z) if (length(z) >= 2) z[2] else NA_character_) |> unlist(use.names = FALSE)
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
compute_logcpm <- function(count_dt) {
  mat <- as.matrix(count_dt[, -1]); rownames(mat) <- count_dt[[1]]; storage.mode(mat) <- "double"
  lib <- colSums(mat); lib[lib <= 0] <- NA_real_
  log2(t(t(mat) / lib) * 1e6 + 1)
}
fit_ci <- function(dt, program_name, cohort, model_type, fraction_covars = character()) {
  dt <- copy(dt)
  dt <- dt[!is.na(score) & !is.na(diagnosis)]
  dt[, diagnosis := factor(as.character(diagnosis), levels = c("Control","ASD"))]
  n_asd <- sum(dt$diagnosis == "ASD")
  n_control <- sum(dt$diagnosis == "Control")
  if (n_asd < 3 || n_control < 3) return(NULL)
  rhs <- c("diagnosis", fraction_covars)
  rhs <- rhs[rhs %in% names(dt)]
  cc <- dt[complete.cases(dt[, ..rhs], score)]
  if (sum(cc$diagnosis == "ASD") < 3 || sum(cc$diagnosis == "Control") < 3) return(NULL)
  fml <- as.formula(paste("score ~", paste(rhs, collapse = " + ")))
  fit <- lm(fml, data = cc)
  sm <- summary(fit)$coefficients
  if (!("diagnosisASD" %in% rownames(sm))) return(NULL)
  beta <- sm["diagnosisASD","Estimate"]; se <- sm["diagnosisASD","Std. Error"]; p <- sm["diagnosisASD","Pr(>|t|)"]
  data.table(
    cohort = cohort, program_name = program_name, model_type = model_type,
    n_samples = nrow(cc), n_asd = sum(cc$diagnosis == "ASD"), n_control = sum(cc$diagnosis == "Control"),
    mean_asd = mean(cc$score[cc$diagnosis == "ASD"]), mean_control = mean(cc$score[cc$diagnosis == "Control"]),
    beta_asd = beta, se = se, ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se,
    p_value = p, formula_used = paste(deparse(fml), collapse = "")
  )
}

prog_file <- file.path(args$base, "tables", "20_midPrenatal_core_programs.tsv")
prog <- fread(prog_file)
prog[, gene_symbol := clean_symbol(gene_symbol)]
prog <- unique(prog[!is.na(program_name) & !is.na(gene_symbol), .(program_name, gene_symbol)])
prog_list <- split(prog$gene_symbol, prog$program_name)
sfari_all <- unique(prog_list[["SFARI_all"]])
top20 <- unique(prog_list[["midPrenatal_SFARI_top20"]])
if (length(sfari_all) == 0 || length(top20) == 0) stop("Could not extract SFARI_all or midPrenatal_SFARI_top20 from core programs")
program_sets <- list(SFARI_all = sfari_all, midPrenatal_SFARI_top20 = top20)

meta_file <- file.path(args$psych_dir, "meta.tsv")
feature_file <- file.path(args$psych_dir, "counts_features.tsv.gz")
matrix_file <- file.path(args$psych_dir, "counts_matrix.mtx.gz")
for (f in c(meta_file, feature_file, matrix_file)) if (!file.exists(f)) stop("Missing required PsychENCODE file: ", f)

log_msg("Reading PsychENCODE reference for NNLS basis")
meta <- fread(meta_file)
stopifnot("Cell_ID" %in% colnames(meta), "annotation" %in% colnames(meta))
broad_levels <- c("EXN","INN","AST","OPC","ODC","MG","END")
meta <- meta[annotation %in% broad_levels]
feat <- fread(cmd = paste("zcat", shQuote(feature_file)), header = FALSE)
gene_symbol <- clean_symbol(if (ncol(feat) >= 2L) feat[[2L]] else feat[[1L]])
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
mtx_gene <- collapse_sparse_by_symbol(mtx, rownames(mtx))
basis_counts <- aggregate_sparse_by_group(mtx_gene, as.character(meta$annotation))
basis_counts <- as.matrix(basis_counts)
lib_basis <- colSums(basis_counts); lib_basis[lib_basis <= 0] <- 1
basis_cpm <- sweep(basis_counts, 2, lib_basis, "/", check.margin = FALSE) * 1e6
basis_logcpm <- safe_log2p1(basis_cpm)

marker_list <- list()
for (cls in broad_levels) {
  this <- basis_logcpm[, cls]
  others <- rowMeans(basis_logcpm[, setdiff(colnames(basis_logcpm), cls), drop = FALSE])
  dt <- data.table(gene_symbol = rownames(basis_logcpm), broad_class = cls, logFC = this - others, avg_expr = this)
  dt <- dt[is.finite(logFC) & logFC > 0]
  setorder(dt, -logFC, -avg_expr)
  if (nrow(dt) > args$markers_per_class) dt <- dt[seq_len(args$markers_per_class)]
  marker_list[[cls]] <- dt
}
marker_dt <- unique(rbindlist(marker_list))
basis <- basis_logcpm[unique(marker_dt$gene_symbol), broad_levels, drop = FALSE]
safe_fwrite(marker_dt, file.path(args$outdir, "tables", "MA_round1B_reference_markers.tsv"))

load_gse102741 <- function(base) {
  rds_file <- file.path(base, "rds", "22_GSE102741_geneSymbol_log2cpm.rds")
  meta_file <- file.path(base, "tables", "21_GSE102741_sample_metadata.tsv")
  obj <- readRDS(rds_file)
  expr <- as.matrix(obj$log2cpm); rownames(expr) <- clean_symbol(rownames(expr))
  meta <- fread(meta_file)
  sample_col <- intersect(c("sample_id","sample","Sample.ID"), names(meta))[1]
  dx_col <- intersect(c("diagnosis","Diagnosis","dx"), names(meta))[1]
  meta <- meta[, .(sample = as.character(get(sample_col)), diagnosis = as.character(get(dx_col)))]
  meta[, diagnosis := fifelse(diagnosis %in% c("ASD","Autism"), "ASD", fifelse(diagnosis %in% c("Control","CTL","CTRL"), "Control", diagnosis))]
  keep <- intersect(colnames(expr), meta$sample)
  expr <- expr[, keep, drop = FALSE]
  meta <- meta[match(keep, sample)]
  list(expr = expr, meta = meta)
}
load_gse64018 <- function(base) {
  counts_file <- file.path(base, "tables", "122_GSE64018_collapsed_counts_by_symbol.tsv.gz")
  map_file <- file.path(base, "meta", "120_GSE64018_explicit_sample_mapping.tsv")
  counts <- fread(cmd = paste("zcat", shQuote(counts_file)))
  expr <- compute_logcpm(counts)
  sm <- fread(map_file)
  norm_col <- intersect(c("file_sample_norm","sample","sample_norm"), names(sm))[1]
  dx_col <- intersect(c("diagnosis","group","dx"), names(sm))[1]
  title_col <- intersect(c("geo_title","title"), names(sm))[1]
  if (is.na(dx_col)) {
    if (!is.na(title_col)) { sm[, diagnosis := fifelse(startsWith(get(title_col), "ASD_"), "ASD", "Control")]; dx_col <- "diagnosis" } else stop("Could not detect diagnosis col for GSE64018")
  }
  sm[, sample := normalize_sample_id(get(norm_col))]
  sm[, diagnosis := fifelse(as.character(get(dx_col)) %in% c("ASD","Autism"), "ASD", fifelse(as.character(get(dx_col)) %in% c("Control","CTL"), "Control", as.character(get(dx_col))))]
  sm <- unique(sm[, .(sample, diagnosis)])
  meta <- data.table(sample = normalize_sample_id(colnames(expr)))
  meta <- merge(meta, sm, by = "sample", all.x = TRUE, sort = FALSE)
  colnames(expr) <- meta$sample
  keep <- !is.na(meta$diagnosis)
  list(expr = expr[, keep, drop = FALSE], meta = meta[keep])
}
load_gandal <- function(rdata_path, gtf = "/gpfs/hpc/home/lijc/lianaoj/reference/gencode.v49.basic.annotation.gtf.gz") {
  gtf_dt <- fread(cmd = paste("zcat", shQuote(gtf), "| awk '$3==\"gene\"'"), sep = "\t", header = FALSE)
  gtf_map <- unique(data.table(
    ensembl_gene_id = clean_symbol(extract_attr(gtf_dt$V9, "gene_id")),
    gene_symbol = clean_symbol(extract_attr(gtf_dt$V9, "gene_name"))
  )[ensembl_gene_id != "" & gene_symbol != ""])
  e <- new.env(parent = emptyenv()); load(rdata_path, envir = e)
  datExpr <- get("datExpr", envir = e); datMeta <- as.data.table(get("datMeta", envir = e))
  ex <- as.matrix(datExpr)
  sample_col <- intersect(c("sample_id","Sample.ID","sample","sampleID"), names(datMeta))[1]
  dx_col <- intersect(c("Diagnosis","diagnosis","dx"), names(datMeta))[1]
  sample_ids <- as.character(datMeta[[sample_col]])
  sample_cols <- intersect(colnames(ex), sample_ids)
  if (length(sample_cols) >= 10L) {
    expr <- ex[, sample_cols, drop = FALSE]; feats <- rownames(expr); meta <- datMeta[match(sample_cols, sample_ids)]; expr_samples <- sample_cols
  } else {
    sample_rows <- intersect(rownames(ex), sample_ids); expr <- t(ex[sample_rows, , drop = FALSE]); feats <- rownames(expr); meta <- datMeta[match(colnames(expr), sample_ids)]; expr_samples <- colnames(expr)
  }
  feats <- clean_symbol(feats)
  sym <- if (mean(grepl("^ENSG", feats)) > 0.5) gtf_map$gene_symbol[match(feats, gtf_map$ensembl_gene_id)] else feats
  keep <- !is.na(sym) & sym != ""
  expr <- expr[keep, , drop = FALSE]
  rownames(expr) <- sym[keep]
  expr <- rowsum(expr, rownames(expr))
  dx_raw <- as.character(meta[[dx_col]])
  diagnosis <- fifelse(dx_raw %in% c("ASD","Autism"), "ASD", fifelse(dx_raw %in% c("Control","CTL","CTRL"), "Control", NA_character_))
  keep2 <- !is.na(diagnosis)
  list(expr = expr[, keep2, drop = FALSE], meta = data.table(sample = expr_samples[keep2], diagnosis = diagnosis[keep2]))
}

deconvolve_cohort <- function(expr, meta, cohort_name) {
  common_genes <- intersect(rownames(basis), rownames(expr))
  if (length(common_genes) < 100) stop("Too few common genes for cohort: ", cohort_name)
  B <- as.matrix(basis[common_genes, broad_levels, drop = FALSE])
  X <- as.matrix(expr[common_genes, , drop = FALSE])
  frac_mat <- matrix(NA_real_, nrow = ncol(X), ncol = ncol(B), dimnames = list(colnames(X), colnames(B)))
  for (j in seq_len(ncol(X))) {
    fit <- nnls::nnls(B, X[, j])
    est <- coef(fit)
    est[!is.finite(est) | est < 0] <- 0
    if (sum(est) <= 0) est[] <- 1 / length(est) else est <- est / sum(est)
    frac_mat[j, ] <- est
  }
  frac_dt <- data.table(sample = rownames(frac_mat), as.data.table(frac_mat))
  frac_dt[, cohort := cohort_name]
  safe_fwrite(frac_dt, file.path(args$outdir, "tables", paste0("MA_round1B_", cohort_name, "_nnls_fractions.tsv")))
  model_dt <- merge(meta, frac_dt, by = "sample")
  fraction_covars <- setdiff(colnames(frac_dt), c("sample","cohort", args$drop_fraction_class))
  res_list <- list()
  for (pn in names(program_sets)) {
    genes <- intersect(program_sets[[pn]], rownames(expr))
    score <- colMeans(expr[genes, , drop = FALSE], na.rm = TRUE)
    tmp <- copy(model_dt)
    tmp[, score := score[match(sample, colnames(expr))]]
    unadj <- fit_ci(tmp, pn, cohort_name, "unadjusted", fraction_covars = character())
    adj <- fit_ci(tmp, pn, cohort_name, "nnls_adjusted", fraction_covars = fraction_covars)
    if (!is.null(unadj)) res_list[[paste(pn, "unadj", sep = "__")]] <- unadj
    if (!is.null(adj)) res_list[[paste(pn, "adj", sep = "__")]] <- adj
  }
  rbindlist(res_list, fill = TRUE)
}

log_msg("Running cohort deconvolution and pre/post CI models")
cohorts <- list(
  GSE102741 = load_gse102741(args$base),
  GSE64018 = load_gse64018(args$base),
  Gandal2022 = load_gandal(args$gandal_rdata)
)
res <- rbindlist(lapply(names(cohorts), function(nm) deconvolve_cohort(cohorts[[nm]]$expr, cohorts[[nm]]$meta, nm)), fill = TRUE)
if (nrow(res) == 0) stop("No pre/post model results generated")

wide <- dcast(
  res[, .(cohort, program_name, model_type, beta_asd, se, ci_low, ci_high, p_value, mean_asd, mean_control, n_samples, n_asd, n_control)],
  cohort + program_name + n_samples + n_asd + n_control ~ model_type,
  value.var = c("beta_asd","se","ci_low","ci_high","p_value","mean_asd","mean_control")
)
wide[, attenuation_adjusted_minus_unadjusted := beta_asd_nnls_adjusted - beta_asd_unadjusted]
wide[, direction_retained := sign(beta_asd_nnls_adjusted) == sign(beta_asd_unadjusted)]
setorder(wide, program_name, cohort)

safe_fwrite(res, file.path(args$outdir, "tables", "MA_round1B_prepost_model_results_long.tsv"))
safe_fwrite(wide, file.path(args$outdir, "tables", "MA_round1B_prepost_model_results_wide.tsv"))

run_summary <- data.table(
  item = c("markers_per_class", "min_cells_per_group", "basis_genes", "drop_fraction_class", "n_result_rows"),
  value = c(args$markers_per_class, args$min_cells_per_group, nrow(basis), args$drop_fraction_class, nrow(res))
)
safe_fwrite(run_summary, file.path(args$outdir, "meta", "MA_round1B_run_summary.tsv"))
log_msg("Completed Round1B pre/post CI export.")
