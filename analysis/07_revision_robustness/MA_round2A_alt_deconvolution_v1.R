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
    nnls_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/reviewer_round3/MA_round1B_prepost_CI",
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
    if (key == "--nnls_dir") out$nnls_dir <- val
    if (key == "--drop_fraction_class") out$drop_fraction_class <- val
    if (key == "--outdir") out$outdir <- val
    i <- i + 2L
  }
  if (is.na(out$outdir)) out$outdir <- file.path(out$base, "reviewer_round4", "MA_round2A_alt_deconvolution")
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "meta"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "MA_round2A_alt_deconvolution.log")

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
    beta_asd = beta, se = se, ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se, p_value = p
  )
}
deconv_positive_ols <- function(B, y) {
  fit <- lm.fit(x = cbind(1, B), y = y)
  coef_raw <- fit$coefficients[-1]
  coef_raw[is.na(coef_raw)] <- 0
  coef_raw[coef_raw < 0] <- 0
  if (sum(coef_raw) <= 0) coef_raw[] <- 1 / length(coef_raw) else coef_raw <- coef_raw / sum(coef_raw)
  coef_raw
}
cor_safe <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

# core programs
prog <- fread(file.path(args$base, "tables", "20_midPrenatal_core_programs.tsv"))
prog[, gene_symbol := clean_symbol(gene_symbol)]
prog <- unique(prog[!is.na(program_name) & !is.na(gene_symbol), .(program_name, gene_symbol)])
prog_list <- split(prog$gene_symbol, prog$program_name)
program_sets <- list(
  SFARI_all = unique(prog_list[["SFARI_all"]]),
  midPrenatal_SFARI_top20 = unique(prog_list[["midPrenatal_SFARI_top20"]])
)

# reference basis from PsychENCODE
meta_file <- file.path(args$psych_dir, "meta.tsv")
feature_file <- file.path(args$psych_dir, "counts_features.tsv.gz")
matrix_file <- file.path(args$psych_dir, "counts_matrix.mtx.gz")
for (f in c(meta_file, feature_file, matrix_file)) if (!file.exists(f)) stop("Missing required PsychENCODE file: ", f)

log_msg("Building alternative deconvolution basis from PsychENCODE")
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
safe_fwrite(marker_dt, file.path(args$outdir, "tables", "MA_round2A_reference_markers.tsv"))

# bulk loaders
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
  expr <- expr[, keep, drop = FALSE]; meta <- meta[match(keep, sample)]
  list(expr = expr, meta = meta)
}
load_gse64018 <- function(base) {
  counts_file <- file.path(base, "tables", "122_GSE64018_collapsed_counts_by_symbol.tsv.gz")
  map_file <- file.path(base, "meta", "120_GSE64018_explicit_sample_mapping.tsv")
  counts <- fread(cmd = paste("zcat", shQuote(counts_file))); expr <- compute_logcpm(counts)
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
  expr <- expr[keep, , drop = FALSE]; rownames(expr) <- sym[keep]; expr <- rowsum(expr, rownames(expr))
  dx_raw <- as.character(meta[[dx_col]])
  diagnosis <- fifelse(dx_raw %in% c("ASD","Autism"), "ASD", fifelse(dx_raw %in% c("Control","CTL","CTRL"), "Control", NA_character_))
  keep2 <- !is.na(diagnosis)
  list(expr = expr[, keep2, drop = FALSE], meta = data.table(sample = expr_samples[keep2], diagnosis = diagnosis[keep2]))
}
cohorts <- list(
  GSE102741 = load_gse102741(args$base),
  GSE64018 = load_gse64018(args$base),
  Gandal2022 = load_gandal(args$gandal_rdata)
)

nnls_frac_paths <- c(
  GSE102741 = file.path(args$nnls_dir, "tables", "MA_round1B_GSE102741_nnls_fractions.tsv"),
  GSE64018 = file.path(args$nnls_dir, "tables", "MA_round1B_GSE64018_nnls_fractions.tsv"),
  Gandal2022 = file.path(args$nnls_dir, "tables", "MA_round1B_Gandal2022_nnls_fractions.tsv")
)
for (f in nnls_frac_paths) if (!file.exists(f)) stop("Missing NNLS fraction file: ", f)
nnls_fracs <- lapply(names(nnls_frac_paths), function(nm) fread(nnls_frac_paths[[nm]]))
names(nnls_fracs) <- names(nnls_frac_paths)

# alt deconvolution + concordance + adjusted effects
frac_conc <- list()
model_res <- list()
method_inventory <- data.table(method = c("NNLS_reference_linear_deconvolution","PositivePart_OLS_reference_linear_deconvolution"),
                               operates_on = c("harmonized bulk logCPM/log-expression", "harmonized bulk logCPM/log-expression"))
safe_fwrite(method_inventory, file.path(args$outdir, "meta", "MA_round2A_method_inventory.tsv"))

for (nm in names(cohorts)) {
  log_msg("Running alternative deconvolution for cohort ", nm)
  expr <- cohorts[[nm]]$expr
  meta <- cohorts[[nm]]$meta
  common_genes <- intersect(rownames(basis), rownames(expr))
  B <- as.matrix(basis[common_genes, broad_levels, drop = FALSE])
  X <- as.matrix(expr[common_genes, , drop = FALSE])
  alt_frac <- matrix(NA_real_, nrow = ncol(X), ncol = ncol(B), dimnames = list(colnames(X), colnames(B)))
  for (j in seq_len(ncol(X))) {
    alt_frac[j, ] <- deconv_positive_ols(B, X[, j])
  }
  alt_frac_dt <- data.table(sample = rownames(alt_frac), as.data.table(alt_frac))
  alt_frac_dt[, cohort := nm]
  safe_fwrite(alt_frac_dt, file.path(args$outdir, "tables", paste0("MA_round2A_", nm, "_altfractions.tsv")))

  nnls_dt <- nnls_fracs[[nm]]
  merged_frac <- merge(nnls_dt, alt_frac_dt, by = c("sample","cohort"), suffixes = c("_nnls","_alt"))
  per_ct <- rbindlist(lapply(broad_levels, function(ct) {
    data.table(
      cohort = nm, cell_type = ct,
      pearson_r = cor_safe(merged_frac[[paste0(ct, "_nnls")]], merged_frac[[paste0(ct, "_alt")]], "pearson"),
      spearman_r = cor_safe(merged_frac[[paste0(ct, "_nnls")]], merged_frac[[paste0(ct, "_alt")]], "spearman"),
      mean_abs_diff = mean(abs(merged_frac[[paste0(ct, "_nnls")]] - merged_frac[[paste0(ct, "_alt")]]), na.rm = TRUE)
    )
  }), fill = TRUE)
  frac_conc[[nm]] <- per_ct

  alt_covars <- setdiff(colnames(alt_frac_dt), c("sample","cohort", args$drop_fraction_class))
  model_dt <- merge(meta, alt_frac_dt, by = "sample")
  rr <- rbindlist(lapply(names(program_sets), function(pn) {
    genes <- intersect(program_sets[[pn]], rownames(expr))
    tmp <- copy(model_dt)
    tmp[, score := colMeans(expr[genes, tmp$sample, drop = FALSE], na.rm = TRUE)]
    fit_ci(tmp, pn, nm, "alt_adjusted", alt_covars)
  }), fill = TRUE)
  model_res[[nm]] <- rr
}

frac_concordance <- rbindlist(frac_conc, fill = TRUE)
safe_fwrite(frac_concordance, file.path(args$outdir, "tables", "MA_round2A_fraction_concordance_by_celltype.tsv"))

alt_models <- rbindlist(model_res, fill = TRUE)
# join with NNLS results from Round1B
round1b_long <- fread(file.path(args$nnls_dir, "tables", "MA_round1B_prepost_model_results_long.tsv"))
nnls_models <- round1b_long[model_type == "nnls_adjusted", .(cohort, program_name, beta_asd_nnls = beta_asd, se_nnls = se, ci_low_nnls = ci_low, ci_high_nnls = ci_high, p_nnls = p_value)]
alt_models2 <- alt_models[, .(cohort, program_name, beta_asd_alt = beta_asd, se_alt = se, ci_low_alt = ci_low, ci_high_alt = ci_high, p_alt = p_value)]
model_cmp <- merge(nnls_models, alt_models2, by = c("cohort","program_name"), all = TRUE)
model_cmp[, delta_alt_minus_nnls := beta_asd_alt - beta_asd_nnls]
model_cmp[, direction_match := sign(beta_asd_alt) == sign(beta_asd_nnls)]
safe_fwrite(alt_models, file.path(args$outdir, "tables", "MA_round2A_alt_adjusted_program_effects.tsv"))
safe_fwrite(model_cmp, file.path(args$outdir, "tables", "MA_round2A_alt_vs_nnls_program_effect_comparison.tsv"))

overall <- frac_concordance[, .(
  mean_pearson_r = mean(pearson_r, na.rm = TRUE),
  mean_spearman_r = mean(spearman_r, na.rm = TRUE),
  mean_abs_diff = mean(mean_abs_diff, na.rm = TRUE)
), by = cohort]
safe_fwrite(overall, file.path(args$outdir, "tables", "MA_round2A_fraction_concordance_overall.tsv"))

run_summary <- data.table(
  item = c("markers_per_class", "basis_genes", "n_fraction_concordance_rows", "n_alt_model_rows"),
  value = c(args$markers_per_class, nrow(basis), nrow(frac_concordance), nrow(alt_models))
)
safe_fwrite(run_summary, file.path(args$outdir, "meta", "MA_round2A_run_summary.tsv"))
log_msg("Completed Round2A alternative deconvolution.")
