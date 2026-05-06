#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOGDIR, "Step08_bulk_integration_and_sensitivity_v4.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOGFILE, append = TRUE)
}

safe_fwrite <- function(x, file) fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")

read_program_genes <- function(program_file, target_program) {
  dt <- fread(program_file)
  nms <- tolower(names(dt))
  pcol <- names(dt)[match("program_name", nms)]
  g_cands <- c("gene_symbol", "symbol", "gene", "standardized_gene_symbol")
  g_hit <- g_cands[g_cands %in% nms]
  if (length(g_hit) == 0) stop("No gene symbol column found in: ", program_file)
  gcol <- names(dt)[match(g_hit[1], nms)]
  out <- unique(toupper(trimws(dt[[gcol]][dt[[pcol]] %in% target_program])))
  out[!is.na(out) & out != ""]
}

read_gse64018_counts <- function(file) {
  if (!file.exists(file)) stop("Missing GSE64018 file: ", file)
  con <- gzfile(file, open = "rt")
  header_line <- readLines(con, n = 1)
  data_line <- readLines(con, n = 1)
  close(con)
  h <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  d <- strsplit(data_line, "\t", fixed = TRUE)[[1]]

  if (length(d) == length(h) + 1) {
    dt <- fread(cmd = paste("zcat", shQuote(file)), header = FALSE, skip = 1)
    setnames(dt, c("gene_symbol", h))
    attr(dt, "header_fix_method") <- "prepended_gene_symbol_to_unlabeled_first_column"
  } else {
    dt <- fread(cmd = paste("zcat", shQuote(file)))
    if (names(dt)[1] == "V1") setnames(dt, 1, "gene_symbol")
    attr(dt, "header_fix_method") <- "direct_fread"
  }
  dt
}

collapse_counts_by_symbol_sum <- function(dt, gene_col) {
  dt2 <- copy(dt)
  dt2[[gene_col]] <- toupper(trimws(dt2[[gene_col]]))
  dt2 <- dt2[!is.na(get(gene_col)) & get(gene_col) != ""]
  num_cols <- setdiff(names(dt2), gene_col)
  for (j in num_cols) dt2[[j]] <- as.numeric(dt2[[j]])
  collapsed <- dt2[, lapply(.SD, sum, na.rm = TRUE), by = gene_col, .SDcols = num_cols]
  setnames(collapsed, gene_col, "gene_symbol")
  collapsed
}

compute_logcpm <- function(count_dt) {
  mat <- as.matrix(count_dt[, -1])
  mode(mat) <- "numeric"
  libsize <- colSums(mat, na.rm = TRUE)
  libsize[libsize <= 0] <- NA_real_
  cpm <- t(t(mat) / libsize) * 1e6
  log2(cpm + 1)
}

calc_auc <- function(x, y) {
  n1 <- length(x); n0 <- length(y)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(c(x, y))
  r1 <- sum(ranks[seq_len(n1)])
  u1 <- r1 - n1 * (n1 + 1) / 2
  as.numeric(u1 / (n1 * n0))
}

run_program_tests_simple <- function(score_dt) {
  score_dt <- copy(score_dt)
  score_dt[, diagnosis := factor(as.character(diagnosis), levels = c("Control", "ASD"))]
  out <- score_dt[, {
    x0 <- score[diagnosis == "Control"]
    x1 <- score[diagnosis == "ASD"]
    n0 <- sum(diagnosis == "Control")
    n1 <- sum(diagnosis == "ASD")

    wil_p <- tryCatch(wilcox.test(x1, x0, exact = FALSE)$p.value, error = function(e) NA_real_)
    lm_fit <- tryCatch(lm(score ~ diagnosis), error = function(e) NULL)
    cf <- tryCatch(summary(lm_fit)$coefficients, error = function(e) NULL)
    lm_beta <- if (!is.null(cf) && "diagnosisASD" %in% rownames(cf)) unname(cf["diagnosisASD", "Estimate"]) else NA_real_
    lm_p    <- if (!is.null(cf) && "diagnosisASD" %in% rownames(cf)) unname(cf["diagnosisASD", "Pr(>|t|)"]) else NA_real_

    list(
      n_samples = uniqueN(sample),
      n_asd = n1,
      n_control = n0,
      mean_score_control = mean(x0, na.rm = TRUE),
      mean_score_asd = mean(x1, na.rm = TRUE),
      delta_asd_minus_control = mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE),
      auc_asd_higher = calc_auc(x1, x0),
      wilcox_p = wil_p,
      lm_beta_asd = lm_beta,
      lm_p = lm_p
    )
  }, by = program_name]

  out[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  out[, fdr_lm := p.adjust(lm_p, method = "BH")]
  out[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(out, fdr_lm, fdr_wilcox, -abs_delta_asd_minus_control)
  out
}

infer_gse64018_dx <- function(sample_names) {
  sn <- as.character(sample_names)
  dx <- rep(NA_character_, length(sn))
  method <- rep("unassigned", length(sn))

  is_asd <- grepl("(^|[^A-Z])ASD([^A-Z]|$)|AUTISM", sn, ignore.case = TRUE, perl = TRUE)
  is_ctl <- grepl("CTL|CTRL|CONTROL", sn, ignore.case = TRUE)

  dx[is_asd & !is_ctl] <- "ASD"
  method[is_asd & !is_ctl] <- "sample_name_regex"
  dx[is_ctl & !is_asd] <- "Control"
  method[is_ctl & !is_asd] <- "sample_name_regex"

  still_na <- is.na(dx)
  if (all(still_na) && length(sn) == 24) {
    dx[seq_len(12)] <- "ASD"
    dx[13:24] <- "Control"
    method[] <- "fallback_file_order_first12ASD_last12Control"
  }

  data.table(sample = sn, diagnosis = dx, infer_method = method)
}

signed_z_from_p <- function(delta, p) {
  if (is.na(delta) || is.na(p) || p <= 0 || p > 1) return(NA_real_)
  sign_delta <- ifelse(delta > 0, 1, ifelse(delta < 0, -1, 0))
  if (sign_delta == 0) return(0)
  sign_delta * qnorm(1 - p / 2)
}

read_step02_counts_from_summary <- function(file) {
  if (!file.exists(file)) stop("Missing Step02 run summary file: ", file)
  dt <- fread(file)
  getv <- function(sec, met) {
    x <- dt[section == sec & metric == met, value]
    if (length(x) == 0) NA_real_ else as.numeric(x[1])
  }
  list(
    n_samples = getv("GSE102741", "n_samples"),
    n_asd = getv("GSE102741", "n_asd"),
    n_control = getv("GSE102741", "n_control")
  )
}

PROGRAM_FILE <- file.path(TABDIR, "20_midPrenatal_core_programs.tsv")
GSE64018_FILE <- file.path(BASE, "GSE64018_countlevel_12asd_12ctl.txt.gz")
STEP02_FILE <- file.path(TABDIR, "24_GSE102741_program_dx_tests.tsv")
STEP02_SUMMARY_FILE <- file.path(METADIR, "24b_Step02_run_summary.tsv")
STEP03_FILE <- file.path(TABDIR, "31_Gandal_program_dx_tests.tsv")
STEP06_FILE <- file.path(TABDIR, "90_crossLayer_program_summary.tsv")
SELECT_FILE <- file.path(TABDIR, "91_final_program_selection.tsv")

log_msg("Reading program definitions: ", PROGRAM_FILE)
programs <- c("SFARI_all", "midPrenatal_SFARI_top20", "midPrenatal_SFARI_top05", "midPrenatal_SFARI_top10")
program_genes <- lapply(programs, function(p) read_program_genes(PROGRAM_FILE, p))
names(program_genes) <- programs

log_msg("Reading GSE64018 bulk counts: ", GSE64018_FILE)
gse64018 <- read_gse64018_counts(GSE64018_FILE)
header_fix_method <- attr(gse64018, "header_fix_method")
gse64018_collapsed <- collapse_counts_by_symbol_sum(gse64018, "gene_symbol")
safe_fwrite(gse64018_collapsed, file.path(TABDIR, "100_GSE64018_collapsed_counts.tsv.gz"))

sample_names <- names(gse64018_collapsed)[-1]
sample_meta <- infer_gse64018_dx(sample_names)
safe_fwrite(sample_meta, file.path(METADIR, "100b_GSE64018_sample_inference.tsv"))
if (any(is.na(sample_meta$diagnosis))) stop("Some GSE64018 sample diagnoses could not be inferred even after fallback.")

expr_log <- compute_logcpm(gse64018_collapsed)
rownames(expr_log) <- gse64018_collapsed$gene_symbol

gse64018_map <- data.table(
  program_name = programs,
  n_genes_program = vapply(program_genes, length, integer(1)),
  n_genes_in_matrix = vapply(program_genes, function(gs) sum(gs %in% rownames(expr_log)), integer(1))
)
safe_fwrite(gse64018_map, file.path(METADIR, "100c_GSE64018_program_mapping_summary.tsv"))

scores_list <- lapply(programs, function(p) {
  genes_use <- intersect(program_genes[[p]], rownames(expr_log))
  vals <- colMeans(expr_log[genes_use, , drop = FALSE], na.rm = TRUE)
  data.table(sample = colnames(expr_log), program_name = p, score = as.numeric(vals))
})
gse64018_scores <- rbindlist(scores_list)
gse64018_scores <- merge(gse64018_scores, sample_meta, by = "sample", all.x = TRUE)
safe_fwrite(gse64018_scores, file.path(TABDIR, "100d_GSE64018_program_scores.tsv.gz"))

gse64018_tests <- run_program_tests_simple(gse64018_scores)
safe_fwrite(gse64018_tests, file.path(TABDIR, "101_GSE64018_program_dx_tests.tsv"))

log_msg("Reading existing GSE102741 summary: ", STEP02_FILE)
gse102741 <- fread(STEP02_FILE)
if (!("n_samples" %in% names(gse102741)) || !("n_asd" %in% names(gse102741)) || !("n_control" %in% names(gse102741))) {
  log_msg("Step02 summary missing n_samples/n_asd/n_control; reading fallback counts from: ", STEP02_SUMMARY_FILE)
  cts <- read_step02_counts_from_summary(STEP02_SUMMARY_FILE)
  gse102741[, n_samples := cts$n_samples]
  gse102741[, n_asd := cts$n_asd]
  gse102741[, n_control := cts$n_control]
}
gse102741_std <- gse102741[, .(
  cohort = "GSE102741",
  program_name,
  n_samples,
  n_asd,
  n_control,
  delta_asd_minus_control,
  lm_p = lm_p,
  fdr_lm = fdr_lm
)]

log_msg("Reading existing Gandal summary: ", STEP03_FILE)
gandal <- fread(STEP03_FILE)
gandal_std <- gandal[, .(
  cohort = "Gandal2022",
  program_name,
  n_samples,
  n_asd,
  n_control,
  delta_asd_minus_control,
  lm_p = lm_region_p,
  fdr_lm = fdr_lm_region
)]

gse64018_std <- gse64018_tests[, .(
  cohort = "GSE64018",
  program_name,
  n_samples,
  n_asd,
  n_control,
  delta_asd_minus_control,
  lm_p,
  fdr_lm
)]

bulk_three <- rbindlist(list(gse102741_std, gandal_std, gse64018_std), use.names = TRUE, fill = TRUE)

dir_tbl <- fread(STEP06_FILE)[, .(program_name, adult_expected_direction)]
bulk_three <- merge(bulk_three, dir_tbl, by = "program_name", all.x = TRUE)
bulk_three[, expected_sign := fifelse(adult_expected_direction == "ASD_lower_than_control", -1,
                               fifelse(adult_expected_direction == "ASD_higher_than_control", 1, NA_real_))]
bulk_three[, signed_z := mapply(signed_z_from_p, delta_asd_minus_control, lm_p)]
bulk_three[, aligned_z := expected_sign * signed_z]
bulk_three[, direction_match := sign(delta_asd_minus_control) == expected_sign]
bulk_three[, support_label := fifelse(direction_match, "direction_matched", "direction_mismatched")]
safe_fwrite(bulk_three, file.path(TABDIR, "102_bulk_threeCohort_program_effects.tsv"))

meta_tbl <- bulk_three[, {
  z <- aligned_z[!is.na(aligned_z)]
  n <- n_samples[!is.na(aligned_z)]
  unweighted_z <- if (length(z) > 0) sum(z) / sqrt(length(z)) else NA_real_
  weighted_z <- if (length(z) > 0) sum(z * sqrt(n)) / sqrt(sum(n)) else NA_real_
  list(
    n_cohorts = .N,
    n_valid_z = length(z),
    n_direction_matched = sum(direction_match, na.rm = TRUE),
    n_direction_mismatched = sum(!direction_match, na.rm = TRUE),
    mean_delta = mean(delta_asd_minus_control, na.rm = TRUE),
    stouffer_z_unweighted = unweighted_z,
    stouffer_p_unweighted_one_sided = ifelse(is.na(unweighted_z), NA_real_, 1 - pnorm(unweighted_z)),
    stouffer_z_weighted = weighted_z,
    stouffer_p_weighted_one_sided = ifelse(is.na(weighted_z), NA_real_, 1 - pnorm(weighted_z))
  )
}, by = .(program_name, adult_expected_direction)]
meta_tbl[, fdr_weighted := p.adjust(stouffer_p_weighted_one_sided, method = "BH")]
meta_tbl[, fdr_unweighted := p.adjust(stouffer_p_unweighted_one_sided, method = "BH")]
setorder(meta_tbl, fdr_weighted, fdr_unweighted, -n_direction_matched)
safe_fwrite(meta_tbl, file.path(TABDIR, "103_bulk_program_meta_summary.tsv"))

sel_tbl <- fread(SELECT_FILE)
anchor_tbl <- merge(sel_tbl[, .(selection_role, program_name, rationale)], meta_tbl, by = "program_name", all.x = TRUE)
anchor_tbl <- merge(anchor_tbl, bulk_three[, .(cohort, program_name, delta_asd_minus_control, lm_p, fdr_lm, direction_match)], by = "program_name", all.x = TRUE)
safe_fwrite(anchor_tbl, file.path(TABDIR, "104_bulk_anchor_programs_for_manuscript.tsv"))

run_summary <- rbindlist(list(
  data.table(section = "GSE64018", metric = "header_fix_method", value = header_fix_method),
  data.table(section = "GSE64018", metric = "n_samples", value = uniqueN(sample_meta$sample)),
  data.table(section = "GSE64018", metric = "n_asd", value = sum(sample_meta$diagnosis == "ASD")),
  data.table(section = "GSE64018", metric = "n_control", value = sum(sample_meta$diagnosis == "Control")),
  data.table(section = "GSE64018", metric = "inference_methods", value = paste(unique(sample_meta$infer_method), collapse = ",")),
  data.table(section = "GSE64018", metric = "best_program_by_lm", value = gse64018_tests[which.min(fdr_lm), program_name]),
  data.table(section = "BulkMeta", metric = "best_program_weighted", value = meta_tbl[which.min(fdr_weighted), program_name]),
  data.table(section = "BulkMeta", metric = "best_program_weighted_p", value = meta_tbl[which.min(fdr_weighted), stouffer_p_weighted_one_sided]),
  data.table(section = "BulkMeta", metric = "best_program_weighted_fdr", value = meta_tbl[which.min(fdr_weighted), fdr_weighted])
))
safe_fwrite(run_summary, file.path(METADIR, "105_Step08_run_summary.tsv"))

log_msg("Step08 v4 completed successfully.")
