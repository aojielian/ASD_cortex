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

LOGFILE <- file.path(LOGDIR, "Step08B_GSE64018_remap_and_reintegration_v3.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOGFILE, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"+|"+$', '', x)
  x <- gsub("^'+|'+$", "", x)
  x <- sub("^(ASD|CTL)_", "", x, ignore.case = FALSE)
  x <- gsub("\\s+", "", x)
  x
}

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
    setnames(dt, c("feature_id", h))
    attr(dt, "header_fix_method") <- "prepended_feature_id_to_unlabeled_first_column"
  } else {
    dt <- fread(cmd = paste("zcat", shQuote(file)))
    if (names(dt)[1] == "V1") setnames(dt, 1, "feature_id")
    attr(dt, "header_fix_method") <- "direct_fread"
  }

  dt
}

strip_ensembl_version <- function(x) sub("\\..*$", "", as.character(x))

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
BRAINSPAN_ROWS <- file.path(BASE, "Gencode_v3c_summarized_to_genes", "rows_metadata.csv")
STEP02_FILE <- file.path(TABDIR, "24_GSE102741_program_dx_tests.tsv")
STEP02_SUMMARY_FILE <- file.path(METADIR, "24b_Step02_run_summary.tsv")
STEP03_FILE <- file.path(TABDIR, "31_Gandal_program_dx_tests.tsv")
STEP06_FILE <- file.path(TABDIR, "90_crossLayer_program_summary.tsv")
SELECT_FILE <- file.path(TABDIR, "91_final_program_selection.tsv")

geo_map <- fread(text = '
gsm\tgeo_title
GSM1562865\tASD_AN02987_ba41-42-22_8.6
GSM1562866\tASD_AN04682_ba41-42-22_8.2
GSM1562867\tASD_UMB5278_ba41-42-22_7.8
GSM1562868\tASD_AN01570_ba41-42-22_7.1
GSM1562869\tASD_AN00493_ba41-42-22_7.3
GSM1562870\tASD_AN12457_ba41-42-22_6
GSM1562871\tASD_AN08166_ba41-42-22_6.4
GSM1562872\tASD_AN08792_ba41-42-22_5.1
GSM1562873\tASD_AN01971_ba41-42-22_2.9
GSM1562874\tASD_AN03632_ba41-42-22_8.1
GSM1562875\tASD_AN08043_ba41-42-22_7.8
GSM1562876\tASD_AN09714_ba41-42-22_7.2
GSM1562877\tCTL_AN17425_ba41-42-22_7.9
GSM1562878\tCTL_AN07444_ba41-42-22_8.2
GSM1562879\tCTL_UMB4590_ba41-42-22_8.3
GSM1562880\tCTL_AN10833_ba41-42-22_1.8
GSM1562881\tCTL_AN14757_ba41-42-22_8
GSM1562882\tCTL_AN19760_ba41-42-22_6.1
GSM1562883\tCTL_AN12137_ba41-42-22_6.4
GSM1562884\tCTL_UMB5079_ba41-42-22_7.9
GSM1562885\tCTL_AN08161_ba41-42-22_8.1
GSM1562886\tCTL_AN08677_ba41-42-22_8
GSM1562887\tCTL_UMB4842_ba41-42-22_8.1
GSM1562888\tCTL_AN13295_ba41-42-22_6.3
')
geo_map[, diagnosis := fifelse(startsWith(geo_title, "ASD_"), "ASD", "Control")]
geo_map[, geo_sample_raw := sub("^(ASD|CTL)_", "", geo_title)]
geo_map[, sample_norm := normalize_sample_id(geo_sample_raw)]

log_msg("Reading program definitions: ", PROGRAM_FILE)
programs <- c("SFARI_all", "midPrenatal_SFARI_top20", "midPrenatal_SFARI_top05", "midPrenatal_SFARI_top10")
program_genes <- lapply(programs, function(p) read_program_genes(PROGRAM_FILE, p))
names(program_genes) <- programs

log_msg("Reading GSE64018 counts: ", GSE64018_FILE)
gse <- read_gse64018_counts(GSE64018_FILE)
header_fix_method <- attr(gse, "header_fix_method")
sample_names_raw <- names(gse)[-1]
sample_names_norm <- normalize_sample_id(sample_names_raw)

file_samples <- data.table(
  file_sample_raw = sample_names_raw,
  sample_norm = sample_names_norm
)

sample_match <- merge(
  file_samples,
  geo_map,
  by = "sample_norm",
  all = TRUE,
  sort = FALSE
)
sample_match[, matched := !is.na(file_sample_raw) & !is.na(gsm)]
safe_fwrite(sample_match, file.path(METADIR, "120_GSE64018_explicit_sample_mapping.tsv"))

if (any(!sample_match$matched)) {
  unmatched_file <- sample_match[is.na(gsm), .(file_sample_raw, sample_norm)]
  unmatched_geo  <- sample_match[is.na(file_sample_raw), .(gsm, geo_title, geo_sample_raw, sample_norm)]
  if (nrow(unmatched_file) > 0) {
    safe_fwrite(unmatched_file, file.path(METADIR, "120b_GSE64018_unmatched_file_samples.tsv"))
  }
  if (nrow(unmatched_geo) > 0) {
    safe_fwrite(unmatched_geo, file.path(METADIR, "120c_GSE64018_unmatched_geo_samples.tsv"))
  }
  stop("Explicit sample mapping did not fully match after normalization. Check 120b/120c files.")
}

log_msg("Reading BrainSpan rows metadata: ", BRAINSPAN_ROWS)
rows <- fread(BRAINSPAN_ROWS)
rows[, ensembl_gene_id_novers := strip_ensembl_version(ensembl_gene_id)]
rows[, gene_symbol_upper := toupper(trimws(as.character(gene_symbol)))]
rowmap <- unique(rows[
  !is.na(ensembl_gene_id_novers) & ensembl_gene_id_novers != "" &
  !is.na(gene_symbol_upper) & gene_symbol_upper != "",
  .(ensembl_gene_id_novers, gene_symbol_upper)
])

gse_map <- data.table(
  feature_id = as.character(gse$feature_id),
  ensembl_gene_id_novers = strip_ensembl_version(gse$feature_id)
)
gse_map <- merge(gse_map, rowmap, by = "ensembl_gene_id_novers", all.x = FALSE, allow.cartesian = FALSE)
gse_map <- unique(gse_map[, .(feature_id, gene_symbol_upper)])
safe_fwrite(gse_map, file.path(TABDIR, "121_GSE64018_ensembl_to_symbol_mapping.tsv.gz"))

gse_dt <- copy(gse)
num_cols <- setdiff(names(gse_dt), "feature_id")
for (j in num_cols) gse_dt[[j]] <- as.numeric(gse_dt[[j]])
gse_mapped <- merge(gse_map, gse_dt, by = "feature_id", all.x = FALSE, allow.cartesian = FALSE)
gse_collapsed <- gse_mapped[, lapply(.SD, sum, na.rm = TRUE), by = gene_symbol_upper, .SDcols = num_cols]
setnames(gse_collapsed, "gene_symbol_upper", "gene_symbol")
safe_fwrite(gse_collapsed, file.path(TABDIR, "122_GSE64018_collapsed_counts_by_symbol.tsv.gz"))

expr_log <- compute_logcpm(gse_collapsed)
rownames(expr_log) <- gse_collapsed$gene_symbol

gse_prog_map <- data.table(
  program_name = programs,
  n_genes_program = vapply(program_genes, length, integer(1)),
  n_genes_in_matrix = vapply(program_genes, function(gs) sum(gs %in% rownames(expr_log)), integer(1))
)
safe_fwrite(gse_prog_map, file.path(METADIR, "123_GSE64018_program_mapping_summary.tsv"))

scores_list <- lapply(programs, function(p) {
  genes_use <- intersect(program_genes[[p]], rownames(expr_log))
  vals <- colMeans(expr_log[genes_use, , drop = FALSE], na.rm = TRUE)
  data.table(
    sample_raw = colnames(expr_log),
    sample_norm = normalize_sample_id(colnames(expr_log)),
    program_name = p,
    score = as.numeric(vals)
  )
})
gse_scores <- rbindlist(scores_list)
gse_scores <- merge(
  gse_scores,
  sample_match[, .(sample_norm, gsm, geo_title, diagnosis)],
  by = "sample_norm",
  all.x = TRUE
)
safe_fwrite(gse_scores, file.path(TABDIR, "124_GSE64018_program_scores.tsv.gz"))

gse_tests <- run_program_tests_simple(gse_scores[, .(sample = sample_raw, program_name, score, diagnosis)])
safe_fwrite(gse_tests, file.path(TABDIR, "125_GSE64018_program_dx_tests.tsv"))

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

gse64018_std <- gse_tests[, .(
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
safe_fwrite(bulk_three, file.path(TABDIR, "126_bulk_threeCohort_program_effects_updated.tsv"))

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
safe_fwrite(meta_tbl, file.path(TABDIR, "127_bulk_program_meta_summary_updated.tsv"))

sel_tbl <- fread(SELECT_FILE)
anchor_tbl <- merge(sel_tbl[, .(selection_role, program_name, rationale)], meta_tbl, by = "program_name", all.x = TRUE)
anchor_tbl <- merge(anchor_tbl, bulk_three[, .(cohort, program_name, delta_asd_minus_control, lm_p, fdr_lm, direction_match)], by = "program_name", all.x = TRUE)
safe_fwrite(anchor_tbl, file.path(TABDIR, "128_bulk_anchor_programs_for_manuscript_updated.tsv"))

run_summary <- rbindlist(list(
  data.table(section = "GSE64018", metric = "header_fix_method", value = header_fix_method),
  data.table(section = "GSE64018", metric = "n_samples", value = uniqueN(sample_match$sample_norm)),
  data.table(section = "GSE64018", metric = "n_asd", value = sum(sample_match$diagnosis == "ASD")),
  data.table(section = "GSE64018", metric = "n_control", value = sum(sample_match$diagnosis == "Control")),
  data.table(section = "GSE64018", metric = "n_feature_ids_mapped", value = uniqueN(gse_map$feature_id)),
  data.table(section = "GSE64018", metric = "n_gene_symbols_collapsed", value = uniqueN(gse_collapsed$gene_symbol)),
  data.table(section = "GSE64018", metric = "best_program_by_lm", value = gse_tests[which.min(fdr_lm), program_name]),
  data.table(section = "BulkMetaUpdated", metric = "best_program_weighted", value = meta_tbl[which.min(fdr_weighted), program_name]),
  data.table(section = "BulkMetaUpdated", metric = "best_program_weighted_p", value = meta_tbl[which.min(fdr_weighted), stouffer_p_weighted_one_sided]),
  data.table(section = "BulkMetaUpdated", metric = "best_program_weighted_fdr", value = meta_tbl[which.min(fdr_weighted), fdr_weighted])
))
safe_fwrite(run_summary, file.path(METADIR, "129_Step08B_run_summary.tsv"))

log_msg("Step08B v3 completed successfully.")
