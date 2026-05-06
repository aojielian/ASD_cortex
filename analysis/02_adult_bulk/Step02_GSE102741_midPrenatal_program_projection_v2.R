#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
RDSDIR <- file.path(OUTDIR, "rds")
LOGDIR <- file.path(OUTDIR, "logs")

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDSDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step02_GSE102741_midPrenatal_program_projection_v2.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

calc_auc <- function(score, is_case) {
  is_case <- as.logical(is_case)
  n1 <- sum(is_case)
  n0 <- sum(!is_case)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  U <- sum(ranks[is_case]) - n1 * (n1 + 1) / 2
  U / (n1 * n0)
}

mean_matched_sample <- function(target_genes, bg_dt) {
  tg <- bg_dt[gene_symbol %in% target_genes]
  out <- character(0)
  for (b in unique(tg$bin)) {
    n_b <- sum(tg$bin == b, na.rm = TRUE)
    pool_b <- bg_dt[!(gene_symbol %in% target_genes) & bin == b, gene_symbol]
    if (length(pool_b) < n_b) {
      pool_b <- bg_dt[!(gene_symbol %in% target_genes), gene_symbol]
    }
    out <- c(out, sample(pool_b, n_b, replace = FALSE))
  }
  unique(out)
}

set.seed(20260319)

sfari_file <- file.path(TABDIR, "01_SFARI_primary_Sand1_standardized.tsv")
stage_z_file <- file.path(TABDIR, "14_BrainSpan_stage_specificity_zscores.tsv.gz")
brainspan_rows_file <- file.path(METADIR, "12c_BrainSpan_rows_metadata_standardized.tsv")
gse_file <- file.path(BASE, "GSE102741_raw_counts_GRCh38.p13_NCBI.tsv.gz")

log_msg("Reading SFARI standardized file: ", sfari_file)
sfari <- fread(sfari_file)
sfari[, gene_symbol := toupper(trimws(gene_symbol))]
sfari <- unique(sfari[gene_symbol != "" & !is.na(gene_symbol)])

log_msg("Reading BrainSpan stage specificity file: ", stage_z_file)
stage_z <- fread(stage_z_file)
stage_z[, gene_symbol := toupper(trimws(gene_symbol))]
if (!"mid_prenatal" %in% names(stage_z)) {
  stop("stage_z file must contain mid_prenatal column.")
}

log_msg("Reading BrainSpan row metadata standardized: ", brainspan_rows_file)
rowmap <- fread(brainspan_rows_file)
rowmap[, gene_symbol := toupper(trimws(gene_symbol))]
rowmap[, entrez_id := trimws(as.character(entrez_id))]
rowmap <- unique(
  rowmap[gene_symbol != "" & !is.na(gene_symbol) & entrez_id != "" & !is.na(entrez_id),
         .(entrez_id, gene_symbol)]
)

core_dt <- copy(stage_z[, .(gene_symbol, mid_prenatal)])
core_dt[, is_sfari := gene_symbol %in% sfari$gene_symbol]
core_dt <- core_dt[!is.na(mid_prenatal)]
core_dt <- core_dt[order(-mid_prenatal)]
core_dt[, rank_desc := seq_len(.N)]
core_dt[, frac_rank := rank_desc / .N]

program_defs <- list(
  midPrenatal_SFARI_top05 = core_dt[is_sfari == TRUE & frac_rank <= 0.05, gene_symbol],
  midPrenatal_SFARI_top10 = core_dt[is_sfari == TRUE & frac_rank <= 0.10, gene_symbol],
  midPrenatal_SFARI_top20 = core_dt[is_sfari == TRUE & frac_rank <= 0.20, gene_symbol],
  SFARI_all = intersect(sfari$gene_symbol, core_dt$gene_symbol)
)

program_tbl <- rbindlist(lapply(names(program_defs), function(nm) {
  gs <- sort(unique(program_defs[[nm]]))
  data.table(program_name = nm, gene_symbol = gs)
}))
safe_fwrite(program_tbl, file.path(TABDIR, "20_midPrenatal_core_programs.tsv"))

program_summary <- program_tbl[, .N, by = program_name]
setnames(program_summary, "N", "n_genes")
safe_fwrite(program_summary, file.path(METADIR, "20b_midPrenatal_core_program_summary.tsv"))
log_msg("Program sizes: ", paste(program_summary$program_name, program_summary$n_genes, sep = "=", collapse = "; "))

log_msg("Reading GSE102741 counts: ", gse_file)
gse <- fread(cmd = paste("zcat", shQuote(gse_file)))
geneid_col <- names(gse)[1]
setnames(gse, geneid_col, "GeneID")
gse[, GeneID := trimws(as.character(GeneID))]

sample_ids <- setdiff(names(gse), "GeneID")
log_msg("GSE102741 samples detected: ", length(sample_ids))

mapped <- merge(gse, rowmap, by.x = "GeneID", by.y = "entrez_id", all.x = FALSE, all.y = FALSE)
log_msg("Mapped GSE rows to gene_symbol: ", nrow(mapped))
if (nrow(mapped) < 10000) {
  stop("Too few GSE rows mapped to gene symbols.")
}

count_cols <- sample_ids
mapped_collapsed <- mapped[, lapply(.SD, sum, na.rm = TRUE), by = gene_symbol, .SDcols = count_cols]
expr_counts <- as.matrix(mapped_collapsed[, ..count_cols])
rownames(expr_counts) <- mapped_collapsed$gene_symbol
mode(expr_counts) <- "numeric"

libsize <- colSums(expr_counts, na.rm = TRUE)
cpm <- t(t(expr_counts) / libsize * 1e6)
log2cpm <- log2(cpm + 1)

saveRDS(list(
  counts = expr_counts,
  log2cpm = log2cpm,
  sample_ids = sample_ids
), file.path(RDSDIR, "22_GSE102741_geneSymbol_log2cpm.rds"))

# NOTE: current grouping is fixed from the prior plan
asd_ids <- paste0("GSM27459", sprintf("%02d", 36:48))
sample_meta <- data.table(sample_id = sample_ids)
sample_meta[, diagnosis := ifelse(sample_id %in% asd_ids, "ASD", "Control")]
sample_meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
sample_meta[, is_asd := diagnosis == "ASD"]
safe_fwrite(sample_meta, file.path(TABDIR, "21_GSE102741_sample_metadata.tsv"))

log_msg("GSE102741 diagnosis counts: ASD=", sum(sample_meta$is_asd), "; Control=", sum(!sample_meta$is_asd))

gene_mean <- rowMeans(log2cpm, na.rm = TRUE)
gene_bg <- data.table(gene_symbol = names(gene_mean), mean_expr = as.numeric(gene_mean))
gene_bg[, bin := cut(mean_expr,
                     breaks = quantile(mean_expr, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
                     include.lowest = TRUE, ordered_result = TRUE)]
gene_bg[, bin := as.character(bin)]

score_long_list <- list()
test_list <- list()

nperm <- 1000

for (pn in names(program_defs)) {
  genes <- intersect(unique(program_defs[[pn]]), rownames(log2cpm))
  if (length(genes) < 5) {
    log_msg("Skipping ", pn, " because too few mapped genes in GSE102741: ", length(genes))
    next
  }

  prog_score <- colMeans(log2cpm[genes, , drop = FALSE], na.rm = TRUE)
  dt_score <- merge(
    sample_meta,
    data.table(sample_id = names(prog_score), score = as.numeric(prog_score)),
    by = "sample_id"
  )
  dt_score[, program_name := pn]

  wt <- wilcox.test(score ~ diagnosis, data = dt_score, exact = FALSE)
  fit <- lm(score ~ diagnosis, data = dt_score)
  coef_tab <- summary(fit)$coefficients
  beta_asd <- coef_tab["diagnosisASD", "Estimate"]
  p_lm <- coef_tab["diagnosisASD", "Pr(>|t|)"]

  mean_case <- dt_score[diagnosis == "ASD", mean(score)]
  mean_ctrl <- dt_score[diagnosis == "Control", mean(score)]
  delta_asd_minus_ctrl <- mean_case - mean_ctrl
  auc <- calc_auc(dt_score$score, dt_score$is_asd)

  perm_vals <- replicate(nperm, {
    gs_perm <- mean_matched_sample(genes, gene_bg)
    mean(colMeans(log2cpm[gs_perm, sample_meta[diagnosis == "ASD", sample_id], drop = FALSE], na.rm = TRUE)) -
      mean(colMeans(log2cpm[gs_perm, sample_meta[diagnosis == "Control", sample_id], drop = FALSE], na.rm = TRUE))
  })

  emp_p_greater <- (sum(perm_vals >= delta_asd_minus_ctrl) + 1) / (length(perm_vals) + 1)
  emp_p_less <- (sum(perm_vals <= delta_asd_minus_ctrl) + 1) / (length(perm_vals) + 1)
  z_perm <- ifelse(sd(perm_vals) > 0, (delta_asd_minus_ctrl - mean(perm_vals)) / sd(perm_vals), NA_real_)

  score_long_list[[pn]] <- dt_score

  test_list[[pn]] <- data.table(
    program_name = pn,
    n_genes_in_program = length(unique(program_defs[[pn]])),
    n_genes_mapped_gse = length(genes),
    mean_score_control = mean_ctrl,
    mean_score_asd = mean_case,
    delta_asd_minus_control = delta_asd_minus_ctrl,
    auc_asd_higher = auc,
    wilcox_p = wt$p.value,
    lm_beta_asd = beta_asd,
    lm_p = p_lm,
    perm_mean_delta = mean(perm_vals),
    perm_sd_delta = sd(perm_vals),
    perm_z = z_perm,
    emp_p_greater = emp_p_greater,
    emp_p_less = emp_p_less
  )
}

score_long <- rbindlist(score_long_list, use.names = TRUE, fill = TRUE)
test_res <- rbindlist(test_list, use.names = TRUE, fill = TRUE)

test_res[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
test_res[, fdr_lm := p.adjust(lm_p, method = "BH")]
test_res[, fdr_emp_greater := p.adjust(emp_p_greater, method = "BH")]
test_res[, fdr_emp_less := p.adjust(emp_p_less, method = "BH")]
test_res[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]

setorder(test_res, fdr_wilcox, wilcox_p, -abs_delta_asd_minus_control)

safe_fwrite(score_long, file.path(TABDIR, "23_GSE102741_program_scores_long.tsv"))
safe_fwrite(test_res, file.path(TABDIR, "24_GSE102741_program_dx_tests.tsv"))

best_program <- test_res[1]

run_summary <- rbindlist(list(
  data.table(section = "Programs", metric = "n_programs_tested", value = nrow(test_res)),
  data.table(section = "Programs", metric = "best_program", value = best_program$program_name),
  data.table(section = "Programs", metric = "best_program_n_genes_mapped", value = best_program$n_genes_mapped_gse),
  data.table(section = "GSE102741", metric = "n_samples", value = nrow(sample_meta)),
  data.table(section = "GSE102741", metric = "n_asd", value = sum(sample_meta$is_asd)),
  data.table(section = "GSE102741", metric = "n_control", value = sum(!sample_meta$is_asd)),
  data.table(section = "BestProgram", metric = "delta_asd_minus_control", value = best_program$delta_asd_minus_control),
  data.table(section = "BestProgram", metric = "auc_asd_higher", value = best_program$auc_asd_higher),
  data.table(section = "BestProgram", metric = "wilcox_p", value = best_program$wilcox_p),
  data.table(section = "BestProgram", metric = "fdr_wilcox", value = best_program$fdr_wilcox),
  data.table(section = "BestProgram", metric = "lm_beta_asd", value = best_program$lm_beta_asd),
  data.table(section = "BestProgram", metric = "lm_p", value = best_program$lm_p),
  data.table(section = "BestProgram", metric = "perm_z", value = best_program$perm_z),
  data.table(section = "BestProgram", metric = "emp_p_greater", value = best_program$emp_p_greater),
  data.table(section = "BestProgram", metric = "emp_p_less", value = best_program$emp_p_less)
))
safe_fwrite(run_summary, file.path(METADIR, "24b_Step02_run_summary.tsv"))

log_msg("Step02 v2 completed successfully.")
log_msg("Best program: ", best_program$program_name)
log_msg("Best program delta ASD-Control: ", best_program$delta_asd_minus_control)
log_msg("Best program Wilcox p: ", best_program$wilcox_p, " ; FDR=", best_program$fdr_wilcox)
log_msg("Best program permutation z: ", best_program$perm_z, " ; emp_p_greater=", best_program$emp_p_greater, " ; emp_p_less=", best_program$emp_p_less)
