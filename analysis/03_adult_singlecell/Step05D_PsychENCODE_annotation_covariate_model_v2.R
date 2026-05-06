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

SCORE_FILE <- file.path(TABDIR, "77_PsychENCODE_cell_program_scores_cleanSubtype.tsv.gz")
META_FILE  <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/meta.tsv"

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step05D_PsychENCODE_annotation_covariate_model_v2.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

empty_dt <- function(cols) {
  out <- data.table(matrix(ncol = length(cols), nrow = 0))
  setnames(out, cols)
  out
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

coef_asd <- function(fit) {
  cf <- summary(fit)$coefficients
  rn <- rownames(cf)
  if ("diagnosisASD" %in% rn) {
    return(list(beta = unname(cf["diagnosisASD", "Estimate"]),
                p = unname(cf["diagnosisASD", "Pr(>|t|)"]),
                coef_name = "diagnosisASD"))
  }
  if ("diagnosisControl" %in% rn) {
    return(list(beta = -unname(cf["diagnosisControl", "Estimate"]),
                p = unname(cf["diagnosisControl", "Pr(>|t|)"]),
                coef_name = "diagnosisControl_flipped"))
  }
  return(list(beta = NA_real_, p = NA_real_, coef_name = NA_character_))
}

for (f in c(SCORE_FILE, META_FILE)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

log_msg("Reading clean score file: ", SCORE_FILE)
score_cell <- fread(SCORE_FILE)
req_score <- c("Cell_ID","score","program_name","annotation","individual_ID")
if (!all(req_score %in% names(score_cell))) {
  stop("Score file must contain: ", paste(req_score, collapse = ", "))
}
score_cell[, Cell_ID := as.character(Cell_ID)]
score_cell[, program_name := as.character(program_name)]
score_cell[, annotation := as.character(annotation)]
score_cell[, individual_ID := as.character(individual_ID)]

# Remove any old covariate columns before merge to avoid .x/.y suffixes
drop_cols <- intersect(
  c("diagnosis","is_asd","Brain_Region","Age","Sex_Chromosome","PMI","RIN"),
  names(score_cell)
)
if (length(drop_cols) > 0) {
  log_msg("Dropping pre-existing columns from score file before merge: ", paste(drop_cols, collapse = ", "))
  score_cell[, (drop_cols) := NULL]
}

log_msg("Reading PsychENCODE meta: ", META_FILE)
meta <- fread(META_FILE)
req_meta <- c("Cell_ID","individual_ID","Diagnosis","Brain_Region","Age","Sex_Chromosome","PMI","RIN")
if (!all(req_meta %in% names(meta))) {
  stop("Meta file must contain: ", paste(req_meta, collapse = ", "))
}
meta[, Cell_ID := as.character(Cell_ID)]
meta[, individual_ID := as.character(individual_ID)]
meta[, Diagnosis := as.character(Diagnosis)]
meta[, Brain_Region := as.character(Brain_Region)]
meta[, Sex_Chromosome := as.character(Sex_Chromosome)]
meta[, Age := suppressWarnings(as.numeric(Age))]
meta[, PMI := suppressWarnings(as.numeric(PMI))]
meta[, RIN := suppressWarnings(as.numeric(RIN))]
meta[, diagnosis := fifelse(Diagnosis == "ASD", "ASD",
                            fifelse(Diagnosis == "CTL", "Control", NA_character_))]
meta <- meta[!is.na(diagnosis)]
meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
meta[, is_asd := diagnosis == "ASD"]

meta_keep <- unique(meta[, .(
  Cell_ID, individual_ID, diagnosis, is_asd, Brain_Region, Age, Sex_Chromosome, PMI, RIN
)])

score_cov <- merge(score_cell, meta_keep, by = c("Cell_ID","individual_ID"), all.x = TRUE)
safe_fwrite(score_cov, file.path(TABDIR, "80_PsychENCODE_cell_program_scores_cleanSubtype_withCov.tsv.gz"))

donor_ann <- score_cov[, .(
  n_cells = uniqueN(Cell_ID),
  mean_score = mean(score, na.rm = TRUE)
), by = .(
  program_name, annotation, individual_ID, diagnosis, is_asd,
  Brain_Region, Age, Sex_Chromosome, PMI, RIN
)]

safe_fwrite(donor_ann, file.path(TABDIR, "80b_PsychENCODE_donorRegion_annotation_program_scores.tsv.gz"))

consistency <- donor_ann[, .(
  n_rows = .N,
  n_regions = uniqueN(Brain_Region),
  n_age = uniqueN(Age),
  n_sex = uniqueN(Sex_Chromosome),
  n_pmi = uniqueN(PMI),
  n_rin = uniqueN(RIN),
  n_dx = uniqueN(diagnosis)
), by = individual_ID]
safe_fwrite(consistency, file.path(METADIR, "80c_PsychENCODE_donor_covariate_consistency.tsv"))

results_list <- list()
programs <- sort(unique(donor_ann$program_name))
annotations <- sort(unique(donor_ann$annotation))

for (pn in programs) {
  for (ann in annotations) {
    dd_all <- donor_ann[program_name == pn & annotation == ann]
    if (nrow(dd_all) < 8) next
    if (sum(dd_all$diagnosis == "ASD") < 3 || sum(dd_all$diagnosis == "Control") < 3) next

    dd_all[, diagnosis := factor(as.character(diagnosis), levels = c("Control", "ASD"))]
    dd_all[, Brain_Region := factor(Brain_Region)]
    dd_all[, Sex_Chromosome := factor(Sex_Chromosome)]

    wt_red <- wilcox.test(mean_score ~ diagnosis, data = dd_all, exact = FALSE)
    fit_red <- lm(mean_score ~ diagnosis, data = dd_all)
    ca_red <- coef_asd(fit_red)

    dd_adj <- dd_all[complete.cases(mean_score, diagnosis, Brain_Region, Age, Sex_Chromosome, PMI, RIN)]
    n_adj <- nrow(dd_adj)
    n_asd_adj <- sum(dd_adj$diagnosis == "ASD")
    n_ctl_adj <- sum(dd_adj$diagnosis == "Control")

    adj_beta <- NA_real_
    adj_p <- NA_real_
    adj_coef_name <- NA_character_
    wt_adj_p <- NA_real_
    auc_adj <- NA_real_
    delta_adj <- NA_real_
    mean_ctl_adj <- NA_real_
    mean_asd_adj <- NA_real_
    n_regions_adj <- if (n_adj > 0) uniqueN(dd_adj$Brain_Region) else NA_integer_

    if (n_adj >= 8 && n_asd_adj >= 3 && n_ctl_adj >= 3) {
      dd_adj[, diagnosis := factor(as.character(diagnosis), levels = c("Control", "ASD"))]
      dd_adj[, Brain_Region := factor(Brain_Region)]
      dd_adj[, Sex_Chromosome := factor(Sex_Chromosome)]

      wt_adj <- wilcox.test(mean_score ~ diagnosis, data = dd_adj, exact = FALSE)
      fit_adj <- lm(mean_score ~ diagnosis + Brain_Region + Age + Sex_Chromosome + PMI + RIN, data = dd_adj)
      ca_adj <- coef_asd(fit_adj)

      wt_adj_p <- wt_adj$p.value
      adj_beta <- ca_adj$beta
      adj_p <- ca_adj$p
      adj_coef_name <- ca_adj$coef_name
      auc_adj <- calc_auc(dd_adj$mean_score, dd_adj$is_asd)
      mean_ctl_adj <- dd_adj[diagnosis == "Control", mean(mean_score)]
      mean_asd_adj <- dd_adj[diagnosis == "ASD", mean(mean_score)]
      delta_adj <- mean_asd_adj - mean_ctl_adj
    }

    results_list[[paste0(pn, "__", ann)]] <- data.table(
      program_name = pn,
      annotation = ann,
      n_donorRegion_total = nrow(dd_all),
      n_asd_total = sum(dd_all$diagnosis == "ASD"),
      n_control_total = sum(dd_all$diagnosis == "Control"),
      mean_score_control_reduced = dd_all[diagnosis == "Control", mean(mean_score)],
      mean_score_asd_reduced = dd_all[diagnosis == "ASD", mean(mean_score)],
      delta_asd_minus_control_reduced = dd_all[diagnosis == "ASD", mean(mean_score)] - dd_all[diagnosis == "Control", mean(mean_score)],
      auc_reduced = calc_auc(dd_all$mean_score, dd_all$is_asd),
      wilcox_p_reduced = wt_red$p.value,
      lm_beta_asd_reduced = ca_red$beta,
      lm_p_reduced = ca_red$p,
      lm_coef_name_reduced = ca_red$coef_name,
      n_complete_adjusted = n_adj,
      n_asd_complete = n_asd_adj,
      n_control_complete = n_ctl_adj,
      n_regions_complete = n_regions_adj,
      mean_score_control_adjustedSet = mean_ctl_adj,
      mean_score_asd_adjustedSet = mean_asd_adj,
      delta_asd_minus_control_adjustedSet = delta_adj,
      auc_adjustedSet = auc_adj,
      wilcox_p_adjustedSet = wt_adj_p,
      lm_beta_asd_adjusted = adj_beta,
      lm_p_adjusted = adj_p,
      lm_coef_name_adjusted = adj_coef_name
    )
  }
}

res <- if (length(results_list) > 0) rbindlist(results_list, use.names = TRUE, fill = TRUE) else empty_dt(c(
  "program_name","annotation","n_donorRegion_total","n_asd_total","n_control_total",
  "mean_score_control_reduced","mean_score_asd_reduced","delta_asd_minus_control_reduced",
  "auc_reduced","wilcox_p_reduced","lm_beta_asd_reduced","lm_p_reduced","lm_coef_name_reduced",
  "n_complete_adjusted","n_asd_complete","n_control_complete","n_regions_complete",
  "mean_score_control_adjustedSet","mean_score_asd_adjustedSet","delta_asd_minus_control_adjustedSet",
  "auc_adjustedSet","wilcox_p_adjustedSet","lm_beta_asd_adjusted","lm_p_adjusted","lm_coef_name_adjusted"
))

if (nrow(res) > 0) {
  res[, fdr_lm_reduced := p.adjust(lm_p_reduced, method = "BH")]
  res[, fdr_lm_adjusted := p.adjust(lm_p_adjusted, method = "BH")]
  res[, fdr_wilcox_reduced := p.adjust(wilcox_p_reduced, method = "BH")]
  res[, fdr_wilcox_adjustedSet := p.adjust(wilcox_p_adjustedSet, method = "BH")]
  res[, abs_delta_adjustedSet := abs(delta_asd_minus_control_adjustedSet)]
  setorder(res, fdr_lm_adjusted, fdr_lm_reduced, -abs_delta_adjustedSet)
}
safe_fwrite(res, file.path(TABDIR, "81_PsychENCODE_annotation_dx_tests_adjusted.tsv"))

best_adj <- if (nrow(res) > 0) res[1] else NULL

run_summary <- rbindlist(list(
  data.table(section = "PsychENCODE", metric = "n_clean_cells_long", value = nrow(score_cov)),
  data.table(section = "PsychENCODE", metric = "n_clean_unique_cells", value = uniqueN(score_cov$Cell_ID)),
  data.table(section = "PsychENCODE", metric = "n_donorRegion_rows", value = nrow(donor_ann)),
  data.table(section = "PsychENCODE", metric = "n_individuals", value = uniqueN(donor_ann$individual_ID)),
  data.table(section = "PsychENCODE", metric = "n_annotations", value = uniqueN(donor_ann$annotation)),
  data.table(section = "PsychENCODE", metric = "n_programs", value = uniqueN(donor_ann$program_name))
), fill = TRUE)

if (!is.null(best_adj)) {
  run_summary <- rbind(run_summary, rbindlist(list(
    data.table(section = "BestAdjusted", metric = "program_name", value = best_adj$program_name),
    data.table(section = "BestAdjusted", metric = "annotation", value = best_adj$annotation),
    data.table(section = "BestAdjusted", metric = "delta_asd_minus_control_adjustedSet", value = best_adj$delta_asd_minus_control_adjustedSet),
    data.table(section = "BestAdjusted", metric = "lm_p_adjusted", value = best_adj$lm_p_adjusted),
    data.table(section = "BestAdjusted", metric = "fdr_lm_adjusted", value = best_adj$fdr_lm_adjusted),
    data.table(section = "BestReduced", metric = "lm_p_reduced", value = best_adj$lm_p_reduced),
    data.table(section = "BestReduced", metric = "fdr_lm_reduced", value = best_adj$fdr_lm_reduced)
  )), fill = TRUE)
}

safe_fwrite(run_summary, file.path(METADIR, "82_Step05D_run_summary.tsv"))

log_msg("Step05D v2 completed successfully.")
