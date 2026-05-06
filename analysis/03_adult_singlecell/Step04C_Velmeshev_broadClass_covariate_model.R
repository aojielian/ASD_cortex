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

SCORE_FILE <- file.path(TABDIR, "51_Velmeshev_cell_program_scores.tsv.gz")
META_FILE  <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/meta.tsv"

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step04C_Velmeshev_broadClass_covariate_model.log")

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

log_msg("Reading score file: ", SCORE_FILE)
score_cell <- fread(SCORE_FILE)
req_score <- c("cell","program_name","score","broad_class","individual")
if (!all(req_score %in% names(score_cell))) {
  stop("Score file must contain: ", paste(req_score, collapse = ", "))
}
score_cell[, cell := as.character(cell)]
score_cell[, program_name := as.character(program_name)]
score_cell[, broad_class := as.character(broad_class)]
score_cell[, individual := as.character(individual)]

# remove possibly stale covariate columns before merge
drop_cols <- intersect(c("diagnosis","is_asd","region","age","sex","PMI","RIN"), names(score_cell))
if (length(drop_cols) > 0) {
  log_msg("Dropping pre-existing columns from score file before merge: ", paste(drop_cols, collapse = ", "))
  score_cell[, (drop_cols) := NULL]
}

log_msg("Reading Velmeshev meta: ", META_FILE)
meta <- fread(META_FILE)
req_meta <- c("cell","individual","region","age","sex","diagnosis","post-mortem interval (hours)","RNA Integrity Number")
if (!all(req_meta %in% names(meta))) {
  stop("Meta file must contain: ", paste(req_meta, collapse = ", "))
}

meta[, cell := as.character(cell)]
meta[, individual := as.character(individual)]
meta[, region := as.character(region)]
meta[, age := suppressWarnings(as.numeric(age))]
meta[, sex := as.character(sex)]
meta[, diagnosis := as.character(diagnosis)]
meta[, PMI := suppressWarnings(as.numeric(`post-mortem interval (hours)`))]
meta[, RIN := suppressWarnings(as.numeric(`RNA Integrity Number`))]
meta <- meta[diagnosis %in% c("ASD", "Control")]
meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
meta[, is_asd := diagnosis == "ASD"]

meta_keep <- unique(meta[, .(cell, individual, diagnosis, is_asd, region, age, sex, PMI, RIN)])

score_cov <- merge(score_cell, meta_keep, by = c("cell","individual"), all.x = TRUE)
safe_fwrite(score_cov, file.path(TABDIR, "56_Velmeshev_cell_program_scores_withCov.tsv.gz"))

donor_broad <- score_cov[, .(
  n_cells = uniqueN(cell),
  mean_score = mean(score, na.rm = TRUE)
), by = .(
  program_name, broad_class, individual, diagnosis, is_asd, region, age, sex, PMI, RIN
)]
safe_fwrite(donor_broad, file.path(TABDIR, "56b_Velmeshev_donorRegion_broadClass_program_scores.tsv.gz"))

consistency <- donor_broad[, .(
  n_rows = .N,
  n_regions = uniqueN(region),
  n_age = uniqueN(age),
  n_sex = uniqueN(sex),
  n_pmi = uniqueN(PMI),
  n_rin = uniqueN(RIN),
  n_dx = uniqueN(diagnosis)
), by = individual]
safe_fwrite(consistency, file.path(METADIR, "56c_Velmeshev_donor_covariate_consistency.tsv"))

fit_reduced_subset <- function(dd) {
  dd <- copy(dd)
  dd <- dd[complete.cases(mean_score, diagnosis, region, age, sex)]
  if (nrow(dd) < 8) {
    return(list(
      dd = dd, formula_used = NA_character_, beta = NA_real_, p = NA_real_,
      coef_name = NA_character_, wilcox_p = NA_real_, auc = NA_real_
    ))
  }
  dd[, diagnosis := factor(as.character(diagnosis), levels = c("Control", "ASD"))]
  dd[, region := factor(region)]
  dd[, sex := factor(sex)]

  vars <- c("diagnosis")
  if (uniqueN(dd$region) > 1) vars <- c(vars, "region")
  if (sum(!is.na(dd$age)) > 1 && uniqueN(dd$age) > 1) vars <- c(vars, "age")
  if (uniqueN(dd$sex) > 1) vars <- c(vars, "sex")

  formula_used <- paste("mean_score ~", paste(vars, collapse = " + "))
  wt <- wilcox.test(mean_score ~ diagnosis, data = dd, exact = FALSE)
  fit <- lm(as.formula(formula_used), data = dd)
  ca <- coef_asd(fit)

  list(
    dd = dd,
    formula_used = formula_used,
    beta = ca$beta,
    p = ca$p,
    coef_name = ca$coef_name,
    wilcox_p = wt$p.value,
    auc = calc_auc(dd$mean_score, dd$is_asd)
  )
}

results_list <- list()
programs <- sort(unique(donor_broad$program_name))
bclasses <- sort(unique(donor_broad$broad_class))

for (pn in programs) {
  for (bc in bclasses) {
    dd_all <- donor_broad[program_name == pn & broad_class == bc]
    if (nrow(dd_all) < 8) next
    if (sum(dd_all$diagnosis == "ASD") < 3 || sum(dd_all$diagnosis == "Control") < 3) next

    fitobj <- fit_reduced_subset(dd_all)
    dd <- fitobj$dd

    n_complete <- nrow(dd)
    n_asd_complete <- sum(dd$diagnosis == "ASD")
    n_ctl_complete <- sum(dd$diagnosis == "Control")

    if (n_complete < 8 || n_asd_complete < 3 || n_ctl_complete < 3) next

    results_list[[paste0(pn, "__", bc)]] <- data.table(
      program_name = pn,
      broad_class = bc,
      n_donorRegion_total = nrow(dd_all),
      n_asd_total = sum(dd_all$diagnosis == "ASD"),
      n_control_total = sum(dd_all$diagnosis == "Control"),
      n_complete_model = n_complete,
      n_asd_complete = n_asd_complete,
      n_control_complete = n_ctl_complete,
      n_regions_complete = uniqueN(dd$region),
      n_sex_levels_complete = uniqueN(dd$sex),
      formula_used = fitobj$formula_used,
      mean_score_control_modelSet = dd[diagnosis == "Control", mean(mean_score)],
      mean_score_asd_modelSet = dd[diagnosis == "ASD", mean(mean_score)],
      delta_asd_minus_control_modelSet = dd[diagnosis == "ASD", mean(mean_score)] - dd[diagnosis == "Control", mean(mean_score)],
      auc_modelSet = fitobj$auc,
      wilcox_p_modelSet = fitobj$wilcox_p,
      lm_beta_asd_reduced = fitobj$beta,
      lm_p_reduced = fitobj$p,
      lm_coef_name_reduced = fitobj$coef_name
    )
  }
}

res <- if (length(results_list) > 0) rbindlist(results_list, use.names = TRUE, fill = TRUE) else empty_dt(c(
  "program_name","broad_class","n_donorRegion_total","n_asd_total","n_control_total",
  "n_complete_model","n_asd_complete","n_control_complete","n_regions_complete","n_sex_levels_complete",
  "formula_used","mean_score_control_modelSet","mean_score_asd_modelSet","delta_asd_minus_control_modelSet",
  "auc_modelSet","wilcox_p_modelSet","lm_beta_asd_reduced","lm_p_reduced","lm_coef_name_reduced"
))

if (nrow(res) > 0) {
  res[, fdr_lm_reduced := p.adjust(lm_p_reduced, method = "BH")]
  res[, fdr_wilcox_modelSet := p.adjust(wilcox_p_modelSet, method = "BH")]
  res[, abs_delta_modelSet := abs(delta_asd_minus_control_modelSet)]
  setorder(res, fdr_lm_reduced, fdr_wilcox_modelSet, -abs_delta_modelSet)
}
safe_fwrite(res, file.path(TABDIR, "57_Velmeshev_broadClass_dx_tests_reduced.tsv"))

best_res <- if (nrow(res) > 0) res[1] else NULL

run_summary <- rbindlist(list(
  data.table(section = "Velmeshev", metric = "n_cells_long", value = nrow(score_cov)),
  data.table(section = "Velmeshev", metric = "n_unique_cells", value = uniqueN(score_cov$cell)),
  data.table(section = "Velmeshev", metric = "n_donorRegion_rows", value = nrow(donor_broad)),
  data.table(section = "Velmeshev", metric = "n_individuals", value = uniqueN(donor_broad$individual)),
  data.table(section = "Velmeshev", metric = "n_broad_classes", value = uniqueN(donor_broad$broad_class)),
  data.table(section = "Velmeshev", metric = "n_programs", value = uniqueN(donor_broad$program_name))
), fill = TRUE)

if (!is.null(best_res)) {
  run_summary <- rbind(run_summary, rbindlist(list(
    data.table(section = "BestReduced", metric = "program_name", value = best_res$program_name),
    data.table(section = "BestReduced", metric = "broad_class", value = best_res$broad_class),
    data.table(section = "BestReduced", metric = "delta_asd_minus_control_modelSet", value = best_res$delta_asd_minus_control_modelSet),
    data.table(section = "BestReduced", metric = "lm_p_reduced", value = best_res$lm_p_reduced),
    data.table(section = "BestReduced", metric = "fdr_lm_reduced", value = best_res$fdr_lm_reduced),
    data.table(section = "BestReduced", metric = "formula_used", value = best_res$formula_used)
  )), fill = TRUE)
}

safe_fwrite(run_summary, file.path(METADIR, "58_Step04C_run_summary.tsv"))

log_msg("Step04C completed successfully.")
