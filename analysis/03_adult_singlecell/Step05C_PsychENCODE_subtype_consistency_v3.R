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

SCORE_FILE <- file.path(TABDIR, "71_PsychENCODE_cell_program_scores.tsv.gz")

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step05C_PsychENCODE_subtype_consistency_v3.log")

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

empty_dt <- function(cols) {
  out <- data.table(matrix(ncol = length(cols), nrow = 0))
  setnames(out, cols)
  out
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

# thresholds
MIN_CELLS_TOTAL <- 100L
MIN_DONORS_TOTAL <- 8L
MIN_ASD_DONORS <- 3L
MIN_CONTROL_DONORS <- 3L
MIN_DOMINANT_PURITY <- 0.80

if (!file.exists(SCORE_FILE)) stop("Missing score file: ", SCORE_FILE)

log_msg("Reading score file: ", SCORE_FILE)
score_cell <- fread(SCORE_FILE)

required_cols <- c("Cell_ID","score","program_name","annotation","Subtype","individual_ID","diagnosis","is_asd")
if (!all(required_cols %in% names(score_cell))) {
  stop("Score file must contain: ", paste(required_cols, collapse = ", "))
}

score_cell[, Cell_ID := as.character(Cell_ID)]
score_cell[, annotation := as.character(annotation)]
score_cell[, Subtype := as.character(Subtype)]
score_cell[, individual_ID := as.character(individual_ID)]
score_cell[, diagnosis := factor(as.character(diagnosis), levels = c("Control", "ASD"))]
score_cell <- score_cell[!is.na(diagnosis)]
score_cell[, is_asd := diagnosis == "ASD"]

# =========================================================
# A. subtype consistency MUST be computed on UNIQUE CELLS,
#    not on long-format cell × program rows
# =========================================================
cell_meta_unique <- unique(
  score_cell[, .(Cell_ID, annotation, Subtype, individual_ID, diagnosis, is_asd)]
)

cross_dt <- cell_meta_unique[, .(
  n_cells = .N,
  n_donors = uniqueN(individual_ID),
  n_asd_donors = uniqueN(individual_ID[diagnosis == "ASD"]),
  n_control_donors = uniqueN(individual_ID[diagnosis == "Control"])
), by = .(annotation, Subtype)]

setorder(cross_dt, Subtype, -n_cells, -n_donors)
safe_fwrite(cross_dt, file.path(TABDIR, "76_PsychENCODE_annotation_subtype_crosstab.tsv"))

subtype_totals <- cell_meta_unique[, .(
  n_cells_total = .N,
  n_donors_total = uniqueN(individual_ID),
  n_asd_donors_total = uniqueN(individual_ID[diagnosis == "ASD"]),
  n_control_donors_total = uniqueN(individual_ID[diagnosis == "Control"])
), by = Subtype]

dominant_ann <- cross_dt[order(Subtype, -n_cells, -n_donors), .SD[1], by = Subtype]
setnames(
  dominant_ann,
  old = c("annotation", "n_cells", "n_donors", "n_asd_donors", "n_control_donors"),
  new = c("dominant_annotation", "dominant_n_cells", "dominant_n_donors", "dominant_n_asd_donors", "dominant_n_control_donors")
)

subtype_summary <- merge(subtype_totals, dominant_ann, by = "Subtype", all.x = TRUE)
subtype_summary[, dominant_purity_cells := dominant_n_cells / pmax(1, n_cells_total)]
subtype_summary[, subtype_keep := (
  n_cells_total >= MIN_CELLS_TOTAL &
  n_donors_total >= MIN_DONORS_TOTAL &
  n_asd_donors_total >= MIN_ASD_DONORS &
  n_control_donors_total >= MIN_CONTROL_DONORS &
  dominant_purity_cells >= MIN_DOMINANT_PURITY
)]
subtype_summary[, keep_reason := fifelse(subtype_keep, "PASS", "FILTERED")]
setorder(subtype_summary, -subtype_keep, -dominant_purity_cells, -n_cells_total)

safe_fwrite(subtype_summary, file.path(METADIR, "76b_PsychENCODE_subtype_consistency_summary.tsv"))
safe_fwrite(subtype_summary[subtype_keep == TRUE], file.path(METADIR, "76c_PsychENCODE_kept_subtypes.tsv"))

# =========================================================
# B. Restrict long-format score table using clean subtype map
# =========================================================
kept <- subtype_summary[subtype_keep == TRUE, .(Subtype, dominant_annotation)]
score_clean <- merge(score_cell, kept, by = "Subtype", all = FALSE)
score_clean <- score_clean[annotation == dominant_annotation]

safe_fwrite(score_clean, file.path(TABDIR, "77_PsychENCODE_cell_program_scores_cleanSubtype.tsv.gz"))

# donor-level subtype program scores
donor_subtype_clean <- if (nrow(score_clean) > 0) {
  score_clean[, .(
    n_cells = uniqueN(Cell_ID),
    mean_score = mean(score, na.rm = TRUE)
  ), by = .(program_name, dominant_annotation, Subtype, individual_ID, diagnosis, is_asd)]
} else empty_dt(c("program_name","dominant_annotation","Subtype","individual_ID","diagnosis","is_asd","n_cells","mean_score"))

safe_fwrite(donor_subtype_clean, file.path(TABDIR, "77b_PsychENCODE_donor_subtype_program_scores_clean.tsv.gz"))

subtype_tests_list <- list()
if (nrow(donor_subtype_clean) > 0) {
  keys <- unique(donor_subtype_clean[, .(program_name, dominant_annotation, Subtype)])
  for (i in seq_len(nrow(keys))) {
    pn <- keys$program_name[i]
    ann <- keys$dominant_annotation[i]
    st <- keys$Subtype[i]
    dd <- donor_subtype_clean[program_name == pn & dominant_annotation == ann & Subtype == st]
    dd[, diagnosis := factor(as.character(diagnosis), levels = c("Control", "ASD"))]
    if (nrow(dd) < MIN_DONORS_TOTAL) next
    if (sum(dd$diagnosis == "ASD") < MIN_ASD_DONORS || sum(dd$diagnosis == "Control") < MIN_CONTROL_DONORS) next

    wt <- wilcox.test(mean_score ~ diagnosis, data = dd, exact = FALSE)
    fit <- lm(mean_score ~ diagnosis, data = dd)
    ca <- coef_asd(fit)

    subtype_tests_list[[paste0(pn, "__", st)]] <- data.table(
      program_name = pn,
      annotation = ann,
      Subtype = st,
      n_donors = nrow(dd),
      n_asd = sum(dd$diagnosis == "ASD"),
      n_control = sum(dd$diagnosis == "Control"),
      mean_score_control = dd[diagnosis == "Control", mean(mean_score)],
      mean_score_asd = dd[diagnosis == "ASD", mean(mean_score)],
      delta_asd_minus_control = dd[diagnosis == "ASD", mean(mean_score)] - dd[diagnosis == "Control", mean(mean_score)],
      auc_asd_higher = calc_auc(dd$mean_score, dd$is_asd),
      wilcox_p = wt$p.value,
      lm_beta_asd = ca$beta,
      lm_p = ca$p,
      lm_coef_name = ca$coef_name
    )
  }
}

subtype_tests_clean <- if (length(subtype_tests_list) > 0) rbindlist(subtype_tests_list, use.names = TRUE, fill = TRUE) else empty_dt(c(
  "program_name","annotation","Subtype","n_donors","n_asd","n_control","mean_score_control","mean_score_asd",
  "delta_asd_minus_control","auc_asd_higher","wilcox_p","lm_beta_asd","lm_p","lm_coef_name"
))

if (nrow(subtype_tests_clean) > 0) {
  subtype_tests_clean[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  subtype_tests_clean[, fdr_lm := p.adjust(lm_p, method = "BH")]
  subtype_tests_clean[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(subtype_tests_clean, fdr_lm, fdr_wilcox, -abs_delta_asd_minus_control)
}

safe_fwrite(subtype_tests_clean, file.path(TABDIR, "78_PsychENCODE_subtype_dx_tests_clean.tsv"))

# donor-level annotation program scores
donor_ann_clean <- if (nrow(score_clean) > 0) {
  score_clean[, .(
    n_cells = uniqueN(Cell_ID),
    mean_score = mean(score, na.rm = TRUE)
  ), by = .(program_name, dominant_annotation, individual_ID, diagnosis, is_asd)]
} else empty_dt(c("program_name","dominant_annotation","individual_ID","diagnosis","is_asd","n_cells","mean_score"))

safe_fwrite(donor_ann_clean, file.path(TABDIR, "78b_PsychENCODE_donor_annotation_program_scores_clean.tsv.gz"))

ann_tests_list <- list()
if (nrow(donor_ann_clean) > 0) {
  keys <- unique(donor_ann_clean[, .(program_name, dominant_annotation)])
  for (i in seq_len(nrow(keys))) {
    pn <- keys$program_name[i]
    ann <- keys$dominant_annotation[i]
    dd <- donor_ann_clean[program_name == pn & dominant_annotation == ann]
    dd[, diagnosis := factor(as.character(diagnosis), levels = c("Control", "ASD"))]
    if (nrow(dd) < MIN_DONORS_TOTAL) next
    if (sum(dd$diagnosis == "ASD") < MIN_ASD_DONORS || sum(dd$diagnosis == "Control") < MIN_CONTROL_DONORS) next

    wt <- wilcox.test(mean_score ~ diagnosis, data = dd, exact = FALSE)
    fit <- lm(mean_score ~ diagnosis, data = dd)
    ca <- coef_asd(fit)

    ann_tests_list[[paste0(pn, "__", ann)]] <- data.table(
      program_name = pn,
      annotation = ann,
      n_donors = nrow(dd),
      n_asd = sum(dd$diagnosis == "ASD"),
      n_control = sum(dd$diagnosis == "Control"),
      mean_score_control = dd[diagnosis == "Control", mean(mean_score)],
      mean_score_asd = dd[diagnosis == "ASD", mean(mean_score)],
      delta_asd_minus_control = dd[diagnosis == "ASD", mean(mean_score)] - dd[diagnosis == "Control", mean(mean_score)],
      auc_asd_higher = calc_auc(dd$mean_score, dd$is_asd),
      wilcox_p = wt$p.value,
      lm_beta_asd = ca$beta,
      lm_p = ca$p,
      lm_coef_name = ca$coef_name
    )
  }
}

ann_tests_clean <- if (length(ann_tests_list) > 0) rbindlist(ann_tests_list, use.names = TRUE, fill = TRUE) else empty_dt(c(
  "program_name","annotation","n_donors","n_asd","n_control","mean_score_control","mean_score_asd",
  "delta_asd_minus_control","auc_asd_higher","wilcox_p","lm_beta_asd","lm_p","lm_coef_name"
))

if (nrow(ann_tests_clean) > 0) {
  ann_tests_clean[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  ann_tests_clean[, fdr_lm := p.adjust(lm_p, method = "BH")]
  ann_tests_clean[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(ann_tests_clean, fdr_lm, fdr_wilcox, -abs_delta_asd_minus_control)
}

safe_fwrite(ann_tests_clean, file.path(TABDIR, "78c_PsychENCODE_annotation_dx_tests_clean.tsv"))

best_subtype <- if (nrow(subtype_tests_clean) > 0) subtype_tests_clean[1] else NULL
best_ann <- if (nrow(ann_tests_clean) > 0) ann_tests_clean[1] else NULL

run_summary <- rbindlist(list(
  data.table(section = "PsychENCODE", metric = "n_subtypes_total", value = uniqueN(cell_meta_unique$Subtype)),
  data.table(section = "PsychENCODE", metric = "n_subtypes_kept", value = nrow(subtype_summary[subtype_keep == TRUE])),
  data.table(section = "PsychENCODE", metric = "min_cells_total", value = MIN_CELLS_TOTAL),
  data.table(section = "PsychENCODE", metric = "min_donors_total", value = MIN_DONORS_TOTAL),
  data.table(section = "PsychENCODE", metric = "min_asd_donors", value = MIN_ASD_DONORS),
  data.table(section = "PsychENCODE", metric = "min_control_donors", value = MIN_CONTROL_DONORS),
  data.table(section = "PsychENCODE", metric = "min_dominant_purity", value = MIN_DOMINANT_PURITY),
  data.table(section = "PsychENCODE", metric = "n_clean_unique_cells", value = uniqueN(score_clean$Cell_ID)),
  data.table(section = "PsychENCODE", metric = "n_clean_long_rows", value = nrow(score_clean)),
  data.table(section = "Programs", metric = "n_programs_tested_clean_subtype", value = uniqueN(subtype_tests_clean$program_name)),
  data.table(section = "Programs", metric = "n_programs_tested_clean_annotation", value = uniqueN(ann_tests_clean$program_name))
), fill = TRUE)

if (!is.null(best_ann)) {
  run_summary <- rbind(run_summary, rbindlist(list(
    data.table(section = "BestAnnotationClean", metric = "program_name", value = best_ann$program_name),
    data.table(section = "BestAnnotationClean", metric = "annotation", value = best_ann$annotation),
    data.table(section = "BestAnnotationClean", metric = "delta_asd_minus_control", value = best_ann$delta_asd_minus_control),
    data.table(section = "BestAnnotationClean", metric = "lm_p", value = best_ann$lm_p),
    data.table(section = "BestAnnotationClean", metric = "fdr_lm", value = best_ann$fdr_lm)
  )), fill = TRUE)
}

if (!is.null(best_subtype)) {
  run_summary <- rbind(run_summary, rbindlist(list(
    data.table(section = "BestSubtypeClean", metric = "program_name", value = best_subtype$program_name),
    data.table(section = "BestSubtypeClean", metric = "Subtype", value = best_subtype$Subtype),
    data.table(section = "BestSubtypeClean", metric = "annotation", value = best_subtype$annotation),
    data.table(section = "BestSubtypeClean", metric = "delta_asd_minus_control", value = best_subtype$delta_asd_minus_control),
    data.table(section = "BestSubtypeClean", metric = "lm_p", value = best_subtype$lm_p),
    data.table(section = "BestSubtypeClean", metric = "fdr_lm", value = best_subtype$fdr_lm)
  )), fill = TRUE)
}

safe_fwrite(run_summary, file.path(METADIR, "79_Step05C_run_summary.tsv"))

log_msg("Step05C v3 completed successfully.")
