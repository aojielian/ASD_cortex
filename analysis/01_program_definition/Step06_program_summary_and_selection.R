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

PROGRAM_SUMMARY_FILE <- file.path(METADIR, "20b_midPrenatal_core_program_summary.tsv")
PROGRAM_GENE_FILE    <- file.path(TABDIR, "20_midPrenatal_core_programs.tsv")
GSE_FILE             <- file.path(TABDIR, "24_GSE102741_program_dx_tests.tsv")
GANDAL_FILE          <- file.path(TABDIR, "31_Gandal_program_dx_tests.tsv")
VEL_FILE             <- file.path(TABDIR, "57_Velmeshev_broadClass_dx_tests_reduced.tsv")
PSY_FILE             <- file.path(TABDIR, "81_PsychENCODE_annotation_dx_tests_adjusted.tsv")
STEP01C_SUMMARY      <- file.path(METADIR, "18_Step01C_run_summary.tsv")

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step06_program_summary_and_selection.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))
}

neglog10 <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  p[p <= 0] <- .Machine$double.xmin
  -log10(p)
}

# -----------------------------
# Validate files
# -----------------------------
need_any_program <- file.exists(PROGRAM_SUMMARY_FILE) || file.exists(PROGRAM_GENE_FILE)
if (!need_any_program) stop("Need either program summary or gene file.")
for (f in c(GSE_FILE, GANDAL_FILE, VEL_FILE, PSY_FILE)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

# -----------------------------
# Program metadata
# -----------------------------
if (file.exists(PROGRAM_SUMMARY_FILE)) {
  prog_meta <- fread(PROGRAM_SUMMARY_FILE)
  setnames(prog_meta, names(prog_meta), tolower(names(prog_meta)))
  if (!all(c("program_name","n_genes") %in% names(prog_meta))) {
    stop("Program summary file must contain program_name and n_genes")
  }
  prog_meta <- prog_meta[, .(program_name, n_genes)]
} else {
  pg <- fread(PROGRAM_GENE_FILE)
  if (!all(c("program_name","gene_symbol") %in% names(pg))) {
    stop("Program gene file must contain program_name and gene_symbol")
  }
  pg[, gene_symbol := toupper(trimws(gene_symbol))]
  prog_meta <- pg[!is.na(gene_symbol) & gene_symbol != "", .(n_genes = uniqueN(gene_symbol)), by = program_name]
}

prog_meta[, program_category := fifelse(
  program_name == "SFARI_all", "broad_sfari",
  fifelse(grepl("top05", program_name, ignore.case = TRUE), "focused_top05",
          fifelse(grepl("top10", program_name, ignore.case = TRUE), "focused_top10",
                  fifelse(grepl("top20", program_name, ignore.case = TRUE), "focused_top20", "other")))
)]
prog_meta[, developmental_concentration_rank := fifelse(
  program_category == "focused_top05", 1L,
  fifelse(program_category == "focused_top10", 2L,
          fifelse(program_category == "focused_top20", 3L,
                  fifelse(program_category == "broad_sfari", 99L, 50L)))
)]

# -----------------------------
# Read layer results
# -----------------------------
gse <- fread(GSE_FILE)
gandal <- fread(GANDAL_FILE)
vel <- fread(VEL_FILE)
psy <- fread(PSY_FILE)

# Step01C global summary (optional)
step01c_note <- NULL
if (file.exists(STEP01C_SUMMARY)) {
  tmp <- fread(STEP01C_SUMMARY)
  if (all(c("section","metric","value") %in% names(tmp))) {
    step01c_note <- tmp
  }
}

# -----------------------------
# Helpers
# -----------------------------
aligned_subset <- function(dt, sign_col, expected_sign, p_col) {
  if (nrow(dt) == 0 || is.na(expected_sign) || expected_sign == 0) return(dt[0])
  dt2 <- dt[sign(get(sign_col)) == expected_sign]
  if (nrow(dt2) == 0) return(dt2)
  setorderv(dt2, cols = p_col, order = 1L, na.last = TRUE)
  dt2
}

pick_first_row <- function(dt) {
  if (nrow(dt) == 0) return(NULL)
  as.list(dt[1])
}

# GSE support p: use aligned direction and best available aligned evidence
get_gse_support <- function(row, expected_sign) {
  if (is.null(row)) return(list(match = NA, support_p = NA_real_, support_type = NA_character_))
  delta <- as.numeric(row$delta_asd_minus_control)
  lm_p <- as.numeric(row$lm_p)
  emp_less <- if ("emp_p_less" %in% names(row)) as.numeric(row$emp_p_less) else NA_real_
  emp_greater <- if ("emp_p_greater" %in% names(row)) as.numeric(row$emp_p_greater) else NA_real_

  match <- !is.na(delta) && sign(delta) == expected_sign

  if (is.na(expected_sign) || expected_sign == 0) {
    return(list(match = NA, support_p = NA_real_, support_type = NA_character_))
  }

  if (expected_sign < 0) {
    cand <- c(
      if (!is.na(delta) && delta < 0) lm_p else NA_real_,
      emp_less
    )
    nm <- c("lm_aligned", "emp_less")
  } else {
    cand <- c(
      if (!is.na(delta) && delta > 0) lm_p else NA_real_,
      emp_greater
    )
    nm <- c("lm_aligned", "emp_greater")
  }

  if (all(is.na(cand))) {
    return(list(match = match, support_p = NA_real_, support_type = NA_character_))
  }
  idx <- which.min(replace(cand, is.na(cand), Inf))
  list(match = match, support_p = cand[idx], support_type = nm[idx])
}

# -----------------------------
# Cross-layer summary by program
# -----------------------------
summary_list <- list()

for (i in seq_len(nrow(prog_meta))) {
  pn <- prog_meta$program_name[i]
  n_genes <- prog_meta$n_genes[i]
  pcat <- prog_meta$program_category[i]
  crank <- prog_meta$developmental_concentration_rank[i]

  gse_row <- gse[program_name == pn]
  gse_row <- if (nrow(gse_row) > 0) gse_row[1] else gse[0]

  gandal_row <- gandal[program_name == pn]
  gandal_row <- if (nrow(gandal_row) > 0) gandal_row[1] else gandal[0]

  if (nrow(gandal_row) == 0) {
    log_msg("Program missing in Gandal table: ", pn)
    next
  }

  gandal_delta <- as.numeric(gandal_row$delta_asd_minus_control)
  expected_sign <- sign(gandal_delta)
  if (is.na(expected_sign) || expected_sign == 0) expected_sign <- NA_integer_

  vel_sub <- vel[program_name == pn]
  vel_aligned <- aligned_subset(vel_sub, "delta_asd_minus_control_modelSet", expected_sign, "lm_p_reduced")
  vel_best <- pick_first_row(vel_aligned)

  psy_sub <- psy[program_name == pn]
  psy_aligned <- aligned_subset(psy_sub, "delta_asd_minus_control_adjustedSet", expected_sign, "lm_p_adjusted")
  psy_best <- pick_first_row(psy_aligned)

  gse_support <- get_gse_support(
    if (nrow(gse_row) > 0) as.list(gse_row[1]) else NULL,
    expected_sign
  )

  score_gandal <- neglog10(as.numeric(gandal_row$lm_region_p))
  score_gse <- ifelse(!is.na(gse_support$support_p), neglog10(gse_support$support_p), 0)
  score_vel <- ifelse(!is.null(vel_best), neglog10(as.numeric(vel_best$lm_p_reduced)), 0)
  score_psy <- ifelse(!is.null(psy_best), neglog10(as.numeric(psy_best$lm_p_adjusted)), 0)

  # heuristic: adult bulk strongest, discovery bulk next, Vel stronger than Psych for current project state
  total_score <- 2.5 * score_gandal + 1.5 * score_gse + 1.0 * score_vel + 0.5 * score_psy

  summary_list[[pn]] <- data.table(
    program_name = pn,
    n_genes = n_genes,
    program_category = pcat,
    developmental_concentration_rank = crank,
    adult_expected_direction = fifelse(
      is.na(expected_sign), "NA",
      fifelse(expected_sign < 0, "ASD_lower_than_control", "ASD_higher_than_control")
    ),

    gse_delta = if (nrow(gse_row) > 0) as.numeric(gse_row$delta_asd_minus_control) else NA_real_,
    gse_lm_p = if (nrow(gse_row) > 0) as.numeric(gse_row$lm_p) else NA_real_,
    gse_fdr_lm = if (nrow(gse_row) > 0 && "fdr_lm" %in% names(gse_row)) as.numeric(gse_row$fdr_lm) else NA_real_,
    gse_support_p = gse_support$support_p,
    gse_support_type = gse_support$support_type,
    gse_direction_match_gandal = gse_support$match,

    gandal_delta = gandal_delta,
    gandal_lm_region_beta = as.numeric(gandal_row$lm_region_beta_asd),
    gandal_lm_region_p = as.numeric(gandal_row$lm_region_p),
    gandal_fdr_lm_region = if ("fdr_lm_region" %in% names(gandal_row)) as.numeric(gandal_row$fdr_lm_region) else NA_real_,

    vel_best_aligned_broad_class = if (!is.null(vel_best)) as.character(vel_best$broad_class) else NA_character_,
    vel_best_aligned_delta = if (!is.null(vel_best)) as.numeric(vel_best$delta_asd_minus_control_modelSet) else NA_real_,
    vel_best_aligned_lm_p = if (!is.null(vel_best)) as.numeric(vel_best$lm_p_reduced) else NA_real_,
    vel_best_aligned_fdr = if (!is.null(vel_best)) as.numeric(vel_best$fdr_lm_reduced) else NA_real_,

    psych_best_aligned_annotation = if (!is.null(psy_best)) as.character(psy_best$annotation) else NA_character_,
    psych_best_aligned_delta = if (!is.null(psy_best)) as.numeric(psy_best$delta_asd_minus_control_adjustedSet) else NA_real_,
    psych_best_aligned_lm_p = if (!is.null(psy_best)) as.numeric(psy_best$lm_p_adjusted) else NA_real_,
    psych_best_aligned_fdr = if (!is.null(psy_best)) as.numeric(psy_best$fdr_lm_adjusted) else NA_real_,

    score_gandal = score_gandal,
    score_gse = score_gse,
    score_vel = score_vel,
    score_psych = score_psy,
    total_evidence_score = total_score
  )
}

summary_dt <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
setorder(summary_dt, -total_evidence_score, developmental_concentration_rank, -n_genes)
safe_fwrite(summary_dt, file.path(TABDIR, "90_crossLayer_program_summary.tsv"))

# -----------------------------
# Final selection
# -----------------------------
# broad anchor: prefer SFARI_all if available; else highest score overall
broad_candidates <- summary_dt[program_category == "broad_sfari"]
if (nrow(broad_candidates) == 0) broad_candidates <- summary_dt
setorder(broad_candidates, -total_evidence_score, -score_gandal, -score_gse)
broad_choice <- broad_candidates[1]

# developmental focused: non-broad with best total score; tie-break smaller concentration rank
focus_candidates <- summary_dt[program_category %in% c("focused_top05","focused_top10","focused_top20")]
setorder(focus_candidates, -total_evidence_score, developmental_concentration_rank, -score_vel, -score_gse)
focus_choice <- focus_candidates[1]

# sensitivity / supplementary: next best focused program excluding chosen one
sens_candidates <- focus_candidates[program_name != focus_choice$program_name]
setorder(sens_candidates, -total_evidence_score, developmental_concentration_rank, -score_vel, -score_gse)
sens_choice <- if (nrow(sens_candidates) > 0) sens_candidates[1] else focus_choice

selection_dt <- rbindlist(list(
  data.table(
    selection_role = "primary_broad_anchor",
    program_name = broad_choice$program_name,
    rationale = sprintf(
      "Best broad program by cross-layer evidence; strongest adult bulk layer support in Gandal (delta=%s, lm_region_p=%s, FDR=%s).",
      fmt_num(broad_choice$gandal_delta, 4),
      fmt_num(broad_choice$gandal_lm_region_p, 3),
      fmt_num(broad_choice$gandal_fdr_lm_region, 3)
    ),
    total_evidence_score = broad_choice$total_evidence_score
  ),
  data.table(
    selection_role = "primary_developmental_focused_anchor",
    program_name = focus_choice$program_name,
    rationale = sprintf(
      "Best focused mid-prenatal program by cross-layer evidence; balances developmental concentration with adult-layer support (GSE support p=%s; best aligned Velmeshev p=%s in %s).",
      fmt_num(focus_choice$gse_support_p, 3),
      fmt_num(focus_choice$vel_best_aligned_lm_p, 3),
      ifelse(is.na(focus_choice$vel_best_aligned_broad_class), "NA", focus_choice$vel_best_aligned_broad_class)
    ),
    total_evidence_score = focus_choice$total_evidence_score
  ),
  data.table(
    selection_role = "supplementary_sensitivity_program",
    program_name = sens_choice$program_name,
    rationale = sprintf(
      "Reserve as sensitivity / supplementary developmental program; useful to show program-size robustness relative to %s.",
      focus_choice$program_name
    ),
    total_evidence_score = sens_choice$total_evidence_score
  )
), use.names = TRUE, fill = TRUE)

safe_fwrite(selection_dt, file.path(TABDIR, "91_final_program_selection.tsv"))

# -----------------------------
# Key message text
# -----------------------------
lines <- c(
  "Cross-layer summary for ASD cortex project",
  "=========================================",
  "",
  sprintf("Primary broad anchor: %s", broad_choice$program_name),
  sprintf("  - Gandal delta = %s ; lm_region_p = %s ; FDR = %s",
          fmt_num(broad_choice$gandal_delta, 4),
          fmt_num(broad_choice$gandal_lm_region_p, 3),
          fmt_num(broad_choice$gandal_fdr_lm_region, 3)),
  sprintf("  - GSE support p = %s ; direction match = %s",
          fmt_num(broad_choice$gse_support_p, 3),
          ifelse(is.na(broad_choice$gse_direction_match_gandal), "NA", as.character(broad_choice$gse_direction_match_gandal))),
  "",
  sprintf("Primary developmental-focused anchor: %s", focus_choice$program_name),
  sprintf("  - Gandal delta = %s ; lm_region_p = %s ; FDR = %s",
          fmt_num(focus_choice$gandal_delta, 4),
          fmt_num(focus_choice$gandal_lm_region_p, 3),
          fmt_num(focus_choice$gandal_fdr_lm_region, 3)),
  sprintf("  - GSE support p = %s ; support type = %s",
          fmt_num(focus_choice$gse_support_p, 3),
          ifelse(is.na(focus_choice$gse_support_type), "NA", focus_choice$gse_support_type)),
  sprintf("  - Best aligned Velmeshev class = %s ; p = %s",
          ifelse(is.na(focus_choice$vel_best_aligned_broad_class), "NA", focus_choice$vel_best_aligned_broad_class),
          fmt_num(focus_choice$vel_best_aligned_lm_p, 3)),
  "",
  sprintf("Suggested supplementary sensitivity program: %s", sens_choice$program_name),
  "",
  "Current interpretation:",
  "  1) Developmental anchoring is the strongest upstream layer.",
  "  2) Adult bulk dysregulation is the strongest adult disease layer.",
  "  3) Single-cell layers are better treated as cell-context / direction-support rather than robust diagnosis-replication layers.",
  ""
)

if (!is.null(step01c_note)) {
  best_stage <- step01c_note[metric == "best_stage", value][1]
  best_stage_auc <- step01c_note[metric == "best_stage_auc", value][1]
  best_stage_fdr <- step01c_note[metric == "best_stage_fdr", value][1]
  lines <- c(lines,
             "Step01C global developmental note:",
             sprintf("  - best_stage = %s ; auc = %s ; fdr = %s",
                     ifelse(length(best_stage) == 0, "NA", best_stage),
                     ifelse(length(best_stage_auc) == 0, "NA", best_stage_auc),
                     ifelse(length(best_stage_fdr) == 0, "NA", best_stage_fdr)))
}

writeLines(lines, file.path(METADIR, "92_program_selection_key_messages.txt"))

# manuscript-ready compact table
manu_dt <- summary_dt[, .(
  program_name,
  n_genes,
  adult_expected_direction,
  gse_support_p,
  gandal_lm_region_p,
  gandal_fdr_lm_region,
  vel_best_aligned_broad_class,
  vel_best_aligned_lm_p,
  psych_best_aligned_annotation,
  psych_best_aligned_lm_p,
  total_evidence_score
)]
safe_fwrite(manu_dt, file.path(TABDIR, "93_program_summary_for_manuscript.tsv"))

run_summary <- rbindlist(list(
  data.table(section = "Programs", metric = "n_programs_summarized", value = nrow(summary_dt)),
  data.table(section = "Selection", metric = "primary_broad_anchor", value = broad_choice$program_name),
  data.table(section = "Selection", metric = "primary_developmental_focused_anchor", value = focus_choice$program_name),
  data.table(section = "Selection", metric = "supplementary_sensitivity_program", value = sens_choice$program_name)
), use.names = TRUE, fill = TRUE)

safe_fwrite(run_summary, file.path(METADIR, "94_Step06_run_summary.tsv"))

log_msg("Step06 completed successfully.")
