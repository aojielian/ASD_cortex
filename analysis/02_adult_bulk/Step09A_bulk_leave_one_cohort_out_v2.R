#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOGDIR, "Step09A_bulk_leave_one_cohort_out_v2.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOGFILE, append = TRUE)
}
safe_fwrite <- function(x, file) fwrite(x, file, sep="\t", quote=FALSE, na="NA")

BULKUPD <- file.path(TABDIR, "126_bulk_threeCohort_program_effects_updated.tsv")
SELUPD <- file.path(TABDIR, "131_final_program_selection_updated.tsv")

bulk <- fread(BULKUPD)
all_cohorts <- sort(unique(bulk$cohort))
pairs <- combn(all_cohorts, 2, simplify = FALSE)

loo_meta <- rbindlist(lapply(pairs, function(cc) {
  x <- bulk[cohort %in% cc]
  x[, {
    z <- aligned_z[!is.na(aligned_z)]
    n <- n_samples[!is.na(aligned_z)]
    uw_z <- if (length(z) > 0) sum(z) / sqrt(length(z)) else NA_real_
    wt_z <- if (length(z) > 0) sum(z * sqrt(n)) / sqrt(sum(n)) else NA_real_
    .(
      cohort_pair = paste(cc, collapse = "__"),
      n_valid_z = length(z),
      n_direction_matched = sum(direction_match, na.rm = TRUE),
      n_direction_mismatched = sum(!direction_match, na.rm = TRUE),
      mean_delta = mean(delta_asd_minus_control, na.rm = TRUE),
      stouffer_z_unweighted = uw_z,
      stouffer_p_unweighted_one_sided = ifelse(is.na(uw_z), NA_real_, 1 - pnorm(uw_z)),
      stouffer_z_weighted = wt_z,
      stouffer_p_weighted_one_sided = ifelse(is.na(wt_z), NA_real_, 1 - pnorm(wt_z))
    )
  }, by = .(program_name, adult_expected_direction)]
}), use.names = TRUE)

loo_meta[, fdr_weighted_withinPair := p.adjust(stouffer_p_weighted_one_sided, method = "BH"), by = cohort_pair]
loo_meta[, fdr_unweighted_withinPair := p.adjust(stouffer_p_unweighted_one_sided, method = "BH"), by = cohort_pair]
setorder(loo_meta, cohort_pair, fdr_weighted_withinPair, fdr_unweighted_withinPair)
safe_fwrite(loo_meta, file.path(TABDIR, "134_bulk_leaveOneOut_meta.tsv"))

sel <- fread(SELUPD)
anchor_summary <- merge(sel, loo_meta, by = "program_name", all.x = TRUE)
role_order <- c("primary_broad_anchor","primary_developmental_context_program","supplementary_sensitivity_program","supplementary_directional_contrast_program")
anchor_summary[, order_role := match(selection_role, role_order)]
anchor_summary[is.na(order_role), order_role := 999L]
setorder(anchor_summary, order_role, cohort_pair, program_name)
anchor_summary[, order_role := NULL]
safe_fwrite(anchor_summary, file.path(TABDIR, "135_bulk_leaveOneOut_anchor_summary.tsv"))

sfari_ok <- loo_meta[program_name == "SFARI_all", all(n_direction_mismatched == 0, na.rm = TRUE)]
run_summary <- rbindlist(list(
  data.table(section="LeaveOneOut", metric="n_pairs", value=length(pairs)),
  data.table(section="LeaveOneOut", metric="cohort_pairs", value=paste(vapply(pairs, paste, collapse="__", character(1)), collapse=",")),
  data.table(section="SFARI_all", metric="all_pairs_direction_matched", value=sfari_ok),
  data.table(section="SFARI_all", metric="best_pair_weighted_p", value=min(loo_meta[program_name=="SFARI_all", stouffer_p_weighted_one_sided], na.rm = TRUE)),
  data.table(section="midPrenatal_SFARI_top20", metric="best_pair_weighted_p", value=min(loo_meta[program_name=="midPrenatal_SFARI_top20", stouffer_p_weighted_one_sided], na.rm = TRUE))
))
safe_fwrite(run_summary, file.path(METADIR, "136_Step09A_run_summary.tsv"))

log_msg("Step09A v2 completed successfully.")
