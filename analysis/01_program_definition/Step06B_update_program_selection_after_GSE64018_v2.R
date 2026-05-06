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

LOGFILE <- file.path(LOGDIR, "Step06B_update_program_selection_after_GSE64018_v2.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOGFILE, append = TRUE)
}
safe_fwrite <- function(x, file) fwrite(x, file, sep="\t", quote=FALSE, na="NA")

CROSSLAYER <- file.path(TABDIR, "90_crossLayer_program_summary.tsv")
METAUPD <- file.path(TABDIR, "127_bulk_program_meta_summary_updated.tsv")
BULKUPD <- file.path(TABDIR, "126_bulk_threeCohort_program_effects_updated.tsv")
PROGFILE <- file.path(TABDIR, "20_midPrenatal_core_programs.tsv")

log_msg("Reading cross-layer summary: ", CROSSLAYER)
cross <- fread(CROSSLAYER)
log_msg("Reading updated bulk meta summary: ", METAUPD)
meta <- fread(METAUPD)
log_msg("Reading updated bulk cohort effects: ", BULKUPD)
bulk <- fread(BULKUPD)
prog <- fread(PROGFILE)

nms <- tolower(names(prog))
pcol <- names(prog)[match("program_name", nms)]
g_candidates <- c("gene_symbol","symbol","gene","standardized_gene_symbol")
g_hit <- g_candidates[g_candidates %in% nms]
if (length(g_hit) == 0) stop("No gene symbol column found in program file.")
gcol <- names(prog)[match(g_hit[1], nms)]

gene_counts <- unique(prog[, .(program_name = get(pcol), gene_symbol = toupper(trimws(get(gcol))))])
gene_counts <- gene_counts[!is.na(gene_symbol) & gene_symbol != "", .(n_genes = uniqueN(gene_symbol)), by = program_name]

summary_dt <- merge(
  cross[, .(program_name, program_category, developmental_concentration_rank, adult_expected_direction)],
  meta,
  by = c("program_name","adult_expected_direction"),
  all = TRUE
)
summary_dt <- merge(summary_dt, gene_counts, by = "program_name", all.x = TRUE)

cohort_wide <- dcast(
  bulk[, .(program_name, cohort, delta_asd_minus_control, lm_p, direction_match)],
  program_name ~ cohort,
  value.var = c("delta_asd_minus_control", "lm_p", "direction_match")
)
summary_dt <- merge(summary_dt, cohort_wide, by = "program_name", all.x = TRUE)

selection <- rbindlist(list(
  data.table(
    selection_role = "primary_broad_anchor",
    program_name = "SFARI_all",
    rationale = paste0(
      "Best broad adult anchor after formal three-cohort bulk integration; 3/3 cohorts direction-matched, weighted Stouffer p=",
      format(summary_dt[program_name=="SFARI_all", stouffer_p_weighted_one_sided], scientific = TRUE, digits = 3),
      ", FDR=",
      format(summary_dt[program_name=="SFARI_all", fdr_weighted], scientific = TRUE, digits = 3), "."
    )
  ),
  data.table(
    selection_role = "primary_developmental_context_program",
    program_name = "midPrenatal_SFARI_top20",
    rationale = paste0(
      "Preferred focused developmental program for the manuscript narrative; strong prenatal concentration but mixed adult bulk behavior (",
      summary_dt[program_name=="midPrenatal_SFARI_top20", n_direction_matched], "/",
      summary_dt[program_name=="midPrenatal_SFARI_top20", n_valid_z],
      " cohorts direction-matched). Treat as developmental-context program rather than adult replicated anchor."
    )
  ),
  data.table(
    selection_role = "supplementary_sensitivity_program",
    program_name = "midPrenatal_SFARI_top05",
    rationale = "Smaller nested prenatal-focused sensitivity set; useful to show robustness / instability across tighter program definitions."
  ),
  data.table(
    selection_role = "supplementary_directional_contrast_program",
    program_name = "midPrenatal_SFARI_top10",
    rationale = "Useful contrast program because adult-direction expectation differs from other prenatal-focused sets and remains less stable across cohorts."
  )
))
safe_fwrite(selection, file.path(TABDIR, "131_final_program_selection_updated.tsv"))

summary_out <- merge(summary_dt, selection[, .(program_name, selection_role)], by = "program_name", all.x = TRUE)
prog_order <- c("SFARI_all","midPrenatal_SFARI_top20","midPrenatal_SFARI_top05","midPrenatal_SFARI_top10")
summary_out[, order_program := match(program_name, prog_order)]
summary_out[is.na(order_program), order_program := 999L]
setorder(summary_out, order_program, program_name)
summary_out[, order_program := NULL]
safe_fwrite(summary_out, file.path(TABDIR, "130_program_summary_for_manuscript_updated.tsv"))

key_lines <- c(
  "Updated cross-layer summary after formal GSE64018 remap and three-cohort bulk integration",
  "==========================================================================",
  "",
  "Primary broad anchor: SFARI_all",
  paste0("  - three-cohort direction match = ",
         summary_dt[program_name=="SFARI_all", n_direction_matched], "/",
         summary_dt[program_name=="SFARI_all", n_valid_z]),
  paste0("  - weighted Stouffer p = ",
         format(summary_dt[program_name=="SFARI_all", stouffer_p_weighted_one_sided], scientific = TRUE, digits = 4),
         " ; FDR = ",
         format(summary_dt[program_name=="SFARI_all", fdr_weighted], scientific = TRUE, digits = 4)),
  paste0("  - cohort deltas: GSE102741=",
         signif(summary_dt[program_name=="SFARI_all", delta_asd_minus_control_GSE102741], 4),
         " ; Gandal2022=",
         signif(summary_dt[program_name=="SFARI_all", delta_asd_minus_control_Gandal2022], 4),
         " ; GSE64018=",
         signif(summary_dt[program_name=="SFARI_all", delta_asd_minus_control_GSE64018], 4)),
  "",
  "Developmental-context anchor: midPrenatal_SFARI_top20",
  paste0("  - developmental concentration rank = ",
         summary_dt[program_name=="midPrenatal_SFARI_top20", developmental_concentration_rank]),
  paste0("  - adult bulk direction match = ",
         summary_dt[program_name=="midPrenatal_SFARI_top20", n_direction_matched], "/",
         summary_dt[program_name=="midPrenatal_SFARI_top20", n_valid_z],
         " ; weighted p = ",
         format(summary_dt[program_name=="midPrenatal_SFARI_top20", stouffer_p_weighted_one_sided], scientific = TRUE, digits = 4)),
  "  - interpretation: strong developmental anchoring, weaker and heterogeneous adult dysregulation.",
  "",
  "Supplementary sets",
  "  - midPrenatal_SFARI_top05: tighter prenatal subset for sensitivity",
  "  - midPrenatal_SFARI_top10: directional-contrast comparator",
  "",
  "Recommended manuscript framing:",
  "  1) strongest result = developmental anchoring in mid-prenatal cortex;",
  "  2) strongest adult disease-layer result = broad SFARI_all downregulation across three bulk cohorts;",
  "  3) focused prenatal subsets should be treated as context-rich but not robust adult replicated anchors."
)
writeLines(key_lines, file.path(METADIR, "132_program_selection_key_messages_updated.txt"))

run_summary <- rbindlist(list(
  data.table(section="Programs", metric="n_programs_summarized", value=summary_dt[, .N]),
  data.table(section="Selection", metric="primary_broad_anchor", value="SFARI_all"),
  data.table(section="Selection", metric="primary_developmental_context_program", value="midPrenatal_SFARI_top20"),
  data.table(section="Selection", metric="supplementary_sensitivity_program", value="midPrenatal_SFARI_top05"),
  data.table(section="Selection", metric="supplementary_directional_contrast_program", value="midPrenatal_SFARI_top10")
))
safe_fwrite(run_summary, file.path(METADIR, "133_Step06B_run_summary.tsv"))

log_msg("Step06B v2 completed successfully.")
