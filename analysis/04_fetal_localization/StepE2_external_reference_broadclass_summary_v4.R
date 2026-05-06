#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

parse_args <- function(x) {
  out <- list(
    base = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex",
    stepd_summary = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/stepD_external_reference_program_scoring/StepD_all_group_program_summary.tsv.gz",
    rna_annotation = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/stepE1_gse162170_rna_annotation/StepE1_GSE162170_RNA_cluster_annotation.tsv",
    outdir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/stepE2_external_reference_broadclass_summary"
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    val <- if (i < length(x)) x[[i + 1L]] else NA_character_
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    out[[key]] <- val
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
outdir <- args$outdir
logdir <- file.path(outdir, "logs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logdir, "StepE2_external_reference_broadclass_summary.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
stop_if_missing <- function(path) if (!file.exists(path)) stop("Missing required file: ", path)

map_nowakowski <- function(x) {
  y <- as.character(x)
  out <- rep("Other", length(y))
  out[grepl("EC|Peric|Endo|Vasc", y, ignore.case = TRUE)] <- "Vascular"
  out[grepl("OPC|Olig", y, ignore.case = TRUE)] <- "OPC_oligo"
  out[grepl("IPC|Cyc", y, ignore.case = TRUE)] <- "IPC_prog"
  out[grepl("^tRG$|^vRG$|RG-div|RG-early|^RG$", y)] <- "RG_like"
  out[grepl("nEN|EN-|IPC-nEN", y)] <- "EN_neurogenic"
  out[grepl("nIN|IN-|MGE|CGE|IN-STR", y)] <- "IN"
  out[grepl("Glyc|AST|MG", y)] <- "Glia_other"
  out
}

map_multiome <- function(x) {
  y <- as.character(x)
  out <- rep("Other", length(y))
  out[grepl("^EC/Peric\\.$|^EC/Peric$", y)] <- "Vascular"
  out[grepl("mGPC/OPC", y)] <- "OPC_oligo"
  out[grepl("^RG$", y)] <- "RG_like"
  out[grepl("nIPC/GluN1|Cyc\\. Prog\\.", y)] <- "IPC_prog"
  out[grepl("^IN", y)] <- "IN"
  out[grepl("^GluN|^SP$", y)] <- "EN_neurogenic"
  out
}

axis_order <- data.table(
  broad_class = c("RG_like", "IPC_prog", "EN_neurogenic", "IN", "OPC_oligo", "Vascular", "Glia_other", "Other"),
  axis_rank = c(1, 2, 3, 4, 5, 6, 7, 8)
)

for (p in c(args$stepd_summary, args$rna_annotation)) stop_if_missing(p)

log_msg("Reading StepD group-level program summary...")
dt <- fread(args$stepd_summary)
log_msg("Reading StepE1 RNA annotation...")
rna_anno <- fread(args$rna_annotation)

# group -> broad class mapping
map_dt <- unique(dt[, .(dataset, group_var, group_id, group_label)])
map_dt[, broad_class := NA_character_]
map_dt[dataset == "Nowakowski_UCSC", broad_class := map_nowakowski(group_label)]
map_dt[dataset == "GSE162170_multiome", broad_class := map_multiome(group_label)]
map_dt[dataset == "GSE162170_RNA", broad_class := rna_anno$final_broad_class[match(group_id, rna_anno$rna_cluster_id)]]
map_dt[dataset == "GSE162170_RNA", fine_label := rna_anno$final_label[match(group_id, rna_anno$rna_cluster_id)]]
map_dt[dataset != "GSE162170_RNA", fine_label := group_label]
map_dt[is.na(broad_class) | broad_class == "", broad_class := "Other"]
map_dt <- merge(map_dt, axis_order, by = "broad_class", all.x = TRUE)
safe_fwrite(map_dt, file.path(outdir, "StepE2_group_to_broad_class.tsv"))

# collapse group summaries into broad-class summaries
merged <- merge(dt, map_dt[, .(dataset, group_var, group_id, broad_class, axis_rank, fine_label)],
                by = c("dataset", "group_var", "group_id"), all.x = TRUE)
merged[is.na(broad_class) | broad_class == "", broad_class := "Other"]
if (!"axis_rank" %in% names(merged)) {
  merged <- merge(merged, axis_order, by = "broad_class", all.x = TRUE)
}
if ("axis_rank.x" %in% names(merged) || "axis_rank.y" %in% names(merged)) {
  if (!("axis_rank.x" %in% names(merged) && "axis_rank.y" %in% names(merged))) stop("axis_rank columns are inconsistent after merge.")
  merged[, axis_rank := fifelse(!is.na(axis_rank.x), axis_rank.x, axis_rank.y)]
  merged[, c("axis_rank.x", "axis_rank.y") := NULL]
}

if (!"axis_rank" %in% names(merged)) stop("axis_rank not found after broad-class mapping merge.")

broad_sum <- merged[, .(
  n_groups = .N,
  total_cells = sum(n_cells, na.rm = TRUE),
  weighted_mean_score = weighted.mean(mean_score, w = pmax(n_cells, 1), na.rm = TRUE),
  weighted_median_of_group_medians = weighted.mean(median_score, w = pmax(n_cells, 1), na.rm = TRUE),
  min_group_rank = min(within_dataset_rank, na.rm = TRUE),
  top_fine_labels = paste(head(fine_label[order(within_dataset_rank)], 5), collapse = ";")
), by = .(dataset, program, broad_class, axis_rank)]
setorder(broad_sum, dataset, program, -weighted_mean_score, axis_rank)
broad_sum[, broad_rank_within_dataset_program := frank(-weighted_mean_score, ties.method = "min"), by = .(dataset, program)]
broad_sum[, broad_mean_z := as.numeric(scale(weighted_mean_score)), by = .(dataset, program)]
safe_fwrite(broad_sum, file.path(outdir, "StepE2_broadclass_program_summary.tsv"))

# top broad classes per dataset/program
broad_top <- broad_sum[order(dataset, program, broad_rank_within_dataset_program)]
broad_top <- broad_top[broad_rank_within_dataset_program <= 5]
safe_fwrite(broad_top, file.path(outdir, "StepE2_top_broadclasses_per_dataset_program.tsv"))

# cross-dataset consensus by program and broad class
consensus <- broad_sum[, .(
  n_datasets = uniqueN(dataset),
  mean_broad_rank = mean(broad_rank_within_dataset_program, na.rm = TRUE),
  mean_broad_z = mean(broad_mean_z, na.rm = TRUE),
  datasets_supporting_top3 = sum(broad_rank_within_dataset_program <= 3, na.rm = TRUE),
  datasets_supporting_top5 = sum(broad_rank_within_dataset_program <= 5, na.rm = TRUE),
  total_cells_across = sum(total_cells, na.rm = TRUE),
  supporting_datasets = paste(dataset[order(broad_rank_within_dataset_program)], collapse = ";")
), by = .(program, broad_class, axis_rank)]
setorder(consensus, program, mean_broad_rank, -datasets_supporting_top3, axis_rank)
safe_fwrite(consensus, file.path(outdir, "StepE2_cross_dataset_program_consensus.tsv"))

# manuscript-oriented summary: developmental shift relative to SFARI_all
wide_rank <- dcast(broad_sum, broad_class + axis_rank + dataset ~ program, value.var = "broad_rank_within_dataset_program")

get_col_or_na <- function(dt, col) {
  if (col %in% names(dt)) dt[[col]] else rep(NA_real_, nrow(dt))
}

wide_rank[, rank_SFARI_all := get_col_or_na(.SD, "score_SFARI_all")]
wide_rank[, rank_top20 := get_col_or_na(.SD, "score_midPrenatal_SFARI_top20")]
wide_rank[, rank_top10 := get_col_or_na(.SD, "score_midPrenatal_SFARI_top10")]
wide_rank[, rank_top05 := get_col_or_na(.SD, "score_midPrenatal_SFARI_top05")]
wide_rank[, delta_top20_vs_SFARI_all := rank_top20 - rank_SFARI_all]
wide_rank[, delta_top10_vs_top20 := rank_top10 - rank_top20]
wide_rank[, delta_top05_vs_top20 := rank_top05 - rank_top20]
safe_fwrite(wide_rank, file.path(outdir, "StepE2_program_rank_shift_vs_SFARI_all.tsv"))

run_summary <- data.table(
  n_group_rows = nrow(dt),
  n_mapped_groups = nrow(map_dt),
  n_broad_summary_rows = nrow(broad_sum),
  n_consensus_rows = nrow(consensus),
  n_rna_clusters_annotated = uniqueN(rna_anno$rna_cluster_id)
)
safe_fwrite(run_summary, file.path(outdir, "StepE2_run_summary.tsv"))

log_msg("Finished StepE2 successfully.")
