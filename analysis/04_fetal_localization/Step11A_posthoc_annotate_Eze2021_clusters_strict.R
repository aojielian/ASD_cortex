#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list(
    base_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    input_results_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/Step11A_fetal_singlecell_localization_final",
    supp_xlsx = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/inputs/Eze2021_Supplementary_Tables.xlsx",
    outdir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/Step11A_fetal_singlecell_localization_final_posthoc_strict"
  )
  if (length(args) == 0) return(res)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args)) stop("Missing value for argument: ", key)
    val <- args[[i + 1]]
    nm <- sub("^--", "", key)
    if (!nm %in% names(res)) stop("Unknown argument: ", nm)
    res[[nm]] <- val
    i <- i + 2
  }
  res
}

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
}

make_clean_names <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

collapse_top_markers <- function(dt, top_n = 5L) {
  dt <- copy(dt)
  dt <- dt[!is.na(cluster) & !is.na(gene)]
  setorderv(dt, c("cluster", "p_val_adj", "p_val", "avg_logfc"), c(1, 1, 1, -1), na.last = TRUE)
  dt[, rank_within_cluster := seq_len(.N), by = cluster]
  dt_top <- dt[rank_within_cluster <= top_n]
  out <- dt_top[, .(
    official_top5_markers = paste(gene, collapse = ","),
    official_marker_n = .N
  ), by = cluster]
  out[]
}

args <- parse_args()
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

log_msg("Starting strict Step11A posthoc Eze2021 annotation from official supplementary tables")
log_msg("input_results_dir=", args$input_results_dir)
log_msg("supp_xlsx=", args$supp_xlsx)

top_states_file <- file.path(args$input_results_dir, "tables", "10_top_developmental_states_by_score.tsv")
summary_file <- file.path(args$input_results_dir, "12_fetal_localization_summary.txt")

stop_if_missing(top_states_file, "Step11A top developmental states file")
stop_if_missing(summary_file, "Step11A summary file")
stop_if_missing(args$supp_xlsx, "Eze2021 supplementary workbook")

top_dt <- fread(top_states_file)
if (!"celltype" %in% colnames(top_dt)) stop("Expected column 'celltype' missing from top states table: ", top_states_file)

ann <- as.data.table(read_excel(args$supp_xlsx, sheet = "5 Cortex Annotations"))
mrk <- as.data.table(read_excel(args$supp_xlsx, sheet = "6 Cortex Clustermarker Genes"))

setnames(ann, old = names(ann), new = make_clean_names(names(ann)))
setnames(mrk, old = names(mrk), new = make_clean_names(names(mrk)))

required_ann <- c("cluster", "cell_type")
required_mrk <- c("cluster", "gene")
if (!all(required_ann %in% names(ann))) stop("Sheet '5 Cortex Annotations' is missing required columns: ", paste(setdiff(required_ann, names(ann)), collapse = ", "))
if (!all(required_mrk %in% names(mrk))) stop("Sheet '6 Cortex Clustermarker Genes' is missing required columns: ", paste(setdiff(required_mrk, names(mrk)), collapse = ", "))

ann <- ann[!is.na(cluster) & !is.na(cell_type)]
ann[, cluster := as.integer(cluster)]
ann[, cell_type := as.character(cell_type)]

official_map <- ann[, .(
  official_cell_type = names(sort(table(cell_type), decreasing = TRUE))[1],
  official_n_cells = .N,
  official_n_unique_cell_types = uniqueN(cell_type)
), by = cluster]
setorder(official_map, cluster)

mrk <- mrk[!is.na(cluster) & !is.na(gene)]
mrk[, cluster := as.integer(cluster)]
mrk[, gene := as.character(gene)]
if ("avg_logfc" %in% names(mrk)) mrk[, avg_logfc := suppressWarnings(as.numeric(avg_logfc))]
if ("p_val_adj" %in% names(mrk)) mrk[, p_val_adj := suppressWarnings(as.numeric(p_val_adj))]
if ("p_val" %in% names(mrk)) mrk[, p_val := suppressWarnings(as.numeric(p_val))]

official_markers <- collapse_top_markers(mrk, top_n = 5L)
setorder(official_markers, cluster)

top_dt[, cluster_id := suppressWarnings(as.integer(celltype))]
annot_dt <- merge(top_dt, official_map, by.x = "cluster_id", by.y = "cluster", all.x = TRUE)
annot_dt <- merge(annot_dt, official_markers, by.x = "cluster_id", by.y = "cluster", all.x = TRUE)

annot_dt[, primary_method := method == "MeanLogExpr"]
annot_dt[, robust_state := !is.na(n_cells) & n_cells >= 5]

meanlog_all <- annot_dt[method == "MeanLogExpr"]
meanlog_robust <- meanlog_all[robust_state == TRUE]

setorder(meanlog_all, program, -mean_score, -n_cells, cluster_id)
setorder(meanlog_robust, program, -mean_score, -n_cells, cluster_id)

fwrite(official_map, file.path(args$outdir, "tables", "20_official_cluster_mapping.tsv"), sep = "\t")
fwrite(official_markers, file.path(args$outdir, "tables", "21_official_cluster_top5_markers.tsv"), sep = "\t")
fwrite(annot_dt, file.path(args$outdir, "tables", "22_all_states_annotated_strict.tsv"), sep = "\t")
fwrite(meanlog_all, file.path(args$outdir, "tables", "23_meanlogexpr_states_annotated_strict.tsv"), sep = "\t")
fwrite(meanlog_robust, file.path(args$outdir, "tables", "24_meanlogexpr_states_annotated_strict_nge5.tsv"), sep = "\t")

top_by_program <- meanlog_robust[, head(.SD, 6), by = program]

lines <- c(
  "Step11A strict posthoc annotation of Eze2021 fetal cluster codes",
  "",
  paste0("Input results directory: ", args$input_results_dir),
  paste0("Official supplementary workbook: ", args$supp_xlsx),
  "",
  "Strict interpretation basis:",
  "- Cluster-to-cell-type mapping was derived directly from Eze et al. 2021 Supplementary Table 5 (sheet: '5 Cortex Annotations').",
  "- Representative cluster markers were derived directly from Eze et al. 2021 Supplementary Table 6 (sheet: '6 Cortex Clustermarker Genes').",
  "- MeanLogExpr is treated as the primary fetal-localization metric; UCell remains supportive only.",
  "- Robust states are defined here as cluster-stage combinations with n_cells >= 5.",
  ""
)

for (pgm in unique(top_by_program$program)) {
  lines <- c(lines, paste0("Program: ", pgm))
  sub <- top_by_program[program == pgm]
  for (i in seq_len(nrow(sub))) {
    row <- sub[i]
    lines <- c(
      lines,
      sprintf("  Cluster %s | official_cell_type=%s | stage=%s | n_cells=%s | mean_score=%.3f | official_top5_markers=%s",
              row$cluster_id,
              ifelse(is.na(row$official_cell_type), "NA", row$official_cell_type),
              ifelse(is.na(row$stage), "NA", as.character(row$stage)),
              ifelse(is.na(row$n_cells), "NA", as.character(row$n_cells)),
              row$mean_score,
              ifelse(is.na(row$official_top5_markers), "NA", row$official_top5_markers))
    )
  }
  lines <- c(lines, "")
}

writeLines(lines, file.path(args$outdir, "25_posthoc_interpretation_strict.txt"))

log_msg("Completed strict Step11A posthoc Eze2021 annotation from official supplementary tables")
