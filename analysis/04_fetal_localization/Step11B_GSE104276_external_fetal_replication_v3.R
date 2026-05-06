#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(
    counts_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/GSE104276_all_pfc_2394_UMI_count_NOERCC.xls.gz",
    readme_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/GSE104276_readme_sample_barcode.xlsx",
    programs_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/20_midPrenatal_core_programs.tsv",
    sfari_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/01_SFARI_primary_Sand1_standardized.tsv",
    outdir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step11B_GSE104276_external_fetal_replication_v3",
    norm_scale_factor = 10000,
    min_cells_per_state = 10
  )
  if (length(args) %% 2 != 0) stop("Arguments must be key-value pairs.")
  if (length(args) > 0) {
    for (i in seq(1, length(args), by = 2)) {
      k <- sub("^--", "", args[i])
      opt[[k]] <- args[i + 1]
    }
  }
  opt$norm_scale_factor <- as.numeric(opt$norm_scale_factor)
  opt$min_cells_per_state <- as.integer(opt$min_cells_per_state)
  opt
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

norm_name <- function(x) {
  x <- tolower(as.character(x))
  gsub("[^a-z0-9]+", "_", x)
}

list_xlsx_sheets <- function(path) {
  if (requireNamespace("readxl", quietly = TRUE)) return(readxl::excel_sheets(path))
  if (requireNamespace("openxlsx", quietly = TRUE)) return(openxlsx::getSheetNames(path))
  stop("Need readxl or openxlsx installed.")
}

read_xlsx_sheet <- function(path, sheet) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    return(as.data.table(readxl::read_xlsx(path, sheet = sheet)))
  }
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    return(as.data.table(openxlsx::read.xlsx(path, sheet = sheet)))
  }
  stop("Need readxl or openxlsx installed.")
}

load_programs <- function(core_path, sfari_path) {
  core <- fread(core_path, sep = "\t", header = TRUE)
  core <- core[program_name %in% c("midPrenatal_SFARI_top20"), .(program_name, gene_symbol)]
  sf <- fread(sfari_path, sep = "\t", header = TRUE)
  nn <- norm_name(names(sf))
  hit <- grep("gene_symbol|symbol|gene", nn, perl = TRUE)
  if (length(hit) == 0) stop("Could not identify gene symbol column in SFARI file.")
  symbol_col <- names(sf)[hit[1]]
  sfari <- data.table(program_name = "SFARI_all", gene_symbol = as.character(sf[[symbol_col]]))
  prog <- rbindlist(list(core, sfari), use.names = TRUE)
  prog[, gene_symbol := toupper(trimws(gene_symbol))]
  prog <- unique(prog[!is.na(gene_symbol) & gene_symbol != ""])
  split(prog$gene_symbol, prog$program_name)
}

map_author_coarse <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("stem", y)] <- "Stem cells"
  out[grepl("gaba", y)] <- "GABAergic neurons"
  out[grepl("neuron", y)] <- "Neurons"
  out[grepl("astro", y)] <- "Astrocytes"
  out[grepl("opc", y)] <- "OPCs"
  out[grepl("micro", y)] <- "Microglia"
  out[grepl("oligo", y)] <- "Oligodendrocytes"
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

infer_sample_id <- function(cell_id) {
  sub("_sc[0-9]+$", "", as.character(cell_id), perl = TRUE)
}

score_marker_set <- function(norm_mat, genes) {
  idx <- which(toupper(rownames(norm_mat)) %in% toupper(genes))
  if (length(idx) == 0) return(rep(NA_real_, ncol(norm_mat)))
  Matrix::colMeans(norm_mat[idx, , drop = FALSE])
}

refine_stem_cells <- function(norm_mat, coarse_vec) {
  rg_markers  <- c("FABP7","HOPX","PAX6","SOX2","VIM","SLC1A3","HES1","PTN")
  ipc_markers <- c("EOMES","PPP1R17","HES6","NHLH1","ASCL1")
  rg_score  <- score_marker_set(norm_mat, rg_markers)
  ipc_score <- score_marker_set(norm_mat, ipc_markers)
  refined <- coarse_vec
  stem_idx <- which(coarse_vec == "Stem cells")
  if (length(stem_idx) > 0) {
    refined[stem_idx] <- ifelse(rg_score[stem_idx] >= ipc_score[stem_idx], "Radial Glial-like", "IPC-like")
  }
  list(refined = refined, rg_score = rg_score, ipc_score = ipc_score)
}

main <- function() {
  opt <- parse_args()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(opt$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

  log_msg("Loading programs")
  programs <- load_programs(opt$programs_file, opt$sfari_file)
  prog_dt <- rbindlist(lapply(names(programs), function(nm) data.table(program_name = nm, gene_symbol = programs[[nm]])))
  fwrite(prog_dt, file.path(opt$outdir, "tables", "00_program_gene_lists.tsv"), sep = "\t")

  log_msg("Reading count matrix")
  dt <- fread(opt$counts_file, sep = "\t", header = TRUE, fill = TRUE, quote = "")
  fwrite(data.table(raw_columns = names(dt)), file.path(opt$outdir, "tables", "01_count_matrix_columns.tsv"), sep = "\t")

  first_col <- names(dt)[1]
  setnames(dt, first_col, "gene_symbol")
  genes <- toupper(as.character(dt$gene_symbol))
  dt[, gene_symbol := NULL]
  cell_ids <- names(dt)

  mat <- as.matrix(dt)
  storage.mode(mat) <- "numeric"
  mat <- Matrix(mat, sparse = TRUE)
  rownames(mat) <- genes
  colnames(mat) <- cell_ids

  log_msg("Reading workbook")
  sheets <- list_xlsx_sheets(opt$readme_file)
  fwrite(data.table(sheet = sheets), file.path(opt$outdir, "tables", "02_sheet_names.tsv"), sep = "\t")
  if (!("SampleInfo" %in% sheets)) {
    stop("Workbook does not contain SampleInfo. Available: ", paste(sheets, collapse = ", "))
  }
  md <- read_xlsx_sheet(opt$readme_file, "SampleInfo")
  fwrite(data.table(sampleinfo_columns = names(md)), file.path(opt$outdir, "tables", "03_sampleinfo_columns.tsv"), sep = "\t")

  setnames(md, names(md)[1], "cell_id")
  nn <- norm_name(names(md))
  celltype_hit <- grep("^cell_types$|cell_types|celltype", nn, perl = TRUE)
  week_hit <- grep("^week$|gest|gw|stage", nn, perl = TRUE)
  if (length(celltype_hit) == 0 || length(week_hit) == 0) {
    stop("Could not identify cell_types or week column in SampleInfo.")
  }
  celltype_col <- names(md)[celltype_hit[1]]
  week_col <- names(md)[week_hit[1]]

  md[, cell_id := as.character(cell_id)]
  md[, cell_types := as.character(get(celltype_col))]
  md[, week := as.character(get(week_col))]
  md[, stage_num := suppressWarnings(as.integer(gsub("[^0-9]", "", week)))]
  md[, stage_label := ifelse(is.na(stage_num), "NA", paste0("GW", sprintf("%02d", stage_num)))]
  md[, sample_id := infer_sample_id(cell_id)]
  md[, author_coarse := map_author_coarse(cell_types)]

  overlap_n <- sum(md$cell_id %in% cell_ids)
  fwrite(data.table(metric = "sampleinfo_cell_overlap", value = overlap_n),
         file.path(opt$outdir, "tables", "04_overlap_diagnostics.tsv"), sep = "\t")
  if (overlap_n == 0) stop("No overlap between SampleInfo cell_id and count matrix cell IDs.")

  md2 <- md[cell_id %in% cell_ids]
  md2[, col_idx := match(cell_id, cell_ids)]
  setorder(md2, col_idx)
  mat <- mat[, md2$col_idx, drop = FALSE]

  log_msg("Normalizing matrix")
  lib <- Matrix::colSums(mat)
  lib[lib == 0] <- 1
  norm <- t(t(mat) / lib) * opt$norm_scale_factor
  norm@x <- log1p(norm@x)

  log_msg("Refining stem cells")
  stem_ref <- refine_stem_cells(norm, md2$author_coarse)
  md2[, refined_class := stem_ref$refined]
  md2[, rg_score := stem_ref$rg_score]
  md2[, ipc_score := stem_ref$ipc_score]

  log_msg("Scoring programs")
  overlap_dt <- rbindlist(lapply(names(programs), function(nm) {
    data.table(program_name = nm,
               genes_requested = length(programs[[nm]]),
               genes_present = sum(rownames(norm) %in% programs[[nm]]))
  }))
  fwrite(overlap_dt, file.path(opt$outdir, "tables", "05_program_overlap.tsv"), sep = "\t")

  for (nm in names(programs)) {
    idx <- which(rownames(norm) %in% programs[[nm]])
    md2[[nm]] <- if (length(idx) == 0) NA_real_ else Matrix::colMeans(norm[idx, , drop = FALSE])
  }
  fwrite(md2, file.path(opt$outdir, "tables", "06_cell_metadata_scored.tsv.gz"), sep = "\t")

  coarse_dt <- md2[, .(
    n_cells = .N,
    SFARI_all = mean(SFARI_all, na.rm = TRUE),
    midPrenatal_SFARI_top20 = mean(midPrenatal_SFARI_top20, na.rm = TRUE)
  ), by = .(sample_id, stage_num, stage_label, author_coarse)]
  fwrite(coarse_dt, file.path(opt$outdir, "tables", "07_sample_authorCoarse_program_summary.tsv"), sep = "\t")

  refined_dt <- md2[, .(
    n_cells = .N,
    SFARI_all = mean(SFARI_all, na.rm = TRUE),
    midPrenatal_SFARI_top20 = mean(midPrenatal_SFARI_top20, na.rm = TRUE),
    mean_rg_score = mean(rg_score, na.rm = TRUE),
    mean_ipc_score = mean(ipc_score, na.rm = TRUE)
  ), by = .(sample_id, stage_num, stage_label, refined_class)]
  fwrite(refined_dt, file.path(opt$outdir, "tables", "08_sample_refinedClass_program_summary.tsv"), sep = "\t")

  refined_long <- melt(
    refined_dt,
    id.vars = c("sample_id","stage_num","stage_label","refined_class","n_cells","mean_rg_score","mean_ipc_score"),
    measure.vars = c("SFARI_all","midPrenatal_SFARI_top20"),
    variable.name = "program_name",
    value.name = "mean_score"
  )
  robust_states <- refined_long[n_cells >= opt$min_cells_per_state]
  fwrite(robust_states, file.path(opt$outdir, "tables", "09_robust_states_long.tsv.gz"), sep = "\t")

  top_states <- robust_states[order(program_name, -mean_score)][, head(.SD, 12), by = program_name]
  fwrite(top_states, file.path(opt$outdir, "tables", "10_top_states_by_program.tsv"), sep = "\t")

  stage_refined <- robust_states[, .(
    total_cells = sum(n_cells),
    weighted_mean_score = weighted.mean(mean_score, w = pmax(n_cells, 1), na.rm = TRUE),
    n_states = .N
  ), by = .(program_name, stage_num, stage_label, refined_class)]
  fwrite(stage_refined, file.path(opt$outdir, "tables", "11_stage_refinedClass_summary.tsv"), sep = "\t")

  stage_conc <- robust_states[, .(
    stage_level_mean_score = weighted.mean(mean_score, w = pmax(n_cells, 1), na.rm = TRUE),
    n_state_groups = .N
  ), by = .(program_name, stage_num, stage_label)]
  fwrite(stage_conc, file.path(opt$outdir, "tables", "12_stage_concentration_summary.tsv"), sep = "\t")

  pA <- ggplot(stage_refined, aes(x = factor(stage_label, levels = unique(stage_label[order(stage_num)])),
                                  y = refined_class, fill = weighted_mean_score)) +
    geom_tile(color = "white") +
    facet_wrap(~ program_name, nrow = 1) +
    theme_bw(base_size = 11) +
    labs(title = "GSE104276 external fetal reference: stage-by-refined-class localization",
         x = "Gestational week", y = "Refined class", fill = "Weighted\nmean score")
  ggsave(file.path(opt$outdir, "plots", "GSE104276_external_replication_A.png"), pA, width = 10, height = 5)

  bar_dt <- copy(top_states)
  bar_dt[, state_label := paste(sample_id, refined_class, sep = " | ")]
  pB <- ggplot(bar_dt, aes(x = mean_score, y = reorder(state_label, mean_score), fill = factor(stage_label))) +
    geom_col() +
    facet_wrap(~ program_name, scales = "free_y", nrow = 1) +
    theme_bw(base_size = 11) +
    labs(title = "Top states in GSE104276 external fetal reference",
         x = "Mean score", y = NULL, fill = "GW")
  ggsave(file.path(opt$outdir, "plots", "GSE104276_external_replication_B.png"), pB, width = 11, height = 6)

  pC <- ggplot(stage_conc, aes(x = stage_num, y = stage_level_mean_score, color = program_name, size = n_state_groups)) +
    geom_line() + geom_point() +
    theme_bw(base_size = 11) +
    labs(title = "Stage concentration in GSE104276 external fetal reference",
         x = "Gestational week", y = "Stage-level mean score", size = "State groups", color = "Program")
  ggsave(file.path(opt$outdir, "plots", "GSE104276_external_replication_C.png"), pC, width = 8, height = 4.5)

  summary_txt <- c(
    "Step11B GSE104276 external fetal replication v3 summary",
    "",
    paste("Counts file:", opt$counts_file),
    paste("Readme file:", opt$readme_file),
    paste("Cells retained after SampleInfo overlap:", ncol(mat)),
    paste("Genes:", nrow(mat)),
    "",
    "Program overlap:",
    capture.output(print(overlap_dt)),
    "",
    "Refined class counts:",
    capture.output(print(md2[, .N, by = refined_class][order(-N)])),
    "",
    "Top states:",
    capture.output(print(top_states[, .(program_name, sample_id, stage_label, refined_class, n_cells, mean_score)]))
  )
  writeLines(summary_txt, file.path(opt$outdir, "13_external_replication_summary.txt"))
  log_msg("DONE. Outputs written to ", opt$outdir)
}

main()
