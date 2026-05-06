#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

parse_args <- function(x){
  out <- list(
    base = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    stage_mean_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/13_BrainSpan_stage_mean_expression_log2p1.tsv.gz",
    zscore_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/14_BrainSpan_stage_specificity_zscores.tsv.gz",
    sfari_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/01_SFARI_primary_Sand1_standardized.tsv",
    frozen_program_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/20_midPrenatal_core_programs.tsv",
    outdir = NA_character_
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    val <- if (i < length(x)) x[[i+1L]] else NA_character_
    if (key == "--base") out$base <- val
    if (key == "--stage_mean_file") out$stage_mean_file <- val
    if (key == "--zscore_file") out$zscore_file <- val
    if (key == "--sfari_file") out$sfari_file <- val
    if (key == "--frozen_program_file") out$frozen_program_file <- val
    if (key == "--outdir") out$outdir <- val
    i <- i + 2L
  }
  if (is.na(out$outdir)) out$outdir <- file.path(out$base, "reviewer_round3", "MA_round1A_midPrenatal_metric")
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "meta"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "MA_round1A_midPrenatal_metric.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
safe_fwrite <- function(x, file) fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x <- sub("\\..*$", "", x)
  x[nchar(x) == 0L] <- NA_character_
  toupper(x)
}
read_any <- function(f) {
  if (!file.exists(f)) stop("Missing required file: ", f)
  if (grepl("\\.gz$", f, ignore.case = TRUE)) return(fread(cmd = paste("zcat", shQuote(f))))
  fread(f)
}
row_sds <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sqrt(rowMeans((mat - mu)^2, na.rm = TRUE))
}

stage_mean <- read_any(args$stage_mean_file)
zs <- read_any(args$zscore_file)
sfari <- read_any(args$sfari_file)
frozen <- read_any(args$frozen_program_file)

if (!("gene_symbol" %in% names(stage_mean))) stop("stage_mean_file must contain gene_symbol")
if (!("gene_symbol" %in% names(zs))) stop("zscore_file must contain gene_symbol")
if (!("gene_symbol" %in% names(sfari))) stop("sfari_file must contain gene_symbol")
if (!all(c("program_name","gene_symbol") %in% names(frozen))) stop("frozen_program_file must contain program_name and gene_symbol")

stage_mean[, gene_symbol := clean_symbol(gene_symbol)]
zs[, gene_symbol := clean_symbol(gene_symbol)]
sfari[, gene_symbol := clean_symbol(gene_symbol)]
frozen[, gene_symbol := clean_symbol(gene_symbol)]

stage_cols <- setdiff(names(stage_mean), "gene_symbol")
stage_cols <- stage_cols[sapply(stage_cols, function(nm) mean(!is.na(suppressWarnings(as.numeric(stage_mean[[nm]])))) > 0.8)]
if (!("mid_prenatal" %in% stage_cols)) stop("Could not find mid_prenatal among numeric stage columns")
for (nm in stage_cols) stage_mean[[nm]] <- as.numeric(stage_mean[[nm]])

stage_mean <- unique(stage_mean[!is.na(gene_symbol), c("gene_symbol", stage_cols), with = FALSE], by = "gene_symbol")
expr_mat <- as.matrix(stage_mean[, ..stage_cols])
mu_g <- rowMeans(expr_mat, na.rm = TRUE)
sd_g <- row_sds(expr_mat)
sd_g[is.na(sd_g) | sd_g == 0] <- NA_real_
mid_prenatal_concentration <- (stage_mean[["mid_prenatal"]] - mu_g) / sd_g

metric_dt <- data.table(
  gene_symbol = stage_mean$gene_symbol,
  mid_prenatal_stage_mean_log2p1 = stage_mean[["mid_prenatal"]],
  across_stage_mean_log2p1 = mu_g,
  across_stage_sd_log2p1 = sd_g,
  mid_prenatal_cortical_concentration = mid_prenatal_concentration
)
metric_dt <- metric_dt[order(-mid_prenatal_cortical_concentration)]
metric_dt[, rank_desc := seq_len(.N)]
metric_dt[, frac_rank := rank_desc / .N]
metric_dt[, is_sfari_eligible := gene_symbol %in% unique(sfari$gene_symbol)]

zs[, mid_prenatal_from_frozen_z := suppressWarnings(as.numeric(mid_prenatal))]
cmp <- merge(metric_dt, zs[, .(gene_symbol, mid_prenatal_from_frozen_z)], by = "gene_symbol", all.x = TRUE)
corr_val <- suppressWarnings(cor(cmp$mid_prenatal_cortical_concentration, cmp$mid_prenatal_from_frozen_z, use = "pairwise.complete.obs", method = "pearson"))

frozen_sets <- split(unique(frozen$gene_symbol), frozen$program_name)
get_frozen <- function(name) if (name %in% names(frozen_sets)) frozen_sets[[name]] else character()
frozen05 <- get_frozen("midPrenatal_SFARI_top05")
frozen10 <- get_frozen("midPrenatal_SFARI_top10")
frozen20 <- get_frozen("midPrenatal_SFARI_top20")

metric_dt[, in_reconstructed_top05 := is_sfari_eligible & frac_rank <= 0.05]
metric_dt[, in_reconstructed_top10 := is_sfari_eligible & frac_rank <= 0.10]
metric_dt[, in_reconstructed_top20 := is_sfari_eligible & frac_rank <= 0.20]
metric_dt[, in_frozen_top05 := gene_symbol %in% frozen05]
metric_dt[, in_frozen_top10 := gene_symbol %in% frozen10]
metric_dt[, in_frozen_top20 := gene_symbol %in% frozen20]

recon05 <- metric_dt[in_reconstructed_top05 == TRUE, gene_symbol]
recon10 <- metric_dt[in_reconstructed_top10 == TRUE, gene_symbol]
recon20 <- metric_dt[in_reconstructed_top20 == TRUE, gene_symbol]

recon_check <- data.table(
  set_name = c("midPrenatal_SFARI_top05", "midPrenatal_SFARI_top10", "midPrenatal_SFARI_top20"),
  frozen_n = c(length(frozen05), length(frozen10), length(frozen20)),
  reconstructed_n = c(length(recon05), length(recon10), length(recon20)),
  overlap_n = c(length(intersect(frozen05, recon05)), length(intersect(frozen10, recon10)), length(intersect(frozen20, recon20))),
  jaccard = c(
    length(intersect(frozen05, recon05)) / length(unique(c(frozen05, recon05))),
    length(intersect(frozen10, recon10)) / length(unique(c(frozen10, recon10))),
    length(intersect(frozen20, recon20)) / length(unique(c(frozen20, recon20)))
  )
)

stage_meta <- data.table(stage_column = stage_cols, included_in_metric = TRUE)

formula_note <- c(
  "mid-prenatal cortical concentration metric (manuscript-facing explicit definition)",
  "",
  "Input:",
  "E_{g,s} = stage-mean BrainSpan cortical log2(expression + 1) value for gene g at developmental stage-group s.",
  paste0("Stage-groups used in the metric: ", paste(stage_cols, collapse = ", "), "."),
  "",
  "Per-gene summary quantities:",
  "mu_g = mean_s E_{g,s}",
  "sigma_g = sd_s E_{g,s}",
  "",
  "Metric:",
  "mid_prenatal_cortical_concentration_g = (E_{g, mid_prenatal} - mu_g) / sigma_g",
  "",
  "Ranking and frozen nested developmental-context programs:",
  "1) Rank all genes by descending mid_prenatal_cortical_concentration_g.",
  "2) Define frac_rank_g = rank_desc_g / N_all_genes.",
  "3) Restrict to genes present in the eligible SFARI pool.",
  "4) Freeze nested subsets as:",
  "   top05: eligible SFARI genes with frac_rank <= 0.05",
  "   top10: eligible SFARI genes with frac_rank <= 0.10",
  "   top20: eligible SFARI genes with frac_rank <= 0.20",
  "",
  "Controls / non-controls:",
  "- This metric is computed on within-gene stage profiles and therefore centers and scales each gene across developmental stage-groups.",
  "- It does not explicitly model gene length or GC content, because it is not a differential-expression bias-correction model; rather, it is a within-gene developmental concentration metric derived from stage-level normalized expression summaries.",
  "",
  paste0("Pearson correlation between recomputed metric and frozen z-score table mid_prenatal column: ", sprintf("%.6f", corr_val))
)
writeLines(formula_note, file.path(args$outdir, "meta", "MA_round1A_midPrenatal_metric_formula_note.txt"))

safe_fwrite(metric_dt, file.path(args$outdir, "tables", "MA_round1A_midPrenatal_metric_gene_table.tsv.gz"))
safe_fwrite(stage_meta, file.path(args$outdir, "tables", "MA_round1A_midPrenatal_metric_stage_columns.tsv"))
safe_fwrite(recon_check, file.path(args$outdir, "tables", "MA_round1A_midPrenatal_metric_reconstruction_check.tsv"))
safe_fwrite(data.table(item = c("n_genes", "n_sfari_eligible", "metric_vs_frozen_mid_prenatal_pearson"),
                       value = c(nrow(metric_dt), sum(metric_dt$is_sfari_eligible), corr_val)),
            file.path(args$outdir, "meta", "MA_round1A_run_summary.tsv"))

log_msg("Completed Round1A metric export.")
