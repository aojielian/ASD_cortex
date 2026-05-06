#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
RDSDIR <- file.path(OUTDIR, "rds")
LOGDIR <- file.path(OUTDIR, "logs")

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDSDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step01C_BrainSpan_stage_specificity_enrichment.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

calc_auc_from_scores <- function(scores, labels) {
  n1 <- sum(labels)
  n0 <- sum(!labels)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  U <- sum(ranks[labels]) - n1 * (n1 + 1) / 2
  U / (n1 * n0)
}

fisher_top_enrichment <- function(is_sfari, is_top) {
  tab <- matrix(c(
    sum(is_sfari & is_top),
    sum(is_sfari & !is_top),
    sum(!is_sfari & is_top),
    sum(!is_sfari & !is_top)
  ), nrow = 2, byrow = TRUE)
  ft <- fisher.test(tab, alternative = "greater")
  list(
    sfari_in_top = tab[1,1],
    sfari_not_top = tab[1,2],
    nonrisk_in_top = tab[2,1],
    nonrisk_not_top = tab[2,2],
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

stage_order <- c(
  "early_prenatal", "mid_prenatal", "late_prenatal",
  "birth", "infancy", "childhood", "adolescence", "adult"
)

sfari_file <- file.path(TABDIR, "01_SFARI_primary_Sand1_standardized.tsv")
brainspan_rds <- file.path(RDSDIR, "12_BrainSpan_cortex_subset.rds")

log_msg("Reading SFARI standardized set: ", sfari_file)
sfari <- fread(sfari_file)
if (!"gene_symbol" %in% names(sfari)) {
  stop("SFARI standardized file must contain gene_symbol column.")
}
sfari[, gene_symbol := toupper(trimws(gene_symbol))]
sfari <- unique(sfari[gene_symbol != "" & !is.na(gene_symbol)])

log_msg("Reading BrainSpan cortex subset RDS: ", brainspan_rds)
obj <- readRDS(brainspan_rds)
if (!all(c("expr", "samples") %in% names(obj))) {
  stop("BrainSpan cortex RDS must contain expr and samples.")
}

expr <- obj$expr
samples <- as.data.table(obj$samples)

if (!is.matrix(expr)) expr <- as.matrix(expr)
mode(expr) <- "numeric"

sample_id_col <- if ("sample_col" %in% names(samples)) {
  "sample_col"
} else if ("expr_col_original" %in% names(samples)) {
  "expr_col_original"
} else {
  stop("Could not identify sample column in obj$samples.")
}

if (!"stage_group" %in% names(samples)) {
  stop("obj$samples must contain stage_group.")
}

samples[, sample_id := as.character(get(sample_id_col))]
samples <- samples[sample_id %in% colnames(expr)]
samples <- samples[stage_group %in% stage_order]

stage_present <- intersect(stage_order, unique(as.character(samples$stage_group)))
samples[, stage_group := factor(as.character(stage_group), levels = stage_present)]

log_msg("Samples retained after stage filtering: ", nrow(samples))
log_msg("Stages present: ", paste(stage_present, collapse = ", "))

if (nrow(samples) < 50) {
  stop("Too few samples retained after stage filtering.")
}

expr <- expr[, samples$sample_id, drop = FALSE]
gene_symbols <- toupper(trimws(rownames(expr)))
keep_gene <- gene_symbols != "" & !is.na(gene_symbols)
expr <- expr[keep_gene, , drop = FALSE]
gene_symbols <- gene_symbols[keep_gene]

log_msg("Collapsing duplicated gene symbols.")
expr_dt <- as.data.table(expr)
expr_dt[, gene_symbol := gene_symbols]
expr_collapsed <- expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol]
expr_use <- as.matrix(expr_collapsed[, -1])
rownames(expr_use) <- expr_collapsed$gene_symbol
mode(expr_use) <- "numeric"

expr_use <- log2(expr_use + 1)

risk_genes <- intersect(sfari$gene_symbol, rownames(expr_use))
nonrisk_genes <- setdiff(rownames(expr_use), risk_genes)

log_msg("Risk genes in BrainSpan cortex matrix: ", length(risk_genes))
if (length(risk_genes) < 50) {
  stop("Too few SFARI genes overlap BrainSpan matrix.")
}

stage_means_list <- list()
for (stg in stage_present) {
  idx <- which(as.character(samples$stage_group) == stg)
  stage_means_list[[stg]] <- rowMeans(expr_use[, idx, drop = FALSE], na.rm = TRUE)
}
stage_means <- as.data.table(stage_means_list)
stage_means[, gene_symbol := rownames(expr_use)]
setcolorder(stage_means, c("gene_symbol", stage_present))
safe_fwrite(stage_means, file.path(TABDIR, "13_BrainSpan_stage_mean_expression_log2p1.tsv.gz"))

mat <- as.matrix(stage_means[, ..stage_present])
row_mu <- rowMeans(mat, na.rm = TRUE)
row_sd <- apply(mat, 1, sd, na.rm = TRUE)
row_sd[row_sd == 0 | is.na(row_sd)] <- 1
zmat <- (mat - row_mu) / row_sd

stage_z <- as.data.table(zmat)
stage_z[, gene_symbol := stage_means$gene_symbol]
setcolorder(stage_z, c("gene_symbol", stage_present))
safe_fwrite(stage_z, file.path(TABDIR, "14_BrainSpan_stage_specificity_zscores.tsv.gz"))

rank_res_list <- list()
top_res_list <- list()

for (stg in stage_present) {
  dt <- data.table(
    gene_symbol = stage_z$gene_symbol,
    score = stage_z[[stg]]
  )
  dt[, is_sfari := gene_symbol %in% risk_genes]
  dt <- dt[!is.na(score)]
  dt[, rank_desc := frank(-score, ties.method = "average")]
  dt <- dt[order(rank_desc)]

  sf <- dt[is_sfari == TRUE, score]
  bg <- dt[is_sfari == FALSE, score]

  wt <- wilcox.test(sf, bg, alternative = "greater", exact = FALSE)
  auc <- calc_auc_from_scores(dt$score, dt$is_sfari)

  rank_res_list[[stg]] <- data.table(
    stage_group = stg,
    n_sfari = sum(dt$is_sfari),
    n_background = sum(!dt$is_sfari),
    mean_z_sfari = mean(sf, na.rm = TRUE),
    mean_z_background = mean(bg, na.rm = TRUE),
    median_z_sfari = median(sf, na.rm = TRUE),
    median_z_background = median(bg, na.rm = TRUE),
    delta_mean_z = mean(sf, na.rm = TRUE) - mean(bg, na.rm = TRUE),
    delta_median_z = median(sf, na.rm = TRUE) - median(bg, na.rm = TRUE),
    auc = auc,
    wilcox_p_greater = wt$p.value
  )

  for (frac in c(0.05, 0.10, 0.20)) {
    top_n <- max(10, floor(nrow(dt) * frac))
    is_top <- seq_len(nrow(dt)) <= top_n
    ft <- fisher_top_enrichment(dt$is_sfari, is_top)

    top_res_list[[paste0(stg, "__", frac)]] <- data.table(
      stage_group = stg,
      top_fraction = frac,
      top_n = top_n,
      sfari_in_top = ft$sfari_in_top,
      sfari_not_top = ft$sfari_not_top,
      nonrisk_in_top = ft$nonrisk_in_top,
      nonrisk_not_top = ft$nonrisk_not_top,
      odds_ratio = ft$odds_ratio,
      fisher_p_greater = ft$p_value
    )
  }
}

rank_res <- rbindlist(rank_res_list)
rank_res[, fdr_wilcox := p.adjust(wilcox_p_greater, method = "BH")]
setorder(rank_res, fdr_wilcox, wilcox_p_greater, -auc)
safe_fwrite(rank_res, file.path(TABDIR, "15_BrainSpan_stage_rank_enrichment.tsv"))

top_res <- rbindlist(top_res_list)
top_res[, fdr_fisher := p.adjust(fisher_p_greater, method = "BH")]
setorder(top_res, fdr_fisher, fisher_p_greater, -odds_ratio)
safe_fwrite(top_res, file.path(TABDIR, "16_BrainSpan_stage_topFraction_enrichment.tsv"))

prenatal_stages <- intersect(c("early_prenatal", "mid_prenatal", "late_prenatal"), stage_present)
postnatal_stages <- intersect(c("infancy", "childhood", "adolescence", "adult"), stage_present)

if (length(prenatal_stages) >= 2 && length(postnatal_stages) >= 2) {
  prenatal_mean <- rowMeans(as.matrix(stage_means[, ..prenatal_stages]), na.rm = TRUE)
  postnatal_mean <- rowMeans(as.matrix(stage_means[, ..postnatal_stages]), na.rm = TRUE)
  bias_dt <- data.table(
    gene_symbol = stage_means$gene_symbol,
    prenatal_mean = prenatal_mean,
    postnatal_mean = postnatal_mean,
    prenatal_bias = prenatal_mean - postnatal_mean
  )
  bias_dt[, is_sfari := gene_symbol %in% risk_genes]
  wt_bias <- wilcox.test(
    bias_dt[is_sfari == TRUE, prenatal_bias],
    bias_dt[is_sfari == FALSE, prenatal_bias],
    alternative = "greater",
    exact = FALSE
  )
  auc_bias <- calc_auc_from_scores(bias_dt$prenatal_bias, bias_dt$is_sfari)

  safe_fwrite(bias_dt, file.path(TABDIR, "17_BrainSpan_prenatal_vs_postnatal_gene_bias.tsv.gz"))

  prenatal_summary <- data.table(
    metric = c("n_prenatal_stages", "n_postnatal_stages", "mean_bias_sfari", "mean_bias_background",
               "median_bias_sfari", "median_bias_background", "delta_mean_bias", "auc", "wilcox_p_greater"),
    value = c(
      length(prenatal_stages),
      length(postnatal_stages),
      mean(bias_dt[is_sfari == TRUE, prenatal_bias], na.rm = TRUE),
      mean(bias_dt[is_sfari == FALSE, prenatal_bias], na.rm = TRUE),
      median(bias_dt[is_sfari == TRUE, prenatal_bias], na.rm = TRUE),
      median(bias_dt[is_sfari == FALSE, prenatal_bias], na.rm = TRUE),
      mean(bias_dt[is_sfari == TRUE, prenatal_bias], na.rm = TRUE) - mean(bias_dt[is_sfari == FALSE, prenatal_bias], na.rm = TRUE),
      auc_bias,
      wt_bias$p.value
    )
  )
} else {
  prenatal_summary <- data.table(
    metric = c("n_prenatal_stages", "n_postnatal_stages", "note"),
    value = c(length(prenatal_stages), length(postnatal_stages), "Not enough stages for prenatal_vs_postnatal comparison")
  )
}
safe_fwrite(prenatal_summary, file.path(METADIR, "17b_BrainSpan_prenatal_vs_postnatal_summary.tsv"))

best_rank <- rank_res[1]
best_top <- top_res[1]

run_summary <- rbindlist(list(
  data.table(section = "BrainSpan", metric = "n_genes_total", value = nrow(stage_means)),
  data.table(section = "BrainSpan", metric = "n_sfari_in_matrix", value = length(risk_genes)),
  data.table(section = "BrainSpan", metric = "stages_present", value = paste(stage_present, collapse = ",")),
  data.table(section = "RankEnrichment", metric = "best_stage", value = best_rank$stage_group),
  data.table(section = "RankEnrichment", metric = "best_stage_auc", value = best_rank$auc),
  data.table(section = "RankEnrichment", metric = "best_stage_wilcox_p", value = best_rank$wilcox_p_greater),
  data.table(section = "RankEnrichment", metric = "best_stage_fdr", value = best_rank$fdr_wilcox),
  data.table(section = "TopFraction", metric = "best_stage", value = best_top$stage_group),
  data.table(section = "TopFraction", metric = "best_top_fraction", value = best_top$top_fraction),
  data.table(section = "TopFraction", metric = "best_odds_ratio", value = best_top$odds_ratio),
  data.table(section = "TopFraction", metric = "best_fisher_p", value = best_top$fisher_p_greater),
  data.table(section = "TopFraction", metric = "best_fdr", value = best_top$fdr_fisher)
))
safe_fwrite(run_summary, file.path(METADIR, "18_Step01C_run_summary.tsv"))

log_msg("Step01C completed successfully.")
log_msg("Best rank-enrichment stage: ", best_rank$stage_group, " ; p=", best_rank$wilcox_p_greater, " ; fdr=", best_rank$fdr_wilcox, " ; auc=", best_rank$auc)
log_msg("Best top-fraction stage: ", best_top$stage_group, " ; frac=", best_top$top_fraction, " ; p=", best_top$fisher_p_greater, " ; fdr=", best_top$fdr_fisher, " ; OR=", best_top$odds_ratio)
