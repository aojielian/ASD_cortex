suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUT  <- file.path(BASE, "step60_MA_strengthening_v3/results/Step60A_BrainSpan_boundary_sensitivity_v1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(OUT, "Step60A_v1_run.log")
sink(log_file, split = TRUE)

cat("Step60A v1: BrainSpan mid-prenatal boundary/window sensitivity\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("BASE:", BASE, "\n")
cat("OUT :", OUT, "\n\n")

safe_fread <- function(path) {
  cat("[READ]", path, "\n")
  if (!file.exists(path)) stop("Missing file: ", path)
  fread(path)
}

safe_write <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, sep = "\t")
  cat("[WRITE]", path, " nrow=", nrow(dt), " ncol=", ncol(dt), "\n")
}

std_gene <- function(x) {
  x <- as.character(x)
  x <- sub("\\..*$", "", x)
  toupper(trimws(x))
}

std_sample <- function(x) {
  x <- as.character(x)
  x <- gsub('"', "", x, fixed = TRUE)
  trimws(x)
}

norm_dx <- function(x) {
  y <- as.character(x)
  yl <- tolower(y)
  out <- rep(NA_character_, length(y))
  out[grepl("asd|autism|case", yl)] <- "ASD"
  out[grepl("control|ctrl|ctl|unaffected|normal|td", yl)] <- "Control"
  out[yl %in% c("2")] <- "ASD"
  out[yl %in% c("1", "0")] <- "Control"
  out
}

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------
stage_z_path <- file.path(BASE, "step01_prepare/tables/14_BrainSpan_stage_specificity_zscores.tsv.gz")
sample_stage_path <- file.path(BASE, "step01_prepare/tables/08_BrainSpan_samples_harmonized.tsv")
stage_only_path <- file.path(BASE, "step01_prepare/tables/10_BrainSpan_cortex_stage_summary.tsv")
sfari_path <- file.path(BASE, "step01_prepare/tables/01_SFARI_primary_Sand1_standardized.tsv")
frozen_path <- file.path(BASE, "reviewer_round1/R1A_sfari_program_overlap_audit/tables/R1A_frozen_program_memberships.tsv")

gse102_expr <- file.path(BASE, "step01_prepare/inputs/GSE102741.standardized_expr.tsv.gz")
gse102_meta <- file.path(BASE, "step01_prepare/inputs/GSE102741.standardized_meta.tsv")
gse640_expr <- file.path(BASE, "step01_prepare/inputs/GSE64018.standardized_expr.tsv.gz")
gse640_meta <- file.path(BASE, "step01_prepare/inputs/GSE64018.standardized_meta.tsv")
gandal_path <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"

input_manifest <- data.table(
  input = c(
    "stage_z_path", "sample_stage_path", "stage_only_path", "sfari_path", "frozen_path",
    "gse102_expr", "gse102_meta", "gse640_expr", "gse640_meta", "gandal_path"
  ),
  path = c(
    stage_z_path, sample_stage_path, stage_only_path, sfari_path, frozen_path,
    gse102_expr, gse102_meta, gse640_expr, gse640_meta, gandal_path
  )
)
input_manifest[, exists := file.exists(path)]
safe_write(input_manifest, file.path(OUT, "00_input_manifest.tsv"))
print(input_manifest)

if (any(!input_manifest$exists)) {
  stop("Some required inputs are missing. See 00_input_manifest.tsv")
}

# ------------------------------------------------------------
# 2. BrainSpan stage/window definitions
# ------------------------------------------------------------
stage_samples <- safe_fread(sample_stage_path)
stage_samples[, is_cortex := as.logical(is_cortex)]

stage_def <- stage_samples[is_cortex == TRUE, .(
  n_cortex_samples = .N,
  n_donors = uniqueN(donor_id),
  age_examples = paste(head(unique(age), 8), collapse = "; "),
  region_groups = paste(sort(unique(region_group[!is.na(region_group)])), collapse = "; ")
), by = stage_group][order(stage_group)]

safe_write(stage_def, file.path(OUT, "01_BrainSpan_cortex_stage_group_summary.tsv"))

boundary_defs <- data.table(
  boundary_id = c(
    "early_prenatal_top20",
    "mid_prenatal_top20_primary",
    "late_prenatal_top20",
    "early_mid_prenatal_top20",
    "mid_late_prenatal_top20",
    "all_prenatal_top20",
    "mid_minus_adjacent_top20"
  ),
  score_formula = c(
    "early_prenatal",
    "mid_prenatal",
    "late_prenatal",
    "mean(early_prenatal, mid_prenatal)",
    "mean(mid_prenatal, late_prenatal)",
    "mean(early_prenatal, mid_prenatal, late_prenatal)",
    "mid_prenatal - mean(early_prenatal, late_prenatal)"
  ),
  interpretation = c(
    "Adjacent earlier prenatal cortical concentration",
    "Primary manuscript mid-prenatal cortical concentration",
    "Adjacent later prenatal cortical concentration",
    "Broader boundary including early-to-mid prenatal window",
    "Broader boundary including mid-to-late prenatal window",
    "Broad prenatal window across early/mid/late prenatal stages",
    "Mid-prenatal specificity relative to adjacent prenatal stages"
  ),
  selection_rule = "SFARI source-pool genes with whole-transcriptome frac_rank <= 0.20"
)
safe_write(boundary_defs, file.path(OUT, "02_boundary_window_definitions.tsv"))

# ------------------------------------------------------------
# 3. Build alternative boundary-ranked gene sets
# ------------------------------------------------------------
z <- safe_fread(stage_z_path)
needed_cols <- c("gene_symbol", "early_prenatal", "mid_prenatal", "late_prenatal")
if (!all(needed_cols %in% names(z))) {
  stop("Stage z-score file missing required columns: ", paste(setdiff(needed_cols, names(z)), collapse = ", "))
}

z[, gene := std_gene(gene_symbol)]

sfari <- safe_fread(sfari_path)
sfari_gene_col <- names(sfari)[grepl("gene|symbol", names(sfari), ignore.case = TRUE)][1]
sfari_genes <- unique(std_gene(sfari[[sfari_gene_col]]))
sfari_genes <- sfari_genes[!is.na(sfari_genes) & sfari_genes != ""]

frozen <- safe_fread(frozen_path)
frozen[, program_name := as.character(program_name)]
frozen[, gene := std_gene(gene_symbol)]

frozen_top05 <- unique(frozen[program_name == "midPrenatal_SFARI_top05", gene])
frozen_top10 <- unique(frozen[program_name == "midPrenatal_SFARI_top10", gene])
frozen_top20 <- unique(frozen[program_name == "midPrenatal_SFARI_top20", gene])
frozen_sfari <- unique(frozen[program_name == "SFARI_all", gene])

# Score definitions
rank_dt <- z[, .(
  gene,
  early_prenatal_top20 = early_prenatal,
  mid_prenatal_top20_primary = mid_prenatal,
  late_prenatal_top20 = late_prenatal,
  early_mid_prenatal_top20 = rowMeans(cbind(early_prenatal, mid_prenatal), na.rm = TRUE),
  mid_late_prenatal_top20 = rowMeans(cbind(mid_prenatal, late_prenatal), na.rm = TRUE),
  all_prenatal_top20 = rowMeans(cbind(early_prenatal, mid_prenatal, late_prenatal), na.rm = TRUE),
  mid_minus_adjacent_top20 = mid_prenatal - rowMeans(cbind(early_prenatal, late_prenatal), na.rm = TRUE)
)]

rank_long <- melt(
  rank_dt,
  id.vars = "gene",
  variable.name = "boundary_id",
  value.name = "boundary_score"
)
rank_long <- rank_long[!is.na(boundary_score)]
rank_long[, is_sfari_source := gene %in% sfari_genes]

# Rank over full gene universe for each boundary
setorder(rank_long, boundary_id, -boundary_score, gene)
rank_long[, rank_desc := seq_len(.N), by = boundary_id]
rank_long[, frac_rank := rank_desc / .N, by = boundary_id]

safe_write(rank_long, file.path(OUT, "03_all_gene_boundary_rankings_long.tsv.gz"))

alt_membership <- rank_long[
  is_sfari_source == TRUE & frac_rank <= 0.20,
  .(boundary_id, gene, boundary_score, rank_desc, frac_rank)
]

safe_write(alt_membership, file.path(OUT, "04_alternative_boundary_top20_gene_sets.tsv"))

alt_size <- alt_membership[, .(
  n_genes = uniqueN(gene),
  mean_boundary_score = mean(boundary_score, na.rm = TRUE),
  median_boundary_score = median(boundary_score, na.rm = TRUE)
), by = boundary_id][order(boundary_id)]

safe_write(alt_size, file.path(OUT, "05_alternative_boundary_gene_set_sizes.tsv"))

overlap_summary <- alt_membership[, {
  g <- unique(gene)
  .(
    n_genes = length(g),
    overlap_frozen_top05 = length(intersect(g, frozen_top05)),
    overlap_frozen_top10 = length(intersect(g, frozen_top10)),
    overlap_frozen_top20 = length(intersect(g, frozen_top20)),
    pct_overlap_frozen_top20 = 100 * length(intersect(g, frozen_top20)) / max(1, length(frozen_top20)),
    jaccard_frozen_top20 = length(intersect(g, frozen_top20)) / length(union(g, frozen_top20)),
    overlap_frozen_SFARI_all = length(intersect(g, frozen_sfari))
  )
}, by = boundary_id][order(-jaccard_frozen_top20)]

safe_write(overlap_summary, file.path(OUT, "06_overlap_with_frozen_midPrenatal_programs.tsv"))

# ------------------------------------------------------------
# 4. Adult bulk scoring helpers
# ------------------------------------------------------------
load_standard_expr <- function(path) {
  dt <- safe_fread(path)
  if (!"gene_symbol" %in% names(dt)) stop("Expected gene_symbol column in ", path)

  genes <- std_gene(dt$gene_symbol)
  sample_cols <- setdiff(names(dt), "gene_symbol")

  for (cc in sample_cols) {
    if (!is.numeric(dt[[cc]])) {
      suppressWarnings(set(dt, j = cc, value = as.numeric(dt[[cc]])))
    }
  }

  mat <- as.matrix(dt[, sample_cols, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  colnames(mat) <- std_sample(sample_cols)

  keep <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep, , drop = FALSE]

  if (anyDuplicated(rownames(mat))) {
    rs <- rowsum(mat, group = rownames(mat), reorder = FALSE)
    n <- as.numeric(table(rownames(mat))[rownames(rs)])
    mat <- rs / n
  }

  # log2 CPM if raw-count-like
  q99 <- suppressWarnings(quantile(mat, 0.99, na.rm = TRUE))
  if (is.finite(q99) && q99 > 50) {
    lib <- colSums(mat, na.rm = TRUE)
    lib[lib == 0] <- NA
    mat <- log2(t(t(mat) / lib * 1e6) + 1)
  }

  mat
}

load_standard_meta <- function(path, sample_names) {
  md <- safe_fread(path)
  nms <- names(md)

  sample_col <- if ("sample" %in% nms) "sample" else if ("sample_id" %in% nms) "sample_id" else nms[1]
  dx_col <- if ("diagnosis" %in% nms) "diagnosis" else {
    tmp <- nms[grepl("diagnosis|dx|condition|group|phenotype", nms, ignore.case = TRUE)]
    if (length(tmp) == 0) stop("Cannot detect diagnosis column in ", path)
    tmp[1]
  }

  md[, sample := std_sample(get(sample_col))]
  md[, diagnosis := norm_dx(get(dx_col))]
  md <- md[sample %in% sample_names]
  md <- md[match(sample_names[sample_names %in% md$sample], sample)]

  md
}

load_mapping <- function() {
  mp_paths <- c(
    file.path(BASE, "step03b_Gandal_covariate_sensitivity/meta/01_BrainSpan_geneSymbol_to_ensembl.tsv"),
    file.path(BASE, "step03c_Gandal_alternative_scoring_sensitivity/meta/01_BrainSpan_geneSymbol_to_ensembl.tsv")
  )
  mp_path <- mp_paths[file.exists(mp_paths)][1]
  if (is.na(mp_path)) return(NULL)

  mp <- fread(mp_path)
  nms <- names(mp)
  sym_col <- nms[grepl("symbol|gene", nms, ignore.case = TRUE)][1]
  ens_col <- nms[grepl("ensembl|ensg", nms, ignore.case = TRUE)][1]

  if (is.na(sym_col) || is.na(ens_col)) return(NULL)

  unique(data.table(
    gene = std_gene(mp[[sym_col]]),
    ensembl = sub("\\..*$", "", as.character(mp[[ens_col]]))
  ))
}

load_gandal_expr_meta <- function(path) {
  e <- new.env()
  load(path, envir = e)

  if (!all(c("datExpr", "datMeta") %in% ls(e))) {
    stop("Gandal RData must contain datExpr and datMeta.")
  }

  mat_raw <- as.matrix(get("datExpr", e))
  storage.mode(mat_raw) <- "numeric"
  md <- as.data.table(get("datMeta", e))

  raw_ids <- rownames(mat_raw)
  raw_clean <- sub("\\..*$", "", raw_ids)
  direct_symbols <- std_gene(raw_ids)
  direct_overlap <- sum(direct_symbols %in% sfari_genes)

  if (direct_overlap >= 50) {
    sym <- direct_symbols
  } else {
    mp <- load_mapping()
    if (is.null(mp)) stop("Missing mapping for Gandal Ensembl IDs.")
    sym <- mp[match(raw_clean, ensembl)]$gene
  }

  keep <- !is.na(sym) & sym != ""
  mat <- mat_raw[keep, , drop = FALSE]
  sym2 <- sym[keep]

  if (anyDuplicated(sym2)) {
    rs <- rowsum(mat, group = sym2, reorder = FALSE)
    n <- as.numeric(table(sym2)[rownames(rs)])
    mat <- rs / n
  } else {
    rownames(mat) <- sym2
  }

  colnames(mat) <- std_sample(colnames(mat))

  # metadata
  if ("sample_id" %in% names(md)) {
    md[, sample := std_sample(sample_id)]
  } else {
    md[, sample := std_sample(rownames(md))]
  }

  dx_col <- if ("Dx" %in% names(md)) "Dx" else {
    tmp <- names(md)[grepl("diagnosis|dx|condition|group|phenotype", names(md), ignore.case = TRUE)]
    if (length(tmp) == 0) stop("Cannot detect Gandal dx column.")
    tmp[1]
  }
  md[, diagnosis := norm_dx(get(dx_col))]

  md <- md[sample %in% colnames(mat)]
  md <- md[match(colnames(mat)[colnames(mat) %in% md$sample], sample)]

  list(mat = mat, meta = md)
}

score_gene_sets <- function(cohort, mat, md, membership) {
  out <- list()

  sample_keep <- intersect(colnames(mat), md$sample)
  mat <- mat[, sample_keep, drop = FALSE]
  md <- md[match(sample_keep, sample)]

  safe_write(md[, .N, by = diagnosis][order(diagnosis)], file.path(OUT, paste0("QC_", cohort, "_diagnosis_counts.tsv")))

  for (bid in unique(membership$boundary_id)) {
    gs <- unique(membership[boundary_id == bid, gene])
    mapped <- intersect(gs, rownames(mat))
    if (length(mapped) < 3) next

    score <- colMeans(mat[mapped, , drop = FALSE], na.rm = TRUE)
    dt <- data.table(
      cohort = cohort,
      sample = names(score),
      boundary_id = bid,
      score = as.numeric(score),
      n_input = length(gs),
      n_mapped = length(mapped)
    )
    dt <- merge(dt, md[, .(sample, diagnosis)], by = "sample", all.x = TRUE)
    dt <- dt[diagnosis %in% c("ASD", "Control")]

    n_asd <- sum(dt$diagnosis == "ASD")
    n_ctrl <- sum(dt$diagnosis == "Control")

    beta <- se <- p <- mean_asd <- mean_ctrl <- NA_real_
    if (n_asd >= 2 && n_ctrl >= 2) {
      dt[, dx_factor := factor(diagnosis, levels = c("Control", "ASD"))]
      fit <- lm(score ~ dx_factor, data = dt)
      co <- summary(fit)$coefficients
      if ("dx_factorASD" %in% rownames(co)) {
        beta <- co["dx_factorASD", "Estimate"]
        se <- co["dx_factorASD", "Std. Error"]
        p <- co["dx_factorASD", "Pr(>|t|)"]
      }
      mean_asd <- mean(dt[diagnosis == "ASD"]$score, na.rm = TRUE)
      mean_ctrl <- mean(dt[diagnosis == "Control"]$score, na.rm = TRUE)
    }

    out[[paste(cohort, bid, sep = "__")]] <- data.table(
      cohort = cohort,
      boundary_id = bid,
      n_input = length(gs),
      n_mapped = length(mapped),
      n_ASD = n_asd,
      n_Control = n_ctrl,
      mean_ASD = mean_asd,
      mean_Control = mean_ctrl,
      beta_ASD_vs_Control = beta,
      se = se,
      ci_low = beta - 1.96 * se,
      ci_high = beta + 1.96 * se,
      p_lm = p,
      direction = fifelse(beta < 0, "ASD_lower", fifelse(beta > 0, "ASD_higher", "NA"))
    )
  }

  rbindlist(out, fill = TRUE)
}

# ------------------------------------------------------------
# 5. Adult bulk effects for alternative boundary gene sets
# ------------------------------------------------------------
bulk_results <- list()

# GSE102741
mat102 <- load_standard_expr(gse102_expr)
md102 <- load_standard_meta(gse102_meta, colnames(mat102))
bulk_results[["GSE102741"]] <- score_gene_sets("GSE102741", mat102, md102, alt_membership)

# GSE64018
mat640 <- load_standard_expr(gse640_expr)
md640 <- load_standard_meta(gse640_meta, colnames(mat640))
bulk_results[["GSE64018"]] <- score_gene_sets("GSE64018", mat640, md640, alt_membership)

# Gandal2022
gd <- load_gandal_expr_meta(gandal_path)
bulk_results[["Gandal2022"]] <- score_gene_sets("Gandal2022", gd$mat, gd$meta, alt_membership)

bulk_dt <- rbindlist(bulk_results, fill = TRUE)
bulk_dt[, FDR_lm := p.adjust(p_lm, method = "BH")]
safe_write(bulk_dt, file.path(OUT, "07_adult_bulk_effects_alternative_boundary_gene_sets.tsv"))

bulk_summary <- bulk_dt[, .(
  n_cohorts = .N,
  n_direction_ASD_lower = sum(beta_ASD_vs_Control < 0, na.rm = TRUE),
  n_direction_ASD_higher = sum(beta_ASD_vs_Control > 0, na.rm = TRUE),
  mean_beta = mean(beta_ASD_vs_Control, na.rm = TRUE),
  min_p = min(p_lm, na.rm = TRUE)
), by = boundary_id][order(-n_direction_ASD_lower, mean_beta)]

safe_write(bulk_summary, file.path(OUT, "08_adult_bulk_direction_summary_by_boundary.tsv"))

# ------------------------------------------------------------
# 6. Figures
# ------------------------------------------------------------
plot_overlap <- merge(
  overlap_summary,
  boundary_defs[, .(boundary_id, interpretation)],
  by = "boundary_id",
  all.x = TRUE
)

plot_overlap[, boundary_id := factor(boundary_id, levels = overlap_summary$boundary_id)]

p1 <- ggplot(plot_overlap, aes(x = boundary_id, y = jaccard_frozen_top20)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("n=", overlap_frozen_top20)), vjust = -0.3, size = 3) +
  theme_bw(base_size = 11) +
  labs(
    title = "Overlap with frozen midPrenatal_SFARI_top20",
    x = NULL,
    y = "Jaccard overlap"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(file.path(OUT, "Figure_Step60A_overlap_with_frozen_top20.pdf"), p1, width = 8, height = 4.2)
ggsave(file.path(OUT, "Figure_Step60A_overlap_with_frozen_top20.png"), p1, width = 8, height = 4.2, dpi = 300)

plot_bulk <- bulk_dt[!is.na(beta_ASD_vs_Control)]
plot_bulk[, boundary_id := factor(boundary_id, levels = boundary_defs$boundary_id)]

p2 <- ggplot(plot_bulk, aes(x = cohort, y = beta_ASD_vs_Control)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 2.6) +
  facet_wrap(~ boundary_id, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 10) +
  labs(
    title = "Adult bulk effects of alternative BrainSpan prenatal-boundary definitions",
    x = NULL,
    y = expression(beta[ASD]~"(ASD vs Control)")
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8.5),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(OUT, "Figure_Step60A_adult_bulk_effects_alternative_boundaries.pdf"), p2, width = 8.5, height = 7.5)
ggsave(file.path(OUT, "Figure_Step60A_adult_bulk_effects_alternative_boundaries.png"), p2, width = 8.5, height = 7.5, dpi = 300)

# Compact manuscript-facing output
compact <- merge(
  bulk_dt,
  overlap_summary[, .(boundary_id, n_genes, overlap_frozen_top20, pct_overlap_frozen_top20, jaccard_frozen_top20)],
  by = "boundary_id",
  all.x = TRUE
)
compact <- merge(
  compact,
  boundary_defs[, .(boundary_id, score_formula, interpretation)],
  by = "boundary_id",
  all.x = TRUE
)

safe_write(compact, file.path(OUT, "09_manuscript_compact_boundary_sensitivity_summary.tsv"))

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("Key outputs:\n")
cat(" - 02_boundary_window_definitions.tsv\n")
cat(" - 04_alternative_boundary_top20_gene_sets.tsv\n")
cat(" - 06_overlap_with_frozen_midPrenatal_programs.tsv\n")
cat(" - 07_adult_bulk_effects_alternative_boundary_gene_sets.tsv\n")
cat(" - 08_adult_bulk_direction_summary_by_boundary.tsv\n")
cat(" - 09_manuscript_compact_boundary_sensitivity_summary.tsv\n")
cat(" - Figure_Step60A_overlap_with_frozen_top20.pdf/png\n")
cat(" - Figure_Step60A_adult_bulk_effects_alternative_boundaries.pdf/png\n")
sink()
