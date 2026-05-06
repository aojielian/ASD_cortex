suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUT  <- file.path(BASE, "step60_MA_strengthening_v3/results/Step60B_detection_filtered_rescoring_v3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(OUT, "Step60B_v3_run.log")
sink(log_file, split = TRUE)

cat("Step60B v3: detection-filtered adult bulk rescoring using standardized inputs\n")
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

load_standard_expr <- function(path) {
  dt <- safe_fread(path)
  if (!"gene_symbol" %in% names(dt)) {
    stop("Expected gene_symbol column in standardized expression file: ", path)
  }

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

# -------------------------
# Program membership
# -------------------------
membership_path <- file.path(BASE, "reviewer_round1/R1A_sfari_program_overlap_audit/tables/R1A_frozen_program_memberships.tsv")
syn_path <- file.path(BASE, "step01_prepare/Step11J_SFARIall_synaptic_submodule_refinement_v1/tables/01_submodule_membership.tsv.gz")

mem0 <- safe_fread(membership_path)
if (!all(c("program_name", "gene_symbol") %in% names(mem0))) {
  stop("Membership file must contain program_name and gene_symbol.")
}

mem <- unique(mem0[, .(
  program = as.character(program_name),
  gene = std_gene(gene_symbol)
)])

if (file.exists(syn_path)) {
  syn <- fread(syn_path)
  if (all(c("submodule", "gene_symbol") %in% names(syn))) {
    syn2 <- syn[grepl("syn", submodule, ignore.case = TRUE),
                .(program = "SFARI_all_synaptic",
                  gene = std_gene(gene_symbol))]
    mem <- unique(rbindlist(list(mem, syn2), fill = TRUE))
  }
}

keep_programs <- c(
  "SFARI_all",
  "SFARI_all_synaptic",
  "midPrenatal_SFARI_top20",
  "midPrenatal_SFARI_top10",
  "midPrenatal_SFARI_top05"
)

mem <- mem[program %in% keep_programs & !is.na(gene) & gene != ""]
safe_write(mem, file.path(OUT, "00_program_membership_used.tsv"))
safe_write(mem[, .(n_genes = uniqueN(gene)), by = program][order(program)],
           file.path(OUT, "00_program_sizes.tsv"))

# -------------------------
# Main scoring function
# -------------------------
score_cohort <- function(cohort, expr_path, meta_path, detection_mean_min = 0, detection_frac_min = 0.50) {
  cat("\n============================================================\n")
  cat("Scoring cohort:", cohort, "\n")
  cat("expr:", expr_path, "\n")
  cat("meta:", meta_path, "\n")
  cat("============================================================\n")

  mat <- load_standard_expr(expr_path)
  md <- load_standard_meta(meta_path, colnames(mat))

  # Detect whether raw counts need CPM log transform
  q99 <- suppressWarnings(quantile(mat, 0.99, na.rm = TRUE))
  if (is.finite(q99) && q99 > 50) {
    lib <- colSums(mat, na.rm = TRUE)
    lib[lib == 0] <- NA
    mat_score <- log2(t(t(mat) / lib * 1e6) + 1)
    transform_used <- "log2_CPM_plus1"
  } else {
    mat_score <- mat
    transform_used <- "as_loaded"
  }

  dx_count <- md[, .N, by = diagnosis][order(diagnosis)]
  safe_write(dx_count, file.path(OUT, paste0("QC_", cohort, "_diagnosis_counts.tsv")))

  sample_keep <- intersect(colnames(mat_score), md$sample)
  mat_score <- mat_score[, sample_keep, drop = FALSE]
  md <- md[match(sample_keep, sample)]

  score_list <- list()
  summary_list <- list()
  gene_qc_list <- list()

  for (prog in unique(mem$program)) {
    genes_frozen <- unique(mem[program == prog, gene])
    genes_mapped <- intersect(genes_frozen, rownames(mat_score))

    if (length(genes_mapped) == 0) {
      summary_list[[paste0(prog, "_nomap")]] <- data.table(
        cohort = cohort,
        program = prog,
        score_mode = "mapped_all",
        n_frozen = length(genes_frozen),
        n_mapped = 0,
        n_detected = 0,
        detection_retention = NA_real_,
        n_scored = 0,
        n_ASD = sum(md$diagnosis == "ASD", na.rm = TRUE),
        n_Control = sum(md$diagnosis == "Control", na.rm = TRUE),
        mean_ASD = NA_real_,
        mean_Control = NA_real_,
        beta_ASD_vs_Control = NA_real_,
        p_lm = NA_real_,
        p_wilcox = NA_real_,
        transform_used = transform_used
      )
      next
    }

    expr_sub <- mat_score[genes_mapped, , drop = FALSE]
    gene_mean <- rowMeans(expr_sub, na.rm = TRUE)
    gene_frac <- rowMeans(expr_sub > detection_mean_min, na.rm = TRUE)
    genes_detected <- names(gene_mean)[gene_mean > detection_mean_min & gene_frac >= detection_frac_min]

    gene_qc_list[[prog]] <- data.table(
      cohort = cohort,
      program = prog,
      gene = genes_mapped,
      mean_expression = as.numeric(gene_mean[genes_mapped]),
      frac_samples_detected = as.numeric(gene_frac[genes_mapped]),
      detected = genes_mapped %in% genes_detected
    )

    for (mode in c("mapped_all", "detected_only")) {
      use_genes <- if (mode == "mapped_all") genes_mapped else genes_detected
      if (length(use_genes) < 3) next

      score <- colMeans(mat_score[use_genes, , drop = FALSE], na.rm = TRUE)

      sc <- data.table(
        cohort = cohort,
        sample = names(score),
        program = prog,
        score_mode = mode,
        score = as.numeric(score),
        n_scored = length(use_genes),
        n_mapped = length(genes_mapped),
        n_detected = length(genes_detected)
      )

      sc <- merge(sc, md[, .(sample, diagnosis)], by = "sample", all.x = TRUE)
      score_list[[paste(prog, mode, sep = "__")]] <- sc

      sc2 <- sc[diagnosis %in% c("ASD", "Control")]
      n_asd <- sum(sc2$diagnosis == "ASD", na.rm = TRUE)
      n_ctrl <- sum(sc2$diagnosis == "Control", na.rm = TRUE)

      mean_asd <- mean(sc2[diagnosis == "ASD"]$score, na.rm = TRUE)
      mean_ctrl <- mean(sc2[diagnosis == "Control"]$score, na.rm = TRUE)

      beta <- p_lm <- p_wilcox <- NA_real_

      if (n_asd >= 2 && n_ctrl >= 2) {
        sc2[, dx_factor := factor(diagnosis, levels = c("Control", "ASD"))]
        fit <- lm(score ~ dx_factor, data = sc2)
        co <- summary(fit)$coefficients
        if ("dx_factorASD" %in% rownames(co)) {
          beta <- co["dx_factorASD", "Estimate"]
          p_lm <- co["dx_factorASD", "Pr(>|t|)"]
        }
        p_wilcox <- tryCatch(
          wilcox.test(score ~ diagnosis, data = sc2)$p.value,
          error = function(e) NA_real_
        )
      }

      summary_list[[paste(prog, mode, sep = "__")]] <- data.table(
        cohort = cohort,
        program = prog,
        score_mode = mode,
        n_frozen = length(genes_frozen),
        n_mapped = length(genes_mapped),
        n_detected = length(genes_detected),
        detection_retention = length(genes_detected) / length(genes_mapped),
        n_scored = length(use_genes),
        n_ASD = n_asd,
        n_Control = n_ctrl,
        mean_ASD = mean_asd,
        mean_Control = mean_ctrl,
        beta_ASD_vs_Control = beta,
        p_lm = p_lm,
        p_wilcox = p_wilcox,
        transform_used = transform_used
      )
    }
  }

  list(
    summary = rbindlist(summary_list, fill = TRUE),
    scores = rbindlist(score_list, fill = TRUE),
    gene_qc = rbindlist(gene_qc_list, fill = TRUE)
  )
}

# -------------------------
# Run cohorts
# -------------------------
inputs <- data.table(
  cohort = c("GSE102741", "GSE64018"),
  expr_path = c(
    file.path(BASE, "step01_prepare/inputs/GSE102741.standardized_expr.tsv.gz"),
    file.path(BASE, "step01_prepare/inputs/GSE64018.standardized_expr.tsv.gz")
  ),
  meta_path = c(
    file.path(BASE, "step01_prepare/inputs/GSE102741.standardized_meta.tsv"),
    file.path(BASE, "step01_prepare/inputs/GSE64018.standardized_meta.tsv")
  )
)

inputs[, expr_exists := file.exists(expr_path)]
inputs[, meta_exists := file.exists(meta_path)]
safe_write(inputs, file.path(OUT, "01_input_manifest.tsv"))

all_summary <- list()
all_scores <- list()
all_gene_qc <- list()

for (i in seq_len(nrow(inputs))) {
  if (!inputs$expr_exists[i] || !inputs$meta_exists[i]) {
    cat("[WARN] Missing input for ", inputs$cohort[i], "\n")
    next
  }

  res <- score_cohort(
    cohort = inputs$cohort[i],
    expr_path = inputs$expr_path[i],
    meta_path = inputs$meta_path[i]
  )

  all_summary[[inputs$cohort[i]]] <- res$summary
  all_scores[[inputs$cohort[i]]] <- res$scores
  all_gene_qc[[inputs$cohort[i]]] <- res$gene_qc
}

summary_dt <- rbindlist(all_summary, fill = TRUE)
scores_dt <- rbindlist(all_scores, fill = TRUE)
gene_qc_dt <- rbindlist(all_gene_qc, fill = TRUE)

if (nrow(summary_dt) > 0) {
  summary_dt[, FDR_lm := p.adjust(p_lm, method = "BH")]
  summary_dt[, direction := fifelse(beta_ASD_vs_Control < 0, "ASD_lower",
                                    fifelse(beta_ASD_vs_Control > 0, "ASD_higher", "NA"))]
}

safe_write(summary_dt, file.path(OUT, "10_detection_filtered_bulk_model_summary_v3.tsv"))
safe_write(scores_dt, file.path(OUT, "11_detection_filtered_bulk_sample_scores_v3.tsv.gz"))
safe_write(gene_qc_dt, file.path(OUT, "12_detection_gene_qc_v3.tsv.gz"))

wide <- dcast(
  summary_dt,
  cohort + program ~ score_mode,
  value.var = c(
    "n_scored",
    "beta_ASD_vs_Control",
    "p_lm",
    "FDR_lm",
    "detection_retention",
    "n_ASD",
    "n_Control"
  )
)
safe_write(wide, file.path(OUT, "13_detection_filter_comparison_wide_v3.tsv"))

# Compact manuscript-facing table
compact <- summary_dt[
  program %in% c("SFARI_all", "SFARI_all_synaptic", "midPrenatal_SFARI_top20"),
  .(
    cohort,
    program,
    score_mode,
    n_mapped,
    n_detected,
    detection_retention,
    beta_ASD_vs_Control,
    p_lm,
    FDR_lm,
    direction
  )
]
safe_write(compact, file.path(OUT, "14_manuscript_compact_detection_filter_summary.tsv"))

# Plot
plot_dt <- summary_dt[
  program %in% c("SFARI_all", "SFARI_all_synaptic", "midPrenatal_SFARI_top20") &
    !is.na(beta_ASD_vs_Control)
]

if (nrow(plot_dt) > 0) {
  p <- ggplot(plot_dt, aes(x = cohort, y = beta_ASD_vs_Control, shape = score_mode)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(size = 3, position = position_dodge(width = 0.45)) +
    facet_wrap(~ program, scales = "free_y") +
    theme_bw(base_size = 11) +
    labs(
      title = "Detection-filtered adult bulk rescoring",
      x = NULL,
      y = expression(beta[ASD]~"(ASD vs Control)"),
      shape = "Scoring set"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  ggsave(file.path(OUT, "Figure_Step60B_detection_filtered_bulk_rescoring_v3.pdf"),
         p, width = 8.5, height = 4.5)
  ggsave(file.path(OUT, "Figure_Step60B_detection_filtered_bulk_rescoring_v3.png"),
         p, width = 8.5, height = 4.5, dpi = 300)
}

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("Outputs in:", OUT, "\n")
sink()
