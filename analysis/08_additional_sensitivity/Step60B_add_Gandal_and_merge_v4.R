suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUT  <- file.path(BASE, "step60_MA_strengthening_v3/results/Step60B_detection_filtered_rescoring_v4_addGandal")
PREV <- file.path(BASE, "step60_MA_strengthening_v3/results/Step60B_detection_filtered_rescoring_v3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(OUT, "Step60B_v4_addGandal_run.log")
sink(log_file, split = TRUE)

cat("Step60B v4: add Gandal2022 detection-filtered rescoring and merge with v3\n")
cat("Started:", as.character(Sys.time()), "\n")
cat("BASE:", BASE, "\n")
cat("OUT :", OUT, "\n\n")

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

pick_existing <- function(paths) {
  x <- paths[file.exists(paths)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

load_mapping <- function() {
  mp_path <- pick_existing(c(
    file.path(BASE, "step03b_Gandal_covariate_sensitivity/meta/01_BrainSpan_geneSymbol_to_ensembl.tsv"),
    file.path(BASE, "step03c_Gandal_alternative_scoring_sensitivity/meta/01_BrainSpan_geneSymbol_to_ensembl.tsv")
  ))

  if (is.na(mp_path)) {
    cat("[WARN] No local BrainSpan symbol-Ensembl mapping found.\n")
    return(NULL)
  }

  cat("[READ mapping]", mp_path, "\n")
  mp <- fread(mp_path)
  nms <- names(mp)

  sym_col <- nms[grepl("symbol|gene", nms, ignore.case = TRUE)][1]
  ens_col <- nms[grepl("ensembl|ensg", nms, ignore.case = TRUE)][1]

  if (is.na(sym_col) || is.na(ens_col)) {
    cat("[WARN] Could not detect symbol/ensembl columns in mapping file.\n")
    print(names(mp))
    return(NULL)
  }

  ans <- unique(data.table(
    gene_symbol = std_gene(mp[[sym_col]]),
    ensembl = sub("\\..*$", "", as.character(mp[[ens_col]]))
  ))
  ans <- ans[!is.na(gene_symbol) & gene_symbol != "" & !is.na(ensembl) & ensembl != ""]
  ans
}

map_gandal_rows_to_symbol <- function(mat, program_genes) {
  raw_ids <- rownames(mat)
  raw_clean <- sub("\\..*$", "", raw_ids)

  direct_symbols <- std_gene(raw_ids)
  direct_overlap <- sum(direct_symbols %in% program_genes)

  cat("Direct symbol overlap with program genes:", direct_overlap, "\n")

  if (direct_overlap >= 50) {
    sym <- direct_symbols
  } else {
    mp <- load_mapping()
    if (is.null(mp)) {
      stop("Cannot map Gandal rownames to gene symbols and direct overlap is too low.")
    }

    mm <- mp[match(raw_clean, ensembl)]
    sym <- mm$gene_symbol

    mapped_overlap <- sum(sym %in% program_genes, na.rm = TRUE)
    cat("Ensembl-mapped overlap with program genes:", mapped_overlap, "\n")

    if (mapped_overlap < 50) {
      stop("Gandal gene mapping failed: mapped overlap too low: ", mapped_overlap)
    }
  }

  keep <- !is.na(sym) & sym != ""
  mat2 <- mat[keep, , drop = FALSE]
  sym2 <- sym[keep]

  if (anyDuplicated(sym2)) {
    rs <- rowsum(mat2, group = sym2, reorder = FALSE)
    n <- as.numeric(table(sym2)[rownames(rs)])
    mat2 <- rs / n
  } else {
    rownames(mat2) <- sym2
  }

  mat2
}

detect_sample_col <- function(md, sample_names) {
  nms <- names(md)
  cand <- nms[grepl("sample|id|subject|donor|specimen|individual", nms, ignore.case = TRUE)]
  if (length(cand) == 0) cand <- nms

  best <- NA_character_
  best_n <- -1
  for (cc in cand) {
    v <- std_sample(md[[cc]])
    n <- sum(v %in% sample_names)
    if (n > best_n) {
      best <- cc
      best_n <- n
    }
  }

  if (best_n <= 0) return(NA_character_)
  best
}

detect_dx_col <- function(md) {
  nms <- names(md)
  preferred <- c("Dx", "dx", "diagnosis", "Diagnosis", "condition", "group", "phenotype")
  hit <- preferred[preferred %in% nms]
  if (length(hit) > 0) return(hit[1])

  cand <- nms[grepl("diagnosis|dx|condition|group|phenotype|disease|disorder|case", nms, ignore.case = TRUE)]
  if (length(cand) > 0) return(cand[1])

  NA_character_
}

score_gandal <- function(mat_score, md, mem, detection_mean_min = 0, detection_frac_min = 0.50) {
  cohort <- "Gandal2022"

  # Gandal normalized datExpr usually does not need CPM transform
  q99 <- suppressWarnings(quantile(mat_score, 0.99, na.rm = TRUE))
  transform_used <- ifelse(is.finite(q99) && q99 > 50, "as_loaded_high_range_warning", "as_loaded")

  sample_names <- std_sample(colnames(mat_score))
  colnames(mat_score) <- sample_names

  scol <- detect_sample_col(md, sample_names)
  if (is.na(scol)) {
    if (!is.null(rownames(md)) && sum(std_sample(rownames(md)) %in% sample_names) > 0) {
      md[, sample := std_sample(rownames(md))]
    } else {
      stop("Could not align Gandal datMeta to datExpr columns.")
    }
  } else {
    md[, sample := std_sample(get(scol))]
  }

  dxcol <- detect_dx_col(md)
  if (is.na(dxcol)) stop("Could not detect Gandal diagnosis column.")

  md[, diagnosis := norm_dx(get(dxcol))]
  md <- md[sample %in% sample_names]
  md <- md[match(sample_names[sample_names %in% md$sample], sample)]

  safe_write(data.table(column_type = c("sample_col", "dx_col"),
                        column_name = c(ifelse(is.na(scol), "rownames(datMeta)", scol), dxcol)),
             file.path(OUT, "QC_Gandal2022_selected_metadata_columns.tsv"))

  safe_write(md[, .N, by = diagnosis][order(diagnosis)],
             file.path(OUT, "QC_Gandal2022_diagnosis_counts.tsv"))

  out_sum <- list()
  out_scores <- list()
  out_gene <- list()

  for (prog in unique(mem$program)) {
    genes_frozen <- unique(mem[program == prog, gene])
    genes_mapped <- intersect(genes_frozen, rownames(mat_score))

    if (length(genes_mapped) == 0) next

    expr_sub <- mat_score[genes_mapped, , drop = FALSE]
    gene_mean <- rowMeans(expr_sub, na.rm = TRUE)
    gene_frac <- rowMeans(expr_sub > detection_mean_min, na.rm = TRUE)
    genes_detected <- names(gene_mean)[gene_mean > detection_mean_min & gene_frac >= detection_frac_min]

    out_gene[[prog]] <- data.table(
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
      out_scores[[paste(prog, mode, sep = "__")]] <- sc

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

      out_sum[[paste(prog, mode, sep = "__")]] <- data.table(
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
    summary = rbindlist(out_sum, fill = TRUE),
    scores = rbindlist(out_scores, fill = TRUE),
    gene_qc = rbindlist(out_gene, fill = TRUE)
  )
}

# -------------------------
# Load membership from v3 if available
# -------------------------
mem_path <- file.path(PREV, "00_program_membership_used.tsv")
if (!file.exists(mem_path)) {
  stop("Missing v3 membership file: ", mem_path)
}
mem <- fread(mem_path)
mem[, gene := std_gene(gene)]

program_genes <- unique(mem$gene)

# -------------------------
# Load Gandal
# -------------------------
gandal_path <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"
if (!file.exists(gandal_path)) {
  stop("Missing Gandal RData: ", gandal_path)
}

e <- new.env()
load(gandal_path, envir = e)

if (!all(c("datExpr", "datMeta") %in% ls(e))) {
  stop("Gandal RData must contain datExpr and datMeta. Found: ", paste(ls(e), collapse = ", "))
}

mat_raw <- as.matrix(get("datExpr", e))
storage.mode(mat_raw) <- "numeric"
md <- as.data.table(get("datMeta", e))

cat("Gandal raw datExpr dim:", paste(dim(mat_raw), collapse = " x "), "\n")
cat("Gandal datMeta dim:", paste(dim(md), collapse = " x "), "\n")

mat_sym <- map_gandal_rows_to_symbol(mat_raw, program_genes)
cat("Gandal symbol-level matrix dim:", paste(dim(mat_sym), collapse = " x "), "\n")

res_g <- score_gandal(mat_sym, md, mem)

safe_write(res_g$summary, file.path(OUT, "20_Gandal2022_detection_filtered_model_summary_v4.tsv"))
safe_write(res_g$scores, file.path(OUT, "21_Gandal2022_detection_filtered_sample_scores_v4.tsv.gz"))
safe_write(res_g$gene_qc, file.path(OUT, "22_Gandal2022_detection_gene_qc_v4.tsv.gz"))

# -------------------------
# Merge with v3 GSE102741/GSE64018 results
# -------------------------
v3_summary_path <- file.path(PREV, "10_detection_filtered_bulk_model_summary_v3.tsv")
if (!file.exists(v3_summary_path)) {
  stop("Missing v3 summary: ", v3_summary_path)
}

v3 <- fread(v3_summary_path)
combined <- rbindlist(list(v3, res_g$summary), fill = TRUE)

combined[, FDR_lm_global := p.adjust(p_lm, method = "BH")]
combined[, direction := fifelse(beta_ASD_vs_Control < 0, "ASD_lower",
                                fifelse(beta_ASD_vs_Control > 0, "ASD_higher", "NA"))]

safe_write(combined, file.path(OUT, "30_combined_three_cohort_detection_filtered_summary_v4.tsv"))

wide <- dcast(
  combined,
  cohort + program ~ score_mode,
  value.var = c(
    "n_scored",
    "beta_ASD_vs_Control",
    "p_lm",
    "FDR_lm_global",
    "detection_retention",
    "n_ASD",
    "n_Control"
  )
)
safe_write(wide, file.path(OUT, "31_combined_three_cohort_detection_filter_wide_v4.tsv"))

compact <- combined[
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
    FDR_lm_global,
    direction
  )
]
safe_write(compact, file.path(OUT, "32_manuscript_compact_three_cohort_detection_summary_v4.tsv"))

plot_dt <- combined[
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
      title = "Detection-filtered adult bulk rescoring across three cohorts",
      x = NULL,
      y = expression(beta[ASD]~"(ASD vs Control)"),
      shape = "Scoring set"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )

  ggsave(file.path(OUT, "Figure_Step60B_three_cohort_detection_filtered_rescoring_v4.pdf"),
         p, width = 9, height = 4.6)
  ggsave(file.path(OUT, "Figure_Step60B_three_cohort_detection_filtered_rescoring_v4.png"),
         p, width = 9, height = 4.6, dpi = 300)
}

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("Outputs in:", OUT, "\n")
sink()
