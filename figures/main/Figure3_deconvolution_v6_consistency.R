#!/usr/bin/env Rscript

# Figure 3: bulk deconvolution and composition-adjusted interpretation
# Standalone, robust version.
# Key fix: Panel C explicitly requires and plots BOTH unadjusted and NNLS-adjusted
# program effects for SFARI_all and midPrenatal_SFARI_top20.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    base_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    outdir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/Figures_main",
    dx_file = NA_character_,
    unadj_file = NA_character_,
    adj_file = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args)) stop("Missing value for ", key)
    val <- args[[i + 1]]
    out[[sub("^--", "", key)]] <- val
    i <- i + 2
  }
  out
}

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(6, 8, 6, 6)
  )

find_first_existing <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

scan_candidate_files <- function(search_dir, patterns, exts = c("tsv","tsv.gz","txt","txt.gz","csv","csv.gz")) {
  if (is.null(search_dir) || is.na(search_dir) || !dir.exists(search_dir)) return(character(0))
  all_files <- list.files(search_dir, full.names = TRUE, recursive = FALSE)
  if (length(all_files) == 0) return(character(0))
  bn <- basename(all_files)
  keep_ext <- Reduce(`|`, lapply(exts, function(ext) {
    grepl(paste0("\\.", gsub("\\.", "\\\\.", ext), "$"), bn, ignore.case = TRUE)
  }))
  cand <- all_files[keep_ext]
  if (length(cand) == 0) return(character(0))
  score <- rep(0L, length(cand))
  for (pat in patterns) score <- score + as.integer(grepl(pat, basename(cand), ignore.case = TRUE, perl = TRUE))
  cand <- cand[score > 0]
  score <- score[score > 0]
  if (length(cand) == 0) return(character(0))
  cand[order(-score, basename(cand))]
}

pick_input_file <- function(opt_value, candidate_paths, label,
                            search_dir = NULL, scan_patterns = NULL, required = TRUE) {
  if (!is.na(opt_value) && nzchar(opt_value)) {
    if (!file.exists(opt_value)) stop(label, " file not found: ", opt_value)
    return(opt_value)
  }
  x <- find_first_existing(candidate_paths)
  if (!is.na(x)) return(x)
  scanned <- character(0)
  if (!is.null(search_dir) && !is.null(scan_patterns)) {
    scanned <- scan_candidate_files(search_dir, scan_patterns)
    if (length(scanned) > 0) return(scanned[1])
  }
  if (!required) return(NA_character_)
  tried <- candidate_paths
  if (length(scanned) > 0) tried <- c(tried, "", "[scanned candidates]", scanned)
  stop("Could not auto-detect ", label, " file.\nTried:\n", paste(tried, collapse = "\n"))
}

read_tsv <- function(path) {
  if (is.na(path) || !nzchar(path)) stop("Tried to read empty path.")
  fread(path, sep = "\t", header = TRUE, fill = TRUE, quote = "", showProgress = FALSE)
}

rename_first <- function(dt, target, candidates) {
  hit <- intersect(candidates, names(dt))
  if (length(hit) > 0 && hit[1] != target) setnames(dt, hit[1], target)
  invisible(dt)
}

coerce_numeric_if_present <- function(dt, cols) {
  for (nm in intersect(cols, names(dt))) {
    dt[, (nm) := suppressWarnings(as.numeric(get(nm)))]
  }
  dt
}

std_program <- function(x) {
  y <- as.character(x)
  y <- gsub("^sfari_all$", "SFARI_all", y, ignore.case = TRUE)
  y <- gsub("^midprenatal_sfari_top20$", "midPrenatal_SFARI_top20", y, ignore.case = TRUE)
  y
}

std_cohort <- function(x) {
  y <- as.character(x)
  y <- gsub("^gse102741$", "GSE102741", y, ignore.case = TRUE)
  y <- gsub("^gse64018$", "GSE64018", y, ignore.case = TRUE)
  y <- gsub("^gandal2022$|^gandal$", "Gandal2022", y, ignore.case = TRUE)
  y
}

make_placeholder_plot <- function(title, subtitle = "Input table unavailable") {
  ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
    geom_text(label = subtitle, size = 4.5) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))
}

normalize_effect_table <- function(dt) {
  x <- copy(dt)
  rename_first(x, "program", c("program","program_name"))
  rename_first(x, "cohort", c("cohort","dataset"))
  rename_first(x, "model", c("model","adjustment","analysis_model"))
  x <- coerce_numeric_if_present(x, c(
    "delta_asd_minus_control", "delta", "asd_minus_control",
    "beta_asd", "lm_beta_asd", "effect",
    "delta_asd_minus_control_adjusted", "delta_asd_minus_control_adjustedSet",
    "lm_beta_asd_adjusted"
  ))
  if ("program" %in% names(x)) x[, program := std_program(program)]
  if ("cohort" %in% names(x)) x[, cohort := std_cohort(cohort)]
  x
}

extract_from_same_table <- function(x, src_label = "") {
  out <- list()

  # explicit paired columns in same table
  if (all(c("program", "cohort") %in% names(x))) {
    raw_cols <- intersect(c("delta_asd_minus_control", "delta", "asd_minus_control", "beta_asd", "lm_beta_asd", "effect"), names(x))
    adj_cols <- intersect(c("delta_asd_minus_control_adjusted", "delta_asd_minus_control_adjustedSet", "lm_beta_asd_adjusted"), names(x))

    if (length(raw_cols) > 0 && length(adj_cols) > 0) {
      out[[1]] <- x[, .(program, cohort, value = as.numeric(get(raw_cols[1])), model = "Unadjusted")]
      out[[2]] <- x[, .(program, cohort, value = as.numeric(get(adj_cols[1])), model = "NNLS-adjusted")]
      ans <- rbindlist(out, fill = TRUE)
      return(ans)
    }
  }

  # long table with model column
  if (all(c("program","cohort","model") %in% names(x))) {
    value_cols <- intersect(c("delta_asd_minus_control","delta","asd_minus_control","beta_asd","lm_beta_asd","effect",
                              "delta_asd_minus_control_adjusted","delta_asd_minus_control_adjustedSet","lm_beta_asd_adjusted"), names(x))
    if (length(value_cols) > 0) {
      ans <- x[, .(program, cohort, value = as.numeric(get(value_cols[1])), model = as.character(model))]
      ans[, model := fifelse(grepl("adjust|nnls|composition", model, ignore.case = TRUE), "NNLS-adjusted", "Unadjusted")]
      return(ans)
    }
  }

  # simple table: infer label from source name
  if (all(c("program","cohort") %in% names(x))) {
    value_cols <- intersect(c("delta_asd_minus_control","delta","asd_minus_control","beta_asd","lm_beta_asd","effect"), names(x))
    if (length(value_cols) > 0) {
      inferred_model <- if (grepl("adjust|nnls|prepost|composition", basename(src_label), ignore.case = TRUE)) {
        "NNLS-adjusted"
      } else {
        "Unadjusted"
      }
      ans <- x[, .(program, cohort, value = as.numeric(get(value_cols[1])), model = inferred_model)]
      return(ans)
    }
  }

  data.table()
}

extract_effects <- function(path) {
  x <- normalize_effect_table(read_tsv(path))
  ans <- extract_from_same_table(x, src_label = path)
  ans <- ans[program %in% c("SFARI_all","midPrenatal_SFARI_top20") &
               cohort %in% c("GSE102741","GSE64018","Gandal2022")]
  ans
}

main <- function() {
  opt <- parse_args()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  tbl_dir <- file.path(opt$base_dir, "tables")

  dx_file <- pick_input_file(
    opt$dx_file,
    c(
      file.path(tbl_dir, "142_bulk_celltype_dx_tests_nnls.tsv"),
      file.path(tbl_dir, "142_bulk_celltype_dx_tests_nnls.tsv.gz"),
      file.path(tbl_dir, "140_bulk_celltype_dx_tests_nnls.tsv"),
      file.path(tbl_dir, "140_bulk_celltype_dx_tests_nnls.tsv.gz")
    ),
    "NNLS cell-type diagnosis"
  )

  # Prefer direct per-cohort unadjusted effects table, not meta summary
  unadj_file <- pick_input_file(
    opt$unadj_file,
    c(
      file.path(tbl_dir, "126_bulk_threeCohort_program_effects_updated.tsv"),
      file.path(tbl_dir, "126_bulk_threeCohort_program_effects_updated.tsv.gz"),
      file.path(tbl_dir, "103_bulk_threeCohort_program_effects.tsv"),
      file.path(tbl_dir, "103_bulk_threeCohort_program_effects.tsv.gz"),
      file.path(tbl_dir, "127_bulk_program_meta_summary_updated.tsv")
    ),
    "unadjusted bulk program effects",
    search_dir = tbl_dir,
    scan_patterns = c("threecohort.*program.*effect", "bulk.*program.*effect", "unadjust"),
    required = TRUE
  )

  adj_file <- pick_input_file(
    opt$adj_file,
    c(
      file.path(tbl_dir, "143_bulk_program_effects_prepost_nnls.tsv"),
      file.path(tbl_dir, "143_bulk_program_effects_prepost_nnls.tsv.gz"),
      file.path(tbl_dir, "141_bulk_program_nnls_adjusted.tsv"),
      file.path(tbl_dir, "141_bulk_program_nnls_adjusted.tsv.gz"),
      file.path(tbl_dir, "128_bulk_program_prepost_nnls_updated.tsv"),
      file.path(tbl_dir, "128_bulk_program_prepost_nnls_updated.tsv.gz")
    ),
    "NNLS-adjusted bulk program effects",
    search_dir = tbl_dir,
    scan_patterns = c("nnls", "adjust", "prepost", "composition"),
    required = TRUE
  )

  dx <- read_tsv(dx_file)
  rename_first(dx, "cell_type", c("cell_type","celltype"))
  rename_first(dx, "delta_asd_minus_control", c("delta_asd_minus_control","delta","asd_minus_control"))
  rename_first(dx, "cohort", c("cohort","dataset"))
  rename_first(dx, "fdr_lm", c("fdr_lm","FDR","fdr"))
  dx <- coerce_numeric_if_present(dx, c("delta_asd_minus_control","fdr_lm"))
  if ("cohort" %in% names(dx)) dx[, cohort := std_cohort(cohort)]

  A <- make_placeholder_plot("Cell-type composition shifts")
  if (all(c("cohort","cell_type","delta_asd_minus_control") %in% names(dx)) && nrow(dx) > 0) {
    dxA <- copy(dx)
    dxA[, cohort := factor(cohort, levels = c("GSE102741","GSE64018","Gandal2022"))]
    cell_levels <- c("END","MG","ODC","OPC","AST","EXN","INN")
    dxA[, cell_type := factor(cell_type, levels = rev(cell_levels))]
    A <- ggplot(dxA, aes(cohort, cell_type, fill = delta_asd_minus_control)) +
      geom_tile(color = "white", linewidth = 0.7) +
      scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#E15759", midpoint = 0,
                           name = "ASD - Control") +
      labs(title = "Cell-type composition shifts", x = NULL, y = NULL) +
      base_theme
  }

  B <- make_placeholder_plot("Top composition shifts")
  if (all(c("cohort","cell_type","delta_asd_minus_control","fdr_lm") %in% names(dx)) && nrow(dx) > 0) {
    Bdt <- copy(dx)
    Bdt[, label := paste(cohort, cell_type, sep = " | ")]
    Bdt <- Bdt[order(-abs(delta_asd_minus_control), fdr_lm)][1:min(.N, 10)]
    B <- ggplot(Bdt, aes(delta_asd_minus_control, reorder(label, delta_asd_minus_control),
                         size = -log10(pmax(fdr_lm, 1e-300)), color = cohort)) +
      geom_point(alpha = 0.95) +
      scale_color_manual(values = c("GSE102741" = "#4C78A8", "GSE64018" = "#F28E2B", "Gandal2022" = "#E15759")) +
      labs(title = "Top composition shifts", x = "ASD - Control", y = NULL, size = expression(-log[10](FDR)), color = "Cohort") +
      base_theme
  }

  raw_long <- extract_effects(unadj_file)
  raw_long[, model := "Unadjusted"]

  adj_long <- extract_effects(adj_file)
  adj_long <- adj_long[model == "NNLS-adjusted"]

  if (nrow(raw_long) == 0) stop("Could not extract unadjusted cohort effects from: ", unadj_file)
  if (nrow(adj_long) == 0) stop("Could not extract NNLS-adjusted cohort effects from: ", adj_file)

  Cadj <- rbind(raw_long, adj_long, fill = TRUE)
  Cadj <- unique(Cadj[, .(program, cohort, value, model)])
  Cadj[, cohort := factor(cohort, levels = c("GSE102741","GSE64018","Gandal2022"))]
  Cadj[, program := factor(program, levels = c("SFARI_all","midPrenatal_SFARI_top20"),
                           labels = c("SFARI_all", "midPrenatal_SFARI_top20"))]
  Cadj[, model := factor(model, levels = c("Unadjusted","NNLS-adjusted"))]

  C <- ggplot(Cadj, aes(x = value, y = cohort, shape = model, color = model)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
    geom_point(size = 3.2) +
    facet_wrap(~program, nrow = 1, scales = "free_x") +
    scale_shape_manual(values = c("Unadjusted" = 16, "NNLS-adjusted" = 17)) +
    scale_color_manual(values = c("Unadjusted" = "#6B6B6B", "NNLS-adjusted" = "#2C7FB8")) +
    labs(title = "Program effects before and after NNLS adjustment",
         x = expression(beta[ASD]), y = "Cohort", shape = "Model", color = "Model") +
    base_theme +
    theme(strip.background = element_rect(fill = "grey95", color = "grey80"),
          strip.text = element_text(face = "bold", size = 11))

  fig <- (A | B) / C +
    plot_layout(heights = c(1, 1.1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

  pdf_file <- file.path(opt$outdir, "Figure3_deconvolution_v6_consistency.pdf")
  png_file <- file.path(opt$outdir, "Figure3_deconvolution_v6_consistency.png")
  ggsave(pdf_file, fig, width = 14, height = 10.5)
  ggsave(png_file, fig, width = 14, height = 10.5, dpi = 300)

  cat("[DONE] Wrote:\n", pdf_file, "\n", png_file, "\n",
      "[INPUTS]\n", dx_file, "\n", unadj_file, "\n", adj_file, "\n")
}

main()
