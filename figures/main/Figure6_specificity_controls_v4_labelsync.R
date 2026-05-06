#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    base_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    outdir   = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/Figures_main"
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

read_tsv_flex <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  fread(path, sep = "\t", header = TRUE, fill = TRUE, quote = "", showProgress = FALSE)
}

norm_name <- function(x) {
  x <- tolower(x)
  gsub("[^a-z0-9]+", "_", x)
}

pick_col <- function(dt, patterns, required = TRUE, label = "column") {
  nms <- names(dt)
  nn  <- norm_name(nms)
  for (pat in patterns) {
    idx <- grep(pat, nn, perl = TRUE)
    if (length(idx) >= 1) return(nms[idx[1]])
  }
  if (required) stop("Could not identify ", label, ". Available columns: ", paste(nms, collapse = ", "))
  NULL
}

standardize_program <- function(x) {
  y <- as.character(x)
  y <- gsub("^midprenatal_sfari_top20$", "midPrenatal_SFARI_top20", y, ignore.case = TRUE)
  y <- gsub("^midprenatal_sfari_top10$", "midPrenatal_SFARI_top10", y, ignore.case = TRUE)
  y <- gsub("^midprenatal_sfari_top05$", "midPrenatal_SFARI_top05", y, ignore.case = TRUE)
  y <- gsub("^sfari_all$", "SFARI_all", y, ignore.case = TRUE)
  y
}

program_short <- function(x) {
  out <- x
  out[x == "SFARI_all"] <- "SFARI_all"
  out[x == "midPrenatal_SFARI_top20"] <- "midPrenatal_SFARI_top20"
  out[x == "midPrenatal_SFARI_top10"] <- "midPrenatal_SFARI_top10"
  out[x == "midPrenatal_SFARI_top05"] <- "midPrenatal_SFARI_top05"
  out
}

standardize_source <- function(x) {
  y <- as.character(x)
  yy <- norm_name(y)
  out <- y
  out[grepl("syn", yy)] <- "SynGO"
  out[grepl("reactome", yy)] <- "Reactome"
  out
}

standardize_control <- function(x) {
  y <- as.character(x)
  yy <- norm_name(y)
  out <- y
  out[grepl("size", yy)] <- "size-matched"
  out[grepl("expression", yy)] <- "expression-matched"
  out
}

metric_info <- function(raw_name) {
  y <- norm_name(raw_name)
  out_metric <- rep(NA_character_, length(y))
  out_rank <- rep(99L, length(y))

  out_metric[grepl("^n_sig_terms_", y)] <- "Sig terms"
  out_rank[grepl("^n_sig_terms_", y)] <- 1L

  out_metric[grepl("^best_neglog10q_", y)] <- "Best q"
  out_rank[grepl("^best_neglog10q_", y)] <- 1L

  out_metric[grepl("^best_q_", y) & is.na(out_metric)] <- "Best q"
  out_rank[grepl("^best_q_", y) & out_rank == 99L] <- 2L

  out_metric[grepl("^focus_n_sig_terms_", y)] <- "Focus sig terms"
  out_rank[grepl("^focus_n_sig_terms_", y)] <- 1L

  out_metric[grepl("^focus_best_neglog10q_", y)] <- "Focus best q"
  out_rank[grepl("^focus_best_neglog10q_", y)] <- 1L

  out_metric[grepl("^focus_best_q_", y) & is.na(out_metric)] <- "Focus best q"
  out_rank[grepl("^focus_best_q_", y) & out_rank == 99L] <- 2L

  data.table(metric = out_metric, metric_rank = out_rank)
}

extract_wide_long <- function(dt, value_patterns, value_name) {
  p_col <- pick_col(dt, c("^(program|program_name)$"), label = "program")
  s_col <- pick_col(dt, c("^source$"), label = "source")
  c_col <- pick_col(dt, c("(control.*type|match.*type|null.*model|matched.*control|control_set)"), label = "control type")

  dt2 <- copy(dt)
  setnames(dt2, c(p_col, s_col, c_col), c("program", "source", "control_type"))

  nms <- names(dt2)
  nn  <- norm_name(nms)
  keep <- rep(FALSE, length(nms))
  for (pat in value_patterns) keep <- keep | grepl(pat, nn, perl = TRUE)
  measure_cols <- nms[keep]
  if (length(measure_cols) == 0) {
    stop("No matching columns found for ", value_name, ". Available columns: ", paste(names(dt), collapse = ", "))
  }

  long <- melt(
    dt2,
    id.vars = c("program", "source", "control_type"),
    measure.vars = measure_cols,
    variable.name = "metric_raw",
    value.name = value_name,
    variable.factor = FALSE
  )

  mi <- metric_info(long$metric_raw)
  long[, metric := mi$metric]
  long[, metric_rank := mi$metric_rank]
  long <- long[!is.na(metric)]
  long[, program := standardize_program(program)]
  long[, source := standardize_source(source)]
  long[, control_type := standardize_control(control_type)]
  long
}

prep_empirical <- function(dt) {
  long <- extract_wide_long(
    dt,
    value_patterns = c("emp_p$", "empirical_p$"),
    value_name = "empirical_p"
  )
  long[, empirical_p := suppressWarnings(as.numeric(empirical_p))]
  setorder(long, program, source, control_type, metric, metric_rank)
  long <- long[, .SD[1], by = .(program, source, control_type, metric)]
  long[, metric_rank := NULL]
  long
}

prep_zscores <- function(dt) {
  long <- extract_wide_long(
    dt,
    value_patterns = c("zscore$", "z_score$", "_z$"),
    value_name = "zscore"
  )
  long[, zscore := suppressWarnings(as.numeric(zscore))]
  setorder(long, program, source, control_type, metric, metric_rank)
  long <- long[, .SD[1], by = .(program, source, control_type, metric)]
  long[, metric_rank := NULL]
  long
}

fmt_num_clean <- function(x, digits = 1) {
  out <- rep("", length(x))
  ok <- !is.na(x)
  vals <- round(x[ok], digits)
  vals[abs(vals) < (0.5 * 10^(-digits))] <- 0
  out[ok] <- sprintf(paste0("%.", digits, "f"), vals)
  out
}

fmt_clip <- function(x, cap = 20, digits = 1) {
  out <- rep("", length(x))
  out[!is.na(x) & x > cap] <- paste0(">", cap)
  out[!is.na(x) & x < -cap] <- paste0("<-", cap)
  keep <- !is.na(x) & x <= cap & x >= -cap
  vals <- round(x[keep], digits)
  vals[abs(vals) < (0.5 * 10^(-digits))] <- 0
  out[keep] <- sprintf(paste0("%.", digits, "f"), vals)
  out
}

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(6, 8, 6, 6)
  )

main <- function() {
  opt <- parse_args()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  tab_dir <- file.path(opt$base_dir, "Step13A_specificity_matched_random_controls", "tables")
  emp_path  <- file.path(tab_dir, "12_empirical_pvalues.tsv")
  z_path    <- file.path(tab_dir, "13_metric_zscores.tsv")

  emp <- prep_empirical(read_tsv_flex(emp_path))
  zdt <- prep_zscores(read_tsv_flex(z_path))

  main_programs <- c("SFARI_all", "midPrenatal_SFARI_top20")
  metric_levels <- c("Sig terms", "Best q", "Focus sig terms", "Focus best q")

  emp <- emp[
    program %in% main_programs &
    source %in% c("SynGO", "Reactome") &
    control_type %in% c("size-matched", "expression-matched") &
    metric %in% metric_levels
  ]
  zdt <- zdt[
    program %in% main_programs &
    source %in% c("SynGO", "Reactome") &
    control_type %in% c("size-matched", "expression-matched") &
    metric %in% metric_levels
  ]

  if (nrow(emp) == 0) {
    stop(
      "No usable empirical rows after filtering.\nPrograms: ", paste(unique(prep_empirical(read_tsv_flex(emp_path))$program), collapse = ", "),
      "\nSources: ", paste(unique(prep_empirical(read_tsv_flex(emp_path))$source), collapse = ", "),
      "\nControls: ", paste(unique(prep_empirical(read_tsv_flex(emp_path))$control_type), collapse = ", "),
      "\nMetrics: ", paste(unique(prep_empirical(read_tsv_flex(emp_path))$metric), collapse = ", ")
    )
  }

  emp[, program_short := factor(program_short(program), levels = c("SFARI_all", "midPrenatal_SFARI_top20"))]
  emp[, source := factor(source, levels = c("SynGO", "Reactome"))]
  emp[, control_type := factor(control_type, levels = c("size-matched", "expression-matched"))]
  emp[, metric := factor(metric, levels = metric_levels)]
  emp[, neglog10_emp := -log10(pmax(empirical_p, 1e-300))]
  emp[!is.finite(neglog10_emp), neglog10_emp := NA_real_]
  emp[, label := fmt_num_clean(neglog10_emp, 1)]

  pA <- ggplot(emp, aes(x = metric, y = program_short, fill = neglog10_emp)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = label), size = 3.2) +
    facet_grid(control_type ~ source) +
    scale_fill_gradientn(
      colors = c("#F7F7F7", "#D9E6F5", "#8CB3DD", "#3B73B9"),
      limits = c(0, 3),
      na.value = "grey95",
      name = expression(-log[10]("empirical p"))
    ) +
    labs(
      title = "Empirical significance against matched-random controls",
      x = NULL,
      y = NULL
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  zdt[, program_short := factor(program_short(program), levels = c("SFARI_all", "midPrenatal_SFARI_top20"))]
  zdt[, source := factor(source, levels = c("SynGO", "Reactome"))]
  zdt[, control_type := factor(control_type, levels = c("size-matched", "expression-matched"))]
  zdt[, metric := factor(metric, levels = metric_levels)]
  zdt[, zscore_cap := pmax(pmin(zscore, 20), -20)]
  zdt[, label := fmt_clip(zscore, cap = 20, digits = 1)]

  pB <- ggplot(zdt, aes(x = metric, y = program_short, fill = zscore_cap)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = label), size = 3.1) +
    facet_grid(control_type ~ source) +
    scale_fill_gradient2(
      low = "#2F6CB3",
      mid = "#F7F7F7",
      high = "#B35806",
      midpoint = 0,
      limits = c(-20, 20),
      oob = scales::squish,
      na.value = "grey95",
      name = "z-score\n(capped)"
    ) +
    labs(
      title = "Standardized deviation from null distributions",
      x = NULL,
      y = NULL
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  merged <- merge(
    emp[, .(program, source, control_type, metric, neglog10_emp)],
    zdt[, .(program, source, control_type, metric, zscore)],
    by = c("program", "source", "control_type", "metric"),
    all = TRUE
  )

  summary_dt <- copy(merged)
  summary_dt[, metric_priority := fifelse(metric == "Focus best q", 1L,
                                   fifelse(metric == "Best q", 2L,
                                   fifelse(metric == "Focus sig terms", 3L,
                                   fifelse(metric == "Sig terms", 4L, 99L))))]
  setorder(summary_dt, program, source, control_type, metric_priority)
  summary_dt <- summary_dt[, .SD[1], by = .(program, source, control_type)]
  summary_dt <- summary_dt[, .(
    neglog10_emp = mean(neglog10_emp, na.rm = TRUE),
    zscore = mean(zscore, na.rm = TRUE),
    metric_used = first(metric)
  ), by = .(program, source)]
  summary_dt[!is.finite(neglog10_emp), neglog10_emp := NA_real_]
  summary_dt[!is.finite(zscore), zscore := NA_real_]
  summary_dt[, program_short := factor(program_short(program), levels = c("SFARI_all", "midPrenatal_SFARI_top20"))]
  summary_dt[, source := factor(source, levels = c("SynGO", "Reactome"))]
  summary_dt[, label := paste0("emp ", fmt_num_clean(neglog10_emp, 1), "\n", "z ", fmt_clip(zscore, cap = 20, digits = 1))]

  pC <- ggplot(summary_dt, aes(x = source, y = program_short, fill = neglog10_emp)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = 3.4, lineheight = 0.95) +
    scale_fill_gradientn(
      colors = c("#F7F7F7", "#D9E6F5", "#8CB3DD", "#3B73B9"),
      limits = c(0, 3),
      na.value = "grey95",
      name = expression(-log[10]("empirical p"))
    ) +
    labs(
      title = "Integrated specificity summary",
      x = NULL,
      y = NULL
    ) +
    base_theme

  fig <- (pA | pB) / pC +
    plot_layout(heights = c(1.05, 0.62)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

  pdf_file <- file.path(opt$outdir, "Figure6_specificity_controls_v4_labelsync.pdf")
  png_file <- file.path(opt$outdir, "Figure6_specificity_controls_v4_labelsync.png")
  note_file <- file.path(opt$outdir, "Figure6_specificity_controls_v4_labelsync_note.txt")

  ggsave(pdf_file, fig, width = 16, height = 11.2)
  ggsave(png_file, fig, width = 16, height = 11.2, dpi = 300)

  writeLines(c(
    "Figure 6 v4 labelsync inputs:",
    emp_path,
    z_path,
    "",
    "Fixes applied:",
    "- -0.0 labels are converted to 0.0",
    "- Panel A color scale capped at 3",
    "- Panel C label style simplified to 'emp x.x / z x.x'",
    "- Panel C title changed to 'Integrated specificity summary'",
    "- Program display labels unified to manuscript-consistent names (SFARI_all, midPrenatal_SFARI_top20)"
  ), note_file)

  cat("[DONE] Wrote:\n", pdf_file, "\n", png_file, "\n", note_file, "\n")
}

main()
