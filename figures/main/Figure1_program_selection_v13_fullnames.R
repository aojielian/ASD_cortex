#!/usr/bin/env Rscript

# Figure 1: program prioritization using the frozen manuscript programs.
# Note: SFARI_all in the manuscript denotes the frozen 386-gene anchor,
# not the broader 435-gene primary SFARI source pool used only for audit steps.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    base_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    outdir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/Figures_main"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args)) stop("Missing value for ", key)
    val <- args[[i + 1]]
    nm <- sub("^--", "", key)
    out[[nm]] <- val
    i <- i + 2
  }
  out
}

short_pair <- function(x) {
  x <- gsub("GSE102741", "102741", x, fixed = TRUE)
  x <- gsub("GSE64018", "64018", x, fixed = TRUE)
  x <- gsub("Gandal2022", "Gandal", x, fixed = TRUE)
  x <- gsub("__", "+", x, fixed = TRUE)
  x
}

pretty_program <- function(x) {
  factor(
    x,
    levels = c("midPrenatal_SFARI_top05", "midPrenatal_SFARI_top10", "midPrenatal_SFARI_top20", "SFARI_all")
  )
}

program_display <- c(
  "midPrenatal_SFARI_top05" = "midPrenatal_SFARI_top05",
  "midPrenatal_SFARI_top10" = "midPrenatal_SFARI_top10",
  "midPrenatal_SFARI_top20" = "midPrenatal_SFARI_top20",
  "SFARI_all" = "SFARI_all"
)

program_colors <- c(
  "midPrenatal_SFARI_top05" = "#999999",
  "midPrenatal_SFARI_top10" = "#7570B3",
  "midPrenatal_SFARI_top20" = "#1B9E77",
  "SFARI_all" = "#D95F02"
)

read_required <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  fread(path)
}

legend_theme <- theme(
  legend.title = element_text(size = 11, face = "bold"),
  legend.text = element_text(size = 10),
  legend.key.height = unit(0.45, "cm"),
  legend.key.width = unit(0.45, "cm"),
  legend.spacing.y = unit(0.1, "cm"),
  legend.box.spacing = unit(0.1, "cm"),
  legend.margin = margin(0, 0, 0, 0)
)

main <- function() {
  opt <- parse_args()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  tbl_dir <- file.path(opt$base_dir, "tables")
  dt130 <- read_required(file.path(tbl_dir, "130_program_summary_for_manuscript_updated.tsv"))
  dt134 <- read_required(file.path(tbl_dir, "134_bulk_leaveOneOut_meta.tsv"))

  wanted <- c("SFARI_all","midPrenatal_SFARI_top20","midPrenatal_SFARI_top10","midPrenatal_SFARI_top05")
  dt130 <- dt130[program_name %in% wanted]
  dt134 <- dt134[program_name %in% wanted]

  req130 <- c("program_name","developmental_concentration_rank","n_genes","fdr_weighted",
              "delta_asd_minus_control_GSE102741","delta_asd_minus_control_GSE64018",
              "delta_asd_minus_control_Gandal2022","mean_delta")
  req134 <- c("program_name","cohort_pair","n_valid_z","n_direction_matched","n_direction_mismatched")
  if (length(setdiff(req130, names(dt130))) > 0) stop("130 missing required columns")
  if (length(setdiff(req134, names(dt134))) > 0) stop("134 missing required columns")

  max_rank <- max(dt130$developmental_concentration_rank, na.rm = TRUE)

  base_theme <- theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.margin = margin(6, 8, 6, 6)
    )

  # A
  a_dt <- copy(dt130)
  a_dt[, prenatal_focus_score := -log10(developmental_concentration_rank / max_rank)]
  a_dt[, program_f := pretty_program(program_name)]

  pA <- ggplot(a_dt, aes(x = prenatal_focus_score, y = program_f, color = program_name)) +
    geom_segment(aes(x = 0, xend = prenatal_focus_score, y = program_f, yend = program_f),
                 linewidth = 0.8, color = "grey78") +
    geom_point(size = 4) +
    scale_color_manual(values = program_colors, guide = "none") +
    labs(title = "Developmental support",
         x = "Prenatal focus score",
         y = NULL) +
    coord_cartesian(xlim = c(-0.05, max(a_dt$prenatal_focus_score) * 1.18)) +
    base_theme

  # B
  b_dt <- copy(dt130)
  b_dt[, prenatal_focus_score := -log10(developmental_concentration_rank / max_rank)]
  b_dt[, adult_support := -log10(pmax(fdr_weighted, 1e-300))]

  pB <- ggplot(b_dt, aes(x = prenatal_focus_score, y = adult_support, size = n_genes, color = program_name)) +
    geom_point(alpha = 0.95, stroke = 0.35) +
    scale_color_manual(
      values = program_colors,
      breaks = names(program_display),
      labels = unname(program_display[names(program_display)]),
      name = "Program"
    ) +
    scale_size_continuous(name = "Gene set size", range = c(5, 13)) +
    labs(title = "Program prioritization",
         x = "Prenatal focus score",
         y = expression("Adult bulk support (" * -log[10] * " FDR)")) +
    coord_cartesian(
      xlim = c(-0.15, max(b_dt$prenatal_focus_score) * 1.12),
      ylim = c(0, max(b_dt$adult_support) * 1.05)
    ) +
    base_theme +
    theme(legend.position = "right") +
    legend_theme

  # C
  c_dt <- melt(
    copy(dt130)[, .(
      program_name,
      GSE102741 = delta_asd_minus_control_GSE102741,
      GSE64018 = delta_asd_minus_control_GSE64018,
      Gandal2022 = delta_asd_minus_control_Gandal2022,
      Meta = mean_delta
    )],
    id.vars = "program_name",
    variable.name = "series",
    value.name = "effect"
  )
  c_dt[, program_f := pretty_program(program_name)]
  series_cols <- c("GSE102741" = "#4C78A8", "GSE64018" = "#F28E2B", "Gandal2022" = "#E15759", "Meta" = "black")
  series_shapes <- c("GSE102741" = 16, "GSE64018" = 16, "Gandal2022" = 16, "Meta" = 18)

  pC <- ggplot(c_dt, aes(x = effect, y = program_f, color = series, shape = series)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
    geom_point(size = 3.3) +
    scale_color_manual(values = series_cols, name = "Series") +
    scale_shape_manual(values = series_shapes, name = "Series") +
    labs(title = "Adult bulk support of frozen programs",
         x = "ASD - Control program effect",
         y = NULL) +
    annotate(
      "text",
      x = min(c_dt$effect, na.rm = TRUE),
      y = 0.6,
      label = "negative = ASD lower",
      hjust = 0,
      size = 3.2,
      color = "grey35"
    ) +
    coord_cartesian(xlim = c(min(c_dt$effect, na.rm = TRUE) - 0.015, max(c_dt$effect, na.rm = TRUE) + 0.015)) +
    base_theme +
    theme(
      legend.position = "right",
      axis.title.x = element_text(margin = margin(t = 6))
    ) +
    legend_theme

  # D
  d_dt <- copy(dt134)
  d_dt[, program_f := pretty_program(program_name)]
  d_dt[, cohort_pair_short := factor(short_pair(cohort_pair),
                                     levels = c("102741+64018", "Gandal+102741", "Gandal+64018"))]
  d_dt[, matched_fraction := n_direction_matched / pmax(n_valid_z, 1)]
  d_dt[, matched_label := paste0(round(100 * matched_fraction), "%")]

  pD <- ggplot(d_dt, aes(x = cohort_pair_short, y = program_f, fill = matched_fraction)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = matched_label), size = 4) +
    scale_fill_gradient(
      low = "#F2F5FA",
      high = "#4A90D6",
      breaks = c(0.5, 1.0),
      labels = c("50%", "100%"),
      name = "Direction matched"
    ) +
    labs(title = "Leave-one-out consistency",
         x = "Cohort pair",
         y = NULL) +
    base_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    legend_theme

  fig <- (pA + pB) / (pC + pD) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

  pdf_file <- file.path(opt$outdir, "Figure1_program_selection_v13_fullnames.pdf")
  png_file <- file.path(opt$outdir, "Figure1_program_selection_v13_fullnames.png")
  ggsave(pdf_file, fig, width = 14, height = 10)
  ggsave(png_file, fig, width = 14, height = 10, dpi = 300)
  cat("[DONE] Wrote:\n", pdf_file, "\n", png_file, "\n")
}

main()
