#!/usr/bin/env Rscript

# Figure 2: adult bulk support and adult neuronal-context localization.
# Harmonization/collapse audits are handled in the legend and supplement, not as extra main panels.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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
    out[[sub("^--", "", key)]] <- val
    i <- i + 2
  }
  out
}

read_required <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  fread(path)
}

programs <- c("SFARI_all", "midPrenatal_SFARI_top20")
program_labels <- c("SFARI_all" = "SFARI_all", "midPrenatal_SFARI_top20" = "midPrenatal_SFARI_top20")
program_colors <- c("SFARI_all" = "#D95F02", "midPrenatal_SFARI_top20" = "#1B9E77")

cohort_colors <- c(
  "GSE102741" = "#4C78A8",
  "GSE64018" = "#F28E2B",
  "Gandal2022" = "#E15759",
  "Meta" = "black"
)
cohort_shapes <- c(
  "GSE102741" = 16,
  "GSE64018" = 16,
  "Gandal2022" = 16,
  "Meta" = 18
)

broad_levels <- c("INN", "EXN", "AST", "OPC", "ODC", "MG", "END")

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(6, 8, 6, 6)
  )

make_loc_dumbbell <- function(dt, dataset_title, xlim_shared, show_y_title = TRUE) {
  x <- copy(dt)[method == "UCell" & program_name %in% programs]
  x <- x[broad_class %in% broad_levels]
  x[, broad_class := factor(broad_class, levels = rev(broad_levels))]
  x[, program_name := factor(program_name, levels = programs)]

  seg <- dcast(
    x[, .(broad_class, program_name, mean_score)],
    broad_class ~ program_name,
    value.var = "mean_score"
  )

  ggplot() +
    geom_segment(
      data = seg,
      aes(
        x = SFARI_all,
        xend = midPrenatal_SFARI_top20,
        y = broad_class,
        yend = broad_class
      ),
      color = "grey80",
      linewidth = 0.8
    ) +
    geom_point(
      data = x,
      aes(x = mean_score, y = broad_class, color = program_name),
      size = 4.2
    ) +
    scale_color_manual(
      values = program_colors,
      breaks = programs,
      labels = unname(program_labels[programs]),
      guide = "none"
    ) +
    coord_cartesian(xlim = c(0, xlim_shared)) +
    labs(
      title = dataset_title,
      x = "Mean UCell score",
      y = if (show_y_title) "Broad class" else NULL
    ) +
    base_theme +
    theme(axis.title.y = if (show_y_title) element_text() else element_blank())
}

main <- function() {
  opt <- parse_args()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  tbl_dir <- file.path(opt$base_dir, "tables")
  dt130 <- read_required(file.path(tbl_dir, "130_program_summary_for_manuscript_updated.tsv"))
  dt31c <- read_required(file.path(tbl_dir, "31c_Gandal_program_leaveOneRegionOut.tsv"))
  dt95p <- read_required(file.path(tbl_dir, "95_PsychENCODE_UCell_AUCell_broadClass_localization.tsv"))
  dt95v <- read_required(file.path(tbl_dir, "95_Velmeshev_UCell_AUCell_broadClass_localization.tsv"))

  dt130 <- dt130[program_name %in% programs]
  dt31c <- dt31c[program_name %in% programs]
  dt95p <- dt95p[program_name %in% programs]
  dt95v <- dt95v[program_name %in% programs]

  a_dt <- melt(
    dt130[, .(
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
  a_dt[, program_f := factor(program_name, levels = rev(programs))]
  a_dt[, series := factor(series, levels = c("GSE102741", "GSE64018", "Gandal2022", "Meta"))]

  pA <- ggplot(a_dt, aes(x = effect, y = program_f, color = series, shape = series)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
    geom_point(size = 3.5) +
    scale_color_manual(values = cohort_colors, name = "Series") +
    scale_shape_manual(values = cohort_shapes, name = "Series") +
    scale_y_discrete(labels = c("midPrenatal_SFARI_top20", "SFARI_all")) +
    labs(
      title = "Adult bulk cohort effects",
      x = "ASD - Control program effect",
      y = NULL
    ) +
    annotate(
      "text",
      x = min(a_dt$effect, na.rm = TRUE),
      y = 0.65,
      label = "negative = ASD lower",
      hjust = 0,
      size = 3.2,
      color = "grey35"
    ) +
    coord_cartesian(
      xlim = c(min(a_dt$effect, na.rm = TRUE) - 0.025, max(a_dt$effect, na.rm = TRUE) + 0.025)
    ) +
    base_theme +
    theme(legend.position = "right")

  b_dt <- copy(dt31c)
  b_dt[, program_f := factor(program_name, levels = programs, labels = c("SFARI_all", "midPrenatal_SFARI_top20"))]
  region_levels <- b_dt[program_name == "SFARI_all"][order(beta_asd), left_out_region]
  b_dt[, left_out_region := factor(left_out_region, levels = region_levels)]
  b_dt[, neglog10_fdr := -log10(pmax(fdr_p_asd, 1e-300))]
  b_dt[, program_name := factor(program_name, levels = programs)]

  pB <- ggplot(b_dt, aes(x = beta_asd, y = left_out_region, color = program_name)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
    geom_segment(
      aes(x = 0, xend = beta_asd, y = left_out_region, yend = left_out_region),
      color = "grey82",
      linewidth = 0.5
    ) +
    geom_point(aes(size = neglog10_fdr), alpha = 0.95) +
    scale_color_manual(
      values = program_colors,
      breaks = programs,
      labels = unname(program_labels[programs]),
      name = "Program"
    ) +
    scale_size_continuous(name = expression(-log[10]*"(FDR)"), range = c(2.5, 6)) +
    facet_wrap(~program_f, nrow = 1, scales = "free_x") +
    labs(
      title = "Gandal2022 leave-one-region-out sensitivity",
      x = expression(beta[ASD]),
      y = "Left-out region"
    ) +
    base_theme +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", size = 11)
    )

  loc_shared_max <- max(
    dt95p[method == "UCell" & program_name %in% programs, mean_score],
    dt95v[method == "UCell" & program_name %in% programs, mean_score],
    na.rm = TRUE
  ) * 1.08

  pC <- make_loc_dumbbell(dt95p, "PsychENCODE localization", loc_shared_max, show_y_title = TRUE)
  pD <- make_loc_dumbbell(dt95v, "Velmeshev localization", loc_shared_max, show_y_title = FALSE)

  fig <- (pA + pB) / (pC + pD) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

  pdf_file <- file.path(opt$outdir, "Figure2_adult_bulk_singlecell_v5_labelsync.pdf")
  png_file <- file.path(opt$outdir, "Figure2_adult_bulk_singlecell_v5_labelsync.png")
  ggsave(pdf_file, fig, width = 15, height = 11)
  ggsave(png_file, fig, width = 15, height = 11, dpi = 300)
  cat("[DONE] Wrote:\n", pdf_file, "\n", png_file, "\n")
}

main()
