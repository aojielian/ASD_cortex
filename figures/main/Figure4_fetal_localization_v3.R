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

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  scale_y_discrete(labels = function(x) gsub(paste0(sep, ".+$"), "", x), ...)
}

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(6, 8, 6, 6)
  )

program_levels <- c("SFARI_all", "midPrenatal_SFARI_top20")
program_labels <- c(
  "SFARI_all" = "SFARI_all",
  "midPrenatal_SFARI_top20" = "midPrenatal_SFARI_top20"
)
program_colors <- c(
  "SFARI_all" = "#D95F02",
  "midPrenatal_SFARI_top20" = "#1B9E77"
)
stage_palette <- c(
  "12" = "#F8766D",
  "13" = "#A3A500",
  "14" = "#00BF7D",
  "19" = "#00B0F6",
  "22" = "#C77CFF"
)
broad_levels <- c("Radial Glial", "IPC", "Neuronal", "Mesenchymal", "Other")

main <- function() {
  opt <- parse_args()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  step11_posthoc <- file.path(opt$base_dir, "Step11A_fetal_singlecell_localization_final_posthoc", "tables")
  step11_final <- file.path(opt$base_dir, "Step11A_fetal_singlecell_localization_final", "tables")

  stage_bc <- read_required(file.path(step11_posthoc, "24_stage_broadclass_summary_meanlog_nge5.tsv"))
  top_states <- read_required(file.path(step11_posthoc, "23_meanlogexpr_states_annotated_nge5.tsv"))
  top_states_raw <- read_required(file.path(step11_final, "10_top_developmental_states_by_score.tsv"))

  stage_bc <- stage_bc[program %in% program_levels]
  top_states <- top_states[program %in% program_levels & primary_method == TRUE & robust_state == TRUE]
  top_states_raw <- top_states_raw[program %in% program_levels & method == "MeanLogExpr"]

  # Restrict the displayed stages to those with substantive main-figure support.
  keep_stages <- c("12", "13", "19", "22")

  stage_bc[, stage := as.character(stage)]
  top_states[, stage := as.character(stage)]
  top_states_raw[, stage := as.character(stage)]

  stage_bc <- stage_bc[stage %in% keep_stages]
  top_states <- top_states[stage %in% keep_stages]
  top_states_raw <- top_states_raw[stage %in% keep_stages]

  stage_levels <- keep_stages
  stage_bc[, stage := factor(stage, levels = stage_levels)]
  top_states[, stage := factor(stage, levels = stage_levels)]
  top_states_raw[, stage := factor(stage, levels = stage_levels)]

  # Panel A
  stage_bc[, broad_cell_type := factor(broad_cell_type, levels = rev(broad_levels))]
  stage_bc[, program_f := factor(program, levels = program_levels, labels = unname(program_labels[program_levels]))]

  grid_dt <- CJ(
    program = program_levels,
    stage = stage_levels,
    broad_cell_type = broad_levels,
    unique = TRUE
  )
  grid_dt[, stage := factor(stage, levels = stage_levels)]
  grid_dt[, broad_cell_type := factor(broad_cell_type, levels = rev(broad_levels))]
  grid_dt <- merge(
    grid_dt,
    stage_bc[, .(program, stage, broad_cell_type, weighted_mean_score)],
    by = c("program", "stage", "broad_cell_type"),
    all.x = TRUE
  )
  grid_dt[, program_f := factor(program, levels = program_levels, labels = unname(program_labels[program_levels]))]

  pA <- ggplot(grid_dt, aes(x = stage, y = broad_cell_type, fill = weighted_mean_score)) +
    geom_tile(color = "white", linewidth = 0.7) +
    facet_wrap(~program_f, nrow = 1) +
    scale_fill_gradientn(
      colors = c("#F7F7F7", "#DCE6F2", "#A8BCDF", "#6E90C7", "#2F64AD"),
      name = "Mean score",
      na.value = "grey97"
    ) +
    labs(
      title = "Coarse fetal developmental-state localization",
      x = "Stage",
      y = "Broad fetal class"
    ) +
    base_theme +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.major = element_blank()
    )

  # Panel B
  top_states_b <- copy(top_states)[order(program, -mean_score)]
  top_states_b <- top_states_b[, head(.SD, 5), by = program]
  top_states_b[, state_label := paste0("C", cluster_id, " | ", refined_label)]
  top_states_b[, program_f := factor(program, levels = program_levels, labels = unname(program_labels[program_levels]))]
  top_states_b[, state_label_f := reorder_within(state_label, mean_score, program_f)]

  stage_breaks_B <- intersect(names(stage_palette), sort(unique(as.character(top_states_b$stage))))

  pB <- ggplot(top_states_b, aes(x = mean_score, y = state_label_f, fill = stage)) +
    geom_col(width = 0.72) +
    facet_wrap(~program_f, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = stage_palette, breaks = stage_breaks_B, drop = TRUE, name = "Stage") +
    scale_y_reordered() +
    labs(
      title = "Top annotated fetal states",
      x = "Mean score",
      y = NULL
    ) +
    base_theme +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", size = 11)
    )

  # Panel C
  stage_summary <- stage_bc[, .(
    stage_weighted_mean = weighted.mean(weighted_mean_score, w = pmax(total_cells_in_states, 1), na.rm = TRUE),
    total_cells = sum(total_cells_in_states, na.rm = TRUE),
    n_state_groups = sum(n_cluster_stage_states, na.rm = TRUE)
  ), by = .(program, stage)]
  stage_summary[, program_f := factor(program, levels = program_levels, labels = unname(program_labels[program_levels]))]
  stage_summary[, stage := factor(as.character(stage), levels = stage_levels)]

  pC <- ggplot(stage_summary, aes(x = stage, y = stage_weighted_mean, color = program_f, size = n_state_groups)) +
    geom_point(position = position_dodge(width = 0.35), alpha = 0.95) +
    scale_color_manual(values = program_colors, name = "Program") +
    scale_size_continuous(name = "Robust\nstate groups", range = c(3.5, 8), breaks = c(1, 2, 3, 4)) +
    labs(
      title = "Stage concentration of fetal localization support",
      x = "Stage",
      y = "Stage-level mean score"
    ) +
    base_theme +
    theme(legend.position = "right")

  fig <- (pA / pB / pC) +
    plot_layout(heights = c(1.05, 1.05, 0.85)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

  pdf_file <- file.path(opt$outdir, "Figure4_fetal_localization_v3.pdf")
  png_file <- file.path(opt$outdir, "Figure4_fetal_localization_v3.png")
  ggsave(pdf_file, fig, width = 14, height = 13)
  ggsave(png_file, fig, width = 14, height = 13, dpi = 300)

  note_file <- file.path(opt$outdir, "Figure4_fetal_localization_v3_note.txt")
  writeLines(c(
    "Figure 4 v3 inputs:",
    file.path(step11_posthoc, "24_stage_broadclass_summary_meanlog_nge5.tsv"),
    file.path(step11_posthoc, "23_meanlogexpr_states_annotated_nge5.tsv"),
    file.path(step11_final, "10_top_developmental_states_by_score.tsv"),
    "",
    "Displayed stages restricted to: 12, 13, 19, 22"
  ), note_file)

  cat("[DONE] Wrote:\n", pdf_file, "\n", png_file, "\n", note_file, "\n")
}

main()
