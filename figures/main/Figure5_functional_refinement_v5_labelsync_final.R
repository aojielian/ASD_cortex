#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(scales)
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

read_tsv_flex <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  fread(path, sep = "\t", header = TRUE, fill = TRUE, quote = "", showProgress = FALSE)
}

coerce_numeric_if_present <- function(dt, cols) {
  for (nm in intersect(cols, names(dt))) {
    dt[, (nm) := suppressWarnings(as.numeric(get(nm)))]
  }
  dt
}

clean_syngo_term <- function(x) {
  x <- gsub(" \\(GO:[0-9]+\\)$", "", x)
  map <- c(
    "presynaptic active zone cytoplasmic component" = "Presynaptic active zone\ncytoplasmic component",
    "integral component of presynaptic active zone membrane" = "Presynaptic active zone\nmembrane",
    "synaptic vesicle docking" = "Synaptic vesicle docking",
    "ligand-gated ion channel activity involved in regulation of presynaptic membrane potential" = "Presynaptic membrane\npotential regulation",
    "postsynaptic density, intracellular component" = "Postsynaptic density\nintracellular component",
    "postsynaptic density" = "Postsynaptic density"
  )
  y <- unname(map[x])
  y[is.na(y)] <- x[is.na(y)]
  y
}

clean_reactome_term <- function(x) {
  map <- c(
    "NEURONAL SYSTEM" = "Neuronal system",
    "PROTEIN PROTEIN INTERACTIONS AT SYNAPSES" = "Synaptic protein interactions",
    "NEUREXINS AND NEUROLIGINS" = "Neurexins and neuroligins",
    "LONG TERM POTENTIATION" = "Long-term potentiation",
    "UNBLOCKING OF NMDA RECEPTORS GLUTAMATE BINDING AND ACTIVATION" = "NMDA receptor activation",
    "CHROMATIN ORGANIZATION" = "Chromatin organization",
    "EPIGENETIC REGULATION OF GENE EXPRESSION" = "Epigenetic regulation",
    "ATP DEPENDENT CHROMATIN REMODELERS" = "ATP-dependent remodelers",
    "FORMATION OF NEURONAL PROGENITOR AND NEURONAL BAF NPBAF AND NBAF" = "Neuronal BAF/npBAF/nBAF",
    "TRANSCRIPTIONAL REGULATION BY RUNX1" = "RUNX1 transcription",
    "RUNX1 INTERACTS WITH CO FACTORS WHOSE PRECISE EFFECT ON RUNX1 TARGETS IS NOT KNOWN" = "RUNX1 cofactors"
  )
  y <- unname(map[x])
  y[is.na(y)] <- x[is.na(y)]
  y
}

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(6, 8, 6, 6)
  )

program_levels <- c("SFARI_all", "midPrenatal_SFARI_top20")
program_short <- c("SFARI_all" = "SFARI_all", "midPrenatal_SFARI_top20" = "midPrenatal_SFARI_top20")
program_short_levels <- c("SFARI_all", "midPrenatal_SFARI_top20")
program_colors_short <- c("SFARI_all" = "#D95F02", "midPrenatal_SFARI_top20" = "#1B9E77")
domain_colors <- c("BP" = "#F8766D", "CC" = "#00BFC4")

main <- function() {
  opt <- parse_args()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  syngo_sig_path <- file.path(opt$base_dir, "Step12A_SynGO_refinement", "tables", "11_significant_enrichment_results.tsv")
  reactome_sig_path <- file.path(opt$base_dir, "Step12B_upstream_regulator_refinement", "tables", "13_reactome_significant_results.tsv")
  focus_path <- file.path(opt$base_dir, "Step12B_upstream_regulator_refinement", "tables", "14_focus_regulatory_terms.tsv")

  syngo_sig <- read_tsv_flex(syngo_sig_path)
  reactome_sig <- read_tsv_flex(reactome_sig_path)
  focus_dt <- read_tsv_flex(focus_path)

  syngo_sig <- coerce_numeric_if_present(syngo_sig, c("overlap_n", "q_value", "neglog10_q"))
  reactome_sig <- coerce_numeric_if_present(reactome_sig, c("overlap_n", "q_value", "neglog10_q"))
  focus_dt <- coerce_numeric_if_present(focus_dt, c("overlap_n", "q_value", "neglog10_q"))

  syngo_a <- syngo_sig[program == "SFARI_all"][order(neglog10_q)]
  syngo_a[, term_short := clean_syngo_term(term_name)]
  syngo_order <- c(
    "Postsynaptic density",
    "Postsynaptic density\nintracellular component",
    "Presynaptic membrane\npotential regulation",
    "Synaptic vesicle docking",
    "Presynaptic active zone\nmembrane",
    "Presynaptic active zone\ncytoplasmic component"
  )
  syngo_a[, term_short := factor(term_short, levels = syngo_order)]

  pA <- ggplot(syngo_a, aes(x = neglog10_q, y = term_short, size = overlap_n, color = domain)) +
    geom_point(alpha = 0.95) +
    scale_color_manual(values = domain_colors, name = "Domain") +
    scale_size_continuous(name = "Overlap", range = c(2.5, 8), breaks = c(5, 10, 15, 20)) +
    labs(
      title = "SynGO-defined synaptic enrichment of the broad anchor",
      x = expression(-log[10](q)),
      y = NULL
    ) +
    base_theme +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.6, "cm")
    ) +
    guides(
      size = guide_legend(order = 1, override.aes = list(alpha = 1)),
      color = guide_legend(order = 2, override.aes = list(size = 4))
    )

  reactome_terms <- c(
    "NEURONAL SYSTEM",
    "PROTEIN PROTEIN INTERACTIONS AT SYNAPSES",
    "NEUREXINS AND NEUROLIGINS",
    "LONG TERM POTENTIATION",
    "UNBLOCKING OF NMDA RECEPTORS GLUTAMATE BINDING AND ACTIVATION",
    "CHROMATIN ORGANIZATION",
    "EPIGENETIC REGULATION OF GENE EXPRESSION",
    "ATP DEPENDENT CHROMATIN REMODELERS",
    "FORMATION OF NEURONAL PROGENITOR AND NEURONAL BAF NPBAF AND NBAF",
    "TRANSCRIPTIONAL REGULATION BY RUNX1",
    "RUNX1 INTERACTS WITH CO FACTORS WHOSE PRECISE EFFECT ON RUNX1 TARGETS IS NOT KNOWN"
  )
  reactome_class <- data.table(
    term_name = reactome_terms,
    functional_axis = c(rep("Synaptic / neuronal", 5), rep("Chromatin / regulatory", 6))
  )

  reactome_b <- reactome_sig[
    program %in% program_levels & term_name %in% reactome_terms,
    .(program, term_name, overlap_n, neglog10_q)
  ]
  grid_b <- CJ(program = program_levels, term_name = reactome_terms, unique = TRUE)
  reactome_b <- merge(grid_b, reactome_b, by = c("program", "term_name"), all.x = TRUE)
  reactome_b <- merge(reactome_b, reactome_class, by = "term_name", all.x = TRUE, sort = FALSE)
  reactome_b[, program_f := factor(unname(program_short[program]), levels = program_short_levels)]
  reactome_b[, functional_axis := factor(functional_axis, levels = c("Synaptic / neuronal", "Chromatin / regulatory"))]
  reactome_b[, term_short := clean_reactome_term(term_name)]
  term_order <- c(
    "Neuronal system",
    "Synaptic protein interactions",
    "Neurexins and neuroligins",
    "Long-term potentiation",
    "NMDA receptor activation",
    "Chromatin organization",
    "Epigenetic regulation",
    "ATP-dependent remodelers",
    "Neuronal BAF/npBAF/nBAF",
    "RUNX1 transcription",
    "RUNX1 cofactors"
  )
  reactome_b[, term_short := factor(term_short, levels = rev(term_order))]
  reactome_b[, overlap_label := fifelse(is.na(overlap_n), "", as.character(as.integer(overlap_n)))]

  pB <- ggplot(reactome_b, aes(x = program_f, y = term_short, fill = neglog10_q)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = overlap_label), size = 3.2) +
    facet_grid(functional_axis ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradientn(
      colors = c("#F7F7F7", "#DCE6F2", "#A8BCDF", "#6E90C7", "#2F64AD"),
      values = rescale(c(0, 4, 7, 10, 14)),
      limits = c(0, 14),
      na.value = "grey95",
      name = expression(-log[10](q)),
      breaks = c(4, 8, 12)
    ) +
    labs(
      title = "Reactome functional contrast",
      x = NULL,
      y = NULL
    ) +
    base_theme +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 10)
    )

  reg_terms <- c(
    "CHROMATIN ORGANIZATION",
    "EPIGENETIC REGULATION OF GENE EXPRESSION",
    "ATP DEPENDENT CHROMATIN REMODELERS",
    "TRANSCRIPTIONAL REGULATION BY RUNX1"
  )

  focus_c <- focus_dt[
    program %in% program_levels & term_name %in% reg_terms,
    .(program, term_name, overlap_n, neglog10_q)
  ]
  if (nrow(focus_c) == 0) {
    focus_c <- reactome_sig[
      program %in% program_levels & term_name %in% reg_terms,
      .(program, term_name, overlap_n, neglog10_q)
    ]
  }
  if (nrow(focus_c) == 0) stop("No focused regulatory terms were recovered.")
  focus_c[, program_f := factor(unname(program_short[program]), levels = program_short_levels)]
  focus_c[, term_short := clean_reactome_term(term_name)]
  focus_order <- c(
    "RUNX1 transcription",
    "ATP-dependent remodelers",
    "Epigenetic regulation",
    "Chromatin organization"
  )
  focus_c[, term_short := factor(term_short, levels = focus_order)]

  pC <- ggplot(focus_c, aes(x = neglog10_q, y = term_short, color = program_f, size = overlap_n)) +
    geom_point(alpha = 0.95, position = position_dodge(width = 0.45)) +
    scale_color_manual(values = program_colors_short, breaks = program_short_levels, drop = FALSE, name = "Program") +
    scale_size_continuous(name = "Overlap", range = c(3.2, 8.5), breaks = c(10, 20, 30)) +
    labs(
      title = "Focused chromatin-regulatory comparison",
      x = expression(-log[10](q)),
      y = NULL
    ) +
    base_theme +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 4)),
      size = guide_legend(order = 2)
    )

  top_row <- pA + pB + plot_layout(widths = c(1.18, 1.0))
  fig <- top_row / pC +
    plot_layout(heights = c(1.02, 0.95)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

  pdf_file <- file.path(opt$outdir, "Figure5_functional_refinement_v5_labelsync_final.pdf")
  png_file <- file.path(opt$outdir, "Figure5_functional_refinement_v5_labelsync_final.png")
  ggsave(pdf_file, fig, width = 15.5, height = 11.5)
  ggsave(png_file, fig, width = 15.5, height = 11.5, dpi = 300)

  note_file <- file.path(opt$outdir, "Figure5_functional_refinement_v5_labelsync_final_note.txt")
  writeLines(c(
    "Figure 5 v5 labelsync final inputs:",
    syngo_sig_path,
    reactome_sig_path,
    focus_path,
    "",
    "Fixes applied:",
    "- Program display labels unified to manuscript-consistent names (SFARI_all, midPrenatal_SFARI_top20).",
    "- Panel C legend now shows both programs.",
    "- Panel C overlap legend reduced to 10 / 20 / 30."
  ), note_file)

  cat("[DONE] Wrote:\n", pdf_file, "\n", png_file, "\n", note_file, "\n")
}

main()
