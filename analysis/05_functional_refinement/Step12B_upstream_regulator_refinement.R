#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(ggplot2)
  library(openxlsx)
})

pick_first_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    base_dir = NULL,
    program_file = NULL,
    requested_programs = NULL,
    q_cutoff = 0.10,
    min_term_size = 5,
    max_term_size = 500,
    species = "Homo sapiens",
    plot_top_n = 12,
    focus_regex = "chromatin|remodel|nucleosome|histone|epigen|transcription|rna polymerase|mediator|acetyl|methyl|deacetyl|splic|mRNA|gene expression|chromosome organization|SWI|BAF|ncor|coregulator|coactivator|corepressor"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key, call. = FALSE)
    if (grepl("=", key, fixed = TRUE)) stop("Use space-separated args, not --key=value. Offending arg: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for ", key, call. = FALSE)
    val <- args[[i + 1]]
    if (startsWith(val, "--")) stop("Missing value for ", key, call. = FALSE)
    nm <- sub("^--", "", key)
    if (!nm %in% names(out)) stop("Unknown argument: ", key, call. = FALSE)
    out[[nm]] <- val
    i <- i + 2
  }
  if (is.null(out$base_dir)) stop("--base_dir is required", call. = FALSE)
  if (is.null(out$program_file)) out$program_file <- file.path(out$base_dir, "tables", "20_midPrenatal_core_programs.tsv")
  if (is.null(out$requested_programs)) out$requested_programs <- "SFARI_all,midPrenatal_SFARI_top20,midPrenatal_SFARI_top10,midPrenatal_SFARI_top05"
  out$q_cutoff <- as.numeric(out$q_cutoff)
  out$min_term_size <- as.integer(out$min_term_size)
  out$max_term_size <- as.integer(out$max_term_size)
  out$plot_top_n <- as.integer(out$plot_top_n)
  out$requested_programs <- trimws(unlist(strsplit(out$requested_programs, ",", fixed = TRUE)))
  out
}

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
}

normalize_gene <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "", x)
  toupper(x)
}

read_programs <- function(program_file, requested_programs) {
  dt <- fread(program_file)
  nms <- names(dt)
  program_col <- pick_first_col(nms, c("program", "program_std", "program_name", "set", "geneset", "gene_set", "geneset_name"))
  gene_col <- pick_first_col(nms, c("gene", "gene_std", "gene_symbol", "symbol", "hgnc_symbol"))
  if (is.na(program_col) || is.na(gene_col)) {
    stop("Could not identify gene/program columns in program_file: ", program_file,
         "; columns=", paste(names(dt), collapse = ", "), call. = FALSE)
  }
  dt <- dt[, .(program = as.character(get(program_col)), gene = as.character(get(gene_col)))]
  dt[, program := trimws(program)]
  dt[, gene_norm := normalize_gene(gene)]
  dt <- dt[program %in% requested_programs & !is.na(gene_norm) & gene_norm != ""]
  dt <- unique(dt[, .(program, gene_norm)])
  if (nrow(dt) == 0) stop("No program genes remained after filtering requested_programs", call. = FALSE)
  dt
}

get_msig_table <- function(species = "Homo sapiens", collection = NULL) {
  # use new interface if available; fall back to deprecated args for older msigdbr
  x <- tryCatch(msigdbr::msigdbr(species = species, collection = collection), error = function(e) NULL)
  if (is.null(x)) {
    x <- tryCatch(msigdbr::msigdbr(species = species, category = collection), error = function(e) NULL)
  }
  if (is.null(x)) stop("Unable to retrieve msigdbr collection ", collection, call. = FALSE)
  as.data.table(x)
}

prep_msig_source <- function(dt, source_name) {
  nms <- names(dt)
  gene_col <- pick_first_col(nms, c("gene_symbol", "human_gene_symbol", "db_gene_symbol"))
  term_col <- pick_first_col(nms, c("gs_name", "geneset", "gs_exact_source", "term_name"))
  sub_col  <- pick_first_col(nms, c("gs_subcollection", "gs_subcollection_name", "gs_cat", "subcategory", "subcollection"))
  if (is.na(gene_col) || is.na(term_col)) {
    stop("Could not identify gene/term columns in msigdbr source ", source_name,
         "; columns=", paste(nms, collapse = ", "), call. = FALSE)
  }
  dt <- dt[, .(
    gene_norm = normalize_gene(get(gene_col)),
    term_id = as.character(get(term_col)),
    subcollection = if (!is.na(sub_col)) as.character(get(sub_col)) else NA_character_
  )]
  dt <- dt[!is.na(gene_norm) & gene_norm != "" & !is.na(term_id) & term_id != ""]
  unique(dt)
}

fetch_sources <- function(species) {
  c3 <- prep_msig_source(get_msig_table(species, collection = "C3"), "C3")
  # keep TFT only
  if ("subcollection" %in% names(c3)) {
    keep <- grepl("TFT", c3$subcollection, ignore.case = TRUE) | grepl("TFT", c3$term_id, ignore.case = TRUE)
    c3 <- c3[keep]
  }
  c3[, source := "TFT"]
  c3[, term_name := term_id]

  c2 <- prep_msig_source(get_msig_table(species, collection = "C2"), "C2")
  if ("subcollection" %in% names(c2)) {
    keep <- grepl("REACTOME", c2$subcollection, ignore.case = TRUE) | grepl("^REACTOME", c2$term_id, ignore.case = TRUE)
    c2 <- c2[keep]
  }
  c2[, source := "REACTOME"]
  c2[, term_name := gsub("^REACTOME_", "", term_id)]
  c2[, term_name := gsub("_", " ", term_name)]

  list(TFT = c3, REACTOME = c2)
}

ora_one_source <- function(program_dt, source_dt, source_name, q_cutoff, min_term_size, max_term_size) {
  universe <- sort(unique(source_dt$gene_norm))
  term_sizes <- source_dt[, .(term_size = uniqueN(gene_norm)), by = .(term_id, term_name)]
  term_keep <- term_sizes[term_size >= min_term_size & term_size <= max_term_size]
  source_use <- source_dt[term_keep, on = c("term_id", "term_name")]
  universe <- sort(unique(source_use$gene_norm))

  results <- list()
  overlap_rows <- list()
  k <- 1L
  for (prog in unique(program_dt$program)) {
    pg <- sort(unique(program_dt[program == prog, gene_norm]))
    pg_map <- intersect(pg, universe)
    n <- length(pg_map)
    if (n == 0L) {
      results[[k]] <- data.table(program = prog, source = source_name, term_id = character(), term_name = character(), universe_size = integer(), program_size = integer(), term_size = integer(), overlap_n = integer(), overlap_genes = character(), gene_ratio = numeric(), bg_ratio = numeric(), fold_enrichment = numeric(), p_value = numeric(), q_value = numeric(), neglog10_q = numeric())
      k <- k + 1L
      next
    }
    by_term <- source_use[, .(term_genes = list(sort(unique(gene_norm))), term_size = uniqueN(gene_norm)), by = .(term_id, term_name)]
    res <- by_term[, {
      overlap <- intersect(pg_map, term_genes[[1]])
      a <- length(overlap)
      b <- n - a
      c <- term_size - a
      d <- length(universe) - a - b - c
      mat <- matrix(c(a, b, c, d), nrow = 2)
      p <- fisher.test(mat, alternative = "greater")$p.value
      .(universe_size = length(universe), program_size = n, term_size = term_size, overlap_n = a,
        overlap_genes = paste(overlap, collapse = ","), gene_ratio = ifelse(n > 0, a / n, NA_real_),
        bg_ratio = term_size / length(universe), fold_enrichment = ifelse((n > 0 && term_size > 0), (a / n) / (term_size / length(universe)), NA_real_), p_value = p)
    }, by = .(term_id, term_name)]
    if (nrow(res) > 0) {
      res[, q_value := p.adjust(p_value, method = "BH")]
      res[, neglog10_q := -log10(pmax(q_value, 1e-300))]
      res[, source := source_name]
      res[, program := prog]
      setcolorder(res, c("program","source","term_id","term_name","universe_size","program_size","term_size","overlap_n","overlap_genes","gene_ratio","bg_ratio","fold_enrichment","p_value","q_value","neglog10_q"))
      ov <- res[overlap_n > 0][, .(program, source, term_id, term_name, gene = unlist(strsplit(overlap_genes, ",", fixed = TRUE))), by = 1:nrow(res[overlap_n > 0])][, rowid := NULL]
      if (nrow(ov) > 0) overlap_rows[[length(overlap_rows) + 1L]] <- ov
      results[[k]] <- res
      k <- k + 1L
    }
  }
  all_res <- rbindlist(results, fill = TRUE)
  sig_res <- all_res[q_value <= q_cutoff & overlap_n > 0][order(program, q_value, -fold_enrichment, term_name)]
  overlap_dt <- if (length(overlap_rows) > 0) rbindlist(overlap_rows, fill = TRUE) else data.table(program=character(), source=character(), term_id=character(), term_name=character(), gene=character())
  list(all = all_res, sig = sig_res, overlap = overlap_dt, universe = universe, term_n = uniqueN(source_use$term_id))
}

write_placeholder_pdf <- function(path, title, lines) {
  pdf(path, width = 10, height = 6)
  plot.new(); title(main = title); text(0.5, 0.6, paste(lines, collapse = "\n"), cex = 1)
  dev.off()
}

plot_dot <- function(dt, path, title, top_n = 12) {
  if (nrow(dt) == 0) {
    write_placeholder_pdf(path, title, c("No significant terms passed the cutoff."))
    return(invisible(NULL))
  }
  dd <- copy(dt)[order(program, q_value, -fold_enrichment)]
  dd <- dd[, head(.SD, top_n), by = program]
  dd[, label := factor(term_name, levels = rev(unique(term_name)))]
  p <- ggplot(dd, aes(x = program, y = label)) +
    geom_point(aes(size = overlap_n, color = neglog10_q)) +
    scale_color_continuous(name = "-log10(q)") +
    scale_size(name = "Overlap") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(path, plot = p, width = 10, height = max(5, 0.3 * nrow(dd) + 2))
}

plot_focus_heatmap <- function(dt, path, title) {
  if (nrow(dt) == 0) {
    write_placeholder_pdf(path, title, c("No focus regulatory/chromatin terms were detected."))
    return(invisible(NULL))
  }
  dd <- copy(dt)[order(source, q_value, -fold_enrichment)]
  dd <- dd[, head(.SD, 12), by = .(program, source)]
  dd[, term_lab := factor(paste0(source, ": ", term_name), levels = rev(unique(paste0(source, ": ", term_name))))]
  p <- ggplot(dd, aes(x = program, y = term_lab, fill = neglog10_q)) +
    geom_tile(color = "white") +
    labs(title = title, x = NULL, y = NULL, fill = "-log10(q)") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(path, plot = p, width = 10, height = max(5, 0.3 * nrow(dd) + 2))
}

main <- function() {
  opt <- parse_args()
  outdir <- file.path(opt$base_dir, "Step12B_upstream_regulator_refinement")
  tdir <- file.path(outdir, "tables")
  fdir <- file.path(outdir, "figures")
  dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  logfile <- file.path(outdir, "00_run.log")
  sink(logfile, split = TRUE)
  on.exit(sink(), add = TRUE)

  log_msg("Starting Step12B upstream regulator refinement")
  log_msg("base_dir=", opt$base_dir)
  log_msg("program_file=", opt$program_file)
  log_msg("requested_programs=", paste(opt$requested_programs, collapse = ", "))

  prog <- read_programs(opt$program_file, opt$requested_programs)
  fwrite(prog, file.path(tdir, "01_program_membership_used.tsv"), sep = "\t")

  sources <- fetch_sources(opt$species)
  source_summary <- data.table(
    source = c("TFT", "REACTOME"),
    n_terms = c(uniqueN(sources$TFT$term_id), uniqueN(sources$REACTOME$term_id)),
    n_genes = c(uniqueN(sources$TFT$gene_norm), uniqueN(sources$REACTOME$gene_norm))
  )

  tft <- ora_one_source(prog, sources$TFT, "TFT", opt$q_cutoff, opt$min_term_size, opt$max_term_size)
  react <- ora_one_source(prog, sources$REACTOME, "REACTOME", opt$q_cutoff, opt$min_term_size, opt$max_term_size)

  fwrite(tft$all, file.path(tdir, "10_tft_all_results.tsv"), sep = "\t")
  fwrite(tft$sig, file.path(tdir, "11_tft_significant_results.tsv"), sep = "\t")
  fwrite(react$all, file.path(tdir, "12_reactome_all_results.tsv"), sep = "\t")
  fwrite(react$sig, file.path(tdir, "13_reactome_significant_results.tsv"), sep = "\t")

  focus_regex <- opt$focus_regex
  focus <- rbindlist(list(
    tft$sig[grepl(focus_regex, term_name, ignore.case = TRUE) | grepl(focus_regex, term_id, ignore.case = TRUE)],
    react$sig[grepl(focus_regex, term_name, ignore.case = TRUE) | grepl(focus_regex, term_id, ignore.case = TRUE)]
  ), fill = TRUE)
  focus <- unique(focus)
  fwrite(focus, file.path(tdir, "14_focus_regulatory_terms.tsv"), sep = "\t")

  domain_summary <- rbindlist(list(
    tft$all[, .(source = "TFT", program, n_terms_tested = .N, n_terms_sig = sum(q_value <= opt$q_cutoff & overlap_n > 0),
                top_term = ifelse(.N > 0, term_name[which.min(q_value)], NA_character_),
                top_q_value = ifelse(.N > 0, min(q_value), NA_real_)), by = program],
    react$all[, .(source = "REACTOME", program, n_terms_tested = .N, n_terms_sig = sum(q_value <= opt$q_cutoff & overlap_n > 0),
                  top_term = ifelse(.N > 0, term_name[which.min(q_value)], NA_character_),
                  top_q_value = ifelse(.N > 0, min(q_value), NA_real_)), by = program]
  ), fill = TRUE)
  fwrite(domain_summary, file.path(tdir, "15_source_summary.tsv"), sep = "\t")

  overlap_summary <- rbindlist(list(
    prog[, .(requested = uniqueN(gene_norm)), by = program],
    data.table(program = unique(prog$program), requested = NA_integer_)
  ))[!duplicated(program)]
  overlap_summary[, TFT_mapped := sapply(program, function(p) uniqueN(intersect(prog[program == p, gene_norm], tft$universe)))]
  overlap_summary[, REACTOME_mapped := sapply(program, function(p) uniqueN(intersect(prog[program == p, gene_norm], react$universe)))]
  fwrite(overlap_summary, file.path(tdir, "00_input_summary.tsv"), sep = "\t")

  plot_dot(tft$sig, file.path(fdir, "21_tft_dotplot.pdf"), "Step12B TFT target-set refinement", opt$plot_top_n)
  plot_dot(react$sig, file.path(fdir, "22_reactome_dotplot.pdf"), "Step12B Reactome refinement", opt$plot_top_n)
  plot_focus_heatmap(focus, file.path(fdir, "23_focus_regulatory_heatmap.pdf"), "Step12B focus regulatory/chromatin modules")

  wb <- createWorkbook()
  addWorksheet(wb, "input_summary"); writeData(wb, "input_summary", overlap_summary)
  addWorksheet(wb, "source_summary"); writeData(wb, "source_summary", source_summary)
  addWorksheet(wb, "tft_sig"); writeData(wb, "tft_sig", tft$sig)
  addWorksheet(wb, "reactome_sig"); writeData(wb, "reactome_sig", react$sig)
  addWorksheet(wb, "focus_terms"); writeData(wb, "focus_terms", focus)
  saveWorkbook(wb, file.path(outdir, "18_upstream_refinement_workbook.xlsx"), overwrite = TRUE)

  summary_file <- file.path(outdir, "17_upstream_summary.txt")
  con <- file(summary_file, open = "wt")
  writeLines(c(
    "Step12B upstream regulator refinement summary",
    "",
    paste0("Program file: ", opt$program_file),
    paste0("Programs tested: ", paste(opt$requested_programs, collapse = ", ")),
    paste0("msigdbr species: ", opt$species),
    paste0("q_cutoff: ", opt$q_cutoff),
    "",
    "MSigDB collection basis:",
    "- TFT source was derived from the C3 regulatory-target collection (TFT subset).",
    "- Pathway source was derived from the C2 canonical-pathway collection (Reactome subset).",
    "",
    "Program overlap with TFT / Reactome universes:"
  ), con)
  for (i in seq_len(nrow(overlap_summary))) {
    rr <- overlap_summary[i]
    writeLines(sprintf("  - %s: requested=%d, TFT_mapped=%d, REACTOME_mapped=%d", rr$program, rr$requested, rr$TFT_mapped, rr$REACTOME_mapped), con)
  }
  writeLines(c("", "Top significant TFT terms by program:"), con)
  for (p in opt$requested_programs) {
    writeLines(sprintf("Program: %s", p), con)
    dd <- tft$sig[program == p][order(q_value, -fold_enrichment)]
    if (nrow(dd) == 0) {
      writeLines("  No TFT terms passed q cutoff.", con)
    } else {
      for (j in seq_len(min(6, nrow(dd)))) {
        rr <- dd[j]
        writeLines(sprintf("  %s | overlap=%d | q=%.3g | FE=%.2f", rr$term_name, rr$overlap_n, rr$q_value, rr$fold_enrichment), con)
      }
    }
    writeLines("", con)
  }
  writeLines(c("Top significant Reactome terms by program:"), con)
  for (p in opt$requested_programs) {
    writeLines(sprintf("Program: %s", p), con)
    dd <- react$sig[program == p][order(q_value, -fold_enrichment)]
    if (nrow(dd) == 0) {
      writeLines("  No Reactome terms passed q cutoff.", con)
    } else {
      for (j in seq_len(min(6, nrow(dd)))) {
        rr <- dd[j]
        writeLines(sprintf("  %s | overlap=%d | q=%.3g | FE=%.2f", rr$term_name, rr$overlap_n, rr$q_value, rr$fold_enrichment), con)
      }
    }
    writeLines("", con)
  }
  close(con)

  log_msg("Completed Step12B upstream regulator refinement")
}

main()
