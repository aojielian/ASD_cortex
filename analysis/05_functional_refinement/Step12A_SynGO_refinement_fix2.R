#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    base_dir = NULL,
    program_file = NULL,
    syngo_file = NULL,
    syngo_format = "auto",
    requested_programs = "SFARI_all,midPrenatal_SFARI_top20,midPrenatal_SFARI_top10,midPrenatal_SFARI_top05",
    min_term_size = 5,
    max_term_size = 500,
    q_cutoff = 0.10,
    plot_top_n = 12
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(args)) stop("Missing value for --", key)
    val <- args[[i + 1L]]
    out[[key]] <- val
    i <- i + 2L
  }
  out$min_term_size <- as.integer(out$min_term_size)
  out$max_term_size <- as.integer(out$max_term_size)
  out$q_cutoff <- as.numeric(out$q_cutoff)
  out$plot_top_n <- as.integer(out$plot_top_n)
  out
}

log_msg <- function(...) {
  txt <- paste0(..., collapse = "")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), txt))
}

pick_first_col <- function(nms, candidates) {
  idx <- match(tolower(candidates), tolower(nms))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_character_)
  nms[idx[1]]
}

read_programs <- function(program_file, requested_programs = NULL) {
  dt <- fread(program_file)
  pcol <- pick_first_col(names(dt), c("program", "program_std", "program_name", "set", "geneset", "gene_set", "geneset_name"))
  gcol <- pick_first_col(names(dt), c("gene", "gene_std", "gene_symbol", "symbol", "hgnc_symbol"))
  if (is.na(pcol) || is.na(gcol)) {
    stop(
      "Could not identify gene/program columns in program_file: ", program_file,
      "; columns=", paste(names(dt), collapse = ", ")
    )
  }
  dt <- dt[!is.na(get(pcol)) & !is.na(get(gcol))]
  dt[, program := as.character(get(pcol))]
  dt[, gene := toupper(as.character(get(gcol)))]
  dt <- unique(dt[, .(program, gene)])
  if (!is.null(requested_programs)) {
    req <- trimws(unlist(strsplit(requested_programs, ",")))
    dt <- dt[program %in% req]
  }
  split(dt$gene, dt$program)
}

read_syngo <- function(syngo_file, syngo_format = "auto") {
  if (!file.exists(syngo_file)) stop("SynGO file not found: ", syngo_file)
  ext <- tolower(tools::file_ext(syngo_file))
  fmt <- tolower(syngo_format)
  if (fmt == "auto") fmt <- ext

  if (fmt %in% c("xlsx", "xls")) {
    sheets <- getSheetNames(syngo_file)
    if (length(sheets) == 0) stop("No sheets found in SynGO file: ", syngo_file)
    # User's official file has a single Sheet1 with these columns.
    found <- NULL
    for (s in sheets) {
      x <- as.data.table(read.xlsx(syngo_file, sheet = s, colNames = TRUE))
      if (nrow(x) == 0) next
      if (all(c("hgnc_symbol", "go_id", "go_name", "go_domain") %in% names(x))) {
        found <- x
        break
      }
    }
    if (is.null(found)) {
      stop("Could not parse any usable SynGO annotation sheet from: ", syngo_file)
    }
    dt <- found[, .(
      gene = toupper(trimws(as.character(hgnc_symbol))),
      term_id = trimws(as.character(go_id)),
      term_name = trimws(as.character(go_name)),
      domain = trimws(as.character(go_domain))
    )]
  } else if (fmt %in% c("tsv", "txt", "csv")) {
    sep <- if (fmt == "csv") "," else "\t"
    x <- fread(syngo_file, sep = sep)
    gcol <- pick_first_col(names(x), c("hgnc_symbol", "gene_symbol", "symbol", "gene"))
    tidcol <- pick_first_col(names(x), c("go_id", "term_id", "id"))
    tncol <- pick_first_col(names(x), c("go_name", "term_name", "name"))
    dcol <- pick_first_col(names(x), c("go_domain", "domain"))
    if (any(is.na(c(gcol, tidcol, tncol, dcol)))) {
      stop("Could not identify required SynGO columns in flat file: ", syngo_file,
           "; columns=", paste(names(x), collapse = ", "))
    }
    dt <- x[, .(
      gene = toupper(trimws(as.character(get(gcol)))),
      term_id = trimws(as.character(get(tidcol))),
      term_name = trimws(as.character(get(tncol))),
      domain = trimws(as.character(get(dcol)))
    )]
  } else {
    stop("Unsupported syngo_format: ", syngo_format)
  }

  dt <- dt[gene != "" & term_id != "" & term_name != "" & domain != ""]
  dt <- unique(dt)
  dt
}

enrich_one_program <- function(program_name, genes, syngo_dt, min_term_size = 5L, max_term_size = 500L) {
  universe <- unique(syngo_dt$gene)
  genes <- unique(toupper(genes))
  prog_genes <- intersect(genes, universe)
  M <- length(universe)
  N <- length(prog_genes)
  if (N == 0) return(data.table())

  term_dt <- syngo_dt[, .(term_genes = list(unique(gene)), term_size = uniqueN(gene), domain = first(domain), term_name = first(term_name)), by = term_id]
  term_dt <- term_dt[term_size >= min_term_size & term_size <= max_term_size]
  if (nrow(term_dt) == 0) return(data.table())

  res <- rbindlist(lapply(seq_len(nrow(term_dt)), function(i) {
    tg <- term_dt$term_genes[[i]]
    overlap <- intersect(prog_genes, tg)
    k <- length(overlap)
    K <- term_dt$term_size[[i]]
    pval <- phyper(q = k - 1, m = K, n = M - K, k = N, lower.tail = FALSE)
    data.table(
      program = program_name,
      term_id = term_dt$term_id[[i]],
      term_name = term_dt$term_name[[i]],
      domain = term_dt$domain[[i]],
      universe_size = M,
      program_size = N,
      term_size = K,
      overlap_n = k,
      overlap_genes = paste(sort(overlap), collapse = ","),
      gene_ratio = if (N > 0) k / N else NA_real_,
      bg_ratio = if (M > 0) K / M else NA_real_,
      fold_enrichment = if (N > 0 && M > 0 && K > 0) (k / N) / (K / M) else NA_real_,
      p_value = pval
    )
  }), fill = TRUE)
  res[, q_value := p.adjust(p_value, method = "BH"), by = program]
  setorder(res, q_value, -fold_enrichment, -overlap_n)
  res
}

safe_filename <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

main <- function() {
  opt <- parse_args()
  if (is.null(opt$base_dir)) stop("--base_dir is required")
  if (is.null(opt$program_file)) stop("--program_file is required")
  if (is.null(opt$syngo_file)) stop("--syngo_file is required")

  outdir <- file.path(opt$base_dir, "Step12A_SynGO_refinement")
  tabdir <- file.path(outdir, "tables")
  figdir <- file.path(outdir, "figures")
  dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

  log_msg("Starting Step12A SynGO refinement (fix2)")
  log_msg("base_dir=", opt$base_dir)
  log_msg("program_file=", opt$program_file)
  log_msg("syngo_file=", opt$syngo_file)

  gene_sets <- read_programs(opt$program_file, opt$requested_programs)
  syngo_dt <- read_syngo(opt$syngo_file, opt$syngo_format)

  input_summary <- data.table(
    metric = c("n_programs", "n_syngo_annotations", "n_syngo_terms", "n_syngo_genes"),
    value = c(length(gene_sets), nrow(syngo_dt), uniqueN(syngo_dt$term_id), uniqueN(syngo_dt$gene))
  )
  fwrite(input_summary, file.path(tabdir, "00_input_summary.tsv"), sep = "\t")

  overlap_tbl <- rbindlist(lapply(names(gene_sets), function(p) {
    gs <- unique(toupper(gene_sets[[p]]))
    present <- intersect(gs, unique(syngo_dt$gene))
    data.table(
      program = p,
      requested_genes = length(gs),
      syngo_mapped_genes = length(present),
      mapping_fraction = round(length(present) / max(1, length(gs)), 4),
      mapped_gene_list = paste(sort(present), collapse = ",")
    )
  }))
  fwrite(overlap_tbl, file.path(tabdir, "01_program_gene_overlap.tsv"), sep = "\t")

  all_res <- rbindlist(lapply(names(gene_sets), function(p) {
    enrich_one_program(p, gene_sets[[p]], syngo_dt, opt$min_term_size, opt$max_term_size)
  }), fill = TRUE)
  if (nrow(all_res) == 0) stop("No enrichment results were produced. Check program/SynGO overlap.")

  all_res[, neglog10_q := -log10(pmax(q_value, 1e-300))]
  fwrite(all_res, file.path(tabdir, "10_all_enrichment_results.tsv"), sep = "\t")

  sig <- all_res[q_value <= opt$q_cutoff & overlap_n >= 1]
  fwrite(sig, file.path(tabdir, "11_significant_enrichment_results.tsv"), sep = "\t")

  lead <- sig[, .(program, term_id, term_name, domain, overlap_n, overlap_genes)]
  fwrite(lead, file.path(tabdir, "12_overlap_genes.tsv"), sep = "\t")

  domain_summary <- all_res[, .(
    n_terms_tested = .N,
    n_terms_sig = sum(q_value <= opt$q_cutoff),
    top_term = term_name[which.min(q_value)],
    top_q_value = min(q_value, na.rm = TRUE)
  ), by = .(program, domain)]
  fwrite(domain_summary, file.path(tabdir, "13_domain_summary.tsv"), sep = "\t")

  wb <- createWorkbook()
  addWorksheet(wb, "input_summary"); writeData(wb, "input_summary", input_summary)
  addWorksheet(wb, "program_overlap"); writeData(wb, "program_overlap", overlap_tbl)
  addWorksheet(wb, "all_results"); writeData(wb, "all_results", all_res)
  addWorksheet(wb, "significant"); writeData(wb, "significant", sig)
  addWorksheet(wb, "domain_summary"); writeData(wb, "domain_summary", domain_summary)
  saveWorkbook(wb, file.path(outdir, "14_SynGO_refinement_workbook.xlsx"), overwrite = TRUE)

  # Plot 1: dotplot of top terms per program
  dot_dt <- copy(all_res)
  dot_dt <- dot_dt[order(program, q_value, -fold_enrichment)]
  dot_dt <- dot_dt[, head(.SD, opt$plot_top_n), by = program]
  dot_dt[, label := paste0(term_name, " [", domain, "]")]
  p1 <- ggplot(dot_dt, aes(x = program, y = reorder(label, neglog10_q), size = overlap_n, color = neglog10_q)) +
    geom_point(alpha = 0.9) +
    labs(x = NULL, y = NULL, title = "SynGO enrichment top terms", color = "-log10(q)", size = "Overlap") +
    theme_bw(base_size = 11)
  ggsave(file.path(figdir, "21_dotplot_top_terms.pdf"), p1, width = 10, height = 8)

  # Plot 2: heatmap of top unique terms
  heat_dt <- copy(sig)
  if (nrow(heat_dt) == 0) heat_dt <- dot_dt
  heat_dt <- heat_dt[order(q_value, -fold_enrichment)]
  top_terms <- unique(heat_dt$term_name)[seq_len(min(opt$plot_top_n, uniqueN(heat_dt$term_name)))]
  heat_dt <- all_res[term_name %in% top_terms]
  p2 <- ggplot(heat_dt, aes(x = program, y = term_name, fill = neglog10_q)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick") +
    labs(x = NULL, y = NULL, title = "SynGO term-program heatmap", fill = "-log10(q)") +
    theme_bw(base_size = 11)
  ggsave(file.path(figdir, "22_heatmap_term_program.pdf"), p2, width = 8, height = 7)

  # Plot 3: domain counts
  dom_plot <- domain_summary[, .(n_sig = sum(n_terms_sig)), by = .(program, domain)]
  p3 <- ggplot(dom_plot, aes(x = program, y = n_sig, fill = domain)) +
    geom_col(position = "stack") +
    labs(x = NULL, y = "# significant terms", title = "SynGO significant terms by domain") +
    theme_bw(base_size = 11)
  ggsave(file.path(figdir, "23_domain_counts.pdf"), p3, width = 8, height = 5)

  # Text summary
  con <- file(file.path(outdir, "12_SynGO_summary.txt"), open = "wt")
  writeLines(c(
    "Step12A SynGO refinement summary",
    "",
    paste0("SynGO file: ", opt$syngo_file),
    paste0("Programs tested: ", paste(names(gene_sets), collapse = ", ")),
    paste0("SynGO unique genes: ", uniqueN(syngo_dt$gene)),
    paste0("SynGO unique terms: ", uniqueN(syngo_dt$term_id)),
    paste0("q_cutoff: ", opt$q_cutoff),
    "",
    "Program overlap with SynGO universe:"
  ), con)
  for (i in seq_len(nrow(overlap_tbl))) {
    writeLines(sprintf("  - %s: requested=%s, mapped=%s, fraction=%s",
                       overlap_tbl$program[i], overlap_tbl$requested_genes[i], overlap_tbl$syngo_mapped_genes[i], overlap_tbl$mapping_fraction[i]), con)
  }
  writeLines("", con)
  writeLines("Top significant SynGO terms by program:", con)
  for (p in names(gene_sets)) {
    writeLines(paste0("Program: ", p), con)
    sub <- sig[program == p][order(q_value, -fold_enrichment)]
    if (nrow(sub) == 0) {
      writeLines("  No terms passed q cutoff.", con)
    } else {
      for (i in seq_len(min(8, nrow(sub)))) {
        writeLines(sprintf("  %s | %s | %s | overlap=%s | q=%.3g | FE=%.2f",
                           sub$term_name[i], sub$term_id[i], sub$domain[i], sub$overlap_n[i], sub$q_value[i], sub$fold_enrichment[i]), con)
      }
    }
    writeLines("", con)
  }
  close(con)

  log_msg("Completed Step12A SynGO refinement (fix2)")
}

main()
