#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
})

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    base_dir = NULL,
    program_file = NULL,
    syngo_file = NULL,
    requested_programs = c("SFARI_all", "midPrenatal_SFARI_top20", "midPrenatal_SFARI_top10", "midPrenatal_SFARI_top05"),
    expression_ref_file = NULL,
    n_perm = 1000L,
    n_bins = 10L,
    q_cutoff = 0.10,
    seed = 123L,
    outdir = NULL
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (grepl("^--[^=]+=", a)) {
      kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
      key <- kv[1]
      val <- paste(kv[-1], collapse = "=")
    } else if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) stop("Missing value for ", a)
      val <- args[[i + 1L]]
      i <- i + 1L
    } else {
      stop("Unexpected argument: ", a)
    }
    if (key %in% c("base_dir", "program_file", "syngo_file", "expression_ref_file", "outdir")) {
      out[[key]] <- val
    } else if (key == "requested_programs") {
      out[[key]] <- trimws(unlist(strsplit(val, ",", fixed = TRUE)))
    } else if (key %in% c("n_perm", "n_bins", "seed")) {
      out[[key]] <- as.integer(val)
    } else if (key == "q_cutoff") {
      out[[key]] <- as.numeric(val)
    } else {
      stop("Unknown argument: --", key)
    }
    i <- i + 1L
  }
  if (is.null(out$base_dir)) stop("--base_dir is required")
  if (is.null(out$program_file)) out$program_file <- file.path(out$base_dir, "tables", "20_midPrenatal_core_programs.tsv")
  if (is.null(out$syngo_file)) out$syngo_file <- file.path(out$base_dir, "inputs", "SynGO_annotations.xlsx")
  if (is.null(out$expression_ref_file)) {
    candidates <- c(
      file.path(out$base_dir, "inputs", "Gandal2022.standardized_expr.tsv.gz"),
      file.path(out$base_dir, "inputs", "GSE102741.standardized_expr.tsv.gz"),
      file.path(out$base_dir, "inputs", "GSE64018.standardized_expr.tsv.gz")
    )
    pick <- candidates[file.exists(candidates)][1]
    out$expression_ref_file <- if (length(pick)) pick else NA_character_
  }
  if (is.null(out$outdir)) out$outdir <- file.path(out$base_dir, "Step13A_specificity_matched_random_controls")
  out
}

ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

norm_gene <- function(x) toupper(trimws(as.character(x)))

pick_first_col <- function(nms, candidates) {
  idx <- match(candidates, nms)
  idx <- idx[!is.na(idx)]
  if (!length(idx)) return(NA_character_)
  nms[idx[1]]
}

read_programs <- function(program_file, requested_programs) {
  dt <- fread(program_file)
  setnames(dt, names(dt), trimws(names(dt)))
  nms <- names(dt)
  program_col <- pick_first_col(nms, c("program", "program_std", "program_name", "set", "geneset", "gene_set", "geneset_name"))
  gene_col <- pick_first_col(nms, c("gene", "gene_std", "gene_symbol", "symbol", "hgnc_symbol"))
  if (is.na(program_col) || is.na(gene_col)) {
    stop("Could not identify gene/program columns in program_file: ", program_file,
         "; columns=", paste(names(dt), collapse = ", "))
  }
  dt <- dt[, .(program = as.character(get(program_col)), gene = norm_gene(get(gene_col)))]
  dt <- dt[program %in% requested_programs & !is.na(gene) & gene != ""]
  unique(dt)
}

read_syngo <- function(syngo_file) {
  if (!file.exists(syngo_file)) stop("SynGO file not found: ", syngo_file)
  ext <- tolower(tools::file_ext(syngo_file))
  dt <- NULL
  if (ext %in% c("xlsx", "xls")) {
    sheets <- openxlsx::getSheetNames(syngo_file)
    dt <- as.data.table(openxlsx::read.xlsx(syngo_file, sheet = sheets[1]))
  } else {
    dt <- fread(syngo_file)
  }
  setnames(dt, names(dt), trimws(names(dt)))
  gene_col <- pick_first_col(names(dt), c("hgnc_symbol", "gene_symbol", "symbol", "gene"))
  term_id_col <- pick_first_col(names(dt), c("go_id", "term_id", "id"))
  term_name_col <- pick_first_col(names(dt), c("go_name", "term_name", "name"))
  domain_col <- pick_first_col(names(dt), c("go_domain", "domain"))
  if (any(is.na(c(gene_col, term_id_col, term_name_col, domain_col)))) {
    stop("Could not parse SynGO annotations from: ", syngo_file,
         "; columns=", paste(names(dt), collapse = ", "))
  }
  out <- unique(dt[, .(
    gene = norm_gene(get(gene_col)),
    term_id = as.character(get(term_id_col)),
    term_name = as.character(get(term_name_col)),
    domain = as.character(get(domain_col))
  )])
  out <- out[!is.na(gene) & gene != "" & !is.na(term_id) & term_id != ""]
  out[, source := "SynGO"]
  out
}

read_expression_reference <- function(expr_file) {
  if (is.na(expr_file) || !nzchar(expr_file) || !file.exists(expr_file)) {
    return(NULL)
  }
  log_msg("Reading expression reference:", expr_file)
  dt <- fread(expr_file)
  setnames(dt, names(dt), trimws(names(dt)))
  gene_col <- pick_first_col(names(dt), c("gene", "gene_symbol", "symbol", "hgnc_symbol"))
  if (is.na(gene_col)) gene_col <- names(dt)[1]
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, gene_col)
  if (!length(num_cols)) stop("No numeric expression columns found in expression reference: ", expr_file)
  dt[, gene := norm_gene(get(gene_col))]
  dt <- dt[!is.na(gene) & gene != ""]
  dt[, mean_expr := rowMeans(.SD, na.rm = TRUE), .SDcols = num_cols]
  unique(dt[, .(gene, mean_expr)])
}

build_source_universes <- function(syngo_dt) {
  reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
  setnames(reactome, names(reactome), trimws(names(reactome)))
  reactome_dt <- unique(reactome[, .(
    gene = norm_gene(gene_symbol),
    term_id = as.character(gs_id),
    term_name = as.character(gs_name),
    domain = "REACTOME"
  )])
  reactome_dt[, source := "Reactome"]

  tft <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "TFT:GTRD"))
  if (!nrow(tft)) {
    tft <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "TFT"))
  }
  setnames(tft, names(tft), trimws(names(tft)))
  tft_dt <- unique(tft[, .(
    gene = norm_gene(gene_symbol),
    term_id = as.character(gs_id),
    term_name = as.character(gs_name),
    domain = "TFT"
  )])
  tft_dt[, source := "TFT"]

  list(SynGO = syngo_dt, Reactome = reactome_dt, TFT = tft_dt)
}

compute_bins <- function(expr_dt, universe_genes, n_bins = 10L) {
  x <- expr_dt[gene %in% universe_genes]
  if (!nrow(x)) return(NULL)
  x <- unique(x, by = "gene")
  qs <- unique(quantile(x$mean_expr, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE, type = 7))
  if (length(qs) <= 2L) {
    x[, expr_bin := 1L]
  } else {
    x[, expr_bin := cut(mean_expr, breaks = qs, include.lowest = TRUE, labels = FALSE)]
  }
  x[, .(gene, mean_expr, expr_bin)]
}

ora_metrics <- function(program_genes, term2gene, q_cutoff = 0.10, focus_pattern = NULL) {
  universe_genes <- sort(unique(term2gene$gene))
  pg <- sort(unique(program_genes[program_genes %in% universe_genes]))
  M <- length(universe_genes)
  N <- length(pg)
  terms <- unique(term2gene[, .(term_id, term_name, domain)])
  term_sets <- split(term2gene$gene, term2gene$term_id)
  term_sizes <- vapply(term_sets, function(v) length(unique(v)), numeric(1))
  term_ids <- names(term_sets)
  overlap_n <- integer(length(term_ids))
  pvals <- numeric(length(term_ids))
  overlaps <- vector("list", length(term_ids))
  for (i in seq_along(term_ids)) {
    tg <- unique(term_sets[[i]])
    k <- sum(pg %in% tg)
    overlap_n[i] <- k
    overlaps[[i]] <- paste(sort(intersect(pg, tg)), collapse = ",")
    pvals[i] <- phyper(q = k - 1L, m = length(tg), n = M - length(tg), k = N, lower.tail = FALSE)
  }
  qvals <- p.adjust(pvals, method = "BH")
  res <- data.table(
    term_id = term_ids,
    term_size = as.integer(term_sizes[term_ids]),
    overlap_n = overlap_n,
    overlap_genes = unlist(overlaps),
    p_value = pvals,
    q_value = qvals
  )
  res <- merge(res, terms, by = "term_id", all.x = TRUE)
  res[, universe_size := M]
  res[, program_size := N]
  res[, gene_ratio := ifelse(program_size > 0, overlap_n / program_size, 0)]
  res[, bg_ratio := ifelse(universe_size > 0, term_size / universe_size, 0)]
  res[, fold_enrichment := ifelse(bg_ratio > 0, gene_ratio / bg_ratio, 0)]
  setorder(res, q_value, -fold_enrichment, -overlap_n)

  out <- list(
    all = res,
    metrics = data.table(
      n_sig_terms = sum(res$q_value < q_cutoff, na.rm = TRUE),
      best_q = ifelse(any(is.finite(res$q_value)), min(res$q_value, na.rm = TRUE), 1),
      best_neglog10q = ifelse(any(is.finite(res$q_value)) && min(res$q_value, na.rm = TRUE) > 0, -log10(min(res$q_value, na.rm = TRUE)), 0),
      best_fold_enrichment = ifelse(nrow(res), res$fold_enrichment[1], 0),
      best_overlap_n = ifelse(nrow(res), res$overlap_n[1], 0)
    )
  )

  if (!is.null(focus_pattern)) {
    focus <- res[grepl(focus_pattern, term_name, ignore.case = TRUE)]
    out$metrics[, focus_n_sig_terms := sum(focus$q_value < q_cutoff, na.rm = TRUE)]
    out$metrics[, focus_best_q := ifelse(nrow(focus), min(focus$q_value, na.rm = TRUE), 1)]
    out$metrics[, focus_best_neglog10q := ifelse(nrow(focus) && min(focus$q_value, na.rm = TRUE) > 0, -log10(min(focus$q_value, na.rm = TRUE)), 0)]
  }
  out
}

sample_expression_matched <- function(obs_genes, expr_bins_dt, universe_genes) {
  u <- expr_bins_dt[gene %in% universe_genes]
  o <- expr_bins_dt[gene %in% obs_genes]
  if (!nrow(u) || !nrow(o)) {
    return(sample(universe_genes, length(obs_genes), replace = FALSE))
  }
  cnt <- o[, .N, by = expr_bin]
  sampled <- character(0)
  for (i in seq_len(nrow(cnt))) {
    b <- cnt$expr_bin[i]
    n_need <- cnt$N[i]
    pool <- setdiff(u[expr_bin == b, gene], sampled)
    if (length(pool) < n_need) {
      pool <- unique(c(pool, setdiff(u$gene, sampled)))
    }
    if (length(pool) < n_need) {
      pool <- unique(c(pool, universe_genes))
    }
    sampled <- c(sampled, sample(pool, n_need, replace = FALSE))
  }
  unique(sampled)
}

empirical_p <- function(rand, obs, larger_better = TRUE) {
  if (!length(rand) || all(is.na(rand))) return(NA_real_)
  if (larger_better) {
    (1 + sum(rand >= obs, na.rm = TRUE)) / (length(rand) + 1)
  } else {
    (1 + sum(rand <= obs, na.rm = TRUE)) / (length(rand) + 1)
  }
}

main <- function() {
  args <- parse_args()
  set.seed(args$seed)

  ensure_dir(args$outdir)
  ensure_dir(file.path(args$outdir, "tables"))
  ensure_dir(file.path(args$outdir, "figures"))

  flog <- file.path(args$outdir, "00_run.log")
  sink(flog, split = TRUE)
  on.exit(sink(), add = TRUE)

  log_msg("Starting Step13A specificity / matched-random control analyses")
  log_msg("base_dir=", args$base_dir)
  log_msg("program_file=", args$program_file)
  log_msg("syngo_file=", args$syngo_file)
  log_msg("expression_ref_file=", args$expression_ref_file)
  log_msg("requested_programs=", paste(args$requested_programs, collapse = ", "))

  programs <- read_programs(args$program_file, args$requested_programs)
  syngo_dt <- read_syngo(args$syngo_file)
  expr_dt <- read_expression_reference(args$expression_ref_file)
  sources <- build_source_universes(syngo_dt)

  focus_patterns <- list(
    SynGO = "presyn|postsyn|synap|active zone|vesicle|neurexin|neuroligin|density",
    Reactome = "chromatin|epigen|remodel|baf|npbaf|nbaf|transcription|runx1|mediator|histone|nucleosome"
  )

  membership <- copy(programs)
  fwrite(membership, file.path(args$outdir, "tables", "01_program_membership_used.tsv"), sep = "\t")

  bg_summary <- rbindlist(lapply(names(sources), function(src) {
    td <- sources[[src]]
    data.table(source = src, universe_size = uniqueN(td$gene), term_size_n = uniqueN(td$term_id))
  }))
  fwrite(bg_summary, file.path(args$outdir, "tables", "02_source_universe_summary.tsv"), sep = "\t")

  if (!is.null(expr_dt)) {
    fwrite(expr_dt, file.path(args$outdir, "tables", "03_expression_reference_summary.tsv"), sep = "\t")
  }

  observed_list <- list()
  all_random <- list()
  emp_list <- list()
  z_list <- list()
  input_rows <- list()

  for (src in c("SynGO", "Reactome")) {
    td <- sources[[src]]
    universe_genes <- sort(unique(td$gene))
    expr_bins <- if (!is.null(expr_dt)) compute_bins(expr_dt, universe_genes, n_bins = args$n_bins) else NULL

    for (prog in args$requested_programs) {
      req_genes <- unique(programs[program == prog, gene])
      obs_genes <- sort(unique(req_genes[req_genes %in% universe_genes]))
      obs <- ora_metrics(obs_genes, td, q_cutoff = args$q_cutoff, focus_pattern = focus_patterns[[src]])
      om <- copy(obs$metrics)
      om[, `:=`(program = prog, source = src, requested_size = length(req_genes), mapped_size = length(obs_genes), control_type = "observed")]
      observed_list[[paste(src, prog, sep = "__")]] <- om
      input_rows[[paste(src, prog, sep = "__")]] <- data.table(program = prog, source = src, requested_size = length(req_genes), mapped_size = length(obs_genes), universe_size = length(universe_genes))

      if (!length(obs_genes)) next

      rand_src <- vector("list", args$n_perm * 2L)
      ridx <- 1L
      for (b in seq_len(args$n_perm)) {
        rand_size <- sample(universe_genes, length(obs_genes), replace = FALSE)
        m1 <- ora_metrics(rand_size, td, q_cutoff = args$q_cutoff, focus_pattern = focus_patterns[[src]])$metrics
        m1[, `:=`(program = prog, source = src, control_type = "size_matched", iter = b)]
        rand_src[[ridx]] <- m1; ridx <- ridx + 1L

        if (!is.null(expr_bins)) {
          rand_expr <- sample_expression_matched(obs_genes, expr_bins, universe_genes)
        } else {
          rand_expr <- sample(universe_genes, length(obs_genes), replace = FALSE)
        }
        m2 <- ora_metrics(rand_expr, td, q_cutoff = args$q_cutoff, focus_pattern = focus_patterns[[src]])$metrics
        m2[, `:=`(program = prog, source = src, control_type = "expression_matched", iter = b)]
        rand_src[[ridx]] <- m2; ridx <- ridx + 1L
      }
      rand_src <- rbindlist(rand_src, fill = TRUE)
      all_random[[paste(src, prog, sep = "__")]] <- rand_src

      metric_names <- setdiff(names(om), c("program", "source", "requested_size", "mapped_size", "control_type"))
      for (ctype in c("size_matched", "expression_matched")) {
        sub <- rand_src[control_type == ctype]
        ep <- data.table(program = prog, source = src, control_type = ctype)
        zs <- data.table(program = prog, source = src, control_type = ctype)
        for (m in metric_names) {
          randv <- sub[[m]]
          obsv <- om[[m]][1]
          larger_better <- !grepl("best_q$|focus_best_q$", m)
          ep[[paste0(m, "_emp_p")]] <- empirical_p(randv, obsv, larger_better = larger_better)
          zs[[paste0(m, "_z")]] <- if (sd(randv, na.rm = TRUE) > 0) (obsv - mean(randv, na.rm = TRUE)) / sd(randv, na.rm = TRUE) else NA_real_
        }
        emp_list[[paste(src, prog, ctype, sep = "__")]] <- ep
        z_list[[paste(src, prog, ctype, sep = "__")]] <- zs
      }
    }
  }

  observed_dt <- rbindlist(observed_list, fill = TRUE)
  random_dt <- rbindlist(all_random, fill = TRUE)
  empirical_dt <- rbindlist(emp_list, fill = TRUE)
  z_dt <- rbindlist(z_list, fill = TRUE)
  input_dt <- rbindlist(input_rows, fill = TRUE)

  fwrite(input_dt, file.path(args$outdir, "tables", "00_input_summary.tsv"), sep = "\t")
  fwrite(observed_dt, file.path(args$outdir, "tables", "10_observed_metrics.tsv"), sep = "\t")
  fwrite(random_dt, file.path(args$outdir, "tables", "11_random_metrics_long.tsv.gz"), sep = "\t")
  fwrite(empirical_dt, file.path(args$outdir, "tables", "12_empirical_pvalues.tsv"), sep = "\t")
  fwrite(z_dt, file.path(args$outdir, "tables", "13_metric_zscores.tsv"), sep = "\t")

  focus_summary <- merge(observed_dt, empirical_dt, by = c("program", "source"), allow.cartesian = TRUE)
  fwrite(focus_summary, file.path(args$outdir, "tables", "14_focus_summary.tsv"), sep = "\t")

  # plots
  pdat <- merge(observed_dt[, .(program, source, n_sig_terms, best_neglog10q, focus_n_sig_terms, focus_best_neglog10q)],
                empirical_dt, by = c("program", "source"), allow.cartesian = TRUE)
  melt_dat <- melt(pdat,
                   id.vars = c("program", "source", "control_type"),
                   measure.vars = patterns("_emp_p$"),
                   variable.name = "metric_idx",
                   value.name = "emp_p")
  metric_names <- names(pdat)[grep("_emp_p$", names(pdat))]
  melt_dat[, metric := sub("_emp_p$", "", metric_names[metric_idx])]
  g1 <- ggplot(melt_dat, aes(x = metric, y = -log10(pmax(emp_p, 1e-6)), color = control_type)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    facet_grid(source ~ program, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Empirical specificity significance", y = "-log10 empirical p", x = "Metric")
  ggsave(file.path(args$outdir, "figures", "21_empirical_significance_dotplot.pdf"), g1, width = 15, height = 8)

  dens_targets <- rbindlist(lapply(c("SynGO", "Reactome"), function(src) {
    metric <- if (src == "SynGO") "best_neglog10q" else "focus_best_neglog10q"
    obs <- observed_dt[source == src, .(program, source, metric_name = metric, observed = get(metric))]
    rnd <- random_dt[source == src, .(program, source, control_type, value = get(metric))]
    merge(rnd, obs, by = c("program", "source"), allow.cartesian = TRUE)
  }), fill = TRUE)
  g2 <- ggplot(dens_targets[is.finite(value)], aes(x = value, color = control_type, fill = control_type)) +
    geom_density(alpha = 0.20) +
    geom_vline(aes(xintercept = observed), linetype = 2, show.legend = FALSE) +
    facet_grid(source ~ program, scales = "free") +
    theme_bw(base_size = 10) +
    labs(title = "Observed metrics against matched-random null distributions", x = "Metric value", y = "Density")
  ggsave(file.path(args$outdir, "figures", "22_null_density_focus_metrics.pdf"), g2, width = 15, height = 8)

  heat_dt <- melt(merge(observed_dt[, .(program, source, n_sig_terms, best_neglog10q, focus_n_sig_terms, focus_best_neglog10q)],
                        z_dt, by = c("program", "source"), allow.cartesian = TRUE),
                  id.vars = c("program", "source", "control_type"),
                  measure.vars = patterns("_z$"), variable.name = "metric_idx", value.name = "zscore")
  z_names <- names(merge(observed_dt[, .(program, source, n_sig_terms, best_neglog10q, focus_n_sig_terms, focus_best_neglog10q)], z_dt, by = c("program", "source"), allow.cartesian = TRUE))[grep("_z$", names(merge(observed_dt[, .(program, source, n_sig_terms, best_neglog10q, focus_n_sig_terms, focus_best_neglog10q)], z_dt, by = c("program", "source"), allow.cartesian = TRUE)))]
  heat_dt[, metric := sub("_z$", "", z_names[metric_idx])]
  g3 <- ggplot(heat_dt, aes(x = metric, y = program, fill = zscore)) +
    geom_tile() +
    facet_grid(source ~ control_type, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Observed-vs-random z-scores", x = "Metric", y = "Program")
  ggsave(file.path(args$outdir, "figures", "23_metric_zscore_heatmap.pdf"), g3, width = 15, height = 8)

  wb <- createWorkbook()
  addWorksheet(wb, "input_summary")
  writeData(wb, "input_summary", input_dt)
  addWorksheet(wb, "observed_metrics")
  writeData(wb, "observed_metrics", observed_dt)
  addWorksheet(wb, "empirical_pvalues")
  writeData(wb, "empirical_pvalues", empirical_dt)
  addWorksheet(wb, "metric_zscores")
  writeData(wb, "metric_zscores", z_dt)
  addWorksheet(wb, "focus_summary")
  writeData(wb, "focus_summary", focus_summary)
  saveWorkbook(wb, file.path(args$outdir, "16_specificity_controls_workbook.xlsx"), overwrite = TRUE)

  lines <- c(
    "Step13A specificity / matched-random control summary",
    "",
    paste0("Program file: ", args$program_file),
    paste0("SynGO file: ", args$syngo_file),
    paste0("Expression reference: ", ifelse(is.null(expr_dt), "none", args$expression_ref_file)),
    paste0("Programs tested: ", paste(args$requested_programs, collapse = ", ")),
    paste0("n_perm: ", args$n_perm),
    paste0("q_cutoff: ", args$q_cutoff),
    "",
    "Observed metrics were computed on source-specific mapped genes.",
    "Empirical p-values compare observed metrics against size-matched and expression-matched random controls.",
    ""
  )
  for (prog in args$requested_programs) {
    lines <- c(lines, paste0("Program: ", prog))
    for (src in c("SynGO", "Reactome")) {
      obs <- observed_dt[program == prog & source == src]
      if (!nrow(obs)) next
      ep <- empirical_dt[program == prog & source == src]
      lines <- c(lines,
                 paste0("  Source: ", src,
                        " | mapped_size=", obs$mapped_size[1],
                        " | n_sig_terms=", obs$n_sig_terms[1],
                        " | best_neglog10q=", sprintf("%.3f", obs$best_neglog10q[1]),
                        " | focus_n_sig_terms=", ifelse(is.na(obs$focus_n_sig_terms[1]), "NA", obs$focus_n_sig_terms[1]),
                        " | focus_best_neglog10q=", ifelse(is.na(obs$focus_best_neglog10q[1]), "NA", sprintf("%.3f", obs$focus_best_neglog10q[1]))))
      for (ctype in c("size_matched", "expression_matched")) {
        e <- ep[control_type == ctype]
        if (!nrow(e)) next
        lines <- c(lines,
                   paste0("    ", ctype,
                          " | n_sig_terms_emp_p=", sprintf("%.4g", e$n_sig_terms_emp_p[1]),
                          " | best_neglog10q_emp_p=", sprintf("%.4g", e$best_neglog10q_emp_p[1]),
                          " | focus_n_sig_terms_emp_p=", ifelse(is.na(e$focus_n_sig_terms_emp_p[1]), "NA", sprintf("%.4g", e$focus_n_sig_terms_emp_p[1])),
                          " | focus_best_neglog10q_emp_p=", ifelse(is.na(e$focus_best_neglog10q_emp_p[1]), "NA", sprintf("%.4g", e$focus_best_neglog10q_emp_p[1]))))
      }
    }
    lines <- c(lines, "")
  }
  writeLines(lines, file.path(args$outdir, "15_specificity_summary.txt"))

  log_msg("Completed Step13A specificity / matched-random control analyses")
}

main()
