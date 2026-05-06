#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !all(is.na(x))) x else y

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    project_root = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex",
    gandal_rdata = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData",
    brainspan_rowmap = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/meta/12c_BrainSpan_rows_metadata_standardized.tsv",
    program_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/20_midPrenatal_core_programs.tsv",
    outdir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step03c_Gandal_alternative_scoring_sensitivity",
    methods = "meanlog2cpm,zmean,ssgsea",
    primary_programs = "SFARI_all,midPrenatal_SFARI_top20"
  )
  if (length(args) == 0) return(out)
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
  out
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

calc_auc <- function(score, is_case) {
  is_case <- as.logical(is_case)
  n1 <- sum(is_case)
  n0 <- sum(!is_case)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  U <- sum(ranks[is_case]) - n1 * (n1 + 1) / 2
  U / (n1 * n0)
}

empty_dt <- function(cols) {
  out <- data.table(matrix(ncol = length(cols), nrow = 0))
  setnames(out, cols)
  out
}

find_first_col <- function(nms, candidates) {
  low <- tolower(nms)
  cand_low <- tolower(candidates)
  idx <- match(cand_low, low)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NULL)
  nms[idx[1]]
}

clean_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "Unknown", "unknown", "Not Reported")] <- NA_character_
  x <- gsub(",", "", x, fixed = TRUE)
  x <- gsub("[^0-9eE+.-]", "", x)
  suppressWarnings(as.numeric(x))
}

clean_factor <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "Unknown", "unknown", "Not Reported")] <- NA_character_
  factor(x)
}

detect_covariates <- function(meta) {
  age_col <- find_first_col(names(meta), c("Age", "age", "AGE", "age_death", "Age.at.Death"))
  pmi_col <- find_first_col(names(meta), c("PMI", "pmi", "postmortem_interval", "PMIhrs", "PMI_hr"))
  rin_col <- find_first_col(names(meta), c("RIN", "rin", "RIN_value"))
  sex_col <- find_first_col(names(meta), c("Sex", "sex", "SEX", "gender", "Gender"))
  batch_col <- find_first_col(names(meta), c("batch", "Batch", "BATCH", "seq_batch", "library_batch", "BatchID"))

  inv <- rbindlist(list(
    data.table(model_term = "age_num", type = "numeric", non_missing = if (!is.null(age_col)) sum(!is.na(clean_numeric(meta[[age_col]]))) else 0L, n_levels = if (!is.null(age_col)) uniqueN(clean_numeric(meta[[age_col]]), na.rm = TRUE) else 0L),
    data.table(model_term = "pmi_num", type = "numeric", non_missing = if (!is.null(pmi_col)) sum(!is.na(clean_numeric(meta[[pmi_col]]))) else 0L, n_levels = if (!is.null(pmi_col)) uniqueN(clean_numeric(meta[[pmi_col]]), na.rm = TRUE) else 0L),
    data.table(model_term = "rin_num", type = "numeric", non_missing = if (!is.null(rin_col)) sum(!is.na(clean_numeric(meta[[rin_col]]))) else 0L, n_levels = if (!is.null(rin_col)) uniqueN(clean_numeric(meta[[rin_col]]), na.rm = TRUE) else 0L),
    data.table(model_term = "sex_factor", type = "factor", non_missing = if (!is.null(sex_col)) sum(!is.na(clean_factor(meta[[sex_col]]))) else 0L, n_levels = if (!is.null(sex_col)) nlevels(droplevels(clean_factor(meta[[sex_col]]))) else 0L),
    data.table(model_term = "batch_factor", type = "factor", non_missing = if (!is.null(batch_col)) sum(!is.na(clean_factor(meta[[batch_col]]))) else 0L, n_levels = if (!is.null(batch_col)) nlevels(droplevels(clean_factor(meta[[batch_col]]))) else 0L)
  ), use.names = TRUE)

  meta2 <- copy(meta)
  if (!is.null(age_col)) meta2[, age_num := clean_numeric(get(age_col))]
  if (!is.null(pmi_col)) meta2[, pmi_num := clean_numeric(get(pmi_col))]
  if (!is.null(rin_col)) meta2[, rin_num := clean_numeric(get(rin_col))]
  if (!is.null(sex_col)) meta2[, sex_factor := droplevels(clean_factor(get(sex_col)))]
  if (!is.null(batch_col)) meta2[, batch_factor := droplevels(clean_factor(get(batch_col)))]

  adj_terms <- character()
  if ("age_num" %in% names(meta2) && sum(!is.na(meta2$age_num)) >= 50 && uniqueN(meta2$age_num, na.rm = TRUE) >= 10) adj_terms <- c(adj_terms, "age_num")
  if ("pmi_num" %in% names(meta2) && sum(!is.na(meta2$pmi_num)) >= 50 && uniqueN(meta2$pmi_num, na.rm = TRUE) >= 10) adj_terms <- c(adj_terms, "pmi_num")
  if ("rin_num" %in% names(meta2) && sum(!is.na(meta2$rin_num)) >= 50 && uniqueN(meta2$rin_num, na.rm = TRUE) >= 10) adj_terms <- c(adj_terms, "rin_num")
  if ("sex_factor" %in% names(meta2) && nlevels(droplevels(meta2$sex_factor)) >= 2) adj_terms <- c(adj_terms, "sex_factor")
  if ("batch_factor" %in% names(meta2) && nlevels(droplevels(meta2$batch_factor)) >= 2) adj_terms <- c(adj_terms, "batch_factor")

  list(
    meta = meta2,
    inventory = inv,
    adj_terms = adj_terms,
    age_col = age_col, pmi_col = pmi_col, rin_col = rin_col, sex_col = sex_col, batch_col = batch_col
  )
}

score_meanlog2cpm <- function(expr, gene_sets) {
  out <- lapply(names(gene_sets), function(pn) {
    gs <- intersect(gene_sets[[pn]], rownames(expr))
    if (length(gs) == 0) return(rep(NA_real_, ncol(expr)))
    colMeans(expr[gs, , drop = FALSE], na.rm = TRUE)
  })
  out <- as.data.table(out)
  setnames(out, names(gene_sets))
  out
}

score_zmean <- function(expr, gene_sets) {
  row_means <- rowMeans(expr, na.rm = TRUE)
  row_sds <- apply(expr, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  zexpr <- (expr - row_means) / row_sds
  out <- lapply(names(gene_sets), function(pn) {
    gs <- intersect(gene_sets[[pn]], rownames(zexpr))
    if (length(gs) == 0) return(rep(NA_real_, ncol(zexpr)))
    colMeans(zexpr[gs, , drop = FALSE], na.rm = TRUE)
  })
  out <- as.data.table(out)
  setnames(out, names(gene_sets))
  out
}

score_ssgsea <- function(expr, gene_sets) {
  # GSVA-free rank-based single-sample enrichment proxy.
  # Returns a normalized mean-rank score in [0,1], robust to scale differences.
  expr <- as.matrix(expr)
  expr <- expr[apply(expr, 1, function(x) any(is.finite(x))), , drop = FALSE]
  G <- nrow(expr)
  if (G < 10) stop("Expression matrix has too few genes for rank-based scoring.")

  ranks <- apply(expr, 2, rank, ties.method = "average", na.last = "keep")
  if (is.null(dim(ranks))) ranks <- matrix(ranks, nrow = G, ncol = ncol(expr))
  rownames(ranks) <- rownames(expr)
  colnames(ranks) <- colnames(expr)

  out <- lapply(names(gene_sets), function(pn) {
    gs <- intersect(gene_sets[[pn]], rownames(ranks))
    k <- length(gs)
    if (k == 0) return(rep(NA_real_, ncol(ranks)))
    mean_rank <- colMeans(ranks[gs, , drop = FALSE], na.rm = TRUE)
    # Normalize so 0 ~ low-ranked set, 1 ~ high-ranked set.
    (mean_rank - 1) / (G - 1)
  })
  out <- as.data.table(out)
  setnames(out, names(gene_sets))
  for (pn in names(gene_sets)) suppressWarnings(set(out, j = pn, value = as.numeric(out[[pn]])))
  out
}

build_score_dt <- function(score_dt, meta_use, score_method) {
  score_dt <- as.data.table(copy(score_dt))
  # force all score columns numeric and avoid duplicate sample_id cbind issues
  for (nm in names(score_dt)) set(score_dt, j = nm, value = as.numeric(score_dt[[nm]]))
  score_dt[, sample_id := meta_use$sample_id]
  setcolorder(score_dt, c("sample_id", setdiff(names(score_dt), "sample_id")))
  dt <- merge(copy(meta_use), score_dt, by = "sample_id", all.x = TRUE, sort = FALSE)
  long <- melt(dt,
               id.vars = setdiff(names(dt), names(score_dt)[names(score_dt) != "sample_id"]),
               measure.vars = names(score_dt)[names(score_dt) != "sample_id"],
               variable.name = "program_name",
               value.name = "score",
               variable.factor = FALSE)
  long[, score := as.numeric(score)]
  long[, score_method := score_method]
  long[]
}

fit_models_for_subset <- function(dt_sub, adj_terms, mapped_n) {
  out <- list()
  dt_sub <- copy(dt_sub)
  dt_sub[, score := as.numeric(score)]
  dt_sub <- dt_sub[!is.na(score)]
  if (uniqueN(dt_sub$diagnosis) < 2 || sum(dt_sub$diagnosis == "ASD") < 2 || sum(dt_sub$diagnosis == "Control") < 2) return(NULL)

  model_specs <- list(
    M0 = c("diagnosis"),
    M1 = c("diagnosis", "region_factor"),
    M2 = c("diagnosis", "region_factor", adj_terms)
  )

  for (mn in names(model_specs)) {
    terms <- model_specs[[mn]]
    terms <- terms[terms %in% names(dt_sub)]
    keep_cols <- unique(c("score", terms))
    dd <- copy(dt_sub[, ..keep_cols])
    dd <- dd[complete.cases(dd)]
    dd[, score := as.numeric(score)]
    if (nrow(dd) < 30) next
    if (!"diagnosis" %in% names(dd) || uniqueN(dd$diagnosis) < 2) next
    if (!is.numeric(dd$score)) next
    form <- as.formula(paste("score ~", paste(terms, collapse = " + ")))
    fit <- lm(form, data = dd)
    co <- summary(fit)$coefficients
    if (!"diagnosisASD" %in% rownames(co)) next
    ci <- suppressWarnings(confint(fit, parm = "diagnosisASD", level = 0.95))
    wt_p <- tryCatch(wilcox.test(score ~ diagnosis, data = dd, exact = FALSE)$p.value, error = function(e) NA_real_)
    out[[mn]] <- data.table(
      model = mn,
      formula = paste(deparse(form, width.cutoff = 500L), collapse = ""),
      beta_asd = unname(co["diagnosisASD", "Estimate"]),
      ci_lo = unname(ci[1]),
      ci_hi = unname(ci[2]),
      p_asd = unname(co["diagnosisASD", "Pr(>|t|)"]),
      aic = AIC(fit),
      n_samples_model = nrow(dd),
      n_genes_mapped_gandal = mapped_n,
      n_samples = nrow(dt_sub),
      n_asd = sum(dt_sub$diagnosis == "ASD"),
      n_control = sum(dt_sub$diagnosis == "Control"),
      n_regions = uniqueN(dt_sub$region),
      mean_score_control = mean(dd[diagnosis == "Control", score], na.rm = TRUE),
      mean_score_asd = mean(dd[diagnosis == "ASD", score], na.rm = TRUE),
      delta_asd_minus_control = mean(dd[diagnosis == "ASD", score], na.rm = TRUE) - mean(dd[diagnosis == "Control", score], na.rm = TRUE),
      auc_asd_higher = calc_auc(dd$score, dd$diagnosis == "ASD"),
      wilcox_p = wt_p
    )
  }
  if (length(out) == 0) return(NULL)
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

make_plot1 <- function(primary_dt, outfile_base) {
  p <- ggplot(primary_dt, aes(x = model, y = beta_asd, color = score_method, group = score_method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 2) +
    geom_line(linewidth = 0.5) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.12) +
    facet_wrap(~ program_name, scales = "free_y") +
    labs(x = NULL, y = "ASD effect (beta)", color = "Scoring") +
    theme_bw(base_size = 11)
  ggsave(paste0(outfile_base, ".pdf"), p, width = 8.5, height = 4.8)
  ggsave(paste0(outfile_base, ".png"), p, width = 8.5, height = 4.8, dpi = 300)
}

make_plot2 <- function(primary_cmp, outfile_base) {
  p <- ggplot(primary_cmp[model != "M0"], aes(x = score_method, y = delta_beta_vs_M0, fill = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    facet_wrap(~ program_name, scales = "free_y") +
    labs(x = "Scoring method", y = expression(Delta*" beta vs M0"), fill = "Model") +
    theme_bw(base_size = 11)
  ggsave(paste0(outfile_base, ".pdf"), p, width = 8.5, height = 4.8)
  ggsave(paste0(outfile_base, ".png"), p, width = 8.5, height = 4.8, dpi = 300)
}

main <- function() {
  opt <- parse_args()

  PROJECT_ROOT <- opt$project_root
  OUTDIR <- opt$outdir
  METADIR <- file.path(OUTDIR, "meta")
  TABDIR <- file.path(OUTDIR, "tables")
  PLOTDIR <- file.path(OUTDIR, "plots")
  LOGDIR <- file.path(OUTDIR, "logs")

  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(PLOTDIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(LOGDIR, "Step03c_Gandal_alternative_scoring_sensitivity_v1.log")
  log_msg <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
    cat(msg, "\n")
    cat(msg, "\n", file = log_file, append = TRUE)
  }

  log_msg("PROJECT_ROOT: ", PROJECT_ROOT)
  log_msg("OUTDIR: ", OUTDIR)
  log_msg("GANDAL_RDATA: ", opt$gandal_rdata)
  log_msg("BRAINSPAN_ROWMAP: ", opt$brainspan_rowmap)
  log_msg("PROGRAM_FILE: ", opt$program_file)
  log_msg("METHODS: ", opt$methods)

  methods <- unique(trimws(strsplit(opt$methods, ",", fixed = TRUE)[[1]]))
  primary_programs <- unique(trimws(strsplit(opt$primary_programs, ",", fixed = TRUE)[[1]]))
  allowed_methods <- c("meanlog2cpm", "zmean", "ssgsea")
  if (!all(methods %in% allowed_methods)) stop("Unsupported methods. Allowed: ", paste(allowed_methods, collapse = ", "))

  if (!file.exists(opt$program_file)) stop("Missing program file: ", opt$program_file)
  if (!file.exists(opt$brainspan_rowmap)) stop("Missing BrainSpan row mapping file: ", opt$brainspan_rowmap)
  if (!file.exists(opt$gandal_rdata)) stop("Missing Gandal RData: ", opt$gandal_rdata)

  program_tbl <- fread(opt$program_file)
  program_tbl[, gene_symbol := toupper(trimws(gene_symbol))]
  program_tbl <- unique(program_tbl[gene_symbol != "" & !is.na(gene_symbol)])
  safe_fwrite(program_tbl[, .N, by = program_name], file.path(METADIR, "00_program_summary.tsv"))

  rowmap <- fread(opt$brainspan_rowmap)
  needed_cols <- c("gene_symbol", "ensembl_gene_id")
  if (!all(needed_cols %in% names(rowmap))) stop("BrainSpan row mapping file must contain columns: gene_symbol, ensembl_gene_id")
  rowmap[, gene_symbol := toupper(trimws(gene_symbol))]
  rowmap[, ensembl_gene_id := toupper(trimws(as.character(ensembl_gene_id)))]
  rowmap[, ensembl_stable := sub("\\..*$", "", ensembl_gene_id)]
  rowmap <- unique(rowmap[gene_symbol != "" & !is.na(gene_symbol) & ensembl_stable != "" & !is.na(ensembl_stable),
                          .(gene_symbol, ensembl_gene_id, ensembl_stable)])
  safe_fwrite(rowmap, file.path(METADIR, "01_BrainSpan_geneSymbol_to_ensembl.tsv"))

  log_msg("Loading Gandal RData ...")
  e <- new.env(parent = emptyenv())
  load(opt$gandal_rdata, envir = e)
  if (!all(c("datExpr", "datMeta") %in% ls(e))) stop("Expected datExpr and datMeta in Gandal RData.")
  datExpr <- get("datExpr", envir = e)
  datMeta <- as.data.table(get("datMeta", envir = e))
  if (!is.matrix(datExpr)) datExpr <- as.matrix(datExpr)
  mode(datExpr) <- "numeric"
  safe_fwrite(data.table(
    object_name = c("datExpr", "datMeta"),
    class = c(paste(class(datExpr), collapse = ";"), paste(class(datMeta), collapse = ";")),
    nrow = c(nrow(datExpr), nrow(datMeta)),
    ncol = c(ncol(datExpr), ncol(datMeta))
  ), file.path(METADIR, "02_RData_object_inventory.tsv"))

  required_meta_cols <- c("sample_id", "Diagnosis", "region")
  if (!all(required_meta_cols %in% names(datMeta))) stop("datMeta must contain columns: sample_id, Diagnosis, region")
  meta <- copy(datMeta)
  meta[, sample_id := as.character(sample_id)]
  meta[, Diagnosis := as.character(Diagnosis)]
  meta[, region := as.character(region)]
  meta[, diagnosis := fifelse(Diagnosis == "ASD", "ASD", fifelse(Diagnosis == "CTL", "Control", NA_character_))]
  meta <- meta[!is.na(diagnosis)]
  meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
  meta[, is_asd := diagnosis == "ASD"]
  meta[, region_factor := factor(region)]

  cov <- detect_covariates(meta)
  meta_use <- cov$meta
  safe_fwrite(cov$inventory, file.path(METADIR, "03_metadata_covariate_inventory.tsv"))
  safe_fwrite(meta_use, file.path(METADIR, "04_metadata_cleaned.tsv"))
  term_dt <- rbindlist(list(
    data.table(model = "M0", term = "diagnosis"),
    data.table(model = "M1", term = c("diagnosis", "region_factor")),
    data.table(model = "M2", term = c("diagnosis", "region_factor", cov$adj_terms))
  ), use.names = TRUE)
  safe_fwrite(unique(term_dt), file.path(METADIR, "05_adjusted_model_terms.tsv"))
  log_msg("Adjusted model additional terms: ", paste(cov$adj_terms, collapse = ", "))

  common_samples <- intersect(colnames(datExpr), meta_use$sample_id)
  if (length(common_samples) < 20) stop("Too few overlapping samples between datExpr and datMeta")
  meta_use <- meta_use[match(common_samples, sample_id)]
  expr_use <- datExpr[, common_samples, drop = FALSE]

  expr_gene_raw <- rownames(expr_use)
  expr_ensembl_stable <- sub("_.*$", "", expr_gene_raw)
  expr_ensembl_stable <- sub("\\..*$", "", expr_ensembl_stable)
  expr_ensembl_stable <- toupper(trimws(expr_ensembl_stable))
  gandal_gene_map <- data.table(gandal_rowname = expr_gene_raw, ensembl_stable = expr_ensembl_stable)
  gandal_gene_map <- merge(gandal_gene_map, rowmap[, .(ensembl_stable, gene_symbol)], by = "ensembl_stable", all.x = TRUE)
  safe_fwrite(gandal_gene_map, file.path(METADIR, "06_Gandal_gene_mapping.tsv"))
  mapped_rows <- !is.na(gandal_gene_map$gene_symbol)
  log_msg("Mapped Gandal rows to gene symbols: ", sum(mapped_rows), " / ", nrow(gandal_gene_map))
  if (sum(mapped_rows) < 10000) stop("Too few Gandal rows mapped to gene symbols.")

  expr_mapped <- expr_use[mapped_rows, , drop = FALSE]
  mapped_symbols <- gandal_gene_map$gene_symbol[mapped_rows]
  expr_dt <- as.data.table(expr_mapped)
  expr_dt[, gene_symbol := mapped_symbols]
  expr_collapsed <- expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol]
  expr_gene_symbol <- as.matrix(expr_collapsed[, -1])
  rownames(expr_gene_symbol) <- expr_collapsed$gene_symbol
  mode(expr_gene_symbol) <- "numeric"
  safe_fwrite(data.table(mapped_rows = sum(mapped_rows), unique_gene_symbols = nrow(expr_gene_symbol)), file.path(METADIR, "07_expression_mapping_summary.tsv"))

  gene_sets <- split(program_tbl$gene_symbol, program_tbl$program_name)
  gene_sets <- lapply(gene_sets, unique)
  program_diag <- rbindlist(lapply(sort(names(gene_sets)), function(pn) {
    gs <- unique(gene_sets[[pn]])
    data.table(program_name = pn, n_genes_program = length(gs), n_genes_in_gandal = length(intersect(gs, rownames(expr_gene_symbol))))
  }))
  safe_fwrite(program_diag, file.path(METADIR, "08_program_mapping_summary.tsv"))

  score_long_list <- list()
  model_res_list <- list()

  for (m in methods) {
    log_msg("Scoring method: ", m)
    if (m == "meanlog2cpm") {
      score_dt <- score_meanlog2cpm(expr_gene_symbol, gene_sets)
    } else if (m == "zmean") {
      score_dt <- score_zmean(expr_gene_symbol, gene_sets)
    } else if (m == "ssgsea") {
      score_dt <- score_ssgsea(expr_gene_symbol, gene_sets)
    } else {
      stop("Unsupported method: ", m)
    }
    score_dt <- as.data.table(score_dt)
    setnames(score_dt, names(gene_sets))
    long <- build_score_dt(score_dt, meta_use, m)
    score_long_list[[m]] <- long

    for (pn in sort(names(gene_sets))) {
      mapped_n <- length(intersect(gene_sets[[pn]], rownames(expr_gene_symbol)))
      if (mapped_n < 5) next
      sub_dt <- long[program_name == pn & !is.na(score)]
      res <- fit_models_for_subset(sub_dt, cov$adj_terms, mapped_n)
      if (is.null(res)) next
      res[, `:=`(program_name = pn, score_method = m)]
      model_res_list[[paste(m, pn, sep = "__")]] <- res
    }
  }

  score_cols <- c("sample_id","Diagnosis","region","diagnosis","is_asd","region_factor","program_name","score","score_method")
  model_cols <- c("program_name","score_method","model","formula","beta_asd","ci_lo","ci_hi","p_asd","aic","n_samples_model",
                  "n_genes_mapped_gandal","n_samples","n_asd","n_control","n_regions","mean_score_control","mean_score_asd",
                  "delta_asd_minus_control","auc_asd_higher","wilcox_p")
  scores_long <- if (length(score_long_list) > 0) rbindlist(score_long_list, use.names = TRUE, fill = TRUE) else empty_dt(score_cols)
  model_dt <- if (length(model_res_list) > 0) rbindlist(model_res_list, use.names = TRUE, fill = TRUE) else empty_dt(model_cols)

  if (nrow(model_dt) > 0) {
    model_dt <- unique(model_dt, by = c("program_name","score_method","model","formula","beta_asd","ci_lo","ci_hi","p_asd","aic","n_samples_model"))
    model_dt[, fdr_p_asd := p.adjust(p_asd, method = "BH"), by = score_method]
    setcolorder(model_dt, c("program_name","score_method","model","formula","beta_asd","ci_lo","ci_hi","p_asd","fdr_p_asd",
                            "aic","n_samples_model","n_genes_mapped_gandal","n_samples","n_asd","n_control","n_regions",
                            "mean_score_control","mean_score_asd","delta_asd_minus_control","auc_asd_higher","wilcox_p"))
  }
  safe_fwrite(scores_long, file.path(TABDIR, "10_Gandal_alt_scoring_scores_long.tsv"))
  safe_fwrite(model_dt, file.path(TABDIR, "11_Gandal_alt_scoring_model_results.tsv"))

  cmp_list <- list()
  if (nrow(model_dt) > 0) {
    for (m in unique(model_dt$score_method)) {
      for (pn in unique(model_dt$program_name)) {
        dd <- copy(model_dt[score_method == m & program_name == pn])
        if (!"M0" %in% dd$model) next
        ref <- dd[model == "M0"][1]
        dd[, ref_model := "M0"]
        dd[, delta_beta_vs_M0 := beta_asd - ref$beta_asd]
        dd[, delta_abs_beta_vs_M0 := abs(beta_asd) - abs(ref$beta_asd)]
        dd[, delta_log10p_vs_M0 := (-log10(p_asd)) - (-log10(ref$p_asd))]
        dd[, delta_aic_vs_M0 := aic - ref$aic]
        cmp_list[[paste(m, pn, sep = "__")]] <- dd
      }
    }
  }
  cmp_dt <- if (length(cmp_list) > 0) rbindlist(cmp_list, use.names = TRUE, fill = TRUE) else empty_dt(c(names(model_cols), "fdr_p_asd","ref_model","delta_beta_vs_M0","delta_abs_beta_vs_M0","delta_log10p_vs_M0","delta_aic_vs_M0"))
  safe_fwrite(cmp_dt, file.path(TABDIR, "12_Gandal_alt_scoring_model_comparison.tsv"))

  primary_dt <- model_dt[program_name %in% primary_programs]
  safe_fwrite(primary_dt, file.path(TABDIR, "13_Gandal_alt_scoring_primary_program_model_results.tsv"))

  if (nrow(primary_dt) > 0) {
    make_plot1(primary_dt, file.path(PLOTDIR, "01_primary_program_alt_scoring_M0_vs_M2"))
  }
  if (nrow(cmp_dt[program_name %in% primary_programs]) > 0) {
    make_plot2(cmp_dt[program_name %in% primary_programs], file.path(PLOTDIR, "02_primary_program_alt_scoring_effect_shift"))
  }

  run_summary <- rbindlist(list(
    data.table(item = "n_samples_used", value = nrow(meta_use)),
    data.table(item = "n_asd", value = sum(meta_use$is_asd)),
    data.table(item = "n_control", value = sum(!meta_use$is_asd)),
    data.table(item = "n_regions", value = uniqueN(meta_use$region)),
    data.table(item = "mapped_rows", value = sum(mapped_rows)),
    data.table(item = "unique_gene_symbols", value = nrow(expr_gene_symbol)),
    data.table(item = "age_column", value = cov$age_col %||% NA_character_),
    data.table(item = "pmi_column", value = cov$pmi_col %||% NA_character_),
    data.table(item = "rin_column", value = cov$rin_col %||% NA_character_),
    data.table(item = "sex_column", value = cov$sex_col %||% NA_character_),
    data.table(item = "batch_column", value = cov$batch_col %||% NA_character_),
    data.table(item = "adjusted_terms", value = paste(cov$adj_terms, collapse = ",")),
    data.table(item = "methods", value = paste(methods, collapse = ",")),
    data.table(item = "n_programs_modeled", value = uniqueN(model_dt$program_name)),
    data.table(item = "n_rows_model_results", value = nrow(model_dt))
  ), use.names = TRUE)
  safe_fwrite(run_summary, file.path(METADIR, "14_run_summary.tsv"))
  log_msg("Step03c completed successfully.")
}

main()
