#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------
# CLI parsing
# -----------------------------
parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))

PROJECT_ROOT <- if (!is.null(argv$project_root)) argv$project_root else "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
GANDAL_RDATA <- if (!is.null(argv$gandal_rdata)) argv$gandal_rdata else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"
BRAINSPAN_ROWMAP <- if (!is.null(argv$brainspan_rowmap)) argv$brainspan_rowmap else file.path(PROJECT_ROOT, "step01_prepare/meta/12c_BrainSpan_rows_metadata_standardized.tsv")
PROGRAM_FILE <- if (!is.null(argv$program_file)) argv$program_file else file.path(PROJECT_ROOT, "step01_prepare/tables/20_midPrenatal_core_programs.tsv")
OUTDIR <- if (!is.null(argv$outdir)) argv$outdir else file.path(PROJECT_ROOT, "step03b_Gandal_covariate_sensitivity")

METADIR <- file.path(OUTDIR, "meta")
TABDIR  <- file.path(OUTDIR, "tables")
PLOTDIR <- file.path(OUTDIR, "plots")
LOGDIR  <- file.path(OUTDIR, "logs")
for (d in c(OUTDIR, METADIR, TABDIR, PLOTDIR, LOGDIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOGFILE <- file.path(LOGDIR, "Step03b_Gandal_covariate_sensitivity_v5.log")
if (file.exists(LOGFILE)) file.remove(LOGFILE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOGFILE, append = TRUE)
}

safe_fwrite <- function(x, file) fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")

log_msg("PROJECT_ROOT:  ", PROJECT_ROOT)
log_msg("OUTDIR:  ", OUTDIR)
log_msg("GANDAL_RDATA:  ", GANDAL_RDATA)
log_msg("BRAINSPAN_ROWMAP:  ", BRAINSPAN_ROWMAP)
log_msg("PROGRAM_FILE:  ", PROGRAM_FILE)

# -----------------------------
# Helpers
# -----------------------------
calc_auc <- function(score, is_case) {
  is_case <- as.logical(is_case)
  n1 <- sum(is_case)
  n0 <- sum(!is_case)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  U <- sum(ranks[is_case]) - n1 * (n1 + 1) / 2
  U / (n1 * n0)
}

clean_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

clean_ensg <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_.*$", "", x)
  x <- sub("\\..*$", "", x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

standardize_dx <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- ifelse(y %in% c("ASD", "AUTISM", "CASE"), "ASD",
                ifelse(y %in% c("CTL", "CTRL", "CONTROL"), "Control", NA_character_))
  out
}

extract_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x <- gsub("[^0-9eE+.-]", "", x)
  suppressWarnings(as.numeric(x))
}

choose_best_numeric <- function(dt, candidates) {
  hits <- intersect(candidates, names(dt))
  if (length(hits) == 0) return(NULL)
  stats <- rbindlist(lapply(hits, function(nm) {
    v <- extract_numeric(dt[[nm]])
    data.table(column = nm,
               non_missing = sum(!is.na(v)),
               unique_non_missing = uniqueN(v[!is.na(v)]))
  }))
  stats <- stats[non_missing > 0 & unique_non_missing > 1]
  if (nrow(stats) == 0) return(NULL)
  setorder(stats, -non_missing, -unique_non_missing)
  stats$column[1]
}

choose_best_factor <- function(dt, candidates, min_non_missing = 10L, min_levels = 2L, max_levels = 50L) {
  hits <- intersect(candidates, names(dt))
  if (length(hits) == 0) return(NULL)
  stats <- rbindlist(lapply(hits, function(nm) {
    v <- as.character(dt[[nm]])
    nn <- sum(!is.na(v) & trimws(v) != "")
    lv <- uniqueN(v[!is.na(v) & trimws(v) != ""])
    data.table(column = nm, non_missing = nn, n_levels = lv)
  }))
  stats <- stats[non_missing >= min_non_missing & n_levels >= min_levels & n_levels <= max_levels]
  if (nrow(stats) == 0) return(NULL)
  setorder(stats, -non_missing, n_levels)
  stats$column[1]
}

confint_or_na <- function(fit, term) {
  out <- tryCatch(confint(fit, parm = term), error = function(e) c(NA_real_, NA_real_))
  as.numeric(out)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

empty_dt <- function(cols) {
  out <- data.table(matrix(ncol = length(cols), nrow = 0))
  setnames(out, cols)
  out
}

# -----------------------------
# Inputs
# -----------------------------
if (!file.exists(GANDAL_RDATA)) stop("Missing GANDAL_RDATA: ", GANDAL_RDATA)
if (!file.exists(BRAINSPAN_ROWMAP)) stop("Missing BRAINSPAN_ROWMAP: ", BRAINSPAN_ROWMAP)
if (!file.exists(PROGRAM_FILE)) stop("Missing PROGRAM_FILE: ", PROGRAM_FILE)

program_tbl <- fread(PROGRAM_FILE)
if (!all(c("program_name", "gene_symbol") %in% names(program_tbl))) {
  stop("PROGRAM_FILE must contain program_name and gene_symbol columns.")
}
program_tbl[, gene_symbol := clean_symbol(gene_symbol)]
program_tbl <- unique(program_tbl[!is.na(gene_symbol) & gene_symbol != ""])

rowmap <- fread(BRAINSPAN_ROWMAP)
needed_cols <- c("gene_symbol", "ensembl_gene_id")
if (!all(needed_cols %in% names(rowmap))) {
  stop("BRAINSPAN_ROWMAP must contain: ", paste(needed_cols, collapse = ", "))
}
rowmap[, gene_symbol := clean_symbol(gene_symbol)]
rowmap[, ensembl_stable := clean_ensg(ensembl_gene_id)]
rowmap <- unique(rowmap[!is.na(gene_symbol) & !is.na(ensembl_stable), .(gene_symbol, ensembl_stable)])

# -----------------------------
# Load Gandal
# -----------------------------
log_msg("Loading Gandal RData ... ")
e <- new.env(parent = emptyenv())
load(GANDAL_RDATA, envir = e)
objs <- ls(e)

if (!("datExpr" %in% objs) || !("datMeta" %in% objs)) {
  stop("Expected datExpr and datMeta in Gandal RData.")
}

datExpr <- get("datExpr", envir = e)
datMeta <- as.data.table(get("datMeta", envir = e))
if (!is.matrix(datExpr)) datExpr <- as.matrix(datExpr)
mode(datExpr) <- "numeric"

# -----------------------------
# Harmonize metadata
# -----------------------------
required_meta_cols <- c("sample_id", "Diagnosis", "region")
if (!all(required_meta_cols %in% names(datMeta))) {
  stop("datMeta must contain columns: ", paste(required_meta_cols, collapse = ", "))
}

meta <- copy(datMeta)
meta[, sample_id := as.character(sample_id)]
meta[, diagnosis := standardize_dx(Diagnosis)]
meta <- meta[!is.na(diagnosis)]
meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
meta[, is_asd := diagnosis == "ASD"]
meta[, region_factor := factor(as.character(region))]

# Identify additional covariates
age_col <- choose_best_numeric(meta, c("age", "Age", "AGE", "age_death", "AOD", "AgeDeath", "Age_at_Death"))
pmi_col <- choose_best_numeric(meta, c("PMI", "pmi", "post_mortem_interval", "PostmortemInterval"))
rin_col <- choose_best_numeric(meta, c("RIN", "rin", "Rin"))
sex_col <- choose_best_factor(meta, c("sex", "Sex", "SEX", "gender", "Gender"), min_levels = 2L, max_levels = 5L)
batch_col <- choose_best_factor(meta, c("batch", "Batch", "BATCH", "seq_batch", "SeqBatch", "plate", "library_batch"), min_levels = 2L, max_levels = 50L)

inv_dt <- data.table(
  covariate_type = c("age", "PMI", "RIN", "sex", "batch"),
  selected_column = c(age_col %||% NA_character_, pmi_col %||% NA_character_, rin_col %||% NA_character_, sex_col %||% NA_character_, batch_col %||% NA_character_)
)

if (!is.null(age_col)) meta[, age_num := extract_numeric(get(age_col))]
if (!is.null(pmi_col)) meta[, pmi_num := extract_numeric(get(pmi_col))]
if (!is.null(rin_col)) meta[, rin_num := extract_numeric(get(rin_col))]
if (!is.null(sex_col)) meta[, sex_factor := factor(as.character(get(sex_col)))]
if (!is.null(batch_col)) meta[, batch_factor := factor(as.character(get(batch_col)))]

cov_inventory <- rbindlist(lapply(c("age_num", "pmi_num", "rin_num", "sex_factor", "batch_factor"), function(nm) {
  if (!(nm %in% names(meta))) return(data.table(model_term = nm, type = NA_character_, non_missing = 0L, n_levels = NA_integer_))
  if (is.numeric(meta[[nm]])) {
    data.table(model_term = nm, type = "numeric", non_missing = sum(!is.na(meta[[nm]])), n_levels = uniqueN(meta[[nm]][!is.na(meta[[nm]])]))
  } else {
    data.table(model_term = nm, type = "factor", non_missing = sum(!is.na(meta[[nm]])), n_levels = nlevels(meta[[nm]]))
  }
}), fill = TRUE)
safe_fwrite(cov_inventory, file.path(METADIR, "03_metadata_covariate_inventory.tsv"))
safe_fwrite(meta, file.path(METADIR, "04_metadata_cleaned.tsv"))

# Determine usable covariates for adjusted model
adj_terms <- c()
if ("age_num" %in% names(meta) && sum(!is.na(meta$age_num)) >= 20 && uniqueN(meta$age_num[!is.na(meta$age_num)]) > 1) adj_terms <- c(adj_terms, "age_num")
if ("pmi_num" %in% names(meta) && sum(!is.na(meta$pmi_num)) >= 20 && uniqueN(meta$pmi_num[!is.na(meta$pmi_num)]) > 1) adj_terms <- c(adj_terms, "pmi_num")
if ("rin_num" %in% names(meta) && sum(!is.na(meta$rin_num)) >= 20 && uniqueN(meta$rin_num[!is.na(meta$rin_num)]) > 1) adj_terms <- c(adj_terms, "rin_num")
if ("sex_factor" %in% names(meta) && nlevels(droplevels(meta$sex_factor)) >= 2) adj_terms <- c(adj_terms, "sex_factor")
if ("batch_factor" %in% names(meta) && nlevels(droplevels(meta$batch_factor)) >= 2) adj_terms <- c(adj_terms, "batch_factor")

term_dt <- rbindlist(list(
  data.table(model = "M0", term = "diagnosis"),
  data.table(model = "M1", term = c("diagnosis", "region_factor")),
  data.table(model = "M2", term = c("diagnosis", "region_factor", adj_terms))
), use.names = TRUE, fill = TRUE)
term_dt <- unique(term_dt)
safe_fwrite(term_dt, file.path(METADIR, "05_adjusted_model_terms.tsv"))
log_msg("Adjusted model additional terms:  ", paste(adj_terms, collapse = ", "))

# -----------------------------
# Expression harmonization
# -----------------------------
if (is.null(rownames(datExpr)) || is.null(colnames(datExpr))) stop("datExpr must have rownames and colnames")
common_samples <- intersect(colnames(datExpr), meta$sample_id)
if (length(common_samples) < 20) stop("Too few overlapping samples: ", length(common_samples))
meta_use <- meta[match(common_samples, sample_id)]
expr_use <- datExpr[, common_samples, drop = FALSE]

expr_ensg <- clean_ensg(rownames(expr_use))
gandal_gene_map <- data.table(gandal_rowname = rownames(expr_use), ensembl_stable = expr_ensg)
gandal_gene_map <- merge(gandal_gene_map, rowmap, by = "ensembl_stable", all.x = TRUE)
safe_fwrite(gandal_gene_map, file.path(METADIR, "06_Gandal_gene_mapping.tsv"))

mapped_rows <- !is.na(gandal_gene_map$gene_symbol)
log_msg("Mapped Gandal rows to gene symbols:  ", sum(mapped_rows), "  /  ", nrow(gandal_gene_map))
if (sum(mapped_rows) < 5000) stop("Too few mapped rows after gene symbol mapping.")

expr_mapped <- expr_use[mapped_rows, , drop = FALSE]
mapped_symbols <- gandal_gene_map$gene_symbol[mapped_rows]
expr_dt <- as.data.table(expr_mapped)
expr_dt[, gene_symbol := mapped_symbols]
expr_collapsed <- expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol]
expr_gene_symbol <- as.matrix(expr_collapsed[, -1])
rownames(expr_gene_symbol) <- expr_collapsed$gene_symbol
mode(expr_gene_symbol) <- "numeric"

# -----------------------------
# Program mapping diagnostics
# -----------------------------
program_diag <- rbindlist(lapply(sort(unique(program_tbl$program_name)), function(pn) {
  gs <- unique(program_tbl[program_name == pn, gene_symbol])
  data.table(
    program_name = pn,
    n_genes_program = length(gs),
    n_genes_in_rowmap = length(intersect(gs, rowmap$gene_symbol)),
    n_genes_in_gandal = length(intersect(gs, rownames(expr_gene_symbol)))
  )
}))
safe_fwrite(program_diag, file.path(METADIR, "08_program_mapping_summary.tsv"))

# -----------------------------
# Fit models
# -----------------------------
score_list <- list()
model_results <- list()
comparison_results <- list()

fit_model <- function(dd, formula_str) {
  fit <- lm(as.formula(formula_str), data = dd)
  sm <- summary(fit)$coefficients
  if (!("diagnosisASD" %in% rownames(sm))) return(NULL)
  ci <- confint_or_na(fit, "diagnosisASD")
  list(
    beta = unname(sm["diagnosisASD", "Estimate"]),
    p = unname(sm["diagnosisASD", "Pr(>|t|)"]),
    ci_lo = ci[1],
    ci_hi = ci[2],
    aic = AIC(fit)
  )
}

for (pn in sort(unique(program_tbl$program_name))) {
  genes <- intersect(program_tbl[program_name == pn, gene_symbol], rownames(expr_gene_symbol))
  if (length(genes) < 5) next

  score <- colMeans(expr_gene_symbol[genes, , drop = FALSE], na.rm = TRUE)
  dd <- data.table(sample_id = common_samples, score = as.numeric(score), program_name = pn)
  dd <- cbind(dd, meta_use[, .(diagnosis, is_asd, region, region_factor, age_num = if ("age_num" %in% names(meta_use)) age_num else NA_real_,
                               pmi_num = if ("pmi_num" %in% names(meta_use)) pmi_num else NA_real_,
                               rin_num = if ("rin_num" %in% names(meta_use)) rin_num else NA_real_,
                               sex_factor = if ("sex_factor" %in% names(meta_use)) sex_factor else factor(NA),
                               batch_factor = if ("batch_factor" %in% names(meta_use)) batch_factor else factor(NA))])
  score_list[[pn]] <- dd

  base_stats <- data.table(
    program_name = pn,
    n_genes_mapped_gandal = length(genes),
    n_samples = nrow(dd),
    n_asd = sum(dd$is_asd),
    n_control = sum(!dd$is_asd),
    n_regions = uniqueN(dd$region),
    mean_score_control = dd[diagnosis == "Control", mean(score, na.rm = TRUE)],
    mean_score_asd = dd[diagnosis == "ASD", mean(score, na.rm = TRUE)],
    delta_asd_minus_control = dd[diagnosis == "ASD", mean(score, na.rm = TRUE)] - dd[diagnosis == "Control", mean(score, na.rm = TRUE)],
    auc_asd_higher = calc_auc(dd$score, dd$is_asd),
    wilcox_p = tryCatch(wilcox.test(score ~ diagnosis, data = dd, exact = FALSE)$p.value, error = function(e) NA_real_)
  )

  model_rows <- list()

  # M0
  m0 <- fit_model(dd, "score ~ diagnosis")
  if (!is.null(m0)) model_rows[["M0"]] <- data.table(model = "M0", beta_asd = m0$beta, ci_lo = m0$ci_lo, ci_hi = m0$ci_hi, p_asd = m0$p, aic = m0$aic)

  # M1
  if (uniqueN(dd$region_factor) >= 2) {
    m1 <- fit_model(dd, "score ~ diagnosis + region_factor")
    if (!is.null(m1)) model_rows[["M1"]] <- data.table(model = "M1", beta_asd = m1$beta, ci_lo = m1$ci_lo, ci_hi = m1$ci_hi, p_asd = m1$p, aic = m1$aic)
  }

  # M2
  if (length(adj_terms) > 0) {
    terms <- c("diagnosis", "region_factor", adj_terms)
    terms <- unique(terms)
    keep <- rep(TRUE, nrow(dd))
    for (tm in adj_terms) {
      keep <- keep & !is.na(dd[[tm]])
    }
    dd2 <- dd[keep]
    # drop unused levels
    if ("sex_factor" %in% names(dd2)) dd2[, sex_factor := droplevels(sex_factor)]
    if ("batch_factor" %in% names(dd2)) dd2[, batch_factor := droplevels(batch_factor)]
    if ("region_factor" %in% names(dd2)) dd2[, region_factor := droplevels(region_factor)]
    if (nrow(dd2) >= 20 && uniqueN(dd2$diagnosis) == 2) {
      formula_m2 <- paste("score ~", paste(terms, collapse = " + "))
      m2 <- tryCatch(fit_model(dd2, formula_m2), error = function(e) NULL)
      if (!is.null(m2)) model_rows[["M2"]] <- data.table(model = "M2", beta_asd = m2$beta, ci_lo = m2$ci_lo, ci_hi = m2$ci_hi, p_asd = m2$p, aic = m2$aic, n_samples_model = nrow(dd2))
    }
  }

  if (length(model_rows) > 0) {
    model_prog <- rbindlist(model_rows, fill = TRUE)
    # add only non-duplicated base fields
    extra_cols <- setdiff(names(base_stats), c("program_name"))
    for (nm in extra_cols) model_prog[[nm]] <- base_stats[[nm]][1]
    model_prog[, program_name := pn]
    setcolorder(model_prog, c("program_name", "model", setdiff(names(model_prog), c("program_name", "model"))))
    model_results[[pn]] <- model_prog

    # comparisons against M0 when available
    if ("M0" %in% model_prog$model) {
      ref_beta <- model_prog[model == "M0", beta_asd][1]
      ref_p <- model_prog[model == "M0", p_asd][1]
      ref_aic <- model_prog[model == "M0", aic][1]
      cmp <- copy(model_prog)
      cmp[, ref_model := "M0"]
      cmp[, delta_beta_vs_M0 := beta_asd - ref_beta]
      cmp[, delta_abs_beta_vs_M0 := abs(beta_asd) - abs(ref_beta)]
      cmp[, delta_log10p_vs_M0 := -log10(p_asd) + log10(ref_p)]
      cmp[, delta_aic_vs_M0 := aic - ref_aic]
      comparison_results[[pn]] <- cmp
    }
  }
}

scores_long <- if (length(score_list) > 0) rbindlist(score_list, fill = TRUE) else empty_dt(c("sample_id", "score", "program_name", "diagnosis"))
model_dt <- if (length(model_results) > 0) rbindlist(model_results, fill = TRUE) else empty_dt(c("program_name", "model"))
cmp_dt <- if (length(comparison_results) > 0) rbindlist(comparison_results, fill = TRUE) else empty_dt(c("program_name", "model"))

if (nrow(model_dt) > 0) {
  model_dt[, fdr_p_asd := p.adjust(p_asd, method = "BH")]
  setorder(model_dt, program_name, model)
}
if (nrow(cmp_dt) > 0) {
  cmp_dt[, fdr_p_asd := p.adjust(p_asd, method = "BH")]
  setorder(cmp_dt, program_name, model)
}

safe_fwrite(scores_long, file.path(TABDIR, "10_Gandal_program_scores_long.tsv"))
safe_fwrite(model_dt, file.path(TABDIR, "11_Gandal_program_model_results.tsv"))
safe_fwrite(cmp_dt, file.path(TABDIR, "12_Gandal_program_model_comparison.tsv"))
if (nrow(model_dt) > 0) safe_fwrite(model_dt[program_name %in% c("SFARI_all", "midPrenatal_SFARI_top20")], file.path(TABDIR, "13_Gandal_primary_program_model_results.tsv"))

# -----------------------------
# Plots
# -----------------------------
primary_dt <- model_dt[program_name %in% c("SFARI_all", "midPrenatal_SFARI_top20")]
if (nrow(primary_dt) > 0) {
  primary_dt[, model := factor(model, levels = c("M0", "M1", "M2"))]
  primary_dt[, program_name := factor(program_name, levels = c("SFARI_all", "midPrenatal_SFARI_top20"))]

  p1 <- ggplot(primary_dt, aes(x = beta_asd, y = model, color = program_name)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
    geom_point(position = position_dodge(width = 0.5), size = 2.2) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", position = position_dodge(width = 0.5), width = 0.2) +
    facet_wrap(~ program_name, scales = "free_x") +
    labs(x = "Diagnosis coefficient (ASD vs Control)", y = "Model") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
  ggsave(file.path(PLOTDIR, "01_primary_program_covariate_sensitivity.png"), p1, width = 8.5, height = 4.8, dpi = 300)
  ggsave(file.path(PLOTDIR, "01_primary_program_covariate_sensitivity.pdf"), p1, width = 8.5, height = 4.8)

  shift_dt <- cmp_dt[program_name %in% c("SFARI_all", "midPrenatal_SFARI_top20")]
  if (nrow(shift_dt) > 0) {
    shift_dt[, model := factor(model, levels = c("M0", "M1", "M2"))]
    shift_dt[, program_name := factor(program_name, levels = c("SFARI_all", "midPrenatal_SFARI_top20"))]
    p2 <- ggplot(shift_dt[model != "M0"], aes(x = model, y = delta_beta_vs_M0, group = program_name, color = program_name)) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey70") +
      geom_point(size = 2.2) +
      geom_line(linewidth = 0.7) +
      labs(x = "Model", y = "Change in diagnosis coefficient vs M0") +
      theme_bw(base_size = 11)
    ggsave(file.path(PLOTDIR, "02_primary_program_model_shift.png"), p2, width = 6.8, height = 4.4, dpi = 300)
    ggsave(file.path(PLOTDIR, "02_primary_program_model_shift.pdf"), p2, width = 6.8, height = 4.4)
  }
}

# -----------------------------
# Run summary
# -----------------------------
run_summary <- rbindlist(list(
  data.table(item = "n_samples_used", value = nrow(meta_use)),
  data.table(item = "n_asd", value = sum(meta_use$is_asd)),
  data.table(item = "n_control", value = sum(!meta_use$is_asd)),
  data.table(item = "n_regions", value = uniqueN(meta_use$region)),
  data.table(item = "mapped_rows", value = sum(mapped_rows)),
  data.table(item = "unique_gene_symbols", value = nrow(expr_gene_symbol)),
  data.table(item = "age_column", value = age_col %||% NA_character_),
  data.table(item = "pmi_column", value = pmi_col %||% NA_character_),
  data.table(item = "rin_column", value = rin_col %||% NA_character_),
  data.table(item = "sex_column", value = sex_col %||% NA_character_),
  data.table(item = "batch_column", value = batch_col %||% NA_character_),
  data.table(item = "adjusted_terms", value = if (length(adj_terms) == 0) "<none>" else paste(adj_terms, collapse = ",")),
  data.table(item = "n_programs_modeled", value = uniqueN(model_dt$program_name)),
  data.table(item = "best_primary_by_M2_p", value = {
    x <- model_dt[program_name %in% c("SFARI_all", "midPrenatal_SFARI_top20") & model == "M2"]
    if (nrow(x) == 0) NA_character_ else x[which.min(p_asd), paste(program_name, signif(beta_asd, 4), signif(p_asd, 4), sep = " | ")]
  })
), fill = TRUE)
safe_fwrite(run_summary, file.path(METADIR, "14_run_summary.tsv"))

log_msg("Step03b Gandal covariate sensitivity v3 completed successfully.")
