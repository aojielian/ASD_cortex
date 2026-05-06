#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")

GANDAL_RDATA <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"
BRAINSPAN_ROWMAP <- file.path(METADIR, "12c_BrainSpan_rows_metadata_standardized.tsv")
PROGRAM_FILE <- file.path(TABDIR, "20_midPrenatal_core_programs.tsv")

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step03_Gandal_bulk_validation_v4.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
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

# -----------------------------
# Inputs
# -----------------------------
if (!file.exists(PROGRAM_FILE)) stop("Missing program file: ", PROGRAM_FILE)
if (!file.exists(BRAINSPAN_ROWMAP)) stop("Missing BrainSpan row mapping file: ", BRAINSPAN_ROWMAP)
if (!file.exists(GANDAL_RDATA)) stop("Missing Gandal RData: ", GANDAL_RDATA)

program_tbl <- fread(PROGRAM_FILE)
program_tbl[, gene_symbol := toupper(trimws(gene_symbol))]
program_tbl <- unique(program_tbl[gene_symbol != "" & !is.na(gene_symbol)])
safe_fwrite(program_tbl[, .N, by = program_name], file.path(METADIR, "30b_Gandal_input_program_summary.tsv"))

rowmap <- fread(BRAINSPAN_ROWMAP)
# required columns from prior Step01A output
needed_cols <- c("gene_symbol", "ensembl_gene_id")
if (!all(needed_cols %in% names(rowmap))) {
  stop("BrainSpan row mapping file must contain columns: ", paste(needed_cols, collapse = ", "))
}
rowmap[, gene_symbol := toupper(trimws(gene_symbol))]
rowmap[, ensembl_gene_id := toupper(trimws(as.character(ensembl_gene_id)))]
rowmap[, ensembl_stable := sub("\\..*$", "", ensembl_gene_id)]
rowmap <- unique(
  rowmap[gene_symbol != "" & !is.na(gene_symbol) &
           ensembl_stable != "" & !is.na(ensembl_stable),
         .(gene_symbol, ensembl_gene_id, ensembl_stable)]
)

safe_fwrite(rowmap, file.path(METADIR, "30c_BrainSpan_geneSymbol_to_ensembl.tsv"))

# -----------------------------
# Load Gandal datExpr/datMeta
# -----------------------------
log_msg("Loading Gandal RData: ", GANDAL_RDATA)
e <- new.env(parent = emptyenv())
load(GANDAL_RDATA, envir = e)

if (!all(c("datExpr", "datMeta") %in% ls(e))) {
  stop("Expected datExpr and datMeta in Gandal RData.")
}

datExpr <- get("datExpr", envir = e)
datMeta <- as.data.table(get("datMeta", envir = e))

if (!is.matrix(datExpr)) datExpr <- as.matrix(datExpr)
mode(datExpr) <- "numeric"

# quick inventory
obj_info <- data.table(
  object_name = c("datExpr", "datMeta"),
  class = c(paste(class(datExpr), collapse = ";"), paste(class(datMeta), collapse = ";")),
  nrow = c(nrow(datExpr), nrow(datMeta)),
  ncol = c(ncol(datExpr), ncol(datMeta))
)
safe_fwrite(obj_info, file.path(METADIR, "30_Gandal_RData_object_inventory.tsv"))

# -----------------------------
# Harmonize metadata
# -----------------------------
required_meta_cols <- c("sample_id", "Diagnosis", "region")
if (!all(required_meta_cols %in% names(datMeta))) {
  stop("datMeta must contain columns: ", paste(required_meta_cols, collapse = ", "))
}

meta <- copy(datMeta)
meta[, sample_id := as.character(sample_id)]
meta[, Diagnosis := as.character(Diagnosis)]
meta[, region := as.character(region)]
meta[, diagnosis := fifelse(Diagnosis == "ASD", "ASD",
                            fifelse(Diagnosis == "CTL", "Control", NA_character_))]
meta <- meta[!is.na(diagnosis)]
meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
meta[, is_asd := diagnosis == "ASD"]

# -----------------------------
# Harmonize expression
# datExpr: genes x samples, rownames like ENSG00000000003.14_1
# -----------------------------
if (is.null(rownames(datExpr)) || is.null(colnames(datExpr))) {
  stop("datExpr must have rownames and colnames.")
}

common_samples <- intersect(colnames(datExpr), meta$sample_id)
if (length(common_samples) < 20) {
  stop("Too few overlapping samples between datExpr colnames and datMeta sample_id: ", length(common_samples))
}

meta_use <- meta[match(common_samples, sample_id)]
expr_use <- datExpr[, common_samples, drop = FALSE]

expr_gene_raw <- rownames(expr_use)
expr_ensembl_stable <- sub("_.*$", "", expr_gene_raw)
expr_ensembl_stable <- sub("\\..*$", "", expr_ensembl_stable)
expr_ensembl_stable <- toupper(trimws(expr_ensembl_stable))

gandal_gene_map <- data.table(
  gandal_rowname = expr_gene_raw,
  ensembl_stable = expr_ensembl_stable
)

gandal_gene_map <- merge(
  gandal_gene_map,
  rowmap[, .(ensembl_stable, gene_symbol)],
  by = "ensembl_stable",
  all.x = TRUE,
  all.y = FALSE
)

safe_fwrite(gandal_gene_map, file.path(METADIR, "30d_Gandal_gene_mapping.tsv"))

mapped_rows <- !is.na(gandal_gene_map$gene_symbol)
log_msg("Gandal rows mapped to gene_symbol: ", sum(mapped_rows), " / ", nrow(gandal_gene_map))

if (sum(mapped_rows) < 10000) {
  stop("Too few Gandal rows mapped to gene_symbol after ENSG parsing.")
}

expr_mapped <- expr_use[mapped_rows, , drop = FALSE]
mapped_symbols <- gandal_gene_map$gene_symbol[mapped_rows]

# collapse duplicate gene symbols by mean
expr_dt <- as.data.table(expr_mapped)
expr_dt[, gene_symbol := mapped_symbols]
expr_collapsed <- expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol]
expr_gene_symbol <- as.matrix(expr_collapsed[, -1])
rownames(expr_gene_symbol) <- expr_collapsed$gene_symbol
mode(expr_gene_symbol) <- "numeric"

log_msg("Collapsed Gandal matrix to unique gene symbols: ", nrow(expr_gene_symbol), " genes x ", ncol(expr_gene_symbol), " samples")

# -----------------------------
# Mapping diagnostics for programs
# -----------------------------
program_diag <- rbindlist(lapply(sort(unique(program_tbl$program_name)), function(pn) {
  gs <- unique(program_tbl[program_name == pn, gene_symbol])
  gs_in_rowmap <- intersect(gs, rowmap$gene_symbol)
  gs_in_gandal <- intersect(gs, rownames(expr_gene_symbol))
  data.table(
    program_name = pn,
    n_genes_program = length(gs),
    n_genes_in_rowmap = length(gs_in_rowmap),
    n_genes_in_gandal = length(gs_in_gandal)
  )
}))
safe_fwrite(program_diag, file.path(METADIR, "30e_Gandal_program_mapping_summary.tsv"))

# -----------------------------
# Score programs
# -----------------------------
score_list <- list()
overall_list <- list()
region_list <- list()
loro_list <- list()

for (pn in sort(unique(program_tbl$program_name))) {
  genes <- intersect(program_tbl[program_name == pn, gene_symbol], rownames(expr_gene_symbol))
  log_msg("Program ", pn, " mapped genes in Gandal: ", length(genes))
  if (length(genes) < 5) next

  score <- colMeans(expr_gene_symbol[genes, , drop = FALSE], na.rm = TRUE)
  dt <- data.table(sample_id = common_samples, score = as.numeric(score))
  dt <- cbind(dt, meta_use[, .(diagnosis, is_asd, region)])
  dt[, program_name := pn]
  dt[, region_factor := factor(region)]

  if (uniqueN(dt$diagnosis) < 2 || sum(dt$diagnosis == "ASD") < 2 || sum(dt$diagnosis == "Control") < 2) next

  score_list[[pn]] <- dt

  wt <- wilcox.test(score ~ diagnosis, data = dt, exact = FALSE)
  auc <- calc_auc(dt$score, dt$is_asd)

  fit_simple <- lm(score ~ diagnosis, data = dt)
  coef_simple <- summary(fit_simple)$coefficients
  beta_simple <- coef_simple["diagnosisASD", "Estimate"]
  p_simple <- coef_simple["diagnosisASD", "Pr(>|t|)"]

  if (uniqueN(dt$region_factor) >= 2) {
    fit_region <- lm(score ~ diagnosis + region_factor, data = dt)
    coef_region <- summary(fit_region)$coefficients
    beta_region <- coef_region["diagnosisASD", "Estimate"]
    p_region <- coef_region["diagnosisASD", "Pr(>|t|)"]
  } else {
    beta_region <- NA_real_
    p_region <- NA_real_
  }

  overall_list[[pn]] <- data.table(
    program_name = pn,
    n_genes_mapped_gandal = length(genes),
    n_samples = nrow(dt),
    n_asd = sum(dt$is_asd),
    n_control = sum(!dt$is_asd),
    n_regions = uniqueN(dt$region),
    mean_score_control = dt[diagnosis == "Control", mean(score)],
    mean_score_asd = dt[diagnosis == "ASD", mean(score)],
    delta_asd_minus_control = dt[diagnosis == "ASD", mean(score)] - dt[diagnosis == "Control", mean(score)],
    auc_asd_higher = auc,
    wilcox_p = wt$p.value,
    lm_beta_asd = beta_simple,
    lm_p = p_simple,
    lm_region_beta_asd = beta_region,
    lm_region_p = p_region
  )

  for (rg in sort(unique(dt$region))) {
    dd <- dt[region == rg]
    if (nrow(dd) < 6) next
    if (sum(dd$diagnosis == "ASD") < 2 || sum(dd$diagnosis == "Control") < 2) next

    wt_rg <- wilcox.test(score ~ diagnosis, data = dd, exact = FALSE)
    region_list[[paste0(pn, "__", rg)]] <- data.table(
      program_name = pn,
      region = rg,
      n_samples = nrow(dd),
      n_asd = sum(dd$diagnosis == "ASD"),
      n_control = sum(dd$diagnosis == "Control"),
      mean_score_control = dd[diagnosis == "Control", mean(score)],
      mean_score_asd = dd[diagnosis == "ASD", mean(score)],
      delta_asd_minus_control = dd[diagnosis == "ASD", mean(score)] - dd[diagnosis == "Control", mean(score)],
      auc_asd_higher = calc_auc(dd$score, dd$is_asd),
      wilcox_p = wt_rg$p.value
    )
  }

  if (uniqueN(dt$region) >= 3) {
    for (rg in sort(unique(dt$region))) {
      dd <- copy(dt[region != rg])
      if (nrow(dd) < 20) next
      if (sum(dd$diagnosis == "ASD") < 2 || sum(dd$diagnosis == "Control") < 2) next
      dd[, region_factor := factor(region)]

      fit_loro <- lm(score ~ diagnosis + region_factor, data = dd)
      coef_loro <- summary(fit_loro)$coefficients

      loro_list[[paste0(pn, "__leaveout__", rg)]] <- data.table(
        program_name = pn,
        left_out_region = rg,
        n_samples = nrow(dd),
        n_regions_remaining = uniqueN(dd$region),
        beta_asd = coef_loro["diagnosisASD", "Estimate"],
        p_asd = coef_loro["diagnosisASD", "Pr(>|t|)"]
      )
    }
  }
}

score_cols <- c("sample_id","score","diagnosis","is_asd","region","program_name","region_factor")
overall_cols <- c("program_name","n_genes_mapped_gandal","n_samples","n_asd","n_control","n_regions",
                  "mean_score_control","mean_score_asd","delta_asd_minus_control","auc_asd_higher",
                  "wilcox_p","lm_beta_asd","lm_p","lm_region_beta_asd","lm_region_p")
region_cols <- c("program_name","region","n_samples","n_asd","n_control","mean_score_control","mean_score_asd",
                 "delta_asd_minus_control","auc_asd_higher","wilcox_p")
loro_cols <- c("program_name","left_out_region","n_samples","n_regions_remaining","beta_asd","p_asd")

scores_long <- if (length(score_list) > 0) rbindlist(score_list, use.names = TRUE, fill = TRUE) else empty_dt(score_cols)
overall_dt <- if (length(overall_list) > 0) rbindlist(overall_list, use.names = TRUE, fill = TRUE) else empty_dt(overall_cols)
region_dt <- if (length(region_list) > 0) rbindlist(region_list, use.names = TRUE, fill = TRUE) else empty_dt(region_cols)
loro_dt <- if (length(loro_list) > 0) rbindlist(loro_list, use.names = TRUE, fill = TRUE) else empty_dt(loro_cols)

if (nrow(overall_dt) > 0) {
  overall_dt[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  overall_dt[, fdr_lm := p.adjust(lm_p, method = "BH")]
  overall_dt[, fdr_lm_region := p.adjust(lm_region_p, method = "BH")]
  overall_dt[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(overall_dt, fdr_lm_region, fdr_lm, fdr_wilcox, -abs_delta_asd_minus_control)
}

if (nrow(region_dt) > 0) {
  region_dt[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  region_dt[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(region_dt, fdr_wilcox, wilcox_p, -abs_delta_asd_minus_control)
}

if (nrow(loro_dt) > 0) {
  loro_dt[, fdr_p_asd := p.adjust(p_asd, method = "BH")]
  loro_dt[, abs_beta_asd := abs(beta_asd)]
  setorder(loro_dt, fdr_p_asd, p_asd, -abs_beta_asd)
}

safe_fwrite(scores_long, file.path(TABDIR, "30_Gandal_program_scores_long.tsv"))
safe_fwrite(overall_dt, file.path(TABDIR, "31_Gandal_program_dx_tests.tsv"))
safe_fwrite(region_dt, file.path(TABDIR, "31b_Gandal_program_region_tests.tsv"))
safe_fwrite(loro_dt, file.path(TABDIR, "31c_Gandal_program_leaveOneRegionOut.tsv"))

if (nrow(overall_dt) > 0) {
  best_program <- overall_dt[1]
  run_summary <- rbindlist(list(
    data.table(section = "Gandal", metric = "metadata_object", value = "datMeta"),
    data.table(section = "Gandal", metric = "expression_object", value = "datExpr"),
    data.table(section = "Gandal", metric = "orientation", value = "genes_by_samples"),
    data.table(section = "Gandal", metric = "n_samples_used", value = nrow(meta_use)),
    data.table(section = "Gandal", metric = "n_asd", value = sum(meta_use$is_asd)),
    data.table(section = "Gandal", metric = "n_control", value = sum(!meta_use$is_asd)),
    data.table(section = "Gandal", metric = "n_regions", value = uniqueN(meta_use$region)),
    data.table(section = "Gandal", metric = "n_rows_mapped_gene_symbol", value = sum(mapped_rows)),
    data.table(section = "Gandal", metric = "n_unique_gene_symbols", value = nrow(expr_gene_symbol)),
    data.table(section = "Programs", metric = "n_programs_tested", value = nrow(overall_dt)),
    data.table(section = "Programs", metric = "best_program", value = best_program$program_name),
    data.table(section = "BestProgram", metric = "n_genes_mapped_gandal", value = best_program$n_genes_mapped_gandal),
    data.table(section = "BestProgram", metric = "delta_asd_minus_control", value = best_program$delta_asd_minus_control),
    data.table(section = "BestProgram", metric = "auc_asd_higher", value = best_program$auc_asd_higher),
    data.table(section = "BestProgram", metric = "wilcox_p", value = best_program$wilcox_p),
    data.table(section = "BestProgram", metric = "fdr_wilcox", value = best_program$fdr_wilcox),
    data.table(section = "BestProgram", metric = "lm_beta_asd", value = best_program$lm_beta_asd),
    data.table(section = "BestProgram", metric = "lm_p", value = best_program$lm_p),
    data.table(section = "BestProgram", metric = "lm_region_beta_asd", value = best_program$lm_region_beta_asd),
    data.table(section = "BestProgram", metric = "lm_region_p", value = best_program$lm_region_p),
    data.table(section = "BestProgram", metric = "fdr_lm_region", value = best_program$fdr_lm_region)
  ))
} else {
  run_summary <- rbindlist(list(
    data.table(section = "Gandal", metric = "metadata_object", value = "datMeta"),
    data.table(section = "Gandal", metric = "expression_object", value = "datExpr"),
    data.table(section = "Gandal", metric = "orientation", value = "genes_by_samples"),
    data.table(section = "Gandal", metric = "n_samples_used", value = nrow(meta_use)),
    data.table(section = "Gandal", metric = "n_asd", value = sum(meta_use$is_asd)),
    data.table(section = "Gandal", metric = "n_control", value = sum(!meta_use$is_asd)),
    data.table(section = "Gandal", metric = "n_regions", value = uniqueN(meta_use$region)),
    data.table(section = "Gandal", metric = "n_rows_mapped_gene_symbol", value = sum(mapped_rows)),
    data.table(section = "Gandal", metric = "n_unique_gene_symbols", value = nrow(expr_gene_symbol)),
    data.table(section = "Programs", metric = "n_programs_tested", value = 0),
    data.table(section = "Programs", metric = "note", value = "No valid program-level tests were generated.")
  ))
}

safe_fwrite(run_summary, file.path(METADIR, "32_Step03_Gandal_run_summary.tsv"))

log_msg("Step03 v4 completed successfully.")
