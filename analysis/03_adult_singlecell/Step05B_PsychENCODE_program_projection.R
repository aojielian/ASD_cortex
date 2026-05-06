#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# -----------------------------
# Output base: keep consistent with ASD_cortex workflow
# -----------------------------
BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
TMPDIR <- file.path(OUTDIR, "tmp_psychencode")

# -----------------------------
# Input data paths: PsychENCODE single-cell
# -----------------------------
PSY_BASE <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024"
META_FILE <- file.path(PSY_BASE, "meta.tsv")
BARCODES_FILE <- file.path(PSY_BASE, "counts_barcodes.tsv.gz")
FEATURES_FILE <- file.path(PSY_BASE, "counts_features.tsv.gz")
MTX_FILE <- file.path(PSY_BASE, "counts_matrix.mtx.gz")

PROGRAM_FILE <- file.path(BASE, "step01_prepare", "tables", "20_midPrenatal_core_programs.tsv")
BRAINSPAN_ROWMAP <- file.path(BASE, "step01_prepare", "meta", "12c_BrainSpan_rows_metadata_standardized.tsv")

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TMPDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step05B_PsychENCODE_program_projection.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

empty_dt <- function(cols) {
  out <- data.table(matrix(ncol = length(cols), nrow = 0))
  setnames(out, cols)
  out
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

# -----------------------------
# Inputs
# -----------------------------
for (f in c(PROGRAM_FILE, BRAINSPAN_ROWMAP, META_FILE, BARCODES_FILE, FEATURES_FILE, MTX_FILE)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

program_tbl <- fread(PROGRAM_FILE)
program_tbl[, gene_symbol := toupper(trimws(gene_symbol))]
program_tbl <- unique(program_tbl[gene_symbol != "" & !is.na(gene_symbol)])
safe_fwrite(program_tbl[, .N, by = program_name], file.path(METADIR, "70b_PsychENCODE_input_program_summary.tsv"))

rowmap <- fread(BRAINSPAN_ROWMAP)
needed_cols <- c("gene_symbol", "ensembl_gene_id")
if (!all(needed_cols %in% names(rowmap))) {
  stop("BrainSpan row mapping file must contain columns: ", paste(needed_cols, collapse = ", "))
}
rowmap[, gene_symbol := toupper(trimws(gene_symbol))]
rowmap[, ensembl_gene_id := toupper(trimws(as.character(ensembl_gene_id)))]
rowmap[, ensembl_stable := sub("\\..*$", "", ensembl_gene_id)]
rowmap <- unique(rowmap[gene_symbol != "" & !is.na(gene_symbol) &
                          ensembl_stable != "" & !is.na(ensembl_stable),
                        .(gene_symbol, ensembl_stable)])
safe_fwrite(rowmap, file.path(METADIR, "70c_BrainSpan_geneSymbol_to_ensembl.tsv"))

meta <- fread(META_FILE)
required_meta <- c("Cell_ID", "individual_ID", "Diagnosis", "annotation", "Subtype")
if (!all(required_meta %in% names(meta))) {
  stop("meta.tsv must contain columns: ", paste(required_meta, collapse = ", "))
}
meta[, Cell_ID := as.character(Cell_ID)]
meta[, individual_ID := as.character(individual_ID)]
meta[, Diagnosis := as.character(Diagnosis)]
meta[, annotation := as.character(annotation)]
meta[, Subtype := as.character(Subtype)]
meta[, Brain_Region := if ("Brain_Region" %in% names(meta)) as.character(Brain_Region) else "NA_region"]
meta[, diagnosis := fifelse(Diagnosis == "ASD", "ASD",
                            fifelse(Diagnosis == "CTL", "Control", NA_character_))]
meta <- meta[!is.na(diagnosis)]
meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
meta[, is_asd := diagnosis == "ASD"]

# -----------------------------
# Read sparse matrix inputs
# -----------------------------
log_msg("Reading barcodes and features")
barcodes <- fread(cmd = paste("zcat", shQuote(BARCODES_FILE)), header = FALSE)
features <- fread(cmd = paste("zcat", shQuote(FEATURES_FILE)), header = FALSE)

log_msg("Reading counts_matrix.mtx.gz (this may take a while)")
mat <- readMM(gzfile(MTX_FILE))
mat <- as(mat, "dgCMatrix")

if (nrow(features) != nrow(mat)) stop("features rows do not match matrix rows.")
if (nrow(barcodes) != ncol(mat)) stop("barcodes rows do not match matrix cols.")

# features: robust handling
feature_id_raw <- as.character(features[[1]])
feature_symbol_raw <- if (ncol(features) >= 2) as.character(features[[2]]) else feature_id_raw

feature_id_upper <- toupper(trimws(feature_id_raw))
feature_symbol_upper <- toupper(trimws(feature_symbol_raw))
feature_ensg_stable <- sub("\\..*$", "", feature_id_upper)

gene_map <- data.table(
  feature_id_raw = feature_id_raw,
  feature_symbol_from_features = feature_symbol_upper,
  ensembl_stable = feature_ensg_stable
)

gene_map <- merge(
  gene_map,
  rowmap[, .(ensembl_stable, gene_symbol_brainspan = gene_symbol)],
  by = "ensembl_stable",
  all.x = TRUE,
  all.y = FALSE
)

gene_map[, gene_symbol_final := fifelse(
  !is.na(feature_symbol_from_features) & feature_symbol_from_features != "" & !grepl("^ENSG", feature_symbol_from_features),
  feature_symbol_from_features,
  gene_symbol_brainspan
)]
gene_map[, gene_symbol_final := toupper(trimws(gene_symbol_final))]
safe_fwrite(gene_map, file.path(METADIR, "70d_PsychENCODE_gene_mapping.tsv"))

colnames(mat) <- as.character(barcodes[[1]])
rownames(mat) <- gene_map$gene_symbol_final

# align cells
common_cells <- intersect(colnames(mat), meta$Cell_ID)
if (length(common_cells) < 100000) {
  stop("Too few overlapping cells between matrix barcodes and meta$Cell_ID: ", length(common_cells))
}
meta_use <- meta[match(common_cells, Cell_ID)]
mat <- mat[, common_cells, drop = FALSE]

# remove unmapped genes
keep_gene <- !is.na(rownames(mat)) & rownames(mat) != ""
mat <- mat[keep_gene, , drop = FALSE]

# collapse duplicated gene symbols by sum
log_msg("Collapsing duplicated gene symbols in sparse matrix")
u_genes <- unique(rownames(mat))
f <- factor(rownames(mat), levels = u_genes)
collapse_mat <- sparseMatrix(
  i = as.integer(f),
  j = seq_along(f),
  x = 1,
  dims = c(length(u_genes), nrow(mat))
)
mat_collapsed <- collapse_mat %*% mat
rownames(mat_collapsed) <- u_genes
colnames(mat_collapsed) <- colnames(mat)
mat <- as(mat_collapsed, "dgCMatrix")
rm(mat_collapsed, collapse_mat); gc()

libsize <- Matrix::colSums(mat)
if (any(libsize == 0)) libsize[libsize == 0] <- 1

log_msg("Cells retained: ", ncol(mat), " ; genes retained: ", nrow(mat))

# -----------------------------
# Program mapping diagnostics
# -----------------------------
program_diag <- rbindlist(lapply(sort(unique(program_tbl$program_name)), function(pn) {
  gs <- unique(program_tbl[program_name == pn, gene_symbol])
  mapped <- intersect(gs, rownames(mat))
  data.table(
    program_name = pn,
    n_genes_program = length(gs),
    n_genes_in_psychencode = length(mapped)
  )
}))
safe_fwrite(program_diag, file.path(METADIR, "70e_PsychENCODE_program_mapping_summary.tsv"))

# -----------------------------
# Score each program at cell level
# -----------------------------
score_list <- list()

for (pn in sort(unique(program_tbl$program_name))) {
  gs <- intersect(program_tbl[program_name == pn, gene_symbol], rownames(mat))
  log_msg("Program ", pn, " mapped genes in PsychENCODE: ", length(gs))
  if (length(gs) < 5) next

  submat <- mat[gs, , drop = FALSE]
  norm_sub <- t(t(submat) / libsize * 1e4)
  score <- Matrix::colMeans(log1p(norm_sub))

  dt <- data.table(
    Cell_ID = colnames(mat),
    score = as.numeric(score),
    program_name = pn
  )
  score_list[[pn]] <- dt
}

score_cell <- if (length(score_list) > 0) rbindlist(score_list, use.names = TRUE, fill = TRUE) else empty_dt(c("Cell_ID","score","program_name"))
score_cell <- merge(
  score_cell,
  meta_use[, .(Cell_ID, annotation, Subtype, individual_ID, Brain_Region, diagnosis, is_asd)],
  by = "Cell_ID",
  all.x = TRUE
)
safe_fwrite(score_cell, file.path(TABDIR, "71_PsychENCODE_cell_program_scores.tsv.gz"))

# -----------------------------
# annotation / subtype localization
# -----------------------------
ann_loc <- if (nrow(score_cell) > 0) {
  score_cell[, .(
    n_cells = .N,
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ), by = .(program_name, annotation)][order(program_name, -mean_score)]
} else empty_dt(c("program_name","annotation","n_cells","mean_score","median_score"))
safe_fwrite(ann_loc, file.path(TABDIR, "72_PsychENCODE_annotation_program_localization.tsv"))

subtype_loc <- if (nrow(score_cell) > 0) {
  score_cell[, .(
    n_cells = .N,
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ), by = .(program_name, annotation, Subtype)][order(program_name, -mean_score)]
} else empty_dt(c("program_name","annotation","Subtype","n_cells","mean_score","median_score"))
safe_fwrite(subtype_loc, file.path(TABDIR, "72b_PsychENCODE_subtype_program_localization.tsv"))

# -----------------------------
# Donor-level aggregation and dx tests
# -----------------------------
donor_subtype <- if (nrow(score_cell) > 0) {
  score_cell[, .(
    n_cells = .N,
    mean_score = mean(score, na.rm = TRUE)
  ), by = .(program_name, annotation, Subtype, individual_ID, diagnosis, is_asd)]
} else empty_dt(c("program_name","annotation","Subtype","individual_ID","diagnosis","is_asd","n_cells","mean_score"))
safe_fwrite(donor_subtype, file.path(TABDIR, "73_PsychENCODE_donor_subtype_program_scores.tsv.gz"))

subtype_tests_list <- list()
if (nrow(donor_subtype) > 0) {
  keys <- unique(donor_subtype[, .(program_name, annotation, Subtype)])
  for (i in seq_len(nrow(keys))) {
    pn <- keys$program_name[i]
    ann <- keys$annotation[i]
    st <- keys$Subtype[i]
    dd <- donor_subtype[program_name == pn & annotation == ann & Subtype == st]
    if (nrow(dd) < 8) next
    if (sum(dd$diagnosis == "ASD") < 2 || sum(dd$diagnosis == "Control") < 2) next

    wt <- wilcox.test(mean_score ~ diagnosis, data = dd, exact = FALSE)
    fit <- lm(mean_score ~ diagnosis, data = dd)
    cf <- summary(fit)$coefficients

    subtype_tests_list[[paste0(pn, "__", st)]] <- data.table(
      program_name = pn,
      annotation = ann,
      Subtype = st,
      n_donors = nrow(dd),
      n_asd = sum(dd$diagnosis == "ASD"),
      n_control = sum(dd$diagnosis == "Control"),
      mean_score_control = dd[diagnosis == "Control", mean(mean_score)],
      mean_score_asd = dd[diagnosis == "ASD", mean(mean_score)],
      delta_asd_minus_control = dd[diagnosis == "ASD", mean(mean_score)] - dd[diagnosis == "Control", mean(mean_score)],
      auc_asd_higher = calc_auc(dd$mean_score, dd$is_asd),
      wilcox_p = wt$p.value,
      lm_beta_asd = cf["diagnosisASD", "Estimate"],
      lm_p = cf["diagnosisASD", "Pr(>|t|)"]
    )
  }
}
subtype_tests <- if (length(subtype_tests_list) > 0) rbindlist(subtype_tests_list, use.names = TRUE, fill = TRUE) else empty_dt(c(
  "program_name","annotation","Subtype","n_donors","n_asd","n_control","mean_score_control","mean_score_asd",
  "delta_asd_minus_control","auc_asd_higher","wilcox_p","lm_beta_asd","lm_p"
))
if (nrow(subtype_tests) > 0) {
  subtype_tests[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  subtype_tests[, fdr_lm := p.adjust(lm_p, method = "BH")]
  subtype_tests[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(subtype_tests, fdr_lm, fdr_wilcox, -abs_delta_asd_minus_control)
}
safe_fwrite(subtype_tests, file.path(TABDIR, "74_PsychENCODE_subtype_dx_tests.tsv"))

donor_ann <- if (nrow(score_cell) > 0) {
  score_cell[, .(
    n_cells = .N,
    mean_score = mean(score, na.rm = TRUE)
  ), by = .(program_name, annotation, individual_ID, diagnosis, is_asd)]
} else empty_dt(c("program_name","annotation","individual_ID","diagnosis","is_asd","n_cells","mean_score"))
safe_fwrite(donor_ann, file.path(TABDIR, "74b_PsychENCODE_donor_annotation_program_scores.tsv.gz"))

ann_tests_list <- list()
if (nrow(donor_ann) > 0) {
  keys <- unique(donor_ann[, .(program_name, annotation)])
  for (i in seq_len(nrow(keys))) {
    pn <- keys$program_name[i]
    ann <- keys$annotation[i]
    dd <- donor_ann[program_name == pn & annotation == ann]
    if (nrow(dd) < 8) next
    if (sum(dd$diagnosis == "ASD") < 2 || sum(dd$diagnosis == "Control") < 2) next

    wt <- wilcox.test(mean_score ~ diagnosis, data = dd, exact = FALSE)
    fit <- lm(mean_score ~ diagnosis, data = dd)
    cf <- summary(fit)$coefficients

    ann_tests_list[[paste0(pn, "__", ann)]] <- data.table(
      program_name = pn,
      annotation = ann,
      n_donors = nrow(dd),
      n_asd = sum(dd$diagnosis == "ASD"),
      n_control = sum(dd$diagnosis == "Control"),
      mean_score_control = dd[diagnosis == "Control", mean(mean_score)],
      mean_score_asd = dd[diagnosis == "ASD", mean(mean_score)],
      delta_asd_minus_control = dd[diagnosis == "ASD", mean(mean_score)] - dd[diagnosis == "Control", mean(mean_score)],
      auc_asd_higher = calc_auc(dd$mean_score, dd$is_asd),
      wilcox_p = wt$p.value,
      lm_beta_asd = cf["diagnosisASD", "Estimate"],
      lm_p = cf["diagnosisASD", "Pr(>|t|)"]
    )
  }
}
ann_tests <- if (length(ann_tests_list) > 0) rbindlist(ann_tests_list, use.names = TRUE, fill = TRUE) else empty_dt(c(
  "program_name","annotation","n_donors","n_asd","n_control","mean_score_control","mean_score_asd",
  "delta_asd_minus_control","auc_asd_higher","wilcox_p","lm_beta_asd","lm_p"
))
if (nrow(ann_tests) > 0) {
  ann_tests[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  ann_tests[, fdr_lm := p.adjust(lm_p, method = "BH")]
  ann_tests[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(ann_tests, fdr_lm, fdr_wilcox, -abs_delta_asd_minus_control)
}
safe_fwrite(ann_tests, file.path(TABDIR, "74c_PsychENCODE_annotation_dx_tests.tsv"))

# -----------------------------
# Run summary
# -----------------------------
best_subtype <- if (nrow(subtype_tests) > 0) subtype_tests[1] else NULL
best_ann <- if (nrow(ann_tests) > 0) ann_tests[1] else NULL

run_summary <- rbindlist(list(
  data.table(section = "PsychENCODE", metric = "n_cells_used", value = ncol(mat)),
  data.table(section = "PsychENCODE", metric = "n_genes_used", value = nrow(mat)),
  data.table(section = "PsychENCODE", metric = "n_asd_cells", value = sum(meta_use$is_asd)),
  data.table(section = "PsychENCODE", metric = "n_control_cells", value = sum(!meta_use$is_asd)),
  data.table(section = "PsychENCODE", metric = "n_individuals", value = uniqueN(meta_use$individual_ID)),
  data.table(section = "PsychENCODE", metric = "n_annotations", value = uniqueN(meta_use$annotation)),
  data.table(section = "PsychENCODE", metric = "n_subtypes", value = uniqueN(meta_use$Subtype)),
  data.table(section = "Programs", metric = "n_programs_tested_annotation", value = uniqueN(ann_tests$program_name)),
  data.table(section = "Programs", metric = "n_programs_tested_subtype", value = uniqueN(subtype_tests$program_name))
), fill = TRUE)

if (!is.null(best_ann)) {
  run_summary <- rbind(run_summary, rbindlist(list(
    data.table(section = "BestAnnotation", metric = "program_name", value = best_ann$program_name),
    data.table(section = "BestAnnotation", metric = "annotation", value = best_ann$annotation),
    data.table(section = "BestAnnotation", metric = "delta_asd_minus_control", value = best_ann$delta_asd_minus_control),
    data.table(section = "BestAnnotation", metric = "lm_p", value = best_ann$lm_p),
    data.table(section = "BestAnnotation", metric = "fdr_lm", value = best_ann$fdr_lm)
  )), fill = TRUE)
}

if (!is.null(best_subtype)) {
  run_summary <- rbind(run_summary, rbindlist(list(
    data.table(section = "BestSubtype", metric = "program_name", value = best_subtype$program_name),
    data.table(section = "BestSubtype", metric = "Subtype", value = best_subtype$Subtype),
    data.table(section = "BestSubtype", metric = "delta_asd_minus_control", value = best_subtype$delta_asd_minus_control),
    data.table(section = "BestSubtype", metric = "lm_p", value = best_subtype$lm_p),
    data.table(section = "BestSubtype", metric = "fdr_lm", value = best_subtype$fdr_lm)
  )), fill = TRUE)
}

safe_fwrite(run_summary, file.path(METADIR, "75_Step05B_run_summary.tsv"))

log_msg("Step05B completed successfully.")
