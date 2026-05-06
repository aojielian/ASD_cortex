#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA"
OUTDIR <- file.path(BASE, "step04_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
TMPDIR <- file.path(OUTDIR, "tmp_rawMatrix")

PROGRAM_FILE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/20_midPrenatal_core_programs.tsv"
BRAINSPAN_ROWMAP <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/meta/12c_BrainSpan_rows_metadata_standardized.tsv"

META_FILE <- file.path(BASE, "meta.tsv")
RAW_ZIP <- file.path(BASE, "rawMatrix.zip")

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TMPDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step04B_Velmeshev_program_projection.log")

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

broad_class_from_cluster <- function(x) {
  y <- as.character(x)
  out <- rep("other", length(y))
  out[grepl("^Microglia$", y, ignore.case = FALSE)] <- "MG"
  out[grepl("^AST-", y)] <- "AST"
  out[grepl("^OPC$", y)] <- "OPC"
  out[grepl("^Oligodendrocytes$", y)] <- "ODC"
  out[grepl("^Endothelial$", y)] <- "END"
  out[grepl("^IN-", y)] <- "INN"
  out[grepl("^L2/3$|^L4$|^L5/6$|^L5/6-CC$|^Neu-", y)] <- "EXN"
  out
}

# -----------------------------
# Inputs
# -----------------------------
for (f in c(PROGRAM_FILE, BRAINSPAN_ROWMAP, META_FILE, RAW_ZIP)) {
  if (!file.exists(f)) stop("Missing required file: ", f)
}

program_tbl <- fread(PROGRAM_FILE)
program_tbl[, gene_symbol := toupper(trimws(gene_symbol))]
program_tbl <- unique(program_tbl[gene_symbol != "" & !is.na(gene_symbol)])

rowmap <- fread(BRAINSPAN_ROWMAP)
rowmap[, gene_symbol := toupper(trimws(gene_symbol))]
rowmap[, ensembl_gene_id := toupper(trimws(as.character(ensembl_gene_id)))]
rowmap[, ensembl_stable := sub("\\..*$", "", ensembl_gene_id)]
rowmap <- unique(rowmap[gene_symbol != "" & !is.na(gene_symbol) &
                          ensembl_stable != "" & !is.na(ensembl_stable),
                        .(gene_symbol, ensembl_stable)])

meta <- fread(META_FILE)
required_meta <- c("cell", "cluster", "sample", "individual", "region", "diagnosis")
if (!all(required_meta %in% names(meta))) {
  stop("meta.tsv must contain columns: ", paste(required_meta, collapse = ", "))
}
meta[, cell := as.character(cell)]
meta[, cluster := as.character(cluster)]
meta[, sample := as.character(sample)]
meta[, individual := as.character(individual)]
meta[, region := as.character(region)]
meta[, diagnosis := as.character(diagnosis)]
meta <- meta[diagnosis %in% c("ASD", "Control")]
meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
meta[, is_asd := diagnosis == "ASD"]
meta[, broad_class := broad_class_from_cluster(cluster)]

# -----------------------------
# Extract rawMatrix.zip
# -----------------------------
log_msg("Listing and extracting rawMatrix.zip")
zip_list <- unzip(RAW_ZIP, list = TRUE)
safe_fwrite(as.data.table(zip_list), file.path(METADIR, "50b_rawMatrix_zip_listing.tsv"))

needed_files <- c("barcodes.tsv", "genes.tsv", "matrix.mtx")
zip_names <- zip_list$Name
if (!all(needed_files %in% zip_names)) {
  stop("rawMatrix.zip must contain: ", paste(needed_files, collapse = ", "))
}
unzip(RAW_ZIP, files = needed_files, exdir = TMPDIR, overwrite = TRUE)

barcodes_file <- file.path(TMPDIR, "barcodes.tsv")
genes_file <- file.path(TMPDIR, "genes.tsv")
mtx_file <- file.path(TMPDIR, "matrix.mtx")

# -----------------------------
# Read sparse matrix
# -----------------------------
log_msg("Reading barcodes and genes")
barcodes <- fread(barcodes_file, header = FALSE)
genes <- fread(genes_file, header = FALSE)

log_msg("Reading sparse matrix.mtx (this may take a while)")
mat <- readMM(mtx_file)
mat <- as(mat, "dgCMatrix")

if (nrow(genes) != nrow(mat)) stop("genes.tsv rows do not match matrix rows.")
if (nrow(barcodes) != ncol(mat)) stop("barcodes.tsv rows do not match matrix cols.")

# interpret genes.tsv
gene_id_raw <- as.character(genes[[1]])
gene_symbol_raw <- if (ncol(genes) >= 2) as.character(genes[[2]]) else gene_id_raw
gene_symbol_upper <- toupper(trimws(gene_symbol_raw))

# if second col is not useful, fall back to BrainSpan ENSG mapping
gene_id_upper <- toupper(trimws(gene_id_raw))
gene_ensg_stable <- sub("\\..*$", "", gene_id_upper)

gene_map <- data.table(
  gene_id_raw = gene_id_raw,
  gene_symbol_from_genes = gene_symbol_upper,
  ensembl_stable = gene_ensg_stable
)

# map via rowmap
gene_map <- merge(
  gene_map,
  rowmap[, .(ensembl_stable, gene_symbol_brainspan = gene_symbol)],
  by = "ensembl_stable",
  all.x = TRUE,
  all.y = FALSE
)

# choose final symbol
gene_map[, gene_symbol_final := fifelse(
  !is.na(gene_symbol_from_genes) & gene_symbol_from_genes != "" & !grepl("^ENSG", gene_symbol_from_genes),
  gene_symbol_from_genes,
  gene_symbol_brainspan
)]
gene_map[, gene_symbol_final := toupper(trimws(gene_symbol_final))]

safe_fwrite(gene_map, file.path(METADIR, "50c_Velmeshev_gene_mapping.tsv"))

# set matrix dimnames
colnames(mat) <- as.character(barcodes[[1]])
rownames(mat) <- gene_map$gene_symbol_final

# align cells
common_cells <- intersect(colnames(mat), meta$cell)
if (length(common_cells) < 50000) {
  stop("Too few overlapping cells between raw matrix barcodes and meta$cell: ", length(common_cells))
}
meta_use <- meta[match(common_cells, cell)]
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

# library-size normalized log1p(CP10K)
libsize <- Matrix::colSums(mat)
if (any(libsize == 0)) {
  libsize[libsize == 0] <- 1
}
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
    n_genes_in_velmeshev = length(mapped)
  )
}))
safe_fwrite(program_diag, file.path(METADIR, "50d_Velmeshev_program_mapping_summary.tsv"))

# -----------------------------
# Score each program at cell level
# -----------------------------
score_list <- list()

for (pn in sort(unique(program_tbl$program_name))) {
  gs <- intersect(program_tbl[program_name == pn, gene_symbol], rownames(mat))
  log_msg("Program ", pn, " mapped genes in Velmeshev: ", length(gs))
  if (length(gs) < 5) next

  submat <- mat[gs, , drop = FALSE]
  norm_sub <- t(t(submat) / libsize * 1e4)
  score <- Matrix::colMeans(log1p(norm_sub))

  dt <- data.table(
    cell = colnames(mat),
    score = as.numeric(score),
    program_name = pn
  )
  score_list[[pn]] <- dt
}

score_cell <- if (length(score_list) > 0) rbindlist(score_list, use.names = TRUE, fill = TRUE) else empty_dt(c("cell","score","program_name"))
score_cell <- merge(score_cell, meta_use[, .(cell, cluster, broad_class, sample, individual, region, diagnosis, is_asd)], by = "cell", all.x = TRUE)
safe_fwrite(score_cell, file.path(TABDIR, "51_Velmeshev_cell_program_scores.tsv.gz"))

# -----------------------------
# Cluster-level localization (all cells)
# -----------------------------
cluster_loc <- if (nrow(score_cell) > 0) {
  score_cell[, .(
    n_cells = .N,
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ), by = .(program_name, cluster, broad_class)][order(program_name, -mean_score)]
} else empty_dt(c("program_name","cluster","broad_class","n_cells","mean_score","median_score"))
safe_fwrite(cluster_loc, file.path(TABDIR, "52_Velmeshev_cluster_program_localization.tsv"))

broad_loc <- if (nrow(score_cell) > 0) {
  score_cell[, .(
    n_cells = .N,
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ), by = .(program_name, broad_class)][order(program_name, -mean_score)]
} else empty_dt(c("program_name","broad_class","n_cells","mean_score","median_score"))
safe_fwrite(broad_loc, file.path(TABDIR, "52b_Velmeshev_broadClass_program_localization.tsv"))

# -----------------------------
# Donor-level aggregation and dx tests
# -----------------------------
donor_cluster <- if (nrow(score_cell) > 0) {
  score_cell[, .(
    n_cells = .N,
    mean_score = mean(score, na.rm = TRUE)
  ), by = .(program_name, cluster, broad_class, individual, diagnosis, is_asd)]
} else empty_dt(c("program_name","cluster","broad_class","individual","diagnosis","is_asd","n_cells","mean_score"))
safe_fwrite(donor_cluster, file.path(TABDIR, "53_Velmeshev_donor_cluster_program_scores.tsv.gz"))

cluster_tests_list <- list()
if (nrow(donor_cluster) > 0) {
  keys <- unique(donor_cluster[, .(program_name, cluster, broad_class)])
  for (i in seq_len(nrow(keys))) {
    pn <- keys$program_name[i]
    cl <- keys$cluster[i]
    bc <- keys$broad_class[i]
    dd <- donor_cluster[program_name == pn & cluster == cl]
    if (nrow(dd) < 6) next
    if (sum(dd$diagnosis == "ASD") < 2 || sum(dd$diagnosis == "Control") < 2) next

    wt <- wilcox.test(mean_score ~ diagnosis, data = dd, exact = FALSE)
    fit <- lm(mean_score ~ diagnosis, data = dd)
    cf <- summary(fit)$coefficients

    cluster_tests_list[[paste0(pn, "__", cl)]] <- data.table(
      program_name = pn,
      cluster = cl,
      broad_class = bc,
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
cluster_tests <- if (length(cluster_tests_list) > 0) rbindlist(cluster_tests_list, use.names = TRUE, fill = TRUE) else empty_dt(c(
  "program_name","cluster","broad_class","n_donors","n_asd","n_control","mean_score_control","mean_score_asd",
  "delta_asd_minus_control","auc_asd_higher","wilcox_p","lm_beta_asd","lm_p"
))
if (nrow(cluster_tests) > 0) {
  cluster_tests[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  cluster_tests[, fdr_lm := p.adjust(lm_p, method = "BH")]
  cluster_tests[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(cluster_tests, fdr_lm, fdr_wilcox, -abs_delta_asd_minus_control)
}
safe_fwrite(cluster_tests, file.path(TABDIR, "54_Velmeshev_cluster_dx_tests.tsv"))

# broad class donor aggregation
donor_broad <- if (nrow(score_cell) > 0) {
  score_cell[, .(
    n_cells = .N,
    mean_score = mean(score, na.rm = TRUE)
  ), by = .(program_name, broad_class, individual, diagnosis, is_asd)]
} else empty_dt(c("program_name","broad_class","individual","diagnosis","is_asd","n_cells","mean_score"))
safe_fwrite(donor_broad, file.path(TABDIR, "54b_Velmeshev_donor_broadClass_program_scores.tsv.gz"))

broad_tests_list <- list()
if (nrow(donor_broad) > 0) {
  keys <- unique(donor_broad[, .(program_name, broad_class)])
  for (i in seq_len(nrow(keys))) {
    pn <- keys$program_name[i]
    bc <- keys$broad_class[i]
    dd <- donor_broad[program_name == pn & broad_class == bc]
    if (nrow(dd) < 6) next
    if (sum(dd$diagnosis == "ASD") < 2 || sum(dd$diagnosis == "Control") < 2) next

    wt <- wilcox.test(mean_score ~ diagnosis, data = dd, exact = FALSE)
    fit <- lm(mean_score ~ diagnosis, data = dd)
    cf <- summary(fit)$coefficients

    broad_tests_list[[paste0(pn, "__", bc)]] <- data.table(
      program_name = pn,
      broad_class = bc,
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
broad_tests <- if (length(broad_tests_list) > 0) rbindlist(broad_tests_list, use.names = TRUE, fill = TRUE) else empty_dt(c(
  "program_name","broad_class","n_donors","n_asd","n_control","mean_score_control","mean_score_asd",
  "delta_asd_minus_control","auc_asd_higher","wilcox_p","lm_beta_asd","lm_p"
))
if (nrow(broad_tests) > 0) {
  broad_tests[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  broad_tests[, fdr_lm := p.adjust(lm_p, method = "BH")]
  broad_tests[, abs_delta_asd_minus_control := abs(delta_asd_minus_control)]
  setorder(broad_tests, fdr_lm, fdr_wilcox, -abs_delta_asd_minus_control)
}
safe_fwrite(broad_tests, file.path(TABDIR, "54c_Velmeshev_broadClass_dx_tests.tsv"))

# -----------------------------
# Run summary
# -----------------------------
best_cluster <- if (nrow(cluster_tests) > 0) cluster_tests[1] else NULL
best_broad <- if (nrow(broad_tests) > 0) broad_tests[1] else NULL

run_summary <- rbindlist(list(
  data.table(section = "Velmeshev", metric = "n_cells_used", value = ncol(mat)),
  data.table(section = "Velmeshev", metric = "n_genes_used", value = nrow(mat)),
  data.table(section = "Velmeshev", metric = "n_asd_cells", value = sum(meta_use$is_asd)),
  data.table(section = "Velmeshev", metric = "n_control_cells", value = sum(!meta_use$is_asd)),
  data.table(section = "Velmeshev", metric = "n_individuals", value = uniqueN(meta_use$individual)),
  data.table(section = "Velmeshev", metric = "n_clusters", value = uniqueN(meta_use$cluster)),
  data.table(section = "Programs", metric = "n_programs_tested_cluster", value = uniqueN(cluster_tests$program_name)),
  data.table(section = "Programs", metric = "n_programs_tested_broadClass", value = uniqueN(broad_tests$program_name))
), fill = TRUE)

if (!is.null(best_cluster)) {
  run_summary <- rbind(run_summary, rbindlist(list(
    data.table(section = "BestCluster", metric = "program_name", value = best_cluster$program_name),
    data.table(section = "BestCluster", metric = "cluster", value = best_cluster$cluster),
    data.table(section = "BestCluster", metric = "delta_asd_minus_control", value = best_cluster$delta_asd_minus_control),
    data.table(section = "BestCluster", metric = "lm_p", value = best_cluster$lm_p),
    data.table(section = "BestCluster", metric = "fdr_lm", value = best_cluster$fdr_lm)
  )), fill = TRUE)
}

if (!is.null(best_broad)) {
  run_summary <- rbind(run_summary, rbindlist(list(
    data.table(section = "BestBroadClass", metric = "program_name", value = best_broad$program_name),
    data.table(section = "BestBroadClass", metric = "broad_class", value = best_broad$broad_class),
    data.table(section = "BestBroadClass", metric = "delta_asd_minus_control", value = best_broad$delta_asd_minus_control),
    data.table(section = "BestBroadClass", metric = "lm_p", value = best_broad$lm_p),
    data.table(section = "BestBroadClass", metric = "fdr_lm", value = best_broad$fdr_lm)
  )), fill = TRUE)
}

safe_fwrite(run_summary, file.path(METADIR, "55_Step04B_run_summary.tsv"))

log_msg("Step04B completed successfully.")
