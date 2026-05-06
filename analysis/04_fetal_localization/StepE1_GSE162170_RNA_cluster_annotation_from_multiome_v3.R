#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

parse_args <- function(x) {
  out <- list(
    base = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex",
    gtf = "/gpfs/hpc/home/lijc/lianaoj/reference/gencode.v49.basic.annotation.gtf.gz",
    rna_counts = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/GSE162170/GSE162170_rna_counts.tsv.gz",
    rna_meta = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/GSE162170/GSE162170_rna_cell_metadata.txt.gz",
    multiome_counts = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/GSE162170/GSE162170_multiome_rna_logcounts.tsv.gz",
    multiome_meta = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/GSE162170/GSE162170_multiome_cell_metadata.txt.gz",
    multiome_cluster_map = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/stepB_external_reference_loading/GSE162170_multiome/GSE162170_multiome_cluster_map.tsv.gz",
    outdir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/stepE1_gse162170_rna_annotation"
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    val <- if (i < length(x)) x[[i + 1L]] else NA_character_
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    out[[key]] <- val
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
outdir <- args$outdir
logdir <- file.path(outdir, "logs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logdir, "StepE1_GSE162170_RNA_cluster_annotation_from_multiome.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
}

read_first_line <- function(path) {
  cmd <- paste("zcat", shQuote(path), "| head -1")
  out <- system(cmd, intern = TRUE)
  if (length(out) == 0L) stop("Could not read first line from: ", path)
  out[[1]]
}

read_matrix_header_cells <- function(path) {
  hdr <- read_first_line(path)
  x <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
  x <- trimws(as.character(x))
  x[nchar(x) > 0]
}

strip_ens_ver <- function(x) sub("\\..*$", "", as.character(x))

parse_gtf_gene_map <- function(gtf) {
  log_msg("Reading gene-level GTF entries from: ", gtf)
  cmd <- paste(
    "zcat", shQuote(gtf),
    "| awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"} $3==\"gene\" {print $1,$9}'"
  )
  dt <- fread(cmd = cmd, sep = "\t", header = FALSE, quote = "")
  if (ncol(dt) < 2L) stop("Failed to parse GTF gene entries from: ", gtf)
  setnames(dt, c("seqname", "attr"))
  dt[, gene_id := sub('.*gene_id "([^"]+)".*', '\\1', attr)]
  dt[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attr)]
  dt[, gene_id_stable := strip_ens_ver(gene_id)]
  dt[, gene_symbol := toupper(trimws(gene_name))]
  dt <- unique(dt[!is.na(gene_id_stable) & gene_id_stable != "" & !is.na(gene_symbol) & gene_symbol != "",
                  .(gene_id_stable, gene_symbol)])
  dt
}

extract_rows_by_gene_ids <- function(path, target_gene_ids, header_cells) {
  target_gene_ids <- unique(as.character(target_gene_ids))
  target_gene_ids <- target_gene_ids[!is.na(target_gene_ids) & target_gene_ids != ""]
  if (length(target_gene_ids) == 0L) stop("No target_gene_ids supplied.")

  tf <- tempfile(fileext = ".txt")
  writeLines(target_gene_ids, tf)
  on.exit(unlink(tf), add = TRUE)

  cmd <- paste0(
    "zcat ", shQuote(path),
    " | tail -n +2 | awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{a[$1]=1; next} ($1 in a){print}' ",
    shQuote(tf), " -"
  )
  dt <- fread(cmd = cmd, sep = "\t", header = FALSE, data.table = TRUE)
  if (nrow(dt) == 0L) stop("No matching genes extracted from: ", path)
  if ((ncol(dt) - 1L) != length(header_cells)) {
    stop("Extracted matrix columns do not match header cells for: ", path,
         "; extracted ncol-1=", ncol(dt) - 1L, "; header cells=", length(header_cells))
  }
  setnames(dt, c("gene_id", header_cells))
  dt
}

as_numeric_matrix <- function(dt, cell_ids) {
  mat <- as.matrix(dt[, ..cell_ids])
  mode(mat) <- "numeric"
  mat
}

dedup_keep_max_mean <- function(mat, symbols) {
  symbols <- toupper(trimws(as.character(symbols)))
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  m <- rowMeans(mat, na.rm = TRUE)
  dd <- data.table(idx = seq_along(symbols), symbol = symbols, mean_expr = m)
  setorder(dd, symbol, -mean_expr, idx)
  keep_idx <- dd[!duplicated(symbol), idx]
  mat2 <- mat[keep_idx, , drop = FALSE]
  rownames(mat2) <- symbols[keep_idx]
  mat2
}

row_z_by_dataset <- function(mat) {
  mu <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1L, sd, na.rm = TRUE)
  sds[is.na(sds) | sds == 0] <- 1
  z <- (mat - mu) / sds
  z[is.na(z)] <- 0
  z
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
}

collapse_cluster_means <- function(mat, clusters) {
  clu <- as.character(clusters)
  lev <- unique(clu)
  out <- sapply(lev, function(k) {
    idx <- which(clu == k)
    if (length(idx) == 1L) mat[, idx] else rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(out))) {
    out <- matrix(out, ncol = 1L, dimnames = list(rownames(mat), lev))
  }
  out
}

score_broad_classes <- function(cluster_z_mat, marker_dt) {
  classes <- unique(marker_dt$broad_class)
  out <- rbindlist(lapply(classes, function(bc) {
    genes <- unique(marker_dt[broad_class == bc, gene_symbol])
    genes <- intersect(genes, rownames(cluster_z_mat))
    if (length(genes) == 0L) return(NULL)

    score_vec <- colMeans(cluster_z_mat[genes, , drop = FALSE], na.rm = TRUE)

    data.table(
      broad_class = rep(bc, length(score_vec)),
      cluster_id = names(score_vec),
      score = as.numeric(score_vec),
      n_marker_genes_used = rep(length(genes), length(score_vec))
    )
  }), fill = TRUE)
  out
}

infer_broad_class_from_multiome_name <- function(x) {
  x <- as.character(x)
  out <- rep("Other", length(x))
  out[grepl("^EC/Peric\\.$|^EC/Peric$", x)] <- "Vascular"
  out[grepl("mGPC/OPC", x)] <- "OPC_oligo"
  out[grepl("^RG$", x)] <- "RG_like"
  out[grepl("nIPC/GluN1|Cyc\\. Prog\\.", x)] <- "IPC_prog"
  out[grepl("^IN", x)] <- "IN"
  out[grepl("^GluN|^SP$", x)] <- "EN_neurogenic"
  out
}

marker_dt <- rbindlist(list(
  data.table(broad_class = "RG_like", gene_symbol = c("VIM","SOX2","PAX6","HES1","HES5","FABP7","HOPX","FAM107A","HMGA2","ID4")),
  data.table(broad_class = "IPC_prog", gene_symbol = c("EOMES","PPP1R17","HMGB2","NEUROD4","ASCL1","MKI67","TOP2A","CENPF","STMN1","PCNA")),
  data.table(broad_class = "EN_neurogenic", gene_symbol = c("DCX","NEUROD2","NEUROD6","SOX11","STMN2","SATB2","BCL11B","TBR1","NHLH1","ELAVL4")),
  data.table(broad_class = "IN", gene_symbol = c("GAD1","GAD2","DLX1","DLX2","DLX5","DLX6","LHX6","RELN","SST","VIP")),
  data.table(broad_class = "OPC_oligo", gene_symbol = c("OLIG1","OLIG2","PDGFRA","CSPG4","SOX10","NKX2-2","PLP1","MBP")),
  data.table(broad_class = "Vascular", gene_symbol = c("CLDN5","FLT1","KDR","EMCN","RGS5","MCAM","PDGFRB","COL4A1","COL4A2"))
), fill = TRUE)
marker_dt[, gene_symbol := toupper(trimws(gene_symbol))]
marker_dt <- unique(marker_dt)

for (p in c(args$gtf, args$rna_counts, args$rna_meta, args$multiome_counts, args$multiome_meta, args$multiome_cluster_map)) stop_if_missing(p)

log_msg("Starting StepE1 RNA-cluster annotation using multiome reference...")
log_msg("outdir = ", outdir)

# GTF gene map + marker mapping
rowmap <- parse_gtf_gene_map(args$gtf)
marker_map <- merge(marker_dt, rowmap, by = "gene_symbol", allow.cartesian = TRUE)
marker_map <- unique(marker_map[, .(broad_class, gene_symbol, gene_id_stable)])
if (nrow(marker_map) == 0L) stop("No marker genes mapped through GTF.")
marker_gene_ids <- unique(marker_map$gene_id_stable)

safe_fwrite(marker_map, file.path(outdir, "StepE1_marker_gene_mapping.tsv"))
safe_fwrite(marker_map[, .(n_mapped = uniqueN(gene_id_stable)), by = broad_class],
            file.path(outdir, "StepE1_marker_gene_mapping_summary.tsv"))

# RNA
log_msg("Reading GSE162170 RNA metadata...")
rna_meta <- fread(args$rna_meta)
rna_meta[, Cell.ID := as.character(Cell.ID)]
rna_meta[, seurat_clusters := as.character(seurat_clusters)]
rna_meta[, DF_classification := as.character(DF_classification)]
rna_primary <- copy(rna_meta)
if ("DF_classification" %in% names(rna_primary)) {
  rna_primary <- rna_primary[DF_classification == "Singlet"]
}
rna_header_cells <- read_matrix_header_cells(args$rna_counts)
rna_cells <- intersect(rna_primary$Cell.ID, rna_header_cells)
rna_primary <- rna_primary[match(rna_cells, Cell.ID)]
log_msg("RNA primary cells retained: ", nrow(rna_primary))

rna_dt <- extract_rows_by_gene_ids(args$rna_counts, marker_gene_ids, rna_header_cells)
rna_dt[, gene_id_stable := strip_ens_ver(gene_id)]
rna_dt <- merge(rna_dt, unique(marker_map[, .(gene_id_stable, gene_symbol)]), by = "gene_id_stable")
rna_mat_counts <- as_numeric_matrix(rna_dt, rna_cells)
libsize <- as.numeric(rna_primary$RNA.Counts)
libsize[!is.finite(libsize) | libsize <= 0] <- NA_real_
if (any(is.na(libsize))) {
  stop("RNA.Counts missing/invalid in RNA metadata; cannot normalize RNA marker matrix.")
}
rna_mat_norm <- log1p(t(t(rna_mat_counts) / libsize) * 1e4)
rna_mat_sym <- dedup_keep_max_mean(rna_mat_norm, rna_dt$gene_symbol)
rna_cluster_avg <- collapse_cluster_means(rna_mat_sym, rna_primary$seurat_clusters)
rna_cluster_z <- row_z_by_dataset(rna_cluster_avg)
rna_broad_scores <- score_broad_classes(rna_cluster_z, marker_dt)
setnames(rna_broad_scores, c("cluster_id", "score"), c("rna_cluster_id", "broad_class_score"))

# Multiome
log_msg("Reading GSE162170 multiome metadata...")
mo_meta <- fread(args$multiome_meta)
mo_meta[, Cell.ID := as.character(Cell.ID)]
mo_meta[, seurat_clusters := as.character(seurat_clusters)]
mo_meta[, DF_classification := as.character(DF_classification)]
cluster_map <- fread(args$multiome_cluster_map)
setnames(cluster_map, names(cluster_map), gsub("\\.", "_", names(cluster_map)))
# expected columns after StepB: Assay Cluster.ID Cluster.Name
assay_col <- names(cluster_map)[tolower(names(cluster_map)) %in% c("assay")]
id_col <- names(cluster_map)[tolower(names(cluster_map)) %in% c("cluster_id")]
name_col <- names(cluster_map)[tolower(names(cluster_map)) %in% c("cluster_name")]
if (length(id_col) != 1L || length(name_col) != 1L) {
  stop("Could not identify Cluster.ID / Cluster.Name columns in multiome cluster map.")
}
rename_old <- c(id_col, name_col)
rename_new <- c("cluster_id", "cluster_name")
if (length(assay_col) == 1L) {
  rename_old <- c(assay_col, rename_old)
  rename_new <- c("assay", rename_new)
}
setnames(cluster_map, rename_old, rename_new)
if (!("assay" %in% names(cluster_map))) cluster_map[, assay := NA_character_]
cluster_map[, assay := as.character(assay)]
cluster_map[, cluster_id := as.character(cluster_id)]
cluster_map[, cluster_name := as.character(cluster_name)]

# keep only the Multiome RNA mapping if assay labels are present; otherwise deduplicate by cluster_id
if (any(!is.na(cluster_map$assay) & cluster_map$assay != "")) {
  cm0 <- copy(cluster_map)
  cluster_map <- cluster_map[tolower(trimws(assay)) == "multiome rna"]
  if (nrow(cluster_map) == 0L) {
    warning("No rows with assay == 'Multiome RNA' found in cluster map; falling back to deduplicated full map.")
    cluster_map <- unique(cm0[, .(cluster_id, cluster_name)])
  }
}
cluster_map <- unique(cluster_map[, .(cluster_id, cluster_name)])
if (anyDuplicated(cluster_map$cluster_id)) {
  dup_ids <- unique(cluster_map$cluster_id[duplicated(cluster_map$cluster_id)])
  stop("cluster_map still has duplicated cluster_id values after filtering: ", paste(head(dup_ids, 10L), collapse = ", "))
}

mo_primary <- copy(mo_meta)
if ("DF_classification" %in% names(mo_primary)) {
  mo_primary <- mo_primary[DF_classification == "Singlet"]
}
mo_primary <- merge(mo_primary, cluster_map, by.x = "seurat_clusters", by.y = "cluster_id", all.x = TRUE)
mo_primary[, cluster_name := fifelse(is.na(cluster_name) | cluster_name == "", seurat_clusters, cluster_name)]
mo_primary[, broad_class_multiome := infer_broad_class_from_multiome_name(cluster_name)]

mo_header_cells <- read_matrix_header_cells(args$multiome_counts)
mo_cells <- intersect(mo_primary$Cell.ID, mo_header_cells)
mo_primary <- mo_primary[match(mo_cells, Cell.ID)]
log_msg("Multiome primary cells retained: ", nrow(mo_primary))

mo_dt <- extract_rows_by_gene_ids(args$multiome_counts, marker_gene_ids, mo_header_cells)
mo_dt[, gene_id_stable := strip_ens_ver(gene_id)]
mo_dt <- merge(mo_dt, unique(marker_map[, .(gene_id_stable, gene_symbol)]), by = "gene_id_stable")
mo_mat <- as_numeric_matrix(mo_dt, mo_cells)
mo_mat_sym <- dedup_keep_max_mean(mo_mat, mo_dt$gene_symbol)
mo_cluster_avg <- collapse_cluster_means(mo_mat_sym, mo_primary$cluster_name)
mo_cluster_z <- row_z_by_dataset(mo_cluster_avg)
mo_broad_scores <- score_broad_classes(mo_cluster_z, marker_dt)
setnames(mo_broad_scores, c("cluster_id", "score"), c("multiome_cluster_name", "broad_class_score"))

# RNA cluster -> multiome cluster correlation on shared marker genes
common_genes <- intersect(rownames(rna_cluster_z), rownames(mo_cluster_z))
if (length(common_genes) < 10L) stop("Too few common marker genes between RNA and multiome cluster matrices.")
log_msg("Common marker genes used for RNA-vs-multiome cluster correlation: ", length(common_genes))

rna_ids <- colnames(rna_cluster_z)
mo_ids  <- colnames(mo_cluster_z)
cor_dt <- rbindlist(lapply(rna_ids, function(rc) {
  x <- rna_cluster_z[common_genes, rc]
  rbindlist(lapply(mo_ids, function(mc) {
    data.table(
      rna_cluster_id = rc,
      multiome_cluster_name = mc,
      correlation = safe_cor(x, mo_cluster_z[common_genes, mc])
    )
  }))
}))

mo_name_to_broad <- unique(mo_primary[, .(multiome_cluster_name = cluster_name, broad_class_multiome)])
cor_dt <- merge(cor_dt, mo_name_to_broad, by = "multiome_cluster_name", all.x = TRUE)
setorder(cor_dt, rna_cluster_id, -correlation, multiome_cluster_name)
cor_dt[, corr_rank := seq_len(.N), by = rna_cluster_id]

best_exact <- cor_dt[corr_rank == 1L, .(
  rna_cluster_id,
  best_multiome_cluster = multiome_cluster_name,
  best_multiome_broad_class = broad_class_multiome,
  best_multiome_cor = correlation
)]
best_exact <- unique(best_exact, by = "rna_cluster_id")
second_exact <- cor_dt[corr_rank == 2L, .(
  rna_cluster_id,
  second_multiome_cluster = multiome_cluster_name,
  second_multiome_cor = correlation
)]
second_exact <- unique(second_exact, by = "rna_cluster_id")

# broad-class assignment from direct marker scores
setorder(rna_broad_scores, rna_cluster_id, -broad_class_score, broad_class)
rna_broad_scores[, broad_rank := seq_len(.N), by = rna_cluster_id]
rb1 <- rna_broad_scores[broad_rank == 1L, .(
  rna_cluster_id,
  best_direct_broad_class = broad_class,
  best_direct_broad_score = broad_class_score,
  n_marker_genes_used
)]
rb1 <- unique(rb1, by = "rna_cluster_id")
rb2 <- rna_broad_scores[broad_rank == 2L, .(
  rna_cluster_id,
  second_direct_broad_class = broad_class,
  second_direct_broad_score = broad_class_score
)]
rb2 <- unique(rb2, by = "rna_cluster_id")

rna_cluster_sizes <- rna_primary[, .N, by = seurat_clusters]
setnames(rna_cluster_sizes, c("seurat_clusters", "N"), c("rna_cluster_id", "n_cells"))
rna_cluster_sizes <- unique(rna_cluster_sizes, by = "rna_cluster_id")

stopifnot(!anyDuplicated(rb1$rna_cluster_id))
stopifnot(!anyDuplicated(rb2$rna_cluster_id))
stopifnot(!anyDuplicated(best_exact$rna_cluster_id))
stopifnot(!anyDuplicated(second_exact$rna_cluster_id))
stopifnot(!anyDuplicated(rna_cluster_sizes$rna_cluster_id))

anno <- merge(rb1, rb2, by = "rna_cluster_id", all = TRUE)
anno <- merge(anno, best_exact, by = "rna_cluster_id", all = TRUE)
anno <- merge(anno, second_exact, by = "rna_cluster_id", all = TRUE)
anno <- merge(anno, rna_cluster_sizes, by = "rna_cluster_id", all = TRUE)
anno[, direct_margin := best_direct_broad_score - second_direct_broad_score]
anno[, exact_margin := best_multiome_cor - second_multiome_cor]
anno[, final_broad_class := fifelse(
  !is.na(best_multiome_broad_class) & !is.na(best_direct_broad_class) & best_multiome_broad_class == best_direct_broad_class,
  best_direct_broad_class,
  fifelse(!is.na(direct_margin) & direct_margin >= 0.15, best_direct_broad_class,
          fifelse(!is.na(best_multiome_broad_class), best_multiome_broad_class, best_direct_broad_class))
)]
anno[, annotation_confidence := fifelse(
  !is.na(best_multiome_broad_class) & !is.na(best_direct_broad_class) & best_multiome_broad_class == best_direct_broad_class &
    !is.na(exact_margin) & exact_margin >= 0.05 & !is.na(direct_margin) & direct_margin >= 0.10,
  "high",
  fifelse(!is.na(best_multiome_broad_class) & !is.na(best_direct_broad_class) & best_multiome_broad_class == best_direct_broad_class,
          "medium",
          fifelse(!is.na(direct_margin) & direct_margin >= 0.15, "medium", "low"))
)]
anno[, final_label := fifelse(
  annotation_confidence == "high" & !is.na(best_multiome_cluster),
  paste0(final_broad_class, "__", best_multiome_cluster),
  final_broad_class
)]
setorder(anno, final_broad_class, -n_cells, rna_cluster_id)

# outputs
safe_fwrite(anno, file.path(outdir, "StepE1_GSE162170_RNA_cluster_annotation.tsv"))
safe_fwrite(rna_broad_scores, file.path(outdir, "StepE1_GSE162170_RNA_cluster_broadclass_marker_scores.tsv"))
safe_fwrite(cor_dt, file.path(outdir, "StepE1_GSE162170_RNA_to_multiome_cluster_correlation.tsv.gz"))
safe_fwrite(mo_name_to_broad, file.path(outdir, "StepE1_GSE162170_multiome_cluster_to_broadclass.tsv"))

run_summary <- data.table(
  rna_primary_cells = nrow(rna_primary),
  multiome_primary_cells = nrow(mo_primary),
  n_marker_symbols = uniqueN(marker_dt$gene_symbol),
  n_marker_gene_ids = uniqueN(marker_map$gene_id_stable),
  n_common_genes_for_correlation = length(common_genes),
  n_rna_clusters = uniqueN(rna_primary$seurat_clusters),
  n_multiome_clusters = uniqueN(mo_primary$cluster_name)
)
safe_fwrite(run_summary, file.path(outdir, "StepE1_run_summary.tsv"))

log_msg("Finished StepE1 successfully.")
