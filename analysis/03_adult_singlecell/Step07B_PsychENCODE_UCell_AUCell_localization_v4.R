#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

patch_matrixStats_useNames <- function() {
  if (!requireNamespace("matrixStats", quietly = TRUE)) return(invisible(FALSE))

  make_patch <- function(orig_fun) {
    force(orig_fun)
    function(...) {
      args <- list(...)
      if (!("useNames" %in% names(args)) || (length(args$useNames) == 1 && isTRUE(is.na(args$useNames)))) {
        args$useNames <- FALSE
      }
      do.call(orig_fun, args)
    }
  }

  patched_any <- FALSE
  for (ns_name in c("matrixStats", "UCell", "AUCell")) {
    if (!isNamespaceLoaded(ns_name) && ns_name != "matrixStats") {
      suppressWarnings(try(loadNamespace(ns_name), silent = TRUE))
    }
    if (!isNamespaceLoaded(ns_name)) next
    for (fun_name in c("colRanks", "rowRanks")) {
      suppressWarnings(try({
        orig_fun <- get(fun_name, envir = asNamespace(ns_name), inherits = FALSE)
        if (is.function(orig_fun)) {
          assignInNamespace(fun_name, make_patch(orig_fun), ns = ns_name)
          patched_any <- TRUE
        }
      }, silent = TRUE))
    }
  }
  invisible(patched_any)
}

if (!requireNamespace("UCell", quietly = TRUE)) stop("Package 'UCell' is not installed.")
if (!requireNamespace("AUCell", quietly = TRUE)) stop("Package 'AUCell' is not installed.")
patched_ok <- patch_matrixStats_useNames()

parse_args <- function(x) {
  out <- list(base = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex", batch_size = NA_integer_)
  i <- 1
  while (i <= length(x)) {
    key <- x[i]
    val <- if (i < length(x)) x[i + 1] else NA_character_
    if (key == "--base") out$base <- val
    if (key == "--batch_size") out$batch_size <- as.integer(val)
    i <- i + 2
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
BASE <- args$base
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

safe_fwrite <- function(x, file) fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")

dedup_keep_max <- function(mat, symbols) {
  symbols <- toupper(trimws(symbols))
  keep_valid <- !is.na(symbols) & symbols != ""
  mat <- mat[keep_valid, , drop = FALSE]
  symbols <- symbols[keep_valid]
  rs <- Matrix::rowSums(mat)
  dt <- data.table(idx = seq_along(symbols), symbol = symbols, rs = rs)
  setorder(dt, symbol, -rs, idx)
  keep_idx <- dt[!duplicated(symbol), idx]
  mat2 <- mat[keep_idx, , drop = FALSE]
  rownames(mat2) <- symbols[keep_idx]
  mat2
}

read_program_genes <- function(program_file, target_program) {
  if (!file.exists(program_file)) stop("Missing program file: ", program_file)
  dt <- fread(program_file)
  nms <- tolower(names(dt))
  if (!("program_name" %in% nms)) stop("program_name column missing in program file.")
  pcol <- names(dt)[match("program_name", nms)]
  g_cands <- c("gene_symbol", "symbol", "gene", "standardized_gene_symbol")
  g_hit <- g_cands[g_cands %in% nms]
  if (length(g_hit) == 0) stop("gene symbol column missing in program file.")
  gcol <- names(dt)[match(g_hit[1], nms)]

  keep_idx <- dt[[pcol]] %in% target_program
  out <- unique(toupper(trimws(dt[[gcol]][keep_idx])))
  out[!is.na(out) & out != ""]
}

score_in_batches <- function(mat, meta, genesets_use, batch_size, dataset_tag, prefix, log_msg, ncores_use = 1L) {
  n_cells <- ncol(mat)
  batches <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / batch_size))
  log_msg("Scoring ", dataset_tag, " with ", length(batches), " batches; batch_size=", batch_size, "; ncores=", ncores_use)

  score_batches <- vector("list", length(batches))

  for (b in seq_along(batches)) {
    idx <- batches[[b]]
    submat <- mat[, idx, drop = FALSE]
    log_msg("Batch ", b, "/", length(batches), ": ", ncol(submat), " cells")

    uc <- UCell::ScoreSignatures_UCell(
      submat,
      features = genesets_use,
      maxRank = min(1500, nrow(submat)),
      ncores = ncores_use
    )
    uc_dt <- as.data.table(uc)
    uc_dt[, cell := rownames(uc)]
    uc_long <- melt(uc_dt, id.vars = "cell", variable.name = "program_name", value.name = "score")
    uc_long[, program_name := sub("_UCell$", "", program_name)]
    uc_long[, method := "UCell"]

    rankings <- AUCell::AUCell_buildRankings(
      submat,
      plotStats = FALSE,
      splitByBlocks = FALSE,
      verbose = FALSE
    )
    auc_obj <- AUCell::AUCell_calcAUC(
      genesets_use,
      rankings,
      aucMaxRank = max(50, ceiling(0.05 * nrow(submat))),
      nCores = ncores_use,
      verbose = FALSE
    )
    auc_mat <- t(as.matrix(AUCell::getAUC(auc_obj)))
    auc_dt <- as.data.table(auc_mat)
    auc_dt[, cell := rownames(auc_mat)]
    auc_long <- melt(auc_dt, id.vars = "cell", variable.name = "program_name", value.name = "score")
    auc_long[, method := "AUCell"]

    score_batches[[b]] <- rbindlist(list(uc_long, auc_long), use.names = TRUE, fill = TRUE)
    rm(submat, uc, uc_dt, uc_long, rankings, auc_obj, auc_mat, auc_dt, auc_long)
    gc()
  }

  scores_long <- rbindlist(score_batches, use.names = TRUE, fill = TRUE)
  score_annot <- merge(scores_long, meta, by = "cell", all.x = TRUE)

  safe_fwrite(score_annot, file.path(TABDIR, paste0(prefix, "_UCell_AUCell_cellScores.tsv.gz")))

  broad_loc <- score_annot[, .(
    n_cells = uniqueN(cell),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ), by = .(method, program_name, broad_class)]
  setorder(broad_loc, method, program_name, -mean_score)
  safe_fwrite(broad_loc, file.path(TABDIR, paste0(prefix, "_UCell_AUCell_broadClass_localization.tsv")))

  sub_loc <- score_annot[, .(
    n_cells = uniqueN(cell),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ), by = .(method, program_name, broad_class, subtype)]
  setorder(sub_loc, method, program_name, -mean_score)
  safe_fwrite(sub_loc, file.path(TABDIR, paste0(prefix, "_UCell_AUCell_subtype_localization.tsv")))

  broad_rank <- copy(broad_loc)
  broad_rank[, rank_by_mean := frank(-mean_score, ties.method = "min"), by = .(method, program_name)]
  setorder(broad_rank, method, program_name, rank_by_mean, broad_class)
  safe_fwrite(broad_rank, file.path(TABDIR, paste0(prefix, "_UCell_AUCell_broadClass_ranked.tsv")))

  donor_broad <- score_annot[, .(
    n_cells = uniqueN(cell),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE)
  ), by = .(method, program_name, broad_class, individual, diagnosis, is_asd)]
  safe_fwrite(donor_broad, file.path(TABDIR, paste0(prefix, "_UCell_AUCell_donor_broadClass_scores.tsv.gz")))

  wide_loc <- dcast(broad_loc, program_name + broad_class ~ method, value.var = "mean_score")
  method_cor <- if (all(c("UCell", "AUCell") %in% names(wide_loc))) {
    wide_loc[, .(
      spearman_rho = suppressWarnings(cor(UCell, AUCell, method = "spearman")),
      pearson_r = suppressWarnings(cor(UCell, AUCell, method = "pearson"))
    ), by = program_name]
  } else {
    data.table(program_name = unique(broad_loc$program_name), spearman_rho = NA_real_, pearson_r = NA_real_)
  }
  safe_fwrite(method_cor, file.path(METADIR, paste0(prefix, "_UCell_AUCell_method_concordance.tsv")))

  list(
    score_annot = score_annot,
    broad_rank = broad_rank,
    method_cor = method_cor,
    n_batches = length(batches)
  )
}


DATASET_TAG <- "PsychENCODE"
PREFIX <- "95_PsychENCODE"
LOGFILE <- file.path(LOGDIR, "Step07_PsychENCODE_UCell_AUCell_localization_v4.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOGFILE, append = TRUE)
}

if (is.na(args$batch_size)) args$batch_size <- 12000L
BATCH_SIZE <- as.integer(args$batch_size)
NCORES_USE <- max(1L, min(8L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))))

DATA_BASE <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024"
META_FILE <- file.path(DATA_BASE, "meta.tsv")
MTX_FILE  <- file.path(DATA_BASE, "counts_matrix.mtx.gz")
FEAT_FILE <- file.path(DATA_BASE, "counts_features.tsv.gz")
BARC_FILE <- file.path(DATA_BASE, "counts_barcodes.tsv.gz")
PROGRAM_FILE <- file.path(TABDIR, "20_midPrenatal_core_programs.tsv")

log_msg("matrixStats useNames compatibility patch applied: ", patched_ok)
log_msg("Reading PsychENCODE meta: ", META_FILE)
meta <- fread(META_FILE)
req <- c("Cell_ID", "annotation", "Subtype", "individual_ID", "Diagnosis")
if (!all(req %in% names(meta))) stop("PsychENCODE meta missing required columns.")
meta[, cell := as.character(Cell_ID)]
meta[, broad_class := as.character(annotation)]
meta[, subtype := as.character(Subtype)]
meta[, individual := as.character(individual_ID)]
meta[, diagnosis := fifelse(Diagnosis == "CTL", "Control",
                            fifelse(Diagnosis == "ASD", "ASD", NA_character_))]
meta[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
meta <- meta[!is.na(broad_class) & !is.na(diagnosis)]
meta[, is_asd := diagnosis == "ASD"]
meta_keep <- meta[, .(cell, broad_class, subtype, individual, diagnosis, is_asd)]

for (f in c(MTX_FILE, FEAT_FILE, BARC_FILE)) if (!file.exists(f)) stop("Missing file: ", f)

fdt <- fread(cmd = paste("zcat", shQuote(FEAT_FILE)), header = FALSE, sep = "\t")
symbols <- if (ncol(fdt) >= 2) fdt[[2]] else fdt[[1]]
barcodes <- fread(cmd = paste("zcat", shQuote(BARC_FILE)), header = FALSE)$V1

log_msg("Reading PsychENCODE sparse matrix: ", MTX_FILE)
con <- gzfile(MTX_FILE, open = "rt")
mat <- readMM(con)
close(con)
mat <- as(mat, "CsparseMatrix")
colnames(mat) <- as.character(barcodes)

log_msg("Deduplicating PsychENCODE gene symbols by max-row-sum representative")
mat <- dedup_keep_max(mat, symbols)

common_cells <- intersect(colnames(mat), meta_keep$cell)
log_msg("PsychENCODE matched cells: ", length(common_cells))
meta_keep <- meta_keep[match(common_cells, cell)]
mat <- mat[, common_cells, drop = FALSE]

genes_sfari <- read_program_genes(PROGRAM_FILE, "SFARI_all")
genes_top20 <- read_program_genes(PROGRAM_FILE, "midPrenatal_SFARI_top20")
log_msg("Raw program sizes before matrix mapping: SFARI_all=", length(genes_sfari), "; midPrenatal_SFARI_top20=", length(genes_top20))

GENESETS <- list(
  SFARI_all = genes_sfari,
  midPrenatal_SFARI_top20 = genes_top20
)

gene_map_summary <- data.table(
  program_name = names(GENESETS),
  n_genes_program = vapply(GENESETS, length, integer(1)),
  n_genes_in_matrix = vapply(GENESETS, function(gs) sum(gs %in% rownames(mat)), integer(1))
)
safe_fwrite(gene_map_summary, file.path(METADIR, paste0(PREFIX, "_geneSet_mapping_summary.tsv")))

GENESETS_USE <- lapply(GENESETS, function(gs) gs[gs %in% rownames(mat)])
if (any(vapply(GENESETS_USE, length, integer(1)) == 0)) stop("At least one gene set has 0 mapped genes.")

res <- score_in_batches(mat, meta_keep, GENESETS_USE, BATCH_SIZE, DATASET_TAG, PREFIX, log_msg, NCORES_USE)

run_summary <- rbindlist(list(
  data.table(section = DATASET_TAG, metric = "n_unique_cells", value = uniqueN(res$score_annot$cell)),
  data.table(section = DATASET_TAG, metric = "n_long_rows_scores", value = nrow(res$score_annot)),
  data.table(section = DATASET_TAG, metric = "n_genes_matrix", value = nrow(mat)),
  data.table(section = DATASET_TAG, metric = "n_broad_classes", value = uniqueN(res$score_annot$broad_class)),
  data.table(section = DATASET_TAG, metric = "n_subtypes", value = uniqueN(res$score_annot$subtype)),
  data.table(section = DATASET_TAG, metric = "batch_size", value = BATCH_SIZE),
  data.table(section = DATASET_TAG, metric = "n_batches", value = res$n_batches)
), fill = TRUE)

for (k in names(GENESETS_USE)) {
  run_summary <- rbind(run_summary, data.table(section = "GeneSet", metric = paste0(k, "_n_genes_in_matrix"), value = length(GENESETS_USE[[k]])), fill = TRUE)
}

best_broad <- res$broad_rank[rank_by_mean == 1]
for (i in seq_len(nrow(best_broad))) {
  run_summary <- rbind(
    run_summary,
    data.table(
      section = paste0("TopBroad_", best_broad$method[i]),
      metric = best_broad$program_name[i],
      value = paste0(best_broad$broad_class[i], " (mean=", formatC(best_broad$mean_score[i], digits = 4, format = "fg"), ")")
    ),
    fill = TRUE
  )
}

safe_fwrite(run_summary, file.path(METADIR, paste0(PREFIX, "_UCell_AUCell_run_summary.tsv")))
log_msg("Step07 scoring completed successfully for ", DATASET_TAG)
