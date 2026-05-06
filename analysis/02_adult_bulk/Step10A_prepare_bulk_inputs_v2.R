#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

option_list <- list(
  make_option("--outdir", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
              help = "Step output directory [default %default]"),
  make_option("--psych_dir", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024",
              help = "PsychENCODE directory containing meta and matrix/feature/barcode files [default %default]"),
  make_option("--psych_meta", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/meta.tsv",
              help = "PsychENCODE meta.tsv path [default %default]"),
  make_option("--gse102741_raw", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/GSE102741_raw_counts_GRCh38.p13_NCBI.tsv.gz",
              help = "GSE102741 raw counts path [default %default]"),
  make_option("--gse64018_raw", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/GSE64018_countlevel_12asd_12ctl.txt.gz",
              help = "GSE64018 raw counts path [default %default]"),
  make_option("--gse64018_mapping", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/meta/120_GSE64018_explicit_sample_mapping.tsv",
              help = "Explicit GSE64018 sample mapping TSV [default %default]"),
  make_option("--gandal_raw_rdata", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/RawData_GeneCounts.RData",
              help = "Gandal raw counts RData [default %default]"),
  make_option("--gandal_meta_rdata", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData",
              help = "Gandal normalized/meta RData [default %default]"),
  make_option("--force_gandal_normalized", action = "store_true", default = FALSE,
              help = "Force use normalized Gandal expression if raw extraction fails [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------
# helpers
# ------------------------------
base_dir <- normalizePath(opt$outdir, mustWork = FALSE)
input_dir <- file.path(base_dir, "inputs")
meta_dir <- file.path(base_dir, "meta")
log_dir  <- file.path(base_dir, "logs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "Step10A_prepare_bulk_inputs_v2.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  if (grepl("\\.gz$", file)) {
    fwrite(x, file = file, sep = "\t", quote = FALSE, compress = "gzip")
  } else {
    fwrite(x, file = file, sep = "\t", quote = FALSE)
  }
}

normalize_id <- function(x) {
  x <- as.character(x)
  x <- gsub('^"|"$', '', x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x
}

normalize_sample_id <- function(x) {
  x <- normalize_id(x)
  x <- gsub("[.]", "-", x)
  x
}

choose_existing <- function(paths) {
  paths <- unique(paths)
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) return(NA_character_)
  ok[1]
}

find_one_file <- function(root, patterns) {
  all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
  hits <- all_files[Reduce(`|`, lapply(patterns, function(p) grepl(p, basename(all_files), ignore.case = TRUE)))]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

looks_like_ensembl <- function(x) mean(grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", x))
looks_like_entrez  <- function(x) mean(grepl("^[0-9]+$", x))
looks_like_symbol  <- function(x) mean(grepl("^[A-Za-z][A-Za-z0-9._-]*$", x))

map_ids_to_symbol <- function(ids) {
  ids <- normalize_id(ids)
  ids[ids %in% c("", "NA", "NaN")] <- NA_character_
  out <- rep(NA_character_, length(ids))
  idx <- !is.na(ids)
  if (!any(idx)) return(out)

  x <- ids[idx]
  if (looks_like_ensembl(x) >= 0.5) {
    keys <- sub("\\..*$", "", x)
    sym <- suppressWarnings(mapIds(org.Hs.eg.db, keys = keys, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
    out[idx] <- unname(sym)
  } else if (looks_like_entrez(x) >= 0.5) {
    sym <- suppressWarnings(mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"))
    out[idx] <- unname(sym)
  } else {
    out[idx] <- x
  }
  out <- normalize_id(out)
  out[out == ""] <- NA_character_
  out
}

collapse_by_gene <- function(dt_or_df, gene_col = "gene_symbol") {
  dt <- as.data.table(dt_or_df)
  stopifnot(gene_col %in% names(dt))
  sample_cols <- setdiff(names(dt), gene_col)
  for (j in sample_cols) set(dt, j = j, value = as.numeric(dt[[j]]))
  dt <- dt[!is.na(get(gene_col)) & get(gene_col) != ""]
  collapsed <- dt[, lapply(.SD, sum, na.rm = TRUE), by = gene_col, .SDcols = sample_cols]
  setorderv(collapsed, gene_col)
  collapsed
}

write_standardized_pair <- function(expr_dt, meta_dt, prefix, input_dir) {
  expr_file <- file.path(input_dir, sprintf("%s.standardized_expr.tsv.gz", prefix))
  meta_file <- file.path(input_dir, sprintf("%s.standardized_meta.tsv", prefix))
  safe_fwrite(expr_dt, expr_file)
  safe_fwrite(meta_dt, meta_file)
  list(expr = expr_file, meta = meta_file)
}

# ------------------------------
# GSE102741
# ------------------------------
prepare_gse102741 <- function(raw_file, input_dir) {
  log_msg("Preparing GSE102741 from: ", raw_file)
  dt <- fread(raw_file)
  first_col <- names(dt)[1]
  setnames(dt, first_col, "feature_id")
  dt[, feature_id := normalize_id(feature_id)]
  dt[, gene_symbol := map_ids_to_symbol(feature_id)]
  expr <- collapse_by_gene(dt[, c("gene_symbol", setdiff(names(dt), c("feature_id", "gene_symbol"))), with = FALSE])

  samples <- setdiff(names(expr), "gene_symbol")
  meta <- data.table(sample = samples)
  meta[, diagnosis := fifelse(sample %in% paste0("GSM", 2745936:2745948), "ASD", "Control")]
  setorderv(meta, "sample")

  out <- write_standardized_pair(expr, meta, "GSE102741", input_dir)
  list(expr_dt = expr, meta_dt = meta, files = out)
}

# ------------------------------
# GSE64018
# ------------------------------
prepare_gse64018 <- function(raw_file, mapping_file, input_dir) {
  log_msg("Preparing GSE64018 from: ", raw_file)
  dt <- fread(raw_file, check.names = FALSE)
  if (!nrow(dt)) stop("GSE64018 raw file is empty.")
  first_name <- names(dt)[1]
  setnames(dt, first_name, "feature_id")
  dt[, feature_id := normalize_id(feature_id)]
  sample_cols <- setdiff(names(dt), "feature_id")
  clean_cols <- normalize_sample_id(sample_cols)
  setnames(dt, old = sample_cols, new = clean_cols)
  dt[, gene_symbol := map_ids_to_symbol(feature_id)]
  expr <- collapse_by_gene(dt[, c("gene_symbol", setdiff(names(dt), c("feature_id", "gene_symbol"))), with = FALSE])

  map_dt <- fread(mapping_file)
  map_dt[, sample_norm := normalize_sample_id(sample_norm)]
  map_dt[, diagnosis := fifelse(toupper(diagnosis) %in% c("CTL", "CTRL", "CONTROL"), "Control", "ASD")]
  meta <- unique(map_dt[, .(sample = sample_norm, diagnosis)])
  meta <- meta[sample %in% setdiff(names(expr), "gene_symbol")]
  setorderv(meta, "sample")

  out <- write_standardized_pair(expr, meta, "GSE64018", input_dir)
  list(expr_dt = expr, meta_dt = meta, files = out)
}

# ------------------------------
# Gandal extraction helpers
# ------------------------------
make_numeric_matrix <- function(obj) {
  if (inherits(obj, "data.table") || is.data.frame(obj)) {
    x <- as.data.frame(obj)
    if (ncol(x) >= 2 && !is.numeric(x[[1]]) && all(vapply(x[, -1, drop = FALSE], is.numeric, logical(1)))) {
      rn <- normalize_id(x[[1]])
      x <- as.matrix(x[, -1, drop = FALSE])
      rownames(x) <- rn
      return(x)
    }
    if (all(vapply(x, is.numeric, logical(1)))) return(as.matrix(x))
    return(NULL)
  }
  if (is.matrix(obj) && is.numeric(obj)) return(obj)
  NULL
}

score_expr_candidate <- function(m) {
  if (is.null(m)) return(-Inf)
  nr <- nrow(m); nc <- ncol(m)
  if (is.null(nr) || is.null(nc)) return(-Inf)
  if (nr < 500 || nc < 20) return(-Inf)
  score <- log10(nr * nc)
  if (nr > nc) score <- score + 2
  if (!is.null(rownames(m))) {
    rn <- head(rownames(m), min(1000, nr))
    score <- score + 3 * looks_like_symbol(rn) + 3 * looks_like_ensembl(rn) + 2 * looks_like_entrez(rn)
  }
  score
}

extract_best_expr_from_rdata <- function(rdata_path) {
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  nms <- ls(e)
  mats <- lapply(nms, function(nm) make_numeric_matrix(get(nm, envir = e)))
  scores <- vapply(mats, score_expr_candidate, numeric(1))
  if (all(!is.finite(scores))) return(NULL)
  best <- which.max(scores)
  m <- mats[[best]]
  nm <- nms[best]

  # orient gene x sample
  if (!is.null(colnames(m)) && looks_like_ensembl(colnames(m)) > looks_like_ensembl(rownames(m)) &&
      ncol(m) > nrow(m)) {
    m <- t(m)
  }
  if (!is.null(colnames(m)) && looks_like_symbol(colnames(m)) > looks_like_symbol(rownames(m)) && ncol(m) > nrow(m)) {
    m <- t(m)
  }
  list(name = nm, matrix = m)
}

extract_best_meta_from_rdata <- function(rdata_path, sample_ids = NULL) {
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  nms <- ls(e)
  dfs <- nms[vapply(nms, function(nm) is.data.frame(get(nm, envir = e)) || inherits(get(nm, envir = e), "data.table"), logical(1))]
  if (length(dfs) == 0) return(NULL)

  best_score <- -Inf
  best_dt <- NULL
  best_name <- NULL
  for (nm in dfs) {
    dt <- as.data.table(get(nm, envir = e))
    sc <- 0
    if (nrow(dt) >= 20) sc <- sc + log10(nrow(dt) + 1)
    dcols <- names(dt)[grepl("diagn|dx|group|phenotype|condition|status", names(dt), ignore.case = TRUE)]
    if (length(dcols)) sc <- sc + 5
    if (!is.null(sample_ids)) {
      for (cc in names(dt)) {
        val <- normalize_sample_id(dt[[cc]])
        sc <- max(sc, sum(val %in% normalize_sample_id(sample_ids)) / max(1, length(sample_ids)) * 8)
      }
      if (!is.null(rownames(dt))) {
        rv <- normalize_sample_id(rownames(dt))
        sc <- max(sc, sum(rv %in% normalize_sample_id(sample_ids)) / max(1, length(sample_ids)) * 8)
      }
    }
    if (sc > best_score) {
      best_score <- sc
      best_dt <- dt
      best_name <- nm
    }
  }
  if (is.null(best_dt)) return(NULL)
  list(name = best_name, meta = best_dt)
}

standardize_gandal_meta <- function(meta_dt, sample_ids) {
  dt <- copy(as.data.table(meta_dt))
  sample_ids_norm <- normalize_sample_id(sample_ids)

  # detect sample column
  best_col <- NULL
  best_match <- -1
  for (cc in names(dt)) {
    vals <- normalize_sample_id(dt[[cc]])
    m <- sum(vals %in% sample_ids_norm)
    if (m > best_match) {
      best_match <- m
      best_col <- cc
    }
  }
  if (!is.null(rownames(dt))) {
    rn_match <- sum(normalize_sample_id(rownames(dt)) %in% sample_ids_norm)
    if (rn_match > best_match) {
      dt[, sample := normalize_sample_id(rownames(meta_dt))]
      best_col <- "sample"
      best_match <- rn_match
    }
  }
  if (is.null(best_col) || best_match == 0) stop("Could not identify Gandal sample column from metadata.")

  dcols <- names(dt)[grepl("diagn|dx|group|phenotype|condition|status", names(dt), ignore.case = TRUE)]
  if (!length(dcols)) stop("Could not identify diagnosis column in Gandal metadata.")
  dcol <- dcols[1]

  dt[, sample := normalize_sample_id(get(best_col))]
  dx_raw <- toupper(normalize_id(dt[[dcol]]))
  dt[, diagnosis := fifelse(dx_raw %in% c("CTL", "CTRL", "CONTROL", "CON"), "Control",
                            fifelse(dx_raw %in% c("ASD", "AUTISM", "CASE"), "ASD", NA_character_))]
  keep_cols <- c("sample", "diagnosis", intersect(c("region", "Region", "age", "Age", "sex", "Sex", "PMI", "pmi", "RIN", "rin"), names(dt)))
  out <- unique(dt[, ..keep_cols])
  out <- out[sample %in% sample_ids_norm & !is.na(diagnosis)]
  setorderv(out, "sample")
  out
}

prepare_gandal <- function(raw_rdata, meta_rdata, input_dir, force_normalized = FALSE) {
  log_msg("Preparing Gandal from raw/meta RData.")
  expr_obj <- tryCatch(extract_best_expr_from_rdata(raw_rdata), error = function(e) NULL)
  if (is.null(expr_obj) && force_normalized) {
    expr_obj <- extract_best_expr_from_rdata(meta_rdata)
  }
  if (is.null(expr_obj)) stop("Could not extract Gandal expression matrix from provided RData files.")
  m <- expr_obj$matrix

  if (is.null(rownames(m))) stop("Gandal expression matrix has no rownames for gene IDs.")
  gene_ids <- normalize_id(rownames(m))
  gene_symbol <- map_ids_to_symbol(gene_ids)
  expr_dt <- as.data.table(m)
  expr_dt[, gene_symbol := gene_symbol]
  setcolorder(expr_dt, c("gene_symbol", setdiff(names(expr_dt), "gene_symbol")))
  old_cols <- setdiff(names(expr_dt), "gene_symbol")
  new_cols <- normalize_sample_id(old_cols)
  setnames(expr_dt, old = old_cols, new = new_cols)
  expr <- collapse_by_gene(expr_dt)

  meta_obj <- extract_best_meta_from_rdata(meta_rdata, sample_ids = setdiff(names(expr), "gene_symbol"))
  if (is.null(meta_obj)) stop("Could not extract Gandal metadata from provided RData.")
  meta <- standardize_gandal_meta(meta_obj$meta, sample_ids = setdiff(names(expr), "gene_symbol"))

  keep_samples <- intersect(setdiff(names(expr), "gene_symbol"), meta$sample)
  expr <- expr[, c("gene_symbol", keep_samples), with = FALSE]
  meta <- meta[sample %in% keep_samples]
  setorderv(meta, "sample")

  out <- write_standardized_pair(expr, meta, "Gandal2022", input_dir)
  list(expr_dt = expr, meta_dt = meta, files = out,
       expr_source = expr_obj$name, meta_source = meta_obj$name)
}

# ------------------------------
# Run
# ------------------------------
log_msg("Starting Step10A bulk input preparation.")
log_msg("PsychENCODE meta path: ", opt$psych_meta)
log_msg("GSE102741 raw path: ", opt$gse102741_raw)
log_msg("GSE64018 raw path: ", opt$gse64018_raw)
log_msg("GSE64018 explicit mapping path: ", opt$gse64018_mapping)
log_msg("Gandal raw RData path: ", opt$gandal_raw_rdata)
log_msg("Gandal meta RData path: ", opt$gandal_meta_rdata)

res102 <- prepare_gse102741(opt$gse102741_raw, input_dir)
res640 <- prepare_gse64018(opt$gse64018_raw, opt$gse64018_mapping, input_dir)
resGan <- prepare_gandal(opt$gandal_raw_rdata, opt$gandal_meta_rdata, input_dir, opt$force_gandal_normalized)

summary_dt <- rbindlist(list(
  data.table(cohort = "GSE102741", expr_file = res102$files$expr, meta_file = res102$files$meta,
             n_genes = nrow(res102$expr_dt), n_samples = nrow(res102$meta_dt)),
  data.table(cohort = "GSE64018", expr_file = res640$files$expr, meta_file = res640$files$meta,
             n_genes = nrow(res640$expr_dt), n_samples = nrow(res640$meta_dt)),
  data.table(cohort = "Gandal2022", expr_file = resGan$files$expr, meta_file = resGan$files$meta,
             n_genes = nrow(resGan$expr_dt), n_samples = nrow(resGan$meta_dt),
             expr_source = resGan$expr_source, meta_source = resGan$meta_source)
), fill = TRUE)

safe_fwrite(summary_dt, file.path(meta_dir, "Step10A_prepare_bulk_inputs_summary.tsv"))

run_summary <- rbindlist(list(
  data.table(section = "GSE102741", metric = c("n_genes", "n_samples", "n_asd", "n_control"),
             value = c(nrow(res102$expr_dt), nrow(res102$meta_dt), sum(res102$meta_dt$diagnosis == "ASD"), sum(res102$meta_dt$diagnosis == "Control"))),
  data.table(section = "GSE64018", metric = c("n_genes", "n_samples", "n_asd", "n_control"),
             value = c(nrow(res640$expr_dt), nrow(res640$meta_dt), sum(res640$meta_dt$diagnosis == "ASD"), sum(res640$meta_dt$diagnosis == "Control"))),
  data.table(section = "Gandal2022", metric = c("n_genes", "n_samples", "n_asd", "n_control", "expr_source", "meta_source"),
             value = c(nrow(resGan$expr_dt), nrow(resGan$meta_dt), sum(resGan$meta_dt$diagnosis == "ASD"), sum(resGan$meta_dt$diagnosis == "Control"), resGan$expr_source, resGan$meta_source))
), fill = TRUE)

safe_fwrite(run_summary, file.path(meta_dir, "144a_Step10A_prepare_run_summary.tsv"))
log_msg("Step10A_prepare_bulk_inputs_v2 completed successfully.")
