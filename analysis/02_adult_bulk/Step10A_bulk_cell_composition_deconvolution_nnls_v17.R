#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(nnls)
})

# ----------------------------- #
# Basic argument parsing
# ----------------------------- #
parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- sub("^--", "", key)
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- x[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

BASE <- if (!is.null(args$base)) args$base else "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare"
PSY_DIR <- if (!is.null(args$psych_dir)) args$psych_dir else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024"
GANDAL_RAW <- if (!is.null(args$gandal_raw)) args$gandal_raw else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/RawData_GeneCounts.RData"
GANDAL_META <- if (!is.null(args$gandal_meta)) args$gandal_meta else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"
GANDAL_ANNO <- if (!is.null(args$gandal_anno)) args$gandal_anno else "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/Gencode_v3c_summarized_to_genes/rows_metadata.csv"
PROGRAM_FILE <- if (!is.null(args$program_file)) args$program_file else "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/20_midPrenatal_core_programs.tsv"
MARKERS_PER_CLASS <- if (!is.null(args$markers_per_class)) as.integer(args$markers_per_class) else 150L

dir.create(BASE, recursive = TRUE, showWarnings = FALSE)
for (d in c(file.path(BASE, "logs"), file.path(BASE, "tables"), file.path(BASE, "meta"), file.path(BASE, "inputs"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOGFILE <- file.path(BASE, "logs", "Step10A_bulk_cell_composition_deconvolution_nnls_v16.log")
if (file.exists(LOGFILE)) file.remove(LOGFILE)

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOGFILE, append = TRUE)
}

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

normalize_diagnosis <- function(x) {
  x <- trimws(as.character(x))
  x_up <- toupper(x)
  out <- ifelse(x_up %in% c("ASD", "AUTISM", "AUTISM SPECTRUM DISORDER"), "ASD",
         ifelse(x_up %in% c("CTL", "CTRL", "CONTROL", "NEUROTYPICAL", "NORMAL"), "Control", NA_character_))
  out
}

read_dt_no_header <- function(path) {
  fread(path, header = FALSE)
}

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

clean_ensembl <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x <- sub("(_[0-9]+)$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

calc_auc <- function(asd_vals, ctl_vals) {
  asd_vals <- as.numeric(asd_vals)
  ctl_vals <- as.numeric(ctl_vals)
  if (length(asd_vals) == 0L || length(ctl_vals) == 0L) return(NA_real_)
  r <- rank(c(asd_vals, ctl_vals), ties.method = "average")
  n1 <- length(asd_vals)
  n0 <- length(ctl_vals)
  u1 <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  as.numeric(u1 / (n1 * n0))
}

align_meta_expr <- function(expr, meta, sample_col = "sample") {
  common <- intersect(colnames(expr), meta[[sample_col]])
  meta2 <- meta[match(common, meta[[sample_col]])]
  expr2 <- expr[, common, drop = FALSE]
  list(expr = expr2, meta = meta2)
}

collapse_matrix_by_symbol <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  grp <- factor(symbols, levels = unique(symbols))
  res <- rowsum(as.matrix(mat), group = grp, reorder = FALSE)
  res
}

collapse_sparse_by_symbol <- function(mtx, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mtx <- mtx[keep, , drop = FALSE]
  symbols <- symbols[keep]
  grp <- factor(symbols, levels = unique(symbols))
  P <- sparseMatrix(
    i = as.integer(grp),
    j = seq_along(grp),
    x = 1,
    dims = c(nlevels(grp), length(grp))
  )
  res <- P %*% mtx
  rownames(res) <- levels(grp)
  res
}

aggregate_sparse_by_group <- function(mtx, groups) {
  keep <- !is.na(groups) & groups != ""
  mtx <- mtx[, keep, drop = FALSE]
  groups <- groups[keep]
  grp <- factor(groups, levels = unique(groups))
  G <- sparseMatrix(
    i = seq_along(grp),
    j = as.integer(grp),
    x = 1,
    dims = c(length(grp), nlevels(grp))
  )
  res <- mtx %*% G
  colnames(res) <- levels(grp)
  res
}

safe_log2p1 <- function(mat) {
  log2(mat + 1)
}

build_reference_markers <- function(psych_dir, out_table, markers_per_class = 150L) {
  log_msg("Reading PsychENCODE meta.")
  meta_file <- file.path(psych_dir, "meta.tsv")
  feature_file <- file.path(psych_dir, "counts_features.tsv.gz")
  barcode_file <- file.path(psych_dir, "counts_barcodes.tsv.gz")
  matrix_file <- file.path(psych_dir, "counts_matrix.mtx.gz")

  meta <- fread(meta_file)
  stopifnot("Cell_ID" %in% colnames(meta), "annotation" %in% colnames(meta))

  broad_levels <- c("EXN", "INN", "AST", "OPC", "ODC", "MG", "END")
  meta <- meta[annotation %in% broad_levels]

  log_msg("Reading PsychENCODE feature/barcode files.")
  feat <- read_dt_no_header(feature_file)
  bc <- read_dt_no_header(barcode_file)

  gene_symbol <- if (ncol(feat) >= 2L) feat[[2L]] else feat[[1L]]
  gene_symbol <- clean_symbol(gene_symbol)

  log_msg("Reading PsychENCODE sparse matrix.")
  mtx <- readMM(matrix_file)
  mtx <- as(mtx, "CsparseMatrix")

  nrow_target <- min(nrow(mtx), length(gene_symbol))
  if (nrow_target < nrow(mtx) || nrow_target < length(gene_symbol)) {
    log_msg("Trimming PsychENCODE features/matrix rows to shared length = ", nrow_target)
    mtx <- mtx[seq_len(nrow_target), , drop = FALSE]
    gene_symbol <- gene_symbol[seq_len(nrow_target)]
  }

  ncol_target <- min(ncol(mtx), nrow(meta))
  if (ncol_target < ncol(mtx) || ncol_target < nrow(meta)) {
    log_msg("Trimming PsychENCODE meta/matrix cols to shared length = ", ncol_target)
    mtx <- mtx[, seq_len(ncol_target), drop = FALSE]
    meta <- meta[seq_len(ncol_target)]
  }

  rownames(mtx) <- gene_symbol
  colnames(mtx) <- meta$Cell_ID

  log_msg("Collapsing duplicated PsychENCODE genes by symbol-sum.")
  mtx_gene <- collapse_sparse_by_symbol(mtx, rownames(mtx))

  meta_use <- copy(meta)
  meta_use$annotation <- factor(meta_use$annotation, levels = broad_levels)

  log_msg("Aggregating counts by broad class using annotation.")
  basis_counts <- aggregate_sparse_by_group(mtx_gene, as.character(meta_use$annotation))
  basis_counts <- as.matrix(basis_counts)

  lib_sizes <- colSums(basis_counts)
  lib_sizes[lib_sizes <= 0] <- 1
  basis_cpm <- sweep(basis_counts, 2, lib_sizes, "/", check.margin = FALSE) * 1e6
  basis_logcpm <- safe_log2p1(basis_cpm)

  marker_list <- list()
  for (cls in broad_levels) {
    this <- basis_logcpm[, cls]
    other_cols <- setdiff(colnames(basis_logcpm), cls)
    others <- rowMeans(basis_logcpm[, other_cols, drop = FALSE])
    dt <- data.table(
      gene_symbol = rownames(basis_logcpm),
      broad_class = cls,
      logFC = as.numeric(this - others),
      avg_expr = as.numeric(this),
      avg_others = as.numeric(others)
    )
    dt <- dt[is.finite(logFC) & !is.na(gene_symbol) & gene_symbol != ""]
    setorder(dt, -logFC, -avg_expr)
    dt <- dt[logFC > 0]
    if (nrow(dt) > markers_per_class) dt <- dt[seq_len(markers_per_class)]
    marker_list[[cls]] <- dt
  }

  marker_dt <- rbindlist(marker_list, use.names = TRUE, fill = TRUE)
  if (nrow(marker_dt) == 0L) {
    stop("No valid reference markers identified from PsychENCODE.")
  }
  write_tsv(marker_dt, out_table)
  basis <- basis_logcpm[unique(marker_dt$gene_symbol), broad_levels, drop = FALSE]
  list(
    basis = basis,
    marker_dt = marker_dt,
    broad_levels = broad_levels
  )
}

select_best_sample_id <- function(meta, sample_names) {
  candidates <- c("sample_id", "Sample_ID", "sample", "Sample")
  best_name <- NULL
  best_overlap <- -1
  for (nm in candidates) {
    if (nm %in% colnames(meta)) {
      ov <- mean(sample_names %in% as.character(meta[[nm]]))
      if (is.finite(ov) && ov > best_overlap) {
        best_overlap <- ov
        best_name <- nm
      }
    }
  }
  row_ov <- mean(sample_names %in% rownames(meta))
  if (is.finite(row_ov) && row_ov > best_overlap) {
    best_overlap <- row_ov
    best_name <- "__rownames__"
  }
  list(name = best_name, overlap = best_overlap)
}

prepare_gandal_standardized <- function(raw_rdata, meta_rdata, anno_csv, out_expr, out_meta, summary_file = NULL) {
  log_msg("Loading Gandal raw counts RData: ", raw_rdata)
  env_raw <- new.env(parent = emptyenv())
  load(raw_rdata, envir = env_raw)

  log_msg("Loading Gandal meta/normalized RData: ", meta_rdata)
  env_meta <- new.env(parent = emptyenv())
  load(meta_rdata, envir = env_meta)

  stopifnot(exists("rsem_gene_counts", envir = env_raw))
  stopifnot(exists("datMeta", envir = env_meta))

  expr <- get("rsem_gene_counts", envir = env_raw)
  meta <- as.data.table(get("datMeta", envir = env_meta))

  sample_names <- colnames(expr)
  pick <- select_best_sample_id(meta, sample_names)
  if (is.null(pick$name)) stop("Could not identify Gandal sample ID column.")
  log_msg("Metadata sample identifiers chosen from ", pick$name, " with best overlap = ", sprintf("%.4f", pick$overlap))

  if (pick$name == "__rownames__") {
    meta[, sample := rownames(get("datMeta", envir = env_meta))]
  } else {
    meta[, sample := as.character(get(pick$name))]
  }

  diag_col <- if ("Diagnosis" %in% colnames(meta)) "Diagnosis" else stop("Diagnosis column not found in Gandal datMeta.")
  meta[, diagnosis := normalize_diagnosis(get(diag_col))]
  meta <- meta[!is.na(diagnosis)]
  meta <- meta[sample %in% sample_names]
  meta <- unique(meta, by = "sample")
  meta <- meta[match(sample_names, sample)]
  meta <- meta[!is.na(sample)]
  expr <- expr[, meta$sample, drop = FALSE]

  anno <- fread(anno_csv)
  stopifnot(all(c("ensembl_gene_id", "gene_symbol") %in% colnames(anno)))
  anno <- unique(anno[, .(ensembl_gene_id = clean_ensembl(ensembl_gene_id),
                          gene_symbol = clean_symbol(gene_symbol))])
  anno <- anno[!is.na(ensembl_gene_id) & !is.na(gene_symbol)]

  ens_clean <- clean_ensembl(rownames(expr))
  idx <- match(ens_clean, anno$ensembl_gene_id)
  gene_symbol <- anno$gene_symbol[idx]

  expr2 <- collapse_matrix_by_symbol(expr, gene_symbol)
  log_msg("Gandal expression after symbol collapse: genes = ", nrow(expr2), " ; samples = ", ncol(expr2), " ; symbol-like rows = ", sum(!is.na(rownames(expr2))))

  expr_dt <- as.data.table(expr2, keep.rownames = "gene_symbol")
  meta_out <- meta[, .(sample, diagnosis)]

  fwrite(expr_dt, out_expr, sep = "\t", quote = FALSE, compress = "gzip")
  fwrite(meta_out, out_meta, sep = "\t", quote = FALSE)

  if (!is.null(summary_file)) {
    sum_dt <- data.table(
      cohort = "Gandal2022",
      expr_file = out_expr,
      meta_file = out_meta,
      n_genes = nrow(expr2),
      n_samples = ncol(expr2),
      expr_source = "rsem_gene_counts",
      meta_source = "datMeta"
    )
    fwrite(sum_dt, summary_file, sep = "\t", quote = FALSE)
  }

  invisible(list(expr = expr2, meta = meta_out))
}

prepare_standard_inputs <- function(base) {
  input_dir <- file.path(base, "inputs")
  meta_dir <- file.path(base, "meta")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

  # GSE102741
  gse102_expr <- file.path(input_dir, "GSE102741.standardized_expr.tsv.gz")
  gse102_meta <- file.path(input_dir, "GSE102741.standardized_meta.tsv")
  if (!file.exists(gse102_expr) || !file.exists(gse102_meta)) {
    stop("Missing GSE102741 standardized inputs: ", gse102_expr, " and/or ", gse102_meta)
  }

  # GSE64018
  gse64018_expr <- file.path(input_dir, "GSE64018.standardized_expr.tsv.gz")
  gse64018_meta <- file.path(input_dir, "GSE64018.standardized_meta.tsv")
  if (!file.exists(gse64018_expr) || !file.exists(gse64018_meta)) {
    stop("Missing GSE64018 standardized inputs: ", gse64018_expr, " and/or ", gse64018_meta)
  }

  # Gandal
  gandal_expr <- file.path(input_dir, "Gandal2022.standardized_expr.tsv.gz")
  gandal_meta <- file.path(input_dir, "Gandal2022.standardized_meta.tsv")
  prepare_gandal_standardized(
    raw_rdata = GANDAL_RAW,
    meta_rdata = GANDAL_META,
    anno_csv = GANDAL_ANNO,
    out_expr = gandal_expr,
    out_meta = gandal_meta
  )

  sum_rows <- list()

  add_summary <- function(cohort, expr_file, meta_file, expr_source = "", meta_source = "") {
    dt_expr <- fread(expr_file, nrows = 2)
    dt_meta <- fread(meta_file)
    n_genes <- as.character(as.integer(fread(cmd = paste("zcat", shQuote(expr_file), "| wc -l"), header = FALSE)[[1]] - 1L))
    if (!grepl("\\.gz$", expr_file)) n_genes <- as.character(as.integer(nrow(fread(expr_file)) ))
    sum_rows[[length(sum_rows) + 1L]] <<- data.table(
      cohort = cohort,
      expr_file = expr_file,
      meta_file = meta_file,
      n_genes = as.integer(n_genes),
      n_samples = ncol(dt_expr) - 1L,
      expr_source = expr_source,
      meta_source = meta_source
    )
  }

  add_summary("GSE102741", gse102_expr, gse102_meta)
  add_summary("GSE64018", gse64018_expr, gse64018_meta)
  add_summary("Gandal2022", gandal_expr, gandal_meta, "rsem_gene_counts", "datMeta")

  prep_dt <- rbindlist(sum_rows, use.names = TRUE)
  fwrite(prep_dt, file.path(meta_dir, "Step10A_prepare_bulk_inputs_summary.tsv"), sep = "\t", quote = FALSE)

  list(
    GSE102741 = list(expr = gse102_expr, meta = gse102_meta),
    GSE64018 = list(expr = gse64018_expr, meta = gse64018_meta),
    Gandal2022 = list(expr = gandal_expr, meta = gandal_meta)
  )
}

read_standardized_expr <- function(path) {
  dt <- fread(path)
  first_col <- colnames(dt)[1L]
  setnames(dt, first_col, "gene_symbol")
  dt[, gene_symbol := clean_symbol(gene_symbol)]
  dt <- dt[!is.na(gene_symbol) & gene_symbol != ""]
  num_cols <- setdiff(colnames(dt), "gene_symbol")
  for (cc in num_cols) set(dt, j = cc, value = as.numeric(dt[[cc]]))
  dt <- dt[, lapply(.SD, sum), by = gene_symbol, .SDcols = num_cols]
  expr <- as.data.frame(dt)
  rownames(expr) <- expr$gene_symbol
  expr$gene_symbol <- NULL
  as.matrix(expr)
}

read_standardized_meta <- function(path) {
  meta <- fread(path)
  stopifnot(all(c("sample", "diagnosis") %in% colnames(meta)))
  meta[, sample := as.character(sample)]
  meta[, diagnosis := normalize_diagnosis(diagnosis)]
  meta <- meta[!is.na(diagnosis)]
  meta
}

estimate_fractions <- function(basis, target) {
  basis <- as.matrix(basis)
  target <- as.matrix(target)
  common <- intersect(rownames(basis), rownames(target))
  basis <- basis[common, , drop = FALSE]
  target <- target[common, , drop = FALSE]

  if (length(common) < 10L) stop("Too few shared marker genes = ", length(common))

  frac <- matrix(NA_real_, nrow = ncol(target), ncol = ncol(basis))
  colnames(frac) <- colnames(basis)
  rownames(frac) <- colnames(target)
  n_fallback <- 0L

  for (i in seq_len(ncol(target))) {
    fit <- nnls::nnls(basis, target[, i])
    x <- coef(fit)
    x[!is.finite(x) | x < 0] <- 0
    if (sum(x) <= 0) {
      x[] <- 1 / length(x)
      n_fallback <- n_fallback + 1L
    } else {
      x <- x / sum(x)
    }
    frac[i, ] <- x
  }

  list(frac = as.data.table(frac, keep.rownames = "sample"), n_fallback = n_fallback, n_shared = length(common))
}

run_fraction_tests <- function(frac_dt, cohort) {
  cell_types <- setdiff(colnames(frac_dt), c("sample", "diagnosis", "cohort"))
  res <- lapply(cell_types, function(ct) {
    x_asd <- frac_dt[diagnosis == "ASD", get(ct)]
    x_ctl <- frac_dt[diagnosis == "Control", get(ct)]
    wilcox_p <- tryCatch(wilcox.test(x_asd, x_ctl)$p.value, error = function(e) NA_real_)
    auc <- calc_auc(x_asd, x_ctl)
    fit <- tryCatch(lm(reformulate("diagnosis", response = ct), data = as.data.frame(frac_dt)), error = function(e) NULL)
    coef_name <- NA_character_
    beta <- NA_real_
    p <- NA_real_
    if (!is.null(fit)) {
      sm <- summary(fit)$coefficients
      coef_name <- rownames(sm)[grep("^diagnosis", rownames(sm))[1]]
      if (!is.na(coef_name)) {
        beta <- sm[coef_name, "Estimate"]
        p <- sm[coef_name, "Pr(>|t|)"]
      }
    }
    data.table(
      cohort = cohort,
      cell_type = ct,
      n_samples = nrow(frac_dt),
      n_asd = sum(frac_dt$diagnosis == "ASD"),
      n_control = sum(frac_dt$diagnosis == "Control"),
      mean_prop_control = mean(x_ctl),
      mean_prop_asd = mean(x_asd),
      delta_asd_minus_control = mean(x_asd) - mean(x_ctl),
      wilcox_p = wilcox_p,
      lm_beta_asd = beta,
      lm_p = p,
      lm_coef_name = coef_name,
      abs_delta = abs(mean(x_asd) - mean(x_ctl))
    )
  })
  out <- rbindlist(res)
  out[, fdr_wilcox := p.adjust(wilcox_p, method = "fdr")]
  out[, fdr_lm := p.adjust(lm_p, method = "fdr")]
  setcolorder(out, c("cohort", "cell_type", "n_samples", "n_asd", "n_control",
                     "mean_prop_control", "mean_prop_asd", "delta_asd_minus_control",
                     "wilcox_p", "lm_beta_asd", "lm_p", "lm_coef_name",
                     "fdr_wilcox", "fdr_lm", "abs_delta"))
  out
}

read_programs <- function(path) {
  dt <- fread(path)
  gene_col <- intersect(c("gene_symbol", "symbol", "Gene", "gene"), colnames(dt))[1]
  prog_col <- intersect(c("program_name", "program", "set_name"), colnames(dt))[1]
  if (is.na(gene_col) || is.na(prog_col)) stop("Program file missing required columns.")
  dt[, gene_symbol := clean_symbol(get(gene_col))]
  dt[, program_name := as.character(get(prog_col))]
  dt <- dt[!is.na(gene_symbol) & gene_symbol != "" & !is.na(program_name) & program_name != ""]
  out <- split(dt$gene_symbol, dt$program_name)
  out <- lapply(out, function(x) unique(x[!is.na(x) & x != ""]))
  out
}

choose_fraction_covariates <- function(frac_dt) {
  covars <- intersect(c("EXN", "AST", "ODC", "MG", "END", "OPC"), colnames(frac_dt))
  covars <- covars[sapply(covars, function(cc) stats::sd(frac_dt[[cc]]) > 0)]
  covars
}

run_program_tests <- function(expr, meta, frac_dt, programs, cohort) {
  common_samples <- Reduce(intersect, list(colnames(expr), meta$sample, frac_dt$sample))
  expr <- expr[, common_samples, drop = FALSE]
  meta <- meta[match(common_samples, meta$sample)]
  frac_dt <- frac_dt[match(common_samples, frac_dt$sample)]

  expr_log <- safe_log2p1(expr)
  covars <- choose_fraction_covariates(frac_dt)
  frac_df <- as.data.frame(frac_dt[, c("sample", covars), with = FALSE])

  out <- lapply(names(programs), function(pn) {
    genes <- intersect(programs[[pn]], rownames(expr_log))
    if (length(genes) == 0L) return(NULL)
    score <- colMeans(expr_log[genes, , drop = FALSE])
    df <- data.frame(
      sample = common_samples,
      diagnosis = factor(meta$diagnosis, levels = c("Control", "ASD")),
      score = as.numeric(score),
      stringsAsFactors = FALSE
    )
    if (length(covars) > 0L) {
      df <- merge(df, frac_df, by = "sample", sort = FALSE)
      df <- df[match(common_samples, df$sample), , drop = FALSE]
    }

    x_asd <- df$score[df$diagnosis == "ASD"]
    x_ctl <- df$score[df$diagnosis == "Control"]
    auc <- calc_auc(x_asd, x_ctl)
    wilcox_p <- tryCatch(wilcox.test(x_asd, x_ctl)$p.value, error = function(e) NA_real_)

    fit_unadj <- lm(score ~ diagnosis, data = df)
    sm_u <- summary(fit_unadj)$coefficients
    coef_u <- grep("^diagnosis", rownames(sm_u), value = TRUE)[1]
    beta_u <- sm_u[coef_u, "Estimate"]
    p_u <- sm_u[coef_u, "Pr(>|t|)"]

    if (length(covars) > 0L) {
      fml <- as.formula(paste("score ~ diagnosis +", paste(covars, collapse = " + ")))
    } else {
      fml <- score ~ diagnosis
    }
    fit_adj <- lm(fml, data = df)
    sm_a <- summary(fit_adj)$coefficients
    coef_a <- grep("^diagnosis", rownames(sm_a), value = TRUE)[1]
    beta_a <- sm_a[coef_a, "Estimate"]
    p_a <- sm_a[coef_a, "Pr(>|t|)"]

    data.table(
      cohort = cohort,
      program_name = pn,
      n_samples = nrow(df),
      n_asd = sum(df$diagnosis == "ASD"),
      n_control = sum(df$diagnosis == "Control"),
      n_genes_in_matrix = length(genes),
      mean_score_control = mean(x_ctl),
      mean_score_asd = mean(x_asd),
      delta_asd_minus_control = mean(x_asd) - mean(x_ctl),
      auc_asd_higher = auc,
      wilcox_p = wilcox_p,
      lm_beta_asd_unadjusted = beta_u,
      lm_p_unadjusted = p_u,
      lm_beta_asd_adjusted = beta_a,
      lm_p_adjusted = p_a,
      lm_coef_name_adjusted = coef_a,
      formula_used = deparse(fml),
      fraction_covariates = paste(covars, collapse = ",")
    )
  })
  out <- rbindlist(out, use.names = TRUE, fill = TRUE)
  out[, fdr_wilcox := p.adjust(wilcox_p, method = "fdr")]
  out[, fdr_lm_unadjusted := p.adjust(lm_p_unadjusted, method = "fdr")]
  out[, fdr_lm_adjusted := p.adjust(lm_p_adjusted, method = "fdr")]
  out
}

main <- function() {
  log_msg("Starting Step10A NNLS deconvolution and composition-adjusted bulk testing (v16).")

  ref <- build_reference_markers(
    psych_dir = PSY_DIR,
    out_table = file.path(BASE, "tables", "140_reference_marker_summary.tsv"),
    markers_per_class = MARKERS_PER_CLASS
  )

  basis <- ref$basis
  programs <- read_programs(PROGRAM_FILE)
  log_msg("Loaded programs: ", paste(sprintf("%s=%d", names(programs), lengths(programs)), collapse = "; "))

  cohort_files <- prepare_standard_inputs(BASE)

  all_frac <- list()
  all_dx <- list()
  all_prog <- list()
  run_summary <- list(
    data.table(section = "Reference", metric = "n_broad_classes", value = length(ref$broad_levels)),
    data.table(section = "Reference", metric = "n_union_markers", value = nrow(ref$marker_dt))
  )

  for (cohort in names(cohort_files)) {
    log_msg("Loading cohort: ", cohort)
    expr <- tryCatch(read_standardized_expr(cohort_files[[cohort]]$expr), error = function(e) e)
    meta <- tryCatch(read_standardized_meta(cohort_files[[cohort]]$meta), error = function(e) e)

    if (inherits(expr, "error") || inherits(meta, "error")) {
      msg <- if (inherits(expr, "error")) expr$message else meta$message
      log_msg("Failed loading ", cohort, " : ", msg)
      run_summary[[length(run_summary) + 1L]] <- data.table(section = cohort, metric = "n_samples", value = 0)
      run_summary[[length(run_summary) + 1L]] <- data.table(section = cohort, metric = "status", value = "load_failed")
      next
    }

    ali <- align_meta_expr(expr, meta, sample_col = "sample")
    expr <- ali$expr
    meta <- ali$meta

    target_genes <- intersect(rownames(expr), rownames(basis))
    log_msg(cohort, ": shared marker genes = ", length(target_genes))

    frac_fit <- tryCatch(
      estimate_fractions(basis = basis[target_genes, , drop = FALSE], target = safe_log2p1(expr[target_genes, , drop = FALSE])),
      error = function(e) e
    )

    if (inherits(frac_fit, "error")) {
      log_msg("Skipping ", cohort, ": ", frac_fit$message)
      run_summary[[length(run_summary) + 1L]] <- data.table(section = cohort, metric = "n_samples", value = ncol(expr))
      run_summary[[length(run_summary) + 1L]] <- data.table(section = cohort, metric = "status", value = frac_fit$message)
      next
    }

    frac_dt <- frac_fit$frac
    frac_dt[, cohort := cohort]
    frac_dt[, diagnosis := meta$diagnosis[match(sample, meta$sample)]]

    dx_dt <- run_fraction_tests(frac_dt, cohort)
    prog_dt <- run_program_tests(expr = expr, meta = meta, frac_dt = frac_dt, programs = programs, cohort = cohort)

    all_frac[[cohort]] <- frac_dt
    all_dx[[cohort]] <- dx_dt
    all_prog[[cohort]] <- prog_dt

    run_summary[[length(run_summary) + 1L]] <- data.table(section = cohort, metric = "n_samples", value = nrow(frac_dt))
    run_summary[[length(run_summary) + 1L]] <- data.table(section = cohort, metric = "status", value = "completed")
    run_summary[[length(run_summary) + 1L]] <- data.table(section = cohort, metric = "n_fallback_uniform", value = frac_fit$n_fallback)
  }

  if (length(all_frac) > 0L) {
    frac_out <- rbindlist(all_frac, use.names = TRUE, fill = TRUE)
    dx_out <- rbindlist(all_dx, use.names = TRUE, fill = TRUE)
    prog_out <- rbindlist(all_prog, use.names = TRUE, fill = TRUE)

    fwrite(frac_out, file.path(BASE, "tables", "141_bulk_nnls_fractions.tsv"), sep = "\t", quote = FALSE)
    fwrite(dx_out, file.path(BASE, "tables", "142_bulk_celltype_dx_tests_nnls.tsv"), sep = "\t", quote = FALSE)
    fwrite(prog_out, file.path(BASE, "tables", "143_bulk_program_dx_tests_adjusted_nnls.tsv"), sep = "\t", quote = FALSE)
  }

  sum_dt <- rbindlist(run_summary, use.names = TRUE, fill = TRUE)
  fwrite(sum_dt, file.path(BASE, "meta", "144_Step10A_run_summary.tsv"), sep = "\t", quote = FALSE)

  log_msg("Step10A_bulk_cell_composition_deconvolution_nnls_v17 completed successfully.")
}

main()
