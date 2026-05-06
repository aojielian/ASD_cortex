suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

logmsg <- function(...) {
  msg <- paste0("[", timestamp(), "] ", paste(..., collapse = ""))
  message(msg)
}

std_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  toupper(x)
}

safe_num_mat <- function(x) {
  if (is.list(x) && !is.data.frame(x) && length(x) == 1L) x <- x[[1]]
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

safe_write <- function(dt, file) {
  fwrite(as.data.table(dt), file = file, sep = "\t", quote = FALSE, na = "NA")
}

is_numeric_matrix_like <- function(x) {
  ok <- TRUE
  if (!(is.matrix(x) || is.data.frame(x))) return(FALSE)
  y <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (is.null(y)) return(FALSE)
  ok <- tryCatch({
    storage.mode(y) <- "numeric"
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}

collect_matrix_candidates <- function(x, path = "obj", depth = 0L, max_depth = 3L) {
  out <- list()
  if (is_numeric_matrix_like(x)) {
    rn <- rownames(x); cn <- colnames(x)
    out[[length(out) + 1L]] <- list(
      path = path,
      object = x,
      nrow = NROW(x),
      ncol = NCOL(x),
      rownames_present = !is.null(rn),
      colnames_present = !is.null(cn)
    )
  }
  if (depth < max_depth && is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    for (i in seq_along(x)) {
      child <- x[[i]]
      child_path <- paste0(path, "$", nms[i])
      out <- c(out, collect_matrix_candidates(child, path = child_path, depth = depth + 1L, max_depth = max_depth))
    }
  }
  out
}

choose_expr_from_rds <- function(obj, syn_genes, outdir = NULL, prefix = "GSE102741") {
  cands <- collect_matrix_candidates(obj, path = "rds")
  if (length(cands) == 0) stop(prefix, ": no numeric matrix-like candidates found in RDS object.")
  inv <- rbindlist(lapply(cands, function(z) {
    rn <- rownames(z$object); cn <- colnames(z$object)
    rn_std <- if (!is.null(rn)) std_gene(rn) else character()
    cn_std <- if (!is.null(cn)) std_gene(cn) else character()
    gene_like_rows <- sum(grepl("^[A-Z0-9._-]+$", rn_std), na.rm = TRUE)
    gene_like_cols <- sum(grepl("^[A-Z0-9._-]+$", cn_std), na.rm = TRUE)
    data.table(
      path = z$path,
      nrow = as.numeric(z$nrow),
      ncol = as.numeric(z$ncol),
      rownames_present = z$rownames_present,
      colnames_present = z$colnames_present,
      row_syn_overlap = sum(unique(rn_std[!is.na(rn_std)]) %in% syn_genes),
      col_syn_overlap = sum(unique(cn_std[!is.na(cn_std)]) %in% syn_genes),
      gene_like_rows = gene_like_rows,
      gene_like_cols = gene_like_cols
    )
  }), fill = TRUE)
  inv[, best_syn_overlap := pmax(row_syn_overlap, col_syn_overlap)]
  inv[, sort_rows := as.numeric(nrow)]
  inv[, sort_cols := as.numeric(ncol)]
  setorderv(inv, c("best_syn_overlap", "sort_rows", "sort_cols", "path"), c(-1L, -1L, -1L, 1L))
  best_path <- inv$path[1]
  best_obj <- cands[[which(vapply(cands, function(z) identical(z$path, best_path), logical(1)))[1]]]$object
  if (!is.null(outdir)) {
    safe_write(inv, file.path(outdir, paste0("00a_", prefix, "_rds_candidate_inventory.tsv")))
  }
  list(expr = safe_num_mat(best_obj), inventory = inv, chosen_path = best_path)
}

choose_sample_column_and_orient <- function(expr, meta) {
  meta_dt <- as.data.table(meta)
  candidate_cols <- names(meta_dt)
  if (length(candidate_cols) == 0) stop("Metadata has no columns.")
  score_dt <- rbindlist(lapply(candidate_cols, function(cc) {
    vals <- as.character(meta_dt[[cc]])
    data.table(
      candidate_col = cc,
      overlap_with_cols = sum(vals %in% colnames(expr), na.rm = TRUE),
      overlap_with_rows = sum(vals %in% rownames(expr), na.rm = TRUE),
      n_unique = uniqueN(vals[!is.na(vals)])
    )
  }))
  setorderv(score_dt, c("overlap_with_cols", "overlap_with_rows", "n_unique", "candidate_col"), c(-1L, -1L, -1L, 1L))
  best <- score_dt[1]
  orientation <- ifelse(best$overlap_with_rows > best$overlap_with_cols, "transpose", "as_is")
  sample_col <- best$candidate_col
  if (orientation == "transpose") expr <- t(expr)
  list(expr = expr, sample_col = sample_col, alignment = score_dt, orientation = orientation)
}

detect_gene_axis <- function(mat, syn_genes) {
  rn <- rownames(mat); cn <- colnames(mat)
  rn_std <- if (!is.null(rn)) std_gene(rn) else character()
  cn_std <- if (!is.null(cn)) std_gene(cn) else character()
  row_overlap <- sum(unique(rn_std[!is.na(rn_std)]) %in% syn_genes)
  col_overlap <- sum(unique(cn_std[!is.na(cn_std)]) %in% syn_genes)
  if (row_overlap >= col_overlap) "rows" else "cols"
}

clean_gandal_ensg <- function(x) {
  x <- as.character(x)
  x <- sub("_[0-9]+$", "", x)          # ENSG....14_1 -> ENSG....14
  x <- sub("\\.[0-9]+$", "", x)        # ENSG....14   -> ENSG....
  x
}

map_ensembl_to_symbol <- function(ensg_ids) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Gandal mapping requires AnnotationDbi and org.Hs.eg.db in the current R environment.")
  }
  mapped <- suppressMessages(suppressWarnings(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(ensg_ids),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )))
  data.table(ensembl = names(mapped), symbol = as.character(mapped))
}

aggregate_by_symbol_mean <- function(mat, symbols_std) {
  keep <- !is.na(symbols_std) & nzchar(symbols_std)
  mat <- mat[keep, , drop = FALSE]
  symbols_std <- symbols_std[keep]
  rs <- rowsum(mat, group = symbols_std, reorder = FALSE)
  cnt <- table(symbols_std)
  rs / as.numeric(cnt[rownames(rs)])
}

extract_dx <- function(meta) {
  cand <- names(meta)[tolower(names(meta)) %in% c("diagnosis","dx","dx_group","asd","group","condition","status")]
  if (length(cand) == 0) {
    cand <- names(meta)[grepl("diagnosis|dx|asd|group|condition|status", names(meta), ignore.case = TRUE)]
  }
  if (length(cand) == 0) stop("Could not detect diagnosis column.")
  dx_col <- cand[1]
  dx_raw <- as.character(meta[[dx_col]])
  dx <- rep(NA_character_, length(dx_raw))
  dx[grepl("^asd$|aut", dx_raw, ignore.case = TRUE)] <- "ASD"
  dx[grepl("^control$|^ctl$|^ctrl$|normal|non[-_ ]?asd", dx_raw, ignore.case = TRUE)] <- "Control"
  # numeric fallback
  suppressWarnings(num <- as.numeric(dx_raw))
  if (any(!is.na(num))) {
    if (all(na.omit(num) %in% c(0,1))) {
      dx[!is.na(num) & num == 1] <- "ASD"
      dx[!is.na(num) & num == 0] <- "Control"
    }
  }
  list(dx = dx, dx_col = dx_col)
}

extract_sample_ids <- function(meta, expr_colnames = NULL) {
  cand <- names(meta)[tolower(names(meta)) %in% c("sample_id","sample","sampleid","iid","id","run","libraryid","rnaseq_id")]
  if (length(cand) == 0) {
    cand <- names(meta)[grepl("sample|iid|library|id$", names(meta), ignore.case = TRUE)]
  }
  if (length(cand) == 0) {
    meta$sample_id_auto <- rownames(meta)
    return(list(sample_id = as.character(meta$sample_id_auto), sample_col = "rownames"))
  }
  if (!is.null(expr_colnames)) {
    ov <- sapply(cand, function(cc) sum(as.character(meta[[cc]]) %in% expr_colnames, na.rm = TRUE))
    best <- cand[which.max(ov)]
  } else {
    best <- cand[1]
  }
  list(sample_id = as.character(meta[[best]]), sample_col = best)
}

build_auto_formula <- function(df) {
  base_terms <- "dx_group"
  cov_patterns <- c("sex","gender","age","rin","pmi","region","batch","site","brain_region")
  covs <- names(df)[tolower(names(df)) %in% cov_patterns]
  if (length(covs) == 0) {
    covs <- names(df)[grepl("sex|gender|age|rin|pmi|region|batch|site", names(df), ignore.case = TRUE)]
  }
  covs <- unique(setdiff(covs, c("sample_id","dx_group","program_score_raw","program_score_z")))
  if (length(covs) == 0) {
    return(reformulate(base_terms, response = "program_score_z"))
  }
  usable <- vapply(covs, function(cc) {
    x <- df[[cc]]
    if (is.null(x)) return(FALSE)
    if (is.list(x) && !is.data.frame(x)) return(FALSE)
    vals <- tryCatch(na.omit(x), error = function(e) x)
    vals <- tryCatch(unique(vals), error = function(e) unique(as.character(vals)))
    length(vals) > 1
  }, logical(1))
  covs <- covs[usable]
  if (length(covs) == 0) {
    return(reformulate(base_terms, response = "program_score_z"))
  }
  covs <- covs[seq_len(min(length(covs), 4))]
  reformulate(c(base_terms, covs), response = "program_score_z")
}

fit_effect <- function(df, formula_obj) {
  df2 <- copy(df)
  df2 <- df2[!is.na(program_score_z) & !is.na(dx_group)]
  if (nrow(df2) < 8) {
    return(data.table(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      p = NA_real_, n_samples = nrow(df2),
                      n_asd = sum(df2$dx_group == "ASD", na.rm = TRUE),
                      n_control = sum(df2$dx_group == "Control", na.rm = TRUE),
                      formula = deparse(formula_obj), status = "too_few_samples"))
  }
  fit <- tryCatch(lm(formula_obj, data = df2), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.table(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      p = NA_real_, n_samples = nrow(df2),
                      n_asd = sum(df2$dx_group == "ASD", na.rm = TRUE),
                      n_control = sum(df2$dx_group == "Control", na.rm = TRUE),
                      formula = deparse(formula_obj), status = "lm_failed"))
  }
  cf <- summary(fit)$coefficients
  rn <- rownames(cf)
  target <- rn[grepl("dx_group", rn)]
  target <- target[grepl("ASD", target)]
  if (length(target) == 0) {
    return(data.table(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      p = NA_real_, n_samples = nrow(df2),
                      n_asd = sum(df2$dx_group == "ASD", na.rm = TRUE),
                      n_control = sum(df2$dx_group == "Control", na.rm = TRUE),
                      formula = deparse(formula_obj), status = "dx_term_missing"))
  }
  tt <- target[1]
  beta <- unname(cf[tt, "Estimate"])
  se <- unname(cf[tt, "Std. Error"])
  p <- unname(cf[tt, grep("Pr\\(", colnames(cf))])
  ci_low <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  data.table(beta = beta, se = se, ci_low = ci_low, ci_high = ci_high,
             p = p, n_samples = nrow(df2),
             n_asd = sum(df2$dx_group == "ASD", na.rm = TRUE),
             n_control = sum(df2$dx_group == "Control", na.rm = TRUE),
             formula = deparse(formula_obj), status = "ok")
}

make_plot <- function(dt, out_pdf, title_txt) {
  if (nrow(dt) == 0) return(invisible(NULL))
  dtp <- copy(dt)
  dtp <- dtp[order(abs_delta_beta)]
  dtp[, gene_f := factor(dropped_gene, levels = dropped_gene)]
  p <- ggplot(dtp, aes(x = abs_delta_beta, y = gene_f)) +
    geom_point(size = 1.3) +
    labs(x = "|Δ beta from full model|", y = NULL, title = title_txt) +
    theme_bw(base_size = 10)
  h <- min(49, max(4.5, 0.12 * nrow(dtp)))
  ggsave(out_pdf, plot = p, width = 7.5, height = h, limitsize = FALSE)
}

prepare_gse102741 <- function(root_dir, syn_genes = NULL, outdir = NULL) {
  expr_path <- file.path(root_dir, "step01_prepare/rds/22_GSE102741_geneSymbol_log2cpm.rds")
  meta_path <- file.path(root_dir, "step01_prepare/tables/21_GSE102741_sample_metadata.tsv")
  obj <- readRDS(expr_path)
  if (is.null(syn_genes)) syn_genes <- character()
  chosen <- if (is.matrix(obj) || is.data.frame(obj)) {
    list(expr = safe_num_mat(obj), inventory = data.table(path = "rds", nrow = NROW(obj), ncol = NCOL(obj)), chosen_path = "rds")
  } else {
    choose_expr_from_rds(obj, syn_genes = syn_genes, outdir = outdir, prefix = "GSE102741")
  }
  expr <- chosen$expr
  meta <- fread(meta_path)
  orient <- choose_sample_column_and_orient(expr, meta)
  expr <- orient$expr
  if (!is.null(outdir)) {
    safe_write(orient$alignment, file.path(outdir, "00d_GSE102741_sample_alignment.tsv"))
    safe_write(data.table(chosen_path = chosen$chosen_path, orientation = orient$orientation, sample_col = orient$sample_col),
               file.path(outdir, "00e_GSE102741_chosen_expr_object.tsv"))
  }
  sid <- extract_sample_ids(meta, expr_colnames = colnames(expr))
  dxi <- extract_dx(meta)
  meta[, sample_id := sid$sample_id]
  meta[, dx_group := dxi$dx]
  list(expr = expr, meta = meta, input_scale = "log2cpm_or_normalized",
       sample_col = sid$sample_col, dx_col = dxi$dx_col, mapping = NULL)
}

prepare_gse64018 <- function(root_dir) {
  expr_path <- file.path(root_dir, "step01_prepare/tables/122_GSE64018_collapsed_counts_by_symbol.tsv.gz")
  meta_path <- file.path(root_dir, "step01_prepare/meta/120_GSE64018_explicit_sample_mapping.tsv")
  dt <- fread(expr_path)
  gene_col <- names(dt)[1]
  genes <- std_gene(dt[[gene_col]])
  expr <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- genes
  meta <- fread(meta_path)
  sid <- extract_sample_ids(meta, expr_colnames = colnames(expr))
  dxi <- extract_dx(meta)
  meta[, sample_id := sid$sample_id]
  meta[, dx_group := dxi$dx]
  # log2cpm-like transform from collapsed counts
  lib <- colSums(expr, na.rm = TRUE)
  expr <- log2(t(t(expr) / pmax(lib, 1) * 1e6) + 1)
  list(expr = expr, meta = meta, input_scale = "log2cpm_from_collapsed_counts",
       sample_col = sid$sample_col, dx_col = dxi$dx_col, mapping = NULL)
}

prepare_gandal <- function(root_dir) {
  rdata_path <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  datExpr <- get("datExpr", envir = env)
  datMeta <- get("datMeta", envir = env)
  datExpr <- safe_num_mat(datExpr)

  ensg_full <- rownames(datExpr)
  ensg_clean <- clean_gandal_ensg(ensg_full)
  map_dt <- map_ensembl_to_symbol(ensg_clean)
  map_dt[, ensembl_std := std_gene(ensembl)]
  map_dt[, symbol_std := std_gene(symbol)]

  idx_dt <- data.table(
    raw_rowname = ensg_full,
    ensembl_clean = ensg_clean,
    ensembl_std = std_gene(ensg_clean)
  )
  idx_dt <- merge(idx_dt, map_dt[, .(ensembl_std, symbol_std)], by = "ensembl_std", all.x = TRUE)
  symbol_std <- idx_dt$symbol_std
  expr_agg <- aggregate_by_symbol_mean(datExpr, symbol_std)

  meta <- as.data.table(datMeta)
  sid <- extract_sample_ids(meta, expr_colnames = colnames(expr_agg))
  dxi <- extract_dx(meta)
  meta[, sample_id := sid$sample_id]
  meta[, dx_group := dxi$dx]

  mapping_sum <- data.table(
    n_rows_datExpr = nrow(datExpr),
    n_unique_ensembl_clean = uniqueN(ensg_clean),
    n_rows_mapped_to_symbol = sum(!is.na(symbol_std)),
    n_unique_symbol = uniqueN(symbol_std[!is.na(symbol_std)]),
    mapped_fraction_rows = mean(!is.na(symbol_std))
  )

  list(expr = expr_agg, meta = meta, input_scale = "normalized_residual_zlike",
       sample_col = sid$sample_col, dx_col = dxi$dx_col, mapping = list(idx = idx_dt, sum = mapping_sum))
}

run_dataset <- function(dataset_name, obj, syn_genes, outdir) {
  expr <- obj$expr
  meta <- copy(obj$meta)
  expr_genes <- std_gene(rownames(expr))
  rownames(expr) <- expr_genes

  common_samples <- intersect(colnames(expr), meta$sample_id)
  expr <- expr[, common_samples, drop = FALSE]
  meta <- meta[match(common_samples, sample_id)]
  keep_dx <- !is.na(meta$dx_group) & meta$dx_group %in% c("ASD", "Control")
  expr <- expr[, keep_dx, drop = FALSE]
  meta <- meta[keep_dx]
  meta[, dx_group := factor(dx_group, levels = c("Control", "ASD"))]

  mapped <- intersect(unique(expr_genes), syn_genes)
  mapped <- mapped[!is.na(mapped)]
  logmsg("[", dataset_name, "] mapped synaptic genes: ", length(mapped), " / ", length(syn_genes))
  if (length(mapped) == 0) {
    return(list(
      summary = data.table(dataset = dataset_name, model = c("diagnosis_only", "auto_adjusted"),
                           input_scale = obj$input_scale,
                           n_program_genes_input = length(syn_genes),
                           n_program_genes_mapped = 0,
                           n_program_genes_topdrop = 0,
                           topN = 5,
                           topN_genes = NA_character_,
                           full_beta = NA_real_, full_p = NA_real_, full_formula = NA_character_,
                           topdrop_beta = NA_real_, topdrop_p = NA_real_, topdrop_delta_beta = NA_real_,
                           topdrop_direction_retained = NA, loo_n_ok = 0L,
                           loo_direction_retained_prop = NA_real_, loo_p_lt_0_05_prop = NA_real_,
                           loo_max_abs_delta_beta = NA_real_, loo_median_abs_delta_beta = NA_real_,
                           strongest_influence_gene = NA_character_),
      infl = data.table(), top_expr = data.table(), effects = data.table(), sample_scores = data.table()
    ))
  }

  expr_sub <- expr[mapped, , drop = FALSE]
  gene_mean <- rowMeans(expr_sub, na.rm = TRUE)
  topN <- min(5L, length(gene_mean))
  top_genes <- names(sort(gene_mean, decreasing = TRUE))[seq_len(topN)]

  top_expr_dt <- data.table(
    gene = names(sort(gene_mean, decreasing = TRUE)),
    mean_expression = as.numeric(sort(gene_mean, decreasing = TRUE)),
    rank_by_mean_expression = seq_along(sort(gene_mean, decreasing = TRUE)),
    dataset = dataset_name
  )

  score_dt <- data.table(sample_id = colnames(expr_sub))
  score_dt[, program_score_raw := colMeans(expr_sub, na.rm = TRUE)]
  score_dt[, program_score_z := as.numeric(scale(program_score_raw))]
  score_dt <- merge(score_dt, meta, by = "sample_id", all.x = TRUE)

  score_dt_top <- copy(score_dt)
  expr_topdrop <- expr_sub[setdiff(rownames(expr_sub), top_genes), , drop = FALSE]
  score_dt_top[, program_score_raw := colMeans(expr_topdrop, na.rm = TRUE)]
  score_dt_top[, program_score_z := as.numeric(scale(program_score_raw))]

  formulas <- list(
    diagnosis_only = as.formula("program_score_z ~ dx_group"),
    auto_adjusted = build_auto_formula(score_dt)
  )

  effects_list <- list()
  summary_rows <- list()
  infl_rows <- list()
  sample_score_rows <- list()

  for (mm in names(formulas)) {
    full_eff <- fit_effect(score_dt, formulas[[mm]])
    top_eff  <- fit_effect(score_dt_top, formulas[[mm]])

    loo_rows <- lapply(rownames(expr_sub), function(g) {
      keep_g <- setdiff(rownames(expr_sub), g)
      dtg <- data.table(sample_id = colnames(expr_sub))
      dtg[, program_score_raw := colMeans(expr_sub[keep_g, , drop = FALSE], na.rm = TRUE)]
      dtg[, program_score_z := as.numeric(scale(program_score_raw))]
      dtg <- merge(dtg, meta, by = "sample_id", all.x = TRUE)
      eff <- fit_effect(dtg, formulas[[mm]])
      eff[, dataset := dataset_name]
      eff[, model := mm]
      eff[, scenario := "leave_one_out"]
      eff[, dropped_gene := g]
      eff[, dropped_genes := g]
      eff
    })
    loo_dt <- rbindlist(loo_rows, fill = TRUE)

    if (nrow(loo_dt) > 0) {
      loo_dt[, full_beta := full_eff$beta[1]]
      loo_dt[, delta_beta := beta - full_beta]
      loo_dt[, abs_delta_beta := abs(delta_beta)]
      infl_dt <- loo_dt[status == "ok"][order(-abs_delta_beta)]
      if (nrow(infl_dt) == 0) infl_dt <- loo_dt[order(-abs_delta_beta)]
    } else {
      infl_dt <- data.table()
    }

    full_dir <- ifelse(is.na(full_eff$beta[1]), NA_character_, ifelse(full_eff$beta[1] < 0, "negative", ifelse(full_eff$beta[1] > 0, "positive", "zero")))
    top_dir  <- ifelse(is.na(top_eff$beta[1]), NA_character_, ifelse(top_eff$beta[1] < 0, "negative", ifelse(top_eff$beta[1] > 0, "positive", "zero")))

    loo_ok <- loo_dt[status == "ok"]
    if (nrow(loo_ok) > 0) {
      loo_dir <- ifelse(loo_ok$beta < 0, "negative", ifelse(loo_ok$beta > 0, "positive", "zero"))
      loo_retained <- mean(loo_dir == full_dir, na.rm = TRUE)
      loo_p_sig <- mean(loo_ok$p < 0.05, na.rm = TRUE)
      loo_max_abs <- max(loo_ok$abs_delta_beta, na.rm = TRUE)
      loo_med_abs <- median(loo_ok$abs_delta_beta, na.rm = TRUE)
      strongest_gene <- loo_ok[which.max(abs_delta_beta)]$dropped_gene[1]
    } else {
      loo_retained <- NA_real_; loo_p_sig <- NA_real_; loo_max_abs <- NA_real_; loo_med_abs <- NA_real_; strongest_gene <- NA_character_
    }

    summary_rows[[mm]] <- data.table(
      dataset = dataset_name,
      model = mm,
      input_scale = obj$input_scale,
      n_program_genes_input = length(syn_genes),
      n_program_genes_mapped = length(mapped),
      n_program_genes_topdrop = nrow(expr_topdrop),
      topN = topN,
      topN_genes = paste(top_genes, collapse = ";"),
      full_beta = full_eff$beta[1],
      full_p = full_eff$p[1],
      full_formula = full_eff$formula[1],
      topdrop_beta = top_eff$beta[1],
      topdrop_p = top_eff$p[1],
      topdrop_delta_beta = top_eff$beta[1] - full_eff$beta[1],
      topdrop_direction_retained = !is.na(full_dir) & !is.na(top_dir) & full_dir == top_dir,
      loo_n_ok = nrow(loo_ok),
      loo_direction_retained_prop = loo_retained,
      loo_p_lt_0_05_prop = loo_p_sig,
      loo_max_abs_delta_beta = loo_max_abs,
      loo_median_abs_delta_beta = loo_med_abs,
      strongest_influence_gene = strongest_gene
    )

    if (nrow(infl_dt) > 0) {
      infl_dt[, dataset := dataset_name]
      infl_dt[, model := mm]
      infl_rows[[mm]] <- infl_dt
      make_plot(
        infl_dt[order(-abs_delta_beta)][1:min(.N, 80)],
        file.path(outdir, "plots", paste0("leave1out_delta_", dataset_name, "_", mm, ".pdf")),
        paste0(dataset_name, " ", mm, " leave-one-out")
      )
    }

    topdrop_plot_dt <- data.table(
      scenario = c("full","top5_dropout"),
      beta = c(full_eff$beta[1], top_eff$beta[1])
    )
    p2 <- ggplot(topdrop_plot_dt, aes(x = scenario, y = beta)) +
      geom_point(size = 2) + geom_line(aes(group = 1)) +
      theme_bw(base_size = 10) +
      labs(title = paste0(dataset_name, " ", mm, " full vs top5-dropout"), y = "ASD effect beta on z-scored program")
    ggsave(file.path(outdir, "plots", paste0("topdrop_", dataset_name, "_", mm, ".pdf")), plot = p2, width = 5.5, height = 4.2)

    sample_score_rows[[mm]] <- rbind(
      data.table(dataset = dataset_name, model = mm, scenario = "full", sample_id = score_dt$sample_id, program_score_raw = score_dt$program_score_raw, program_score_z = score_dt$program_score_z),
      data.table(dataset = dataset_name, model = mm, scenario = "top5_dropout", sample_id = score_dt_top$sample_id, program_score_raw = score_dt_top$program_score_raw, program_score_z = score_dt_top$program_score_z)
    )
  }

  list(
    summary = rbindlist(summary_rows, fill = TRUE),
    infl = rbindlist(infl_rows, fill = TRUE),
    top_expr = top_expr_dt,
    effects = NULL,
    sample_scores = rbindlist(sample_score_rows, fill = TRUE)
  )
}

main <- function() {
  root_dir <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
  outdir <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3A_synaptic_leave1out_v6")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

  syn_inventory <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3A_synaptic_leave1out/01_synaptic_gene_inventory.tsv")
  syn_dt <- fread(syn_inventory)
  gene_col <- intersect(c("gene","std_gene_id","gene_symbol","symbol"), names(syn_dt))
  if (length(gene_col) == 0) stop("Cannot find gene column in synaptic inventory.")
  syn_genes <- unique(na.omit(std_gene(syn_dt[[gene_col[1]]])))

  safe_write(data.table(
    time = timestamp(),
    root_dir = root_dir,
    outdir = outdir,
    syn_inventory = syn_inventory,
    n_synaptic_genes = length(syn_genes)
  ), file.path(outdir, "00_run_info.txt"))

  safe_write(data.table(gene = syn_genes), file.path(outdir, "01_synaptic_gene_inventory.tsv"))

  logmsg("Preparing datasets ...")
  ds1 <- prepare_gse102741(root_dir, syn_genes = syn_genes, outdir = outdir)
  ds2 <- prepare_gse64018(root_dir)
  ds3 <- prepare_gandal(root_dir)

  if (!is.null(ds3$mapping)) {
    safe_write(ds3$mapping$sum, file.path(outdir, "00b_gandal_mapping_summary.tsv"))
    idx <- ds3$mapping$idx
    idx[, mapped_to_symbol := !is.na(symbol_std)]
    safe_write(idx[1:min(.N, 5000)], file.path(outdir, "00c_gandal_mapping_examples.tsv"))
  }

  logmsg("Running leave-one-out sensitivity for GSE102741")
  r1 <- run_dataset("GSE102741", ds1, syn_genes, outdir)
  logmsg("Running leave-one-out sensitivity for GSE64018")
  r2 <- run_dataset("GSE64018", ds2, syn_genes, outdir)
  logmsg("Running leave-one-out sensitivity for Gandal2022")
  r3 <- run_dataset("Gandal2022", ds3, syn_genes, outdir)

  all_summary <- rbindlist(list(r1$summary, r2$summary, r3$summary), fill = TRUE)
  all_infl <- rbindlist(list(r1$infl, r2$infl, r3$infl), fill = TRUE)
  all_top  <- rbindlist(list(r1$top_expr, r2$top_expr, r3$top_expr), fill = TRUE)
  all_scores <- rbindlist(list(r1$sample_scores, r2$sample_scores, r3$sample_scores), fill = TRUE)

  safe_write(all_top, file.path(outdir, "02_top_expression_ranked_genes.tsv"))
  safe_write(all_infl, file.path(outdir, "05_synaptic_most_influential_genes.tsv"))
  safe_write(all_scores, file.path(outdir, "06_synaptic_full_vs_topdrop_sample_scores.tsv.gz"))
  safe_write(all_summary, file.path(outdir, "04_synaptic_leave1out_summary.tsv"))

  logmsg("Done. Main summary written to: ", file.path(outdir, "04_synaptic_leave1out_summary.tsv"))
}

main()
