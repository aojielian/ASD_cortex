
suppressPackageStartupMessages({
  library(data.table)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

logmsg <- function(...) {
  message("[", timestamp(), "] ", paste(..., collapse = ""))
}

safe_write <- function(dt, file) {
  fwrite(as.data.table(dt), file = file, sep = "\t", quote = FALSE, na = "NA")
}

std_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  toupper(x)
}

choose_gene_col <- function(dt) {
  cand <- intersect(names(dt), c("gene","Gene","symbol","SYMBOL","gene_symbol","std_gene_id","Gene.Symbol"))
  if (length(cand) > 0) return(cand[1])
  cand2 <- names(dt)[grepl("gene|symbol", names(dt), ignore.case = TRUE)]
  if (length(cand2) > 0) return(cand2[1])
  names(dt)[1]
}

extract_program_from_table <- function(program_path, target_program) {
  dt <- fread(program_path)
  target_program <- as.character(target_program)

  if (target_program %in% names(dt)) {
    genes <- unique(na.omit(std_gene(dt[[target_program]])))
    genes <- genes[nzchar(genes)]
    if (length(genes) > 0) {
      return(list(genes = genes, source = program_path, method = paste0("wide_column:", target_program)))
    }
  }

  prog_cols <- names(dt)[grepl("program|set|module|signature|name", names(dt), ignore.case = TRUE)]
  gene_cols <- intersect(names(dt), c("gene","Gene","symbol","SYMBOL","gene_symbol","std_gene_id","Gene.Symbol"))
  if (length(prog_cols) > 0 && length(gene_cols) > 0) {
    for (pc in prog_cols) {
      hit <- dt[as.character(get(pc)) == target_program]
      if (nrow(hit) > 0) {
        for (gc in gene_cols) {
          genes <- unique(na.omit(std_gene(hit[[gc]])))
          genes <- genes[nzchar(genes)]
          if (length(genes) > 0) {
            return(list(genes = genes, source = program_path, method = paste0("long_format:", pc, "->", gc)))
          }
        }
      }
    }
  }

  hit_rows <- rep(FALSE, nrow(dt))
  for (cc in names(dt)) {
    hit_rows <- hit_rows | (as.character(dt[[cc]]) == target_program)
  }
  if (any(hit_rows)) {
    for (gc in names(dt)) {
      vals <- unique(na.omit(std_gene(dt[[gc]][hit_rows])))
      vals <- vals[nzchar(vals)]
      vals <- vals[grepl("^[A-Z0-9._-]+$", vals)]
      if (length(vals) >= 10) {
        return(list(genes = vals, source = program_path, method = paste0("fallback_rows:", gc)))
      }
    }
  }

  stop("Could not extract program ", target_program, " from ", program_path)
}

is_numeric_matrix_like <- function(x) {
  if (!(is.matrix(x) || is.data.frame(x))) return(FALSE)
  y <- tryCatch(as.matrix(x), error = function(e) NULL)
  if (is.null(y)) return(FALSE)
  ok <- tryCatch({
    storage.mode(y) <- "numeric"
    TRUE
  }, error = function(e) FALSE)
  isTRUE(ok)
}

safe_num_mat <- function(x) {
  if (is.list(x) && !is.data.frame(x) && length(x) == 1L) x <- x[[1]]
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

collect_matrix_candidates <- function(x, path = "obj", depth = 0L, max_depth = 3L) {
  out <- list()
  if (is_numeric_matrix_like(x)) {
    out[[length(out) + 1L]] <- list(
      path = path,
      object = x,
      nrow = NROW(x),
      ncol = NCOL(x),
      rownames_present = !is.null(rownames(x)),
      colnames_present = !is.null(colnames(x))
    )
  }
  if (depth < max_depth && is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    for (i in seq_along(x)) {
      out <- c(out, collect_matrix_candidates(x[[i]], path = paste0(path, "$", nms[i]), depth = depth + 1L, max_depth = max_depth))
    }
  }
  out
}

choose_expr_from_rds <- function(obj, target_genes, outdir = NULL, prefix = "GSE102741") {
  cands <- collect_matrix_candidates(obj, path = "rds")
  if (length(cands) == 0) stop(prefix, ": no numeric matrix-like candidates found in RDS object.")
  inv <- rbindlist(lapply(cands, function(z) {
    rn <- rownames(z$object); cn <- colnames(z$object)
    rn_std <- if (!is.null(rn)) std_gene(rn) else character()
    cn_std <- if (!is.null(cn)) std_gene(cn) else character()
    data.table(
      path = z$path,
      nrow = as.numeric(z$nrow),
      ncol = as.numeric(z$ncol),
      rownames_present = z$rownames_present,
      colnames_present = z$colnames_present,
      row_target_overlap = sum(unique(rn_std[!is.na(rn_std)]) %in% target_genes),
      col_target_overlap = sum(unique(cn_std[!is.na(cn_std)]) %in% target_genes)
    )
  }), fill = TRUE)
  inv[, best_overlap := pmax(row_target_overlap, col_target_overlap)]
  setorderv(inv, c("best_overlap", "nrow", "ncol", "path"), c(-1L, -1L, -1L, 1L))
  best_path <- inv$path[1]
  best_idx <- which(vapply(cands, function(z) identical(z$path, best_path), logical(1)))[1]
  best_obj <- cands[[best_idx]]$object
  if (!is.null(outdir)) {
    safe_write(inv, file.path(outdir, paste0("08_", prefix, "_rds_candidate_inventory.tsv")))
  }
  list(expr = safe_num_mat(best_obj), chosen_path = best_path, inventory = inv)
}

extract_dx <- function(meta) {
  cand <- names(meta)[tolower(names(meta)) %in% c("diagnosis","dx","dx_group","asd","group","condition","status")]
  if (length(cand) == 0) cand <- names(meta)[grepl("diagnosis|dx|asd|group|condition|status", names(meta), ignore.case = TRUE)]
  if (length(cand) == 0) stop("Could not detect diagnosis column.")
  dx_col <- cand[1]
  dx_raw <- as.character(meta[[dx_col]])
  dx <- rep(NA_character_, length(dx_raw))
  dx[grepl("^asd$|aut", dx_raw, ignore.case = TRUE)] <- "ASD"
  dx[grepl("^control$|^ctl$|^ctrl$|normal|non[-_ ]?asd", dx_raw, ignore.case = TRUE)] <- "Control"
  suppressWarnings(num <- as.numeric(dx_raw))
  if (any(!is.na(num)) && all(na.omit(num) %in% c(0,1))) {
    dx[!is.na(num) & num == 1] <- "ASD"
    dx[!is.na(num) & num == 0] <- "Control"
  }
  list(dx = dx, dx_col = dx_col)
}

extract_sample_ids <- function(meta, expr_colnames = NULL) {
  cand <- names(meta)[tolower(names(meta)) %in% c("sample_id","sample","sampleid","iid","id","run","libraryid","rnaseq_id")]
  if (length(cand) == 0) cand <- names(meta)[grepl("sample|iid|library|id$", names(meta), ignore.case = TRUE)]
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
  if (orientation == "transpose") expr <- t(expr)
  list(expr = expr, sample_col = best$candidate_col, orientation = orientation, alignment = score_dt)
}

clean_gandal_ensg <- function(x) {
  x <- as.character(x)
  x <- sub("_[0-9]+$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

map_ensembl_to_symbol <- function(ensg_ids) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("This script requires AnnotationDbi and org.Hs.eg.db for Gandal ENSG mapping.")
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

detect_covariates <- function(meta) {
  covs <- names(meta)[grepl("region|sex|gender|age|rin|pmi|batch|site", names(meta), ignore.case = TRUE)]
  covs <- unique(setdiff(covs, c("sample_id","dx_group")))
  covs
}

pick_best_covariates <- function(meta, max_covs = 4L) {
  priority <- c("region","sex","age","rin","pmi","batch","site")
  covs <- detect_covariates(meta)
  if (length(covs) == 0) return(character())
  usable <- vapply(covs, function(cc) {
    x <- meta[[cc]]
    if (is.null(x)) return(FALSE)
    if (is.list(x) && !is.data.frame(x)) return(FALSE)
    vals <- tryCatch(unique(na.omit(x)), error = function(e) unique(as.character(x)))
    length(vals) > 1
  }, logical(1))
  covs <- covs[usable]
  if (length(covs) == 0) return(character())
  score <- vapply(covs, function(cc) {
    nm <- tolower(cc)
    hit <- match(TRUE, vapply(priority, function(pp) grepl(pp, nm, fixed = TRUE), logical(1)))
    ifelse(is.na(hit), 100, hit)
  }, numeric(1))
  covs <- covs[order(score, covs)]
  covs[seq_len(min(length(covs), max_covs))]
}

build_formula <- function(covs) {
  if (length(covs) == 0) return(as.formula("program_score_z ~ dx_group"))
  reformulate(c("dx_group", covs), response = "program_score_z")
}

fit_model <- function(df, formula_obj) {
  dd <- copy(df)
  keep_terms <- unique(c("sample_id", "dx_group", all.vars(formula_obj)))
  keep_terms <- keep_terms[keep_terms %in% names(dd)]
  dd <- dd[, ..keep_terms]
  dd <- na.omit(dd)
  if (nrow(dd) < 8) {
    return(data.table(
      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
      n_samples = nrow(dd), n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE), n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
      formula = deparse(formula_obj), status = "too_few_complete_cases"
    ))
  }
  fit <- tryCatch(lm(formula_obj, data = dd), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.table(
      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
      n_samples = nrow(dd), n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE), n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
      formula = deparse(formula_obj), status = "lm_failed"
    ))
  }
  cf <- summary(fit)$coefficients
  rn <- rownames(cf)
  target <- rn[grepl("dx_group", rn) & grepl("ASD", rn)]
  if (length(target) == 0) {
    return(data.table(
      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
      n_samples = nrow(dd), n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE), n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
      formula = deparse(formula_obj), status = "dx_term_missing"
    ))
  }
  tt <- target[1]
  beta <- unname(cf[tt, "Estimate"])
  se <- unname(cf[tt, "Std. Error"])
  p <- unname(cf[tt, grep("Pr\\(", colnames(cf))])
  data.table(
    beta = beta, se = se, ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se, p = p,
    n_samples = nrow(dd), n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE), n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
    formula = deparse(formula_obj), status = "ok"
  )
}

make_pairwise_overlap <- function(dt_mapped_list) {
  programs <- names(dt_mapped_list)
  out <- list()
  idx <- 1L
  for (i in seq_along(programs)) {
    for (j in i:length(programs)) {
      p1 <- programs[i]; p2 <- programs[j]
      g1 <- sort(unique(dt_mapped_list[[p1]]))
      g2 <- sort(unique(dt_mapped_list[[p2]]))
      inter <- intersect(g1, g2)
      uni <- union(g1, g2)
      out[[idx]] <- data.table(
        program1 = p1,
        program2 = p2,
        n1 = length(g1),
        n2 = length(g2),
        n_intersection = length(inter),
        n_union = length(uni),
        jaccard = ifelse(length(uni) > 0, length(inter) / length(uni), NA_real_),
        identical_sets = identical(g1, g2),
        p1_subset_p2 = all(g1 %in% g2),
        p2_subset_p1 = all(g2 %in% g1)
      )
      idx <- idx + 1L
    }
  }
  rbindlist(out, fill = TRUE)
}

score_program <- function(expr, genes) {
  mapped <- sort(intersect(rownames(expr), genes))
  raw <- colMeans(expr[mapped, , drop = FALSE], na.rm = TRUE)
  z <- as.numeric(scale(raw))
  data.table(sample_id = colnames(expr), program_score_raw = raw, program_score_z = z, n_genes_mapped = length(mapped))
}

prepare_gse102741 <- function(root_dir, target_genes, outdir) {
  expr_path <- file.path(root_dir, "step01_prepare/rds/22_GSE102741_geneSymbol_log2cpm.rds")
  meta_path <- file.path(root_dir, "step01_prepare/tables/21_GSE102741_sample_metadata.tsv")
  obj <- readRDS(expr_path)
  chosen <- if (is.matrix(obj) || is.data.frame(obj)) {
    list(expr = safe_num_mat(obj), chosen_path = "rds")
  } else {
    choose_expr_from_rds(obj, target_genes = target_genes, outdir = outdir, prefix = "GSE102741")
  }
  expr <- chosen$expr
  meta <- fread(meta_path)
  orient <- choose_sample_column_and_orient(expr, meta)
  expr <- orient$expr
  safe_write(orient$alignment, file.path(outdir, "09_GSE102741_sample_alignment.tsv"))
  safe_write(data.table(chosen_path = chosen$chosen_path, orientation = orient$orientation, sample_col = orient$sample_col),
             file.path(outdir, "10_GSE102741_chosen_expr_object.tsv"))
  sid <- extract_sample_ids(meta, expr_colnames = colnames(expr))
  dxi <- extract_dx(meta)
  meta[, sample_id := sid$sample_id]
  meta[, dx_group := dxi$dx]
  meta <- meta[sample_id %in% colnames(expr)]
  expr <- expr[, meta$sample_id, drop = FALSE]
  rownames(expr) <- std_gene(rownames(expr))
  list(dataset = "GSE102741", expr = expr, meta = meta, raw_gene_id_type = "gene_symbol_like",
       mapping_strategy = "direct_standardized_symbol_from_chosen_RDS_object", chosen_object = chosen$chosen_path)
}

prepare_gse64018 <- function(root_dir) {
  expr_path <- file.path(root_dir, "step01_prepare/tables/122_GSE64018_collapsed_counts_by_symbol.tsv.gz")
  meta_path <- file.path(root_dir, "step01_prepare/meta/120_GSE64018_explicit_sample_mapping.tsv")
  expr_dt <- fread(expr_path)
  gene_col <- names(expr_dt)[1]
  genes <- std_gene(expr_dt[[gene_col]])
  expr <- as.matrix(expr_dt[, -1, with = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- genes
  lib <- colSums(expr, na.rm = TRUE)
  expr <- log2(t(t(expr) / pmax(lib, 1) * 1e6) + 1)
  meta <- fread(meta_path)
  sid <- extract_sample_ids(meta, expr_colnames = colnames(expr))
  dxi <- extract_dx(meta)
  meta[, sample_id := sid$sample_id]
  meta[, dx_group := dxi$dx]
  meta <- meta[sample_id %in% colnames(expr)]
  expr <- expr[, meta$sample_id, drop = FALSE]
  list(dataset = "GSE64018", expr = expr, meta = meta, raw_gene_id_type = "gene_symbol_like",
       mapping_strategy = "direct_standardized_symbol_from_collapsed_counts_first_column", chosen_object = gene_col)
}

prepare_gandal <- function() {
  rdata_path <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  datExpr <- get("datExpr", envir = env)
  datMeta <- as.data.table(get("datMeta", envir = env))
  datExpr <- safe_num_mat(datExpr)
  raw_ids <- rownames(datExpr)
  ensg_clean <- clean_gandal_ensg(raw_ids)
  map_dt <- map_ensembl_to_symbol(ensg_clean)
  map_dt[, ensembl_std := std_gene(ensembl)]
  map_dt[, symbol_std := std_gene(symbol)]
  idx_dt <- data.table(raw_gene_id = raw_ids, ensembl_clean = ensg_clean, ensembl_std = std_gene(ensg_clean))
  idx_dt <- merge(idx_dt, map_dt[, .(ensembl_std, symbol_std)], by = "ensembl_std", all.x = TRUE)
  expr <- aggregate_by_symbol_mean(datExpr, idx_dt$symbol_std)
  sid <- extract_sample_ids(datMeta, expr_colnames = colnames(expr))
  dxi <- extract_dx(datMeta)
  datMeta[, sample_id := sid$sample_id]
  datMeta[, dx_group := dxi$dx]
  datMeta <- datMeta[sample_id %in% colnames(expr)]
  expr <- expr[, datMeta$sample_id, drop = FALSE]
  list(dataset = "Gandal2022", expr = expr, meta = datMeta, raw_gene_id_type = "ENSG.version_suffix",
       mapping_strategy = "clean_ENSG_suffix_and_version_then_ENSEMBL_to_SYMBOL_via_org.Hs.eg.db", chosen_object = "datExpr",
       mapping_index = idx_dt)
}

main <- function() {
  root_dir <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
  outdir <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3D_best_available_covariate_models")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  sfari_path <- file.path(root_dir, "step01_prepare/tables/01_SFARI_primary_Sand1_standardized.tsv")
  program_path <- file.path(root_dir, "step01_prepare/tables/20_midPrenatal_core_programs.tsv")
  syn_inventory_path <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3A_synaptic_leave1out_v6/01_synaptic_gene_inventory.tsv")

  sfari_dt <- fread(sfari_path)
  sfari_col <- choose_gene_col(sfari_dt)
  sfari_all <- unique(na.omit(std_gene(sfari_dt[[sfari_col]])))
  sfari_all <- sfari_all[nzchar(sfari_all)]

  syn_dt <- fread(syn_inventory_path)
  syn_col <- choose_gene_col(syn_dt)
  sfari_syn <- unique(na.omit(std_gene(syn_dt[[syn_col]])))
  sfari_syn <- sfari_syn[nzchar(sfari_syn)]

  top20_obj <- extract_program_from_table(program_path, "midPrenatal_SFARI_top20")
  top20 <- unique(na.omit(std_gene(top20_obj$genes)))
  top20 <- top20[nzchar(top20)]

  program_inventory <- rbindlist(list(
    data.table(program = "SFARI_all", gene = sfari_all, source = sfari_path, method = paste0("column:", sfari_col)),
    data.table(program = "SFARI_all_synaptic", gene = sfari_syn, source = syn_inventory_path, method = paste0("column:", syn_col)),
    data.table(program = "midPrenatal_SFARI_top20", gene = top20, source = top20_obj$source, method = top20_obj$method)
  ), fill = TRUE)
  safe_write(program_inventory, file.path(outdir, "01_program_inventory.tsv"))

  target_genes <- unique(program_inventory$gene)

  logmsg("Preparing datasets ...")
  ds_list <- list(
    prepare_gse102741(root_dir, target_genes = target_genes, outdir = outdir),
    prepare_gse64018(root_dir),
    prepare_gandal()
  )

  if (!is.null(ds_list[[3]]$mapping_index)) {
    safe_write(ds_list[[3]]$mapping_index, file.path(outdir, "11_Gandal_row_mapping_table.tsv.gz"))
  }

  gene_space_rows <- list()
  meta_inv_rows <- list()
  formula_rows <- list()
  model_rows <- list()
  score_rows <- list()
  mapped_gene_rows <- list()
  reviewer_rows <- list()
  overlap_rows <- list()

  for (ds in ds_list) {
    dataset <- ds$dataset
    expr <- ds$expr
    meta <- copy(ds$meta)
    meta[, dx_group := factor(dx_group, levels = c("Control","ASD"))]

    expr_genes <- unique(rownames(expr))
    gene_space_rows[[dataset]] <- data.table(
      dataset = dataset,
      raw_gene_id_type = ds$raw_gene_id_type,
      mapping_strategy = ds$mapping_strategy,
      chosen_object = ds$chosen_object,
      n_dataset_genes = length(expr_genes),
      n_samples = ncol(expr)
    )

    meta_inv_rows[[dataset]] <- rbindlist(lapply(names(meta), function(cc) {
      x <- meta[[cc]]
      data.table(
        dataset = dataset,
        variable = cc,
        class = paste(class(x), collapse = "|"),
        n_nonmissing = sum(!is.na(x)),
        n_unique = uniqueN(x[!is.na(x)]),
        is_candidate_covariate = grepl("region|sex|gender|age|rin|pmi|batch|site", cc, ignore.case = TRUE)
      )
    }), fill = TRUE)

    best_covs <- pick_best_covariates(meta)
    formula_rows[[dataset]] <- data.table(
      dataset = dataset,
      model = c("diagnosis_only", "best_available"),
      covariates = c("", paste(best_covs, collapse = ";")),
      formula = c("program_score_z ~ dx_group", deparse(build_formula(best_covs))),
      fallback_to_diagnosis_only = c(FALSE, length(best_covs) == 0)
    )

    mapped_list_ds <- list()

    for (pg in unique(program_inventory$program)) {
      pg_genes <- unique(program_inventory[program == pg, gene])
      mapped <- sort(intersect(expr_genes, pg_genes))
      dropped <- sort(setdiff(pg_genes, expr_genes))
      mapped_list_ds[[pg]] <- mapped

      sc <- score_program(expr, pg_genes)
      sc[, dataset := dataset]
      sc[, program := pg]
      sc <- merge(sc, meta, by = "sample_id", all.x = TRUE)
      score_rows[[paste(dataset, pg, sep = "::")]] <- sc

      diag_formula <- as.formula("program_score_z ~ dx_group")
      best_formula <- build_formula(best_covs)
      fits <- list(
        diagnosis_only = fit_model(sc, diag_formula),
        best_available = fit_model(sc, best_formula)
      )

      overlap_syn <- length(intersect(mapped, mapped_list_ds[["SFARI_all_synaptic"]]))
      jacc_syn <- ifelse(length(union(mapped, mapped_list_ds[["SFARI_all_synaptic"]])) > 0,
                         length(intersect(mapped, mapped_list_ds[["SFARI_all_synaptic"]])) / length(union(mapped, mapped_list_ds[["SFARI_all_synaptic"]])),
                         NA_real_)
      collapse_flag <- identical(sort(mapped), sort(mapped_list_ds[["SFARI_all_synaptic"]]))

      for (mm in names(fits)) {
        rr <- copy(fits[[mm]])
        rr[, dataset := dataset]
        rr[, program := pg]
        rr[, model := mm]
        rr[, n_program_input := length(pg_genes)]
        rr[, n_program_mapped := length(mapped)]
        rr[, n_program_dropped := length(dropped)]
        rr[, mapped_fraction_of_program := ifelse(length(pg_genes) > 0, length(mapped) / length(pg_genes), NA_real_)]
        rr[, overlap_with_synaptic := overlap_syn]
        rr[, jaccard_with_synaptic := jacc_syn]
        rr[, mapped_set_identical_to_SFARI_all_synaptic := collapse_flag]
        rr[, note := fifelse(pg == "SFARI_all" & collapse_flag, "Mapped set collapsed to SFARI_all_synaptic in this dataset",
                             fifelse(pg == "SFARI_all_synaptic", "Reference synaptic mapped set",
                                     fifelse(pg == "midPrenatal_SFARI_top20", "Developmental-context program", "")))]
        model_rows[[paste(dataset, pg, mm, sep = "::")]] <- rr
      }

      if (length(mapped) > 0) {
        mapped_gene_rows[[paste(dataset, pg, sep = "::")]] <- data.table(dataset = dataset, program = pg, gene = mapped)
      } else {
        mapped_gene_rows[[paste(dataset, pg, sep = "::")]] <- data.table(dataset = character(), program = character(), gene = character())
      }
    }

    overlap_ds <- make_pairwise_overlap(mapped_list_ds)
    overlap_ds[, dataset := dataset]
    overlap_rows[[dataset]] <- overlap_ds
  }

  gene_space_dt <- rbindlist(gene_space_rows, fill = TRUE)
  meta_inventory_dt <- rbindlist(meta_inv_rows, fill = TRUE)
  formula_dt <- rbindlist(formula_rows, fill = TRUE)
  model_long_dt <- rbindlist(model_rows, fill = TRUE)
  score_dt <- rbindlist(score_rows, fill = TRUE)
  mapped_gene_dt <- rbindlist(mapped_gene_rows, fill = TRUE)
  overlap_dt <- rbindlist(overlap_rows, fill = TRUE)

  # wide/freeze summary
  diag_dt <- copy(model_long_dt[model == "diagnosis_only"])
  setnames(diag_dt,
           old = c("beta","se","ci_low","ci_high","p","formula","status"),
           new = c("diagnosis_beta","diagnosis_se","diagnosis_ci_low","diagnosis_ci_high","diagnosis_p","diagnosis_formula","diagnosis_status"))
  best_dt <- copy(model_long_dt[model == "best_available"])
  setnames(best_dt,
           old = c("beta","se","ci_low","ci_high","p","formula","status"),
           new = c("best_beta","best_se","best_ci_low","best_ci_high","best_p","best_formula","best_status"))

  keep_diag <- c("dataset","program","n_program_input","n_program_mapped","n_program_dropped","mapped_fraction_of_program",
                 "overlap_with_synaptic","jaccard_with_synaptic","mapped_set_identical_to_SFARI_all_synaptic","note",
                 "diagnosis_beta","diagnosis_se","diagnosis_ci_low","diagnosis_ci_high","diagnosis_p","diagnosis_formula","diagnosis_status")
  keep_best <- c("dataset","program","best_beta","best_se","best_ci_low","best_ci_high","best_p","best_formula","best_status")
  freeze_dt <- merge(diag_dt[, ..keep_diag], best_dt[, ..keep_best], by = c("dataset","program"), all = TRUE)
  freeze_dt <- merge(freeze_dt, gene_space_dt[, .(dataset, n_dataset_genes)], by = "dataset", all.x = TRUE)
  freeze_dt[, direction_diagnosis := fifelse(is.na(diagnosis_beta), NA_character_, fifelse(diagnosis_beta > 0, "positive", fifelse(diagnosis_beta < 0, "negative", "zero")))]
  freeze_dt[, direction_best := fifelse(is.na(best_beta), NA_character_, fifelse(best_beta > 0, "positive", fifelse(best_beta < 0, "negative", "zero")))]
  freeze_dt[, direction_retained_best_vs_diagnosis := !is.na(direction_diagnosis) & !is.na(direction_best) & direction_diagnosis == direction_best]
  freeze_dt <- merge(freeze_dt, formula_dt[model == "best_available", .(dataset, best_available_covariates = covariates, best_available_fallback_to_diagnosis_only = fallback_to_diagnosis_only)], by = "dataset", all.x = TRUE)

  reviewer_dt <- copy(freeze_dt)
  reviewer_dt[, interpretation_note := fifelse(
    dataset %in% c("GSE102741","GSE64018") & program == "SFARI_all" & mapped_set_identical_to_SFARI_all_synaptic,
    "Input broad anchor collapsed to mapped synaptic set in this dataset",
    fifelse(dataset == "Gandal2022" & program == "SFARI_all",
            "Mapped broad anchor retained beyond synaptic core in Gandal2022",
            fifelse(program == "midPrenatal_SFARI_top20",
                    "Developmental-context comparator remained distinct from synaptic set",
                    note)))]
  setorderv(reviewer_dt, c("dataset","program"), c(1L,1L))

  # outputs
  safe_write(data.table(
    time = timestamp(),
    root_dir = root_dir,
    outdir = outdir,
    sfari_path = sfari_path,
    program_path = program_path,
    syn_inventory_path = syn_inventory_path
  ), file.path(outdir, "00_run_info.txt"))

  safe_write(gene_space_dt, file.path(outdir, "02_dataset_gene_space_summary.tsv"))
  safe_write(meta_inventory_dt, file.path(outdir, "03_metadata_inventory.tsv"))
  safe_write(formula_dt, file.path(outdir, "04_formula_inventory.tsv"))
  safe_write(model_long_dt, file.path(outdir, "05_model_results_long.tsv"))
  safe_write(freeze_dt, file.path(outdir, "06_covariate_freeze_summary.tsv"))
  safe_write(reviewer_dt, file.path(outdir, "07_reviewer_facing_covariate_summary.tsv"))
  safe_write(mapped_gene_dt, file.path(outdir, "12_mapped_program_genes.tsv.gz"))
  safe_write(overlap_dt, file.path(outdir, "13_pairwise_program_overlap_by_dataset.tsv"))
  safe_write(score_dt, file.path(outdir, "14_sample_score_table.tsv.gz"))

  logmsg("Done. Reviewer-facing covariate summary written to: ", file.path(outdir, "07_reviewer_facing_covariate_summary.tsv"))
}

main()
