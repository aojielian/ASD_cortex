suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

logmsg <- function(...) {
  message("[", timestamp(), "] ", paste(..., collapse = ""))
}

std_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  toupper(x)
}

safe_write <- function(dt, file) {
  fwrite(as.data.table(dt), file = file, sep = "\t", quote = FALSE, na = "NA")
}

extract_first_existing <- function(nms, patterns) {
  nms[which(tolower(nms) %in% tolower(patterns))[1]]
}

choose_gene_col <- function(dt) {
  cand <- intersect(names(dt), c("gene", "Gene", "symbol", "SYMBOL", "gene_symbol", "std_gene_id", "Gene.Symbol"))
  if (length(cand) > 0) return(cand[1])
  cand2 <- names(dt)[grepl("gene|symbol", names(dt), ignore.case = TRUE)]
  if (length(cand2) > 0) return(cand2[1])
  names(dt)[1]
}

extract_program_from_table <- function(program_path, target_program) {
  dt <- fread(program_path)
  target_program <- as.character(target_program)

  # wide format: target program is a column name
  if (target_program %in% names(dt)) {
    genes <- unique(na.omit(std_gene(dt[[target_program]])))
    genes <- genes[nzchar(genes)]
    if (length(genes) > 0) {
      return(list(genes = genes, source = program_path, method = paste0("wide_column:", target_program)))
    }
  }

  # long format: one column contains program name, another contains genes
  prog_cols <- names(dt)[grepl("program|set|module|signature|name", names(dt), ignore.case = TRUE)]
  gene_cols <- intersect(names(dt), c("gene", "Gene", "symbol", "SYMBOL", "gene_symbol", "std_gene_id", "Gene.Symbol"))
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

  # fallback: find rows containing target program anywhere, then take gene-like columns
  hit_rows <- rep(FALSE, nrow(dt))
  for (cc in names(dt)) {
    hit_rows <- hit_rows | (as.character(dt[[cc]]) == target_program)
  }
  if (any(hit_rows)) {
    candidate_cols <- setdiff(names(dt), names(dt)[sapply(names(dt), function(cc) any(as.character(dt[[cc]]) == target_program))])
    if (length(candidate_cols) == 0) candidate_cols <- names(dt)
    for (gc in candidate_cols) {
      vals <- unique(na.omit(std_gene(dt[[gc]][hit_rows])))
      vals <- vals[nzchar(vals)]
      vals <- vals[grepl("^[A-Z0-9._-]+$", vals)]
      if (length(vals) >= 10) {
        return(list(genes = vals, source = program_path, method = paste0("fallback_rows:", gc)))
      }
    }
  }

  stop("Could not extract program: ", target_program, " from ", program_path)
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

build_best_formula <- function(df) {
  priority <- c("region","sex","age","rin","pmi","batch","site")
  covs <- names(df)[grepl("region|sex|gender|age|rin|pmi|batch|site", names(df), ignore.case = TRUE)]
  covs <- unique(setdiff(covs, c("sample_id","dx_group","program","program_score_raw","program_score_z")))
  if (length(covs) == 0) return(as.formula("program_score_z ~ dx_group"))
  score <- vapply(covs, function(cc) {
    x <- df[[cc]]
    if (is.list(x) && !is.data.frame(x)) return(0)
    if (length(unique(na.omit(x))) <= 1) return(0)
    nm <- tolower(cc)
    rank_match <- match(TRUE, vapply(priority, function(pp) grepl(pp, nm, fixed = TRUE), logical(1)))
    ifelse(is.na(rank_match), 100, rank_match)
  }, numeric(1))
  covs <- covs[order(score, covs)]
  covs <- covs[seq_len(min(length(covs), 4))]
  reformulate(c("dx_group", covs), response = "program_score_z")
}

fit_model <- function(df, formula_obj) {
  dd <- copy(df)
  dd <- dd[!is.na(program_score_z) & !is.na(dx_group)]
  if (nrow(dd) < 8) {
    return(list(fit = NULL, stat = data.table(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                                             p = NA_real_, n_samples = nrow(dd),
                                             n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE),
                                             n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
                                             formula = deparse(formula_obj), status = "too_few_samples")))
  }
  keep_terms <- unique(c("sample_id", "dx_group", all.vars(formula_obj)))
  keep_terms <- keep_terms[keep_terms %in% names(dd)]
  dd <- dd[, ..keep_terms]
  dd <- na.omit(dd)
  if (nrow(dd) < 8) {
    return(list(fit = NULL, stat = data.table(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                                             p = NA_real_, n_samples = nrow(dd),
                                             n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE),
                                             n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
                                             formula = deparse(formula_obj), status = "too_few_complete_cases")))
  }
  fit <- tryCatch(lm(formula_obj, data = dd), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(fit = NULL, stat = data.table(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                                             p = NA_real_, n_samples = nrow(dd),
                                             n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE),
                                             n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
                                             formula = deparse(formula_obj), status = "lm_failed")))
  }
  cf <- summary(fit)$coefficients
  rn <- rownames(cf)
  target <- rn[grepl("dx_group", rn) & grepl("ASD", rn)]
  if (length(target) == 0) {
    return(list(fit = fit, stat = data.table(beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                                            p = NA_real_, n_samples = nrow(dd),
                                            n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE),
                                            n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
                                            formula = deparse(formula_obj), status = "dx_term_missing")))
  }
  tt <- target[1]
  beta <- unname(cf[tt, "Estimate"])
  se <- unname(cf[tt, "Std. Error"])
  p <- unname(cf[tt, grep("Pr\\(", colnames(cf))])
  stat <- data.table(beta = beta, se = se, ci_low = beta - 1.96*se, ci_high = beta + 1.96*se,
                     p = p, n_samples = nrow(dd),
                     n_asd = sum(dd$dx_group == "ASD", na.rm = TRUE),
                     n_control = sum(dd$dx_group == "Control", na.rm = TRUE),
                     formula = deparse(formula_obj), status = "ok")
  list(fit = fit, stat = stat, model_df = dd, dx_term = tt)
}

score_program <- function(expr_log2cpm, genes) {
  gg <- intersect(rownames(expr_log2cpm), genes)
  raw <- colMeans(expr_log2cpm[gg, , drop = FALSE], na.rm = TRUE)
  z <- as.numeric(scale(raw))
  data.table(sample_id = colnames(expr_log2cpm), n_genes_mapped = length(gg), program_score_raw = raw, program_score_z = z)
}

main <- function() {
  root_dir <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
  outdir <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3B_GSE64018_divergence")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

  expr_path <- file.path(root_dir, "step01_prepare/tables/122_GSE64018_collapsed_counts_by_symbol.tsv.gz")
  meta_path <- file.path(root_dir, "step01_prepare/meta/120_GSE64018_explicit_sample_mapping.tsv")
  sfari_path <- file.path(root_dir, "step01_prepare/tables/01_SFARI_primary_Sand1_standardized.tsv")
  program_path <- file.path(root_dir, "step01_prepare/tables/20_midPrenatal_core_programs.tsv")
  syn_inventory <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3A_synaptic_leave1out_v6/01_synaptic_gene_inventory.tsv")

  stopifnot(file.exists(expr_path), file.exists(meta_path), file.exists(sfari_path), file.exists(program_path), file.exists(syn_inventory))

  sfari_dt <- fread(sfari_path)
  sfari_col <- choose_gene_col(sfari_dt)
  sfari_all <- unique(na.omit(std_gene(sfari_dt[[sfari_col]])))
  sfari_all <- sfari_all[nzchar(sfari_all)]

  syn_dt <- fread(syn_inventory)
  syn_col <- choose_gene_col(syn_dt)
  sfari_syn <- unique(na.omit(std_gene(syn_dt[[syn_col]])))
  sfari_syn <- sfari_syn[nzchar(sfari_syn)]

  top20 <- extract_program_from_table(program_path, "midPrenatal_SFARI_top20")$genes

  program_inventory <- rbindlist(list(
    data.table(program = "SFARI_all", gene = sfari_all),
    data.table(program = "midPrenatal_SFARI_top20", gene = top20),
    data.table(program = "SFARI_all_synaptic", gene = sfari_syn)
  ))
  safe_write(program_inventory, file.path(outdir, "01_program_inventory.tsv"))

  counts_dt <- fread(expr_path)
  gene_col <- names(counts_dt)[1]
  genes <- std_gene(counts_dt[[gene_col]])
  expr <- as.matrix(counts_dt[, -1, with = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- genes

  libsize <- colSums(expr, na.rm = TRUE)
  expr_log2cpm <- log2(t(t(expr) / pmax(libsize, 1) * 1e6) + 1)

  meta <- fread(meta_path)
  sid <- extract_sample_ids(meta, expr_colnames = colnames(expr_log2cpm))
  dxi <- extract_dx(meta)
  meta[, sample_id := sid$sample_id]
  meta[, dx_group := dxi$dx]
  meta[, dx_group := factor(dx_group, levels = c("Control", "ASD"))]

  align_dt <- data.table(
    sample_id = meta$sample_id,
    in_expr = meta$sample_id %in% colnames(expr_log2cpm),
    dx_group = meta$dx_group
  )
  safe_write(align_dt, file.path(outdir, "02_sample_alignment.tsv"))

  meta_inv <- rbindlist(lapply(names(meta), function(cc) {
    x <- meta[[cc]]
    data.table(
      variable = cc,
      class = paste(class(x), collapse = "|"),
      n_nonmissing = sum(!is.na(x)),
      n_unique = uniqueN(x[!is.na(x)]),
      is_numeric = is.numeric(x) || suppressWarnings(all(!is.na(as.numeric(as.character(na.omit(x)))))),
      is_candidate_covariate = grepl("region|sex|gender|age|rin|pmi|batch|site", cc, ignore.case = TRUE)
    )
  }), fill = TRUE)
  safe_write(meta_inv, file.path(outdir, "03_metadata_inventory.tsv"))

  common_samples <- intersect(colnames(expr_log2cpm), meta$sample_id)
  expr_log2cpm <- expr_log2cpm[, common_samples, drop = FALSE]
  meta <- meta[match(common_samples, sample_id)]
  meta <- meta[!is.na(dx_group) & dx_group %in% c("Control", "ASD")]

  score_rows <- list()
  model_rows <- list()
  loo_rows <- list()
  infl_rows <- list()
  meta_assoc_rows <- list()

  for (pg in unique(program_inventory$program)) {
    genes_pg <- program_inventory[program == pg, unique(gene)]
    sc <- score_program(expr_log2cpm, genes_pg)
    sc[, program := pg]
    sc <- merge(sc, meta, by = "sample_id", all.x = TRUE)
    score_rows[[pg]] <- sc

    form_diag <- as.formula("program_score_z ~ dx_group")
    form_best <- build_best_formula(sc)
    forms <- list(diagnosis_only = form_diag, best_available = form_best)

    for (mm in names(forms)) {
      fm <- fit_model(sc, forms[[mm]])
      stat <- copy(fm$stat)
      stat[, program := pg]
      stat[, model := mm]
      stat[, n_genes_input := length(genes_pg)]
      stat[, n_genes_mapped := unique(sc$n_genes_mapped)[1]]
      model_rows[[paste(pg, mm, sep = "::")]] <- stat

      if (!is.null(fm$fit) && stat$status[1] == "ok") {
        mdf <- copy(fm$model_df)
        fit <- fm$fit
        infl <- lm.influence(fit, do.coef = TRUE)
        cooks <- cooks.distance(fit)

        dx_term <- fm$dx_term
        dfb <- tryCatch(dfbetas(fit), error = function(e) NULL)
        if (!is.null(dfb)) {
          if (is.matrix(dfb)) {
            dx_dfb <- dfb[, dx_term]
          } else {
            dx_dfb <- as.numeric(dfb)
          }
        } else {
          dx_dfb <- rep(NA_real_, nrow(mdf))
        }

        inf_dt <- data.table(
          sample_id = mdf$sample_id,
          program = pg,
          model = mm,
          dx_group = mdf$dx_group,
          cooks_distance = as.numeric(cooks),
          abs_dfbeta_dx = abs(as.numeric(dx_dfb)),
          dfbeta_dx = as.numeric(dx_dfb),
          fitted = as.numeric(fitted(fit)),
          residual = as.numeric(resid(fit))
        )
        infl_rows[[paste(pg, mm, sep = "::")]] <- inf_dt

        # leave-one-sample-out
        full_beta <- stat$beta[1]
        loo_dt <- rbindlist(lapply(seq_len(nrow(mdf)), function(i) {
          dd <- mdf[-i]
          fit_i <- tryCatch(lm(forms[[mm]], data = dd), error = function(e) NULL)
          if (is.null(fit_i)) {
            data.table(sample_id = mdf$sample_id[i], beta = NA_real_, p = NA_real_, status = "lm_failed")
          } else {
            cf <- summary(fit_i)$coefficients
            rn <- rownames(cf)
            target <- rn[grepl("dx_group", rn) & grepl("ASD", rn)]
            if (length(target) == 0) {
              data.table(sample_id = mdf$sample_id[i], beta = NA_real_, p = NA_real_, status = "dx_term_missing")
            } else {
              tt <- target[1]
              data.table(sample_id = mdf$sample_id[i],
                         beta = unname(cf[tt, "Estimate"]),
                         p = unname(cf[tt, grep("Pr\\(", colnames(cf))]),
                         status = "ok")
            }
          }
        }), fill = TRUE)
        loo_dt[, program := pg]
        loo_dt[, model := mm]
        loo_dt[, full_beta := full_beta]
        loo_dt[, delta_beta := beta - full_beta]
        loo_dt[, abs_delta_beta := abs(delta_beta)]
        loo_dt <- merge(loo_dt, unique(mdf[, .(sample_id, dx_group)]), by = "sample_id", all.x = TRUE)
        loo_rows[[paste(pg, mm, sep = "::")]] <- loo_dt

        # plots
        p1 <- ggplot(mdf, aes(x = dx_group, y = program_score_z)) +
          geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitter(width = 0.08), size = 2) +
          theme_bw(base_size = 10) +
          labs(title = paste0("GSE64018 ", pg, " score by diagnosis"), x = NULL, y = "Program score (z)")
        ggsave(file.path(outdir, "plots", paste0("score_box_", pg, "_", mm, ".pdf")), p1, width = 5.2, height = 4.2)

        p2 <- ggplot(loo_dt[status == "ok"], aes(x = reorder(sample_id, abs_delta_beta), y = abs_delta_beta, color = dx_group)) +
          geom_point(size = 2) + coord_flip() + theme_bw(base_size = 9) +
          labs(title = paste0("Leave-one-sample-out influence: ", pg, " / ", mm), x = NULL, y = "|Δ beta|")
        h <- min(16, max(5, 0.22 * nrow(loo_dt)))
        ggsave(file.path(outdir, "plots", paste0("leave1sample_", pg, "_", mm, ".pdf")), p2, width = 7.2, height = h, limitsize = FALSE)

        topinf <- inf_dt[order(-abs_dfbeta_dx)][1:min(.N, 12)]
        p3 <- ggplot(topinf, aes(x = reorder(sample_id, abs_dfbeta_dx), y = abs_dfbeta_dx, fill = dx_group)) +
          geom_col() + coord_flip() + theme_bw(base_size = 9) +
          labs(title = paste0("Top DFBETA samples: ", pg, " / ", mm), x = NULL, y = "|DFBETA(dx)|")
        h2 <- min(10, max(4.5, 0.35 * nrow(topinf)))
        ggsave(file.path(outdir, "plots", paste0("dfbeta_top_", pg, "_", mm, ".pdf")), p3, width = 6.5, height = h2, limitsize = FALSE)
      }

      # metadata univariable screen
      covars <- names(sc)[grepl("region|sex|gender|age|rin|pmi|batch|site", names(sc), ignore.case = TRUE)]
      covars <- setdiff(unique(covars), c("dx_group", "sample_id", "program", "program_score_raw", "program_score_z"))
      if (length(covars) > 0) {
        for (cv in covars) {
          x <- sc[[cv]]
          if (is.list(x) && !is.data.frame(x)) next
          dfu <- sc[, .(program_score_z, var = get(cv))]
          dfu <- dfu[!is.na(program_score_z) & !is.na(var)]
          if (nrow(dfu) < 6 || uniqueN(dfu$var) < 2) next
          if (is.numeric(dfu$var) || suppressWarnings(all(!is.na(as.numeric(as.character(dfu$var)))))) {
            dfu[, var_num := as.numeric(as.character(var))]
            fit_u <- tryCatch(lm(program_score_z ~ var_num, data = dfu), error = function(e) NULL)
            if (!is.null(fit_u)) {
              cf <- summary(fit_u)$coefficients
              rr <- data.table(program = pg, model = mm, variable = cv, var_type = "numeric",
                               beta = unname(cf["var_num", "Estimate"]),
                               p = unname(cf["var_num", grep("Pr\\(", colnames(cf))]),
                               n = nrow(dfu), n_unique = uniqueN(dfu$var_num))
              meta_assoc_rows[[paste(pg, mm, cv, sep = "::")]] <- rr
            }
          } else {
            fit_u <- tryCatch(anova(lm(program_score_z ~ as.factor(var), data = dfu)), error = function(e) NULL)
            if (!is.null(fit_u) && nrow(fit_u) >= 1) {
              rr <- data.table(program = pg, model = mm, variable = cv, var_type = "categorical",
                               beta = NA_real_,
                               p = fit_u$`Pr(>F)`[1],
                               n = nrow(dfu), n_unique = uniqueN(dfu$var))
              meta_assoc_rows[[paste(pg, mm, cv, sep = "::")]] <- rr
            }
          }
        }
      }
    }
  }

  score_dt_all <- rbindlist(score_rows, fill = TRUE)
  model_dt <- rbindlist(model_rows, fill = TRUE)
  loo_dt_all <- rbindlist(loo_rows, fill = TRUE)
  infl_dt_all <- rbindlist(infl_rows, fill = TRUE)
  meta_assoc_dt <- rbindlist(meta_assoc_rows, fill = TRUE)

  top_samples_loo <- loo_dt_all[status == "ok"][order(program, model, -abs_delta_beta)]
  top_samples_loo <- top_samples_loo[, head(.SD, 10), by = .(program, model)]

  top_samples_dfb <- infl_dt_all[order(program, model, -abs_dfbeta_dx)]
  top_samples_dfb <- top_samples_dfb[, head(.SD, 10), by = .(program, model)]

  summary_dt <- merge(
    model_dt,
    loo_dt_all[status == "ok", .(
      loo_n_ok = .N,
      loo_direction_retained_prop = mean(sign(beta) == sign(first(full_beta)), na.rm = TRUE),
      loo_max_abs_delta_beta = max(abs_delta_beta, na.rm = TRUE),
      loo_median_abs_delta_beta = median(abs_delta_beta, na.rm = TRUE),
      strongest_loo_sample = sample_id[which.max(abs_delta_beta)]
    ), by = .(program, model)],
    by = c("program", "model"),
    all.x = TRUE
  )

  safe_write(data.table(
    time = timestamp(),
    root_dir = root_dir,
    outdir = outdir,
    expr_path = expr_path,
    meta_path = meta_path,
    sfari_path = sfari_path,
    program_path = program_path,
    syn_inventory = syn_inventory,
    n_samples_expr = ncol(expr_log2cpm),
    n_samples_meta = nrow(meta),
    sample_col = sid$sample_col,
    dx_col = dxi$dx_col
  ), file.path(outdir, "00_run_info.txt"))

  safe_write(score_dt_all, file.path(outdir, "04_program_score_sample_table.tsv.gz"))
  safe_write(summary_dt, file.path(outdir, "05_model_summary.tsv"))
  safe_write(loo_dt_all, file.path(outdir, "06_leave1sampleout.tsv.gz"))
  safe_write(infl_dt_all, file.path(outdir, "07_influence_table.tsv.gz"))
  safe_write(top_samples_loo, file.path(outdir, "08_top_leave1sample_influential_samples.tsv"))
  safe_write(top_samples_dfb, file.path(outdir, "09_top_dfbeta_samples.tsv"))
  safe_write(meta_assoc_dt, file.path(outdir, "10_metadata_univariable_associations.tsv"))

  logmsg("Done. Main summary written to: ", file.path(outdir, "05_model_summary.tsv"))
}

main()
