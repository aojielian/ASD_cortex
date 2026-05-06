#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- list(
    outdir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step10C_adult_sc_donor_pseudobulk_validation_raw",
    programs_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/20_midPrenatal_core_programs.tsv",
    sfari_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/01_SFARI_primary_Sand1_standardized.tsv",
    psych_meta = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/meta.tsv",
    psych_features = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/counts_features.tsv.gz",
    psych_barcodes = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/counts_barcodes.tsv.gz",
    psych_matrix = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/counts_matrix.mtx.gz",
    vel_zip = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/rawMatrix.zip",
    min_cells_per_donor_class = 20,
    norm_scale_factor = 10000
  )
  if (length(args) %% 2 != 0) stop("Arguments must be key-value pairs.")
  if (length(args) > 0) {
    for (i in seq(1, length(args), by = 2)) {
      k <- sub("^--", "", args[i])
      v <- args[i + 1]
      opt[[k]] <- v
    }
  }
  opt$min_cells_per_donor_class <- as.integer(opt$min_cells_per_donor_class)
  opt$norm_scale_factor <- as.numeric(opt$norm_scale_factor)
  opt
}

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

norm_name <- function(x) {
  x <- tolower(x)
  gsub("[^a-z0-9]+", "_", x)
}

pick_col <- function(nms, patterns, label) {
  nn <- norm_name(nms)
  for (pat in patterns) {
    idx <- grep(pat, nn, perl = TRUE)
    if (length(idx) > 0) return(nms[idx[1]])
  }
  stop("Could not identify ", label, ". Available columns: ", paste(nms, collapse = ", "))
}

standardize_dx <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("^asd$|autism|case|patient", y)] <- "ASD"
  out[grepl("^control$|^ctl$|^ctrl$|healthy|typical|neurotyp", y)] <- "Control"
  out
}

map_psych_broad_class <- function(x) {
  y <- tolower(as.character(x))
  out <- rep(NA_character_, length(y))
  out[grepl("(^|[^a-z])(exn|excit|glut|pyram|projection|excitatory)($|[^a-z])", y)] <- "EXN"
  out[grepl("(^|[^a-z])(inn|inhib|gaba|interneuron|inhibitory)($|[^a-z])", y)] <- "INN"
  out[grepl("astro|ast", y)] <- "AST"
  out[grepl("opc", y)] <- "OPC"
  out[grepl("oligo|odc", y)] <- "ODC"
  out[grepl("microgl|mg", y)] <- "MG"
  out[grepl("endo|end|vascular|pericyte|vlmc", y)] <- "END"
  out
}

map_vel_broad_class <- function(cluster_vec) {
  y <- tolower(as.character(cluster_vec))
  out <- rep(NA_character_, length(y))

  # glia / non-neuronal
  out[grepl("astro|^ast-|ast-fb|ast-pp", y)] <- "AST"
  out[grepl("^opc$|\\bopc\\b", y)] <- "OPC"
  out[grepl("oligo", y)] <- "ODC"
  out[grepl("micro|mg\\b|macroph", y)] <- "MG"
  out[grepl("endo|pericyte|vascular|vlmc", y)] <- "END"

  # inhibitory neuronal patterns
  inn_pat <- paste(
    c("vip", "sst", "pvalb", "pv\\b", "lamp5", "reln", "sncg", "gad", "gaba",
      "inh", "interneuron", "lhx6", "l2_3_int", "l4_int", "l5_6_int", "dlx", "cck", "in-sv2c", "sv2c"),
    collapse = "|"
  )
  out[is.na(out) & grepl(inn_pat, y)] <- "INN"

  # excitatory neuronal patterns
  exn_pat <- paste(
    c("^l2", "^l3", "^l4", "^l5", "^l6", "^neu-", "excit", "glut", "cux", "satb2", "tbr1",
      "fezf2", "foxp2", "themis", "it\\b", "et\\b", "ct\\b", "cpn", "callosal"),
    collapse = "|"
  )
  out[is.na(out) & grepl(exn_pat, y)] <- "EXN"

  out
}

read_features <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t", fill = TRUE)
  if (ncol(dt) >= 2) {
    gene_symbol <- as.character(dt[[2]])
    gene_id <- as.character(dt[[1]])
  } else {
    gene_symbol <- as.character(dt[[1]])
    gene_id <- as.character(dt[[1]])
  }
  gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <- gene_id[is.na(gene_symbol) | gene_symbol == ""]
  data.table(gene_id = gene_id, gene_symbol = gene_symbol)
}

read_barcodes <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t")
  as.character(dt[[1]])
}

subset_unique_genes <- function(mat, gene_symbols) {
  keep <- !duplicated(gene_symbols) & !is.na(gene_symbols) & gene_symbols != ""
  mat2 <- mat[keep, , drop = FALSE]
  rownames(mat2) <- toupper(gene_symbols[keep])
  mat2
}

load_programs <- function(core_path, sfari_path) {
  core <- fread(core_path, sep = "\t", header = TRUE)
  core <- core[program_name %in% c("midPrenatal_SFARI_top20"), .(program_name, gene_symbol)]
  sf <- fread(sfari_path, sep = "\t", header = TRUE)
  nn <- norm_name(names(sf))
  hit <- grep("gene_symbol|symbol|gene", nn, perl = TRUE)
  if (length(hit) == 0) stop("Could not identify gene symbol column in SFARI file.")
  symbol_col <- names(sf)[hit[1]]
  sfari <- data.table(program_name = "SFARI_all", gene_symbol = as.character(sf[[symbol_col]]))
  prog <- rbindlist(list(core, sfari), use.names = TRUE)
  prog[, gene_symbol := toupper(trimws(gene_symbol))]
  prog <- unique(prog[!is.na(gene_symbol) & gene_symbol != ""])
  split(prog$gene_symbol, prog$program_name)
}

read_psych_dataset <- function(meta_path, features_path, barcodes_path, matrix_path, keep_classes = c("EXN", "INN")) {
  log_msg("Reading PsychENCODE metadata: ", meta_path)
  meta <- fread(meta_path, sep = "\t", header = TRUE, fill = TRUE, quote = "")
  cell_col <- pick_col(names(meta), c("^(barcode|barcodes|cell|cell_id|cellid)$", "(barcode|cell).*id$", "^x$"), "PsychENCODE cell column")
  donor_col <- pick_col(names(meta), c("^(donor|donor_id|subject|subject_id|individual|individual_id|sample_id)$", "(donor|subject|individual)"), "PsychENCODE donor column")
  dx_col <- pick_col(names(meta), c("^(diagnosis|dx|group|condition|phenotype|asd)$", "(diagnosis|case|control|asd)"), "PsychENCODE diagnosis column")
  annot_col <- pick_col(names(meta), c("^(annotation|cell_type|celltype|broad_class|broad_cell_type|class|cluster_annotation)$", "(annotation|cell_type|celltype|broad|cluster)"), "PsychENCODE annotation column")
  setnames(meta, c(cell_col, donor_col, dx_col, annot_col), c("cell_id", "donor_id", "diagnosis_raw", "annotation_raw"))
  meta[, diagnosis := standardize_dx(diagnosis_raw)]
  meta[, broad_class := map_psych_broad_class(annotation_raw)]
  meta <- meta[!is.na(diagnosis) & !is.na(broad_class)]
  meta <- meta[broad_class %in% keep_classes]
  log_msg("PsychENCODE retained metadata rows after filtering: ", nrow(meta))

  feats <- read_features(features_path)
  bcs <- read_barcodes(barcodes_path)
  mat <- readMM(matrix_path)
  mat <- as(mat, "CsparseMatrix")

  bc_dt <- data.table(cell_id = bcs, col_idx = seq_along(bcs))
  meta2 <- merge(meta, bc_dt, by = "cell_id", all.x = FALSE, all.y = FALSE)
  if (nrow(meta2) == 0) stop("PsychENCODE: no barcode overlap between metadata and matrix.")
  setorder(meta2, col_idx)
  mat <- mat[, meta2$col_idx, drop = FALSE]
  mat <- subset_unique_genes(mat, feats$gene_symbol)

  list(dataset = "PsychENCODE", counts = mat,
       meta = meta2[, .(cell_id, donor_id, diagnosis, annotation_raw, broad_class)])
}

read_velmeshev_dataset <- function(zip_path, keep_classes = c("EXN", "INN")) {
  if (!file.exists(zip_path)) stop("Velmeshev zip not found: ", zip_path)
  tmp_dir <- file.path(tempdir(), paste0("vel_zip_", as.integer(Sys.time())))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  log_msg("Unzipping Velmeshev raw matrix to ", tmp_dir)
  utils::unzip(zip_path, files = c("meta.txt", "barcodes.tsv", "genes.tsv", "matrix.mtx"), exdir = tmp_dir)

  meta <- fread(file.path(tmp_dir, "meta.txt"), sep = "\t", header = TRUE, fill = TRUE, quote = "")
  setnames(meta, c("cell", "individual", "diagnosis", "cluster"), c("cell_id", "donor_id", "diagnosis_raw", "annotation_raw"))
  meta[, diagnosis := standardize_dx(diagnosis_raw)]
  meta[, broad_class := map_vel_broad_class(annotation_raw)]

  # save mapping summary
  map_sum <- meta[, .N, by = .(annotation_raw, broad_class)][order(annotation_raw)]
  fwrite(map_sum, file.path(dirname(zip_path), "Velmeshev_cluster_to_broadClass_mapping_check.tsv"), sep = "\t")

  meta <- meta[!is.na(diagnosis) & !is.na(broad_class)]
  meta <- meta[broad_class %in% keep_classes]
  log_msg("Velmeshev retained metadata rows after filtering: ", nrow(meta))
  if (nrow(meta) == 0) stop("Velmeshev: no cells retained after broad-class mapping. Check mapping summary file.")

  bcs <- read_barcodes(file.path(tmp_dir, "barcodes.tsv"))
  feats <- read_features(file.path(tmp_dir, "genes.tsv"))
  mat <- readMM(file.path(tmp_dir, "matrix.mtx"))
  mat <- as(mat, "CsparseMatrix")

  bc_dt <- data.table(cell_id = bcs, col_idx = seq_along(bcs))
  meta2 <- merge(meta, bc_dt, by = "cell_id", all.x = FALSE, all.y = FALSE)
  if (nrow(meta2) == 0) stop("Velmeshev: no barcode overlap between metadata and matrix after mapping.")
  setorder(meta2, col_idx)
  mat <- mat[, meta2$col_idx, drop = FALSE]
  mat <- subset_unique_genes(mat, feats$gene_symbol)

  list(dataset = "Velmeshev", counts = mat,
       meta = meta2[, .(cell_id, donor_id, diagnosis, annotation_raw, broad_class)])
}

cell_program_scores <- function(counts, programs, scale_factor = 10000) {
  lib <- Matrix::colSums(counts)
  lib[lib == 0] <- 1
  norm <- t(t(counts) / lib) * scale_factor
  norm@x <- log1p(norm@x)
  genes <- toupper(rownames(norm))
  score_list <- list()
  map_dt <- data.table()
  for (nm in names(programs)) {
    idx <- which(genes %in% programs[[nm]])
    score_list[[nm]] <- if (length(idx) == 0) rep(NA_real_, ncol(norm)) else Matrix::colMeans(norm[idx, , drop = FALSE])
    map_dt <- rbind(map_dt, data.table(program_name = nm, genes_requested = length(programs[[nm]]), genes_present = length(idx)))
  }
  list(scores = as.data.table(score_list), map = map_dt)
}

aggregate_sparse_by_group <- function(counts, groups) {
  f <- factor(groups)
  mm <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
  agg <- counts %*% mm
  colnames(agg) <- levels(f)
  agg
}

pseudobulk_program_scores <- function(counts, meta, programs, min_cells = 20, scale_factor = 10000) {
  group_dt <- copy(meta)
  group_dt[, group_id := paste(dataset, broad_class, donor_id, diagnosis, sep = "||")]
  cell_n <- group_dt[, .N, by = group_id]
  keep_groups <- cell_n[N >= min_cells, group_id]
  keep_idx <- which(group_dt$group_id %in% keep_groups)
  if (length(keep_idx) == 0) return(NULL)

  agg_counts <- aggregate_sparse_by_group(counts[, keep_idx, drop = FALSE], group_dt$group_id[keep_idx])
  pb_meta <- unique(group_dt[group_id %in% colnames(agg_counts), .(dataset, broad_class, donor_id, diagnosis, group_id)])
  pb_meta[, order_idx := match(group_id, colnames(agg_counts))]
  pb_meta <- pb_meta[!is.na(order_idx)]
  setorder(pb_meta, order_idx)
  agg_counts <- agg_counts[, pb_meta$group_id, drop = FALSE]
  pb_meta[, order_idx := NULL]

  lib <- Matrix::colSums(agg_counts)
  lib[lib == 0] <- 1
  norm <- t(t(agg_counts) / lib) * scale_factor
  norm@x <- log1p(norm@x)

  genes <- toupper(rownames(norm))
  out <- copy(pb_meta)
  for (nm in names(programs)) {
    idx <- which(genes %in% programs[[nm]])
    out[[nm]] <- if (length(idx) == 0) NA_real_ else Matrix::colMeans(norm[idx, , drop = FALSE])
  }
  out[, n_cells := cell_n[match(group_id, cell_n$group_id), N]]
  out
}

donor_mean_scores <- function(score_dt, meta, dataset_name, min_cells = 20) {
  dt <- cbind(copy(meta), score_dt)
  dt[, dataset := dataset_name]
  long <- melt(dt,
               id.vars = c("dataset", "cell_id", "donor_id", "diagnosis", "annotation_raw", "broad_class"),
               variable.name = "program_name",
               value.name = "cell_score")
  donor <- long[, .(
    n_cells = .N,
    mean_score = mean(cell_score, na.rm = TRUE),
    median_score = median(cell_score, na.rm = TRUE)
  ), by = .(dataset, donor_id, diagnosis, broad_class, program_name)]
  donor[n_cells >= min_cells]
}

fit_models <- function(dt, value_col) {
  out <- list()
  for (ds in unique(dt$dataset)) {
    for (bc in unique(dt$broad_class)) {
      for (pg in unique(dt$program_name)) {
        sub <- dt[dataset == ds & broad_class == bc & program_name == pg & !is.na(get(value_col))]
        if (nrow(sub) < 6 || length(unique(sub$diagnosis)) < 2) next
        sub[, diagnosis := factor(diagnosis, levels = c("Control", "ASD"))]
        w_p <- tryCatch(wilcox.test(sub[diagnosis == "ASD", get(value_col)],
                                    sub[diagnosis == "Control", get(value_col)])$p.value,
                        error = function(e) NA_real_)
        delta <- mean(sub[diagnosis == "ASD", get(value_col)], na.rm = TRUE) -
          mean(sub[diagnosis == "Control", get(value_col)], na.rm = TRUE)
        fit <- tryCatch(lm(as.formula(paste(value_col, "~ diagnosis")), data = sub), error = function(e) NULL)
        beta <- p_lm <- ci_lo <- ci_hi <- NA_real_
        if (!is.null(fit)) {
          sm <- summary(fit)$coefficients
          rn <- rownames(sm)
          idx <- grep("^diagnosisASD$", rn)
          if (length(idx) == 1) {
            beta <- sm[idx, "Estimate"]
            p_lm <- sm[idx, "Pr(>|t|)"]
            ci <- tryCatch(confint(fit)[idx, ], error = function(e) c(NA_real_, NA_real_))
            ci_lo <- ci[1]
            ci_hi <- ci[2]
          }
        }
        out[[length(out) + 1]] <- data.table(
          dataset = ds,
          broad_class = bc,
          program_name = pg,
          metric = value_col,
          n_donors = nrow(sub),
          n_asd = sum(sub$diagnosis == "ASD"),
          n_control = sum(sub$diagnosis == "Control"),
          mean_control = mean(sub[diagnosis == "Control", get(value_col)], na.rm = TRUE),
          mean_asd = mean(sub[diagnosis == "ASD", get(value_col)], na.rm = TRUE),
          delta_asd_minus_control = delta,
          wilcox_p = w_p,
          lm_beta_asd = beta,
          lm_p = p_lm,
          ci_lo = ci_lo,
          ci_hi = ci_hi
        )
      }
    }
  }
  if (length(out) == 0) return(data.table())
  rbindlist(out, fill = TRUE)
}

stouffer_meta <- function(stat_dt) {
  if (nrow(stat_dt) == 0) return(data.table())
  dt <- copy(stat_dt)
  dt[, aligned_z := qnorm(pmax(pmin(1 - lm_p/2, 1 - 1e-15), 1e-15)) * sign(lm_beta_asd)]
  dt[, non_missing := !is.na(aligned_z)]
  dt[non_missing == FALSE, aligned_z := 0]
  dt[, .(
    n_datasets = sum(non_missing),
    mean_beta = mean(lm_beta_asd, na.rm = TRUE),
    direction_matched = sum(sign(lm_beta_asd) == sign(mean(lm_beta_asd, na.rm = TRUE)), na.rm = TRUE),
    stouffer_z = sum(aligned_z, na.rm = TRUE) / sqrt(sum(non_missing)),
    stouffer_p_two_sided = 2 * pnorm(-abs(sum(aligned_z, na.rm = TRUE) / sqrt(sum(non_missing))))
  ), by = .(broad_class, program_name, metric)]
}

plot_donor_box <- function(dt, value_col, outfile) {
  p <- ggplot(dt, aes(x = diagnosis, y = .data[[value_col]], fill = diagnosis)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.8) +
    geom_jitter(width = 0.12, size = 1.3, alpha = 0.75) +
    facet_grid(broad_class ~ program_name + dataset, scales = "free_y") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = "grey80")) +
    labs(x = NULL, y = value_col, title = paste("Donor-level", value_col))
  ggsave(outfile, p, width = 14, height = 7)
}

plot_forest <- function(dt, outfile, title_text) {
  if (nrow(dt) == 0) return(invisible(NULL))
  d2 <- copy(dt)
  d2[, label := paste(dataset, broad_class, program_name, sep = " | ")]
  p <- ggplot(d2, aes(x = lm_beta_asd, y = reorder(label, lm_beta_asd), xmin = ci_lo, xmax = ci_hi, color = dataset)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_errorbar(aes(y = reorder(label, lm_beta_asd), xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.18) +
    geom_point(size = 2.2) +
    theme_bw(base_size = 11) +
    labs(x = "ASD effect (linear model beta)", y = NULL, title = title_text)
  ggsave(outfile, p, width = 12, height = max(4, 0.25 * nrow(d2)))
}

main <- function() {
  opt <- parse_args()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(opt$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

  programs <- load_programs(opt$programs_file, opt$sfari_file)
  psych <- read_psych_dataset(opt$psych_meta, opt$psych_features, opt$psych_barcodes, opt$psych_matrix)
  vel <- read_velmeshev_dataset(opt$vel_zip)

  datasets <- list(psych, vel)
  donor_score_all <- list()
  pb_all <- list()
  map_all <- list()

  for (obj in datasets) {
    log_msg("Scoring cell-level programs for ", obj$dataset)
    cs <- cell_program_scores(obj$counts, programs, scale_factor = opt$norm_scale_factor)
    map_dt <- copy(cs$map)
    map_dt[, dataset := obj$dataset]
    map_all[[obj$dataset]] <- map_dt

    donor_dt <- donor_mean_scores(cs$scores, obj$meta, obj$dataset, min_cells = opt$min_cells_per_donor_class)
    donor_score_all[[obj$dataset]] <- donor_dt

    log_msg("Aggregating pseudobulk for ", obj$dataset)
    meta_pb <- copy(obj$meta)
    meta_pb[, dataset := obj$dataset]
    pb_dt <- pseudobulk_program_scores(obj$counts, meta_pb, programs,
                                       min_cells = opt$min_cells_per_donor_class,
                                       scale_factor = opt$norm_scale_factor)
    if (!is.null(pb_dt)) {
      pb_long <- melt(pb_dt,
                      id.vars = c("dataset", "broad_class", "donor_id", "diagnosis", "group_id", "n_cells"),
                      measure.vars = names(programs),
                      variable.name = "program_name",
                      value.name = "pseudobulk_score")
      pb_all[[obj$dataset]] <- pb_long
    }
  }

  donor_scores <- rbindlist(donor_score_all, fill = TRUE)
  pb_scores <- rbindlist(pb_all, fill = TRUE)
  map_dt <- rbindlist(map_all, fill = TRUE)

  fwrite(map_dt, file.path(opt$outdir, "tables", "01_program_gene_overlap.tsv"), sep = "\t")
  fwrite(donor_scores, file.path(opt$outdir, "tables", "02_donor_mean_cell_scores.tsv.gz"), sep = "\t")
  fwrite(pb_scores, file.path(opt$outdir, "tables", "03_donor_pseudobulk_scores.tsv.gz"), sep = "\t")

  donor_stats <- fit_models(donor_scores, "mean_score")
  pb_stats <- fit_models(pb_scores, "pseudobulk_score")
  fwrite(donor_stats, file.path(opt$outdir, "tables", "04_donor_mean_score_case_control_stats.tsv"), sep = "\t")
  fwrite(pb_stats, file.path(opt$outdir, "tables", "05_donor_pseudobulk_case_control_stats.tsv"), sep = "\t")

  donor_meta <- stouffer_meta(donor_stats)
  pb_meta <- stouffer_meta(pb_stats)
  fwrite(donor_meta, file.path(opt$outdir, "tables", "06_donor_mean_score_meta_summary.tsv"), sep = "\t")
  fwrite(pb_meta, file.path(opt$outdir, "tables", "07_donor_pseudobulk_meta_summary.tsv"), sep = "\t")

  plot_donor_box(donor_scores, "mean_score", file.path(opt$outdir, "plots", "01_donor_mean_score_boxplots.png"))
  plot_donor_box(pb_scores, "pseudobulk_score", file.path(opt$outdir, "plots", "02_donor_pseudobulk_score_boxplots.png"))
  plot_forest(donor_stats, file.path(opt$outdir, "plots", "03_donor_mean_score_forest.png"), "Donor-level mean cell-score effects")
  plot_forest(pb_stats, file.path(opt$outdir, "plots", "04_donor_pseudobulk_forest.png"), "Donor-level pseudobulk score effects")

  run_summary <- data.table(
    item = c("datasets", "programs", "min_cells_per_donor_class", "n_donor_score_rows", "n_pseudobulk_rows"),
    value = c(paste(vapply(datasets, `[[`, "", "dataset"), collapse = ", "),
              paste(names(programs), collapse = ", "),
              as.character(opt$min_cells_per_donor_class),
              as.character(nrow(donor_scores)),
              as.character(nrow(pb_scores)))
  )
  fwrite(run_summary, file.path(opt$outdir, "tables", "08_run_summary.tsv"), sep = "\t")
  log_msg("DONE. Outputs written to: ", opt$outdir)
}

main()
