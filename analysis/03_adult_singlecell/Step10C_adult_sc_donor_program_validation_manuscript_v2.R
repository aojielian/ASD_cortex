#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
PROJECT_ROOT <- if (!is.null(args$project_root)) args$project_root else "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- if (!is.null(args$outdir)) args$outdir else file.path(PROJECT_ROOT, "step10C_adult_sc_donor_program_validation_manuscript")
PROGRAMS_FILE <- if (!is.null(args$programs_file)) args$programs_file else file.path(PROJECT_ROOT, "step01_prepare", "tables", "20_midPrenatal_core_programs.tsv")
SFARI_FILE <- if (!is.null(args$sfari_file)) args$sfari_file else file.path(PROJECT_ROOT, "step01_prepare", "tables", "01_SFARI_primary_Sand1_standardized.tsv")
PSYCH_META <- if (!is.null(args$psych_meta)) args$psych_meta else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/meta.tsv"
PSYCH_FEATURES <- if (!is.null(args$psych_features)) args$psych_features else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/counts_features.tsv.gz"
PSYCH_BARCODES <- if (!is.null(args$psych_barcodes)) args$psych_barcodes else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/counts_barcodes.tsv.gz"
PSYCH_MATRIX <- if (!is.null(args$psych_matrix)) args$psych_matrix else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/counts_matrix.mtx.gz"
VEL_ZIP <- if (!is.null(args$vel_zip)) args$vel_zip else "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/rawMatrix.zip"
MIN_CELLS <- if (!is.null(args$min_cells_per_donor_class)) as.integer(args$min_cells_per_donor_class) else 20L
SCALE_FACTOR <- if (!is.null(args$norm_scale_factor)) as.numeric(args$norm_scale_factor) else 10000
KEEP_CLASSES <- if (!is.null(args$keep_classes)) trimws(strsplit(args$keep_classes, ",", fixed = TRUE)[[1]]) else c("EXN","INN","AST","OPC","ODC","MG","END")
PRIMARY_CLASSES <- if (!is.null(args$primary_classes)) trimws(strsplit(args$primary_classes, ",", fixed = TRUE)[[1]]) else c("EXN","INN")

TABDIR <- file.path(OUTDIR, "tables")
PLOTDIR <- file.path(OUTDIR, "plots")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOGDIR, "Step10C_adult_sc_donor_program_validation_manuscript_v2.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOGFILE, append = TRUE)
}

safe_fwrite <- function(x, file) {
  if (grepl("\\.gz$", file, ignore.case = TRUE)) {
    fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
  } else {
    fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
  }
}

norm_name <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))

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
  out[grepl("endo|vascular|pericyte|vlmc", y)] <- "END"
  out
}

map_vel_broad_class <- function(cluster_vec) {
  y <- tolower(as.character(cluster_vec))
  out <- rep(NA_character_, length(y))
  out[grepl("astro|^ast-|ast-fb|ast-pp", y)] <- "AST"
  out[grepl("^opc$|\\bopc\\b", y)] <- "OPC"
  out[grepl("oligo", y)] <- "ODC"
  out[grepl("micro|mg\\b|macroph", y)] <- "MG"
  out[grepl("endo|pericyte|vascular|vlmc", y)] <- "END"
  inn_pat <- paste(c("vip", "sst", "pvalb", "pv\\b", "lamp5", "reln", "sncg", "gad", "gaba",
                     "inh", "interneuron", "lhx6", "l2_3_int", "l4_int", "l5_6_int", "dlx", "cck", "in-sv2c", "sv2c"), collapse = "|")
  out[is.na(out) & grepl(inn_pat, y)] <- "INN"
  exn_pat <- paste(c("^l2", "^l3", "^l4", "^l5", "^l6", "^neu-", "excit", "glut", "cux", "satb2", "tbr1",
                     "fezf2", "foxp2", "themis", "it\\b", "et\\b", "ct\\b", "cpn", "callosal"), collapse = "|")
  out[is.na(out) & grepl(exn_pat, y)] <- "EXN"
  out
}

read_features <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t", fill = TRUE)
  if (ncol(dt) >= 2L) {
    gene_symbol <- as.character(dt[[2L]])
    gene_id <- as.character(dt[[1L]])
  } else {
    gene_symbol <- as.character(dt[[1L]])
    gene_id <- as.character(dt[[1L]])
  }
  gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <- gene_id[is.na(gene_symbol) | gene_symbol == ""]
  data.table(gene_id = gene_id, gene_symbol = gene_symbol)
}

read_barcodes <- function(path) {
  dt <- fread(path, header = FALSE, sep = "\t")
  as.character(dt[[1L]])
}

subset_unique_genes <- function(mat, gene_symbols) {
  keep <- !duplicated(gene_symbols) & !is.na(gene_symbols) & gene_symbols != ""
  mat2 <- mat[keep, , drop = FALSE]
  rownames(mat2) <- toupper(gene_symbols[keep])
  mat2
}

load_programs <- function(core_path, sfari_path) {
  core <- fread(core_path)
  core <- core[program_name %in% c("midPrenatal_SFARI_top20"), .(program_name, gene_symbol)]
  sf <- fread(sfari_path)
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

read_psych_dataset <- function() {
  log_msg("Reading PsychENCODE metadata: ", PSYCH_META)
  meta <- fread(PSYCH_META, sep = "\t", header = TRUE, fill = TRUE, quote = "")
  cell_col <- pick_col(names(meta), c("^(barcode|barcodes|cell|cell_id|cellid)$", "(barcode|cell).*id$", "^x$"), "PsychENCODE cell column")
  donor_col <- pick_col(names(meta), c("^(donor|donor_id|subject|subject_id|individual|individual_id|sample_id)$", "(donor|subject|individual)"), "PsychENCODE donor column")
  dx_col <- pick_col(names(meta), c("^(diagnosis|dx|group|condition|phenotype|asd)$", "(diagnosis|case|control|asd)"), "PsychENCODE diagnosis column")
  annot_col <- pick_col(names(meta), c("^(annotation|cell_type|celltype|broad_class|broad_cell_type|class|cluster_annotation)$", "(annotation|cell_type|celltype|broad|cluster)"), "PsychENCODE annotation column")
  setnames(meta, c(cell_col, donor_col, dx_col, annot_col), c("cell_id", "donor_id", "diagnosis_raw", "annotation_raw"))
  meta[, diagnosis := standardize_dx(diagnosis_raw)]
  meta[, broad_class := map_psych_broad_class(annotation_raw)]
  meta <- meta[!is.na(diagnosis) & !is.na(broad_class)]
  meta <- meta[broad_class %in% KEEP_CLASSES]
  log_msg("PsychENCODE retained metadata rows after class filtering: ", nrow(meta))

  feats <- read_features(PSYCH_FEATURES)
  bcs <- read_barcodes(PSYCH_BARCODES)
  mat <- readMM(PSYCH_MATRIX)
  mat <- as(mat, "CsparseMatrix")

  bc_dt <- data.table(cell_id = bcs, col_idx = seq_along(bcs))
  meta2 <- merge(meta, bc_dt, by = "cell_id", all.x = FALSE, all.y = FALSE)
  if (nrow(meta2) == 0) stop("PsychENCODE: no barcode overlap between metadata and matrix.")
  setorder(meta2, col_idx)
  mat <- mat[, meta2$col_idx, drop = FALSE]
  mat <- subset_unique_genes(mat, feats$gene_symbol)

  list(dataset = "PsychENCODE",
       counts = mat,
       meta = meta2[, .(cell_id, donor_id, diagnosis, annotation_raw, broad_class)])
}

read_velmeshev_dataset <- function() {
  if (!file.exists(VEL_ZIP)) stop("Velmeshev zip not found: ", VEL_ZIP)
  tmp_dir <- file.path(tempdir(), paste0("vel_zip_", as.integer(Sys.time())))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  log_msg("Unzipping Velmeshev raw matrix to ", tmp_dir)
  utils::unzip(VEL_ZIP, files = c("meta.txt", "barcodes.tsv", "genes.tsv", "matrix.mtx"), exdir = tmp_dir)

  meta <- fread(file.path(tmp_dir, "meta.txt"), sep = "\t", header = TRUE, fill = TRUE, quote = "")
  setnames(meta, c("cell", "individual", "diagnosis", "cluster"), c("cell_id", "donor_id", "diagnosis_raw", "annotation_raw"))
  meta[, diagnosis := standardize_dx(diagnosis_raw)]
  meta[, broad_class := map_vel_broad_class(annotation_raw)]
  map_sum <- meta[, .N, by = .(annotation_raw, broad_class)][order(annotation_raw)]
  safe_fwrite(map_sum, file.path(TABDIR, "00_velmeshev_cluster_to_broadClass_mapping.tsv"))
  meta <- meta[!is.na(diagnosis) & !is.na(broad_class)]
  meta <- meta[broad_class %in% KEEP_CLASSES]
  log_msg("Velmeshev retained metadata rows after class filtering: ", nrow(meta))
  if (nrow(meta) == 0) stop("Velmeshev: no cells retained after broad-class mapping.")

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

  list(dataset = "Velmeshev",
       counts = mat,
       meta = meta2[, .(cell_id, donor_id, diagnosis, annotation_raw, broad_class)])
}

cell_program_scores <- function(counts, programs) {
  lib <- Matrix::colSums(counts)
  lib[lib == 0] <- 1
  norm <- t(t(counts) / lib) * SCALE_FACTOR
  norm@x <- log1p(norm@x)
  genes <- toupper(rownames(norm))
  score_list <- list()
  map_dt <- list()
  for (nm in names(programs)) {
    idx <- which(genes %in% programs[[nm]])
    score_list[[nm]] <- if (length(idx) == 0) rep(NA_real_, ncol(norm)) else Matrix::colMeans(norm[idx, , drop = FALSE])
    map_dt[[nm]] <- data.table(program_name = nm, genes_requested = length(programs[[nm]]), genes_present = length(idx))
  }
  list(scores = as.data.table(score_list), map = rbindlist(map_dt, use.names = TRUE, fill = TRUE))
}

aggregate_sparse_by_group <- function(counts, groups) {
  f <- factor(groups)
  mm <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)))
  agg <- counts %*% mm
  colnames(agg) <- levels(f)
  agg
}

pseudobulk_program_scores <- function(counts, meta, programs) {
  group_dt <- copy(meta)
  group_dt[, group_id := paste(dataset, broad_class, donor_id, diagnosis, sep = "||")]

  group_summary <- group_dt[, .(
    dataset = first(dataset),
    broad_class = first(broad_class),
    donor_id = first(donor_id),
    diagnosis = first(diagnosis),
    n_cells = .N
  ), by = group_id]

  group_summary <- group_summary[n_cells >= MIN_CELLS]
  if (nrow(group_summary) == 0) return(NULL)

  keep_idx <- which(group_dt$group_id %in% group_summary$group_id)
  if (length(keep_idx) == 0) return(NULL)

  agg_counts <- aggregate_sparse_by_group(counts[, keep_idx, drop = FALSE], group_dt$group_id[keep_idx])
  if (ncol(agg_counts) == 0) return(NULL)

  keep_groups <- intersect(colnames(agg_counts), group_summary$group_id)
  if (length(keep_groups) == 0) return(NULL)

  agg_counts <- agg_counts[, keep_groups, drop = FALSE]
  group_summary <- group_summary[match(keep_groups, group_id)]

  stopifnot(ncol(agg_counts) == nrow(group_summary))
  stopifnot(all(colnames(agg_counts) == group_summary$group_id))

  lib <- Matrix::colSums(agg_counts)
  lib[lib == 0] <- 1
  norm <- t(t(agg_counts) / lib) * SCALE_FACTOR
  norm@x <- log1p(norm@x)
  genes <- toupper(rownames(norm))
  out <- copy(group_summary)
  for (nm in names(programs)) {
    idx <- which(genes %in% programs[[nm]])
    out[[nm]] <- if (length(idx) == 0) NA_real_ else Matrix::colMeans(norm[idx, , drop = FALSE])
  }
  setcolorder(out, c("dataset", "broad_class", "donor_id", "diagnosis", "group_id", "n_cells", names(programs)))
  out[]
}

donor_mean_scores <- function(score_dt, meta, dataset_name) {
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
  donor[n_cells >= MIN_CELLS]
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
                                    sub[diagnosis == "Control", get(value_col)],
                                    exact = FALSE)$p.value,
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
        out[[length(out) + 1L]] <- data.table(
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
  ans <- rbindlist(out, fill = TRUE)
  ans[, fdr_lm := p.adjust(lm_p, method = "BH")]
  ans[, fdr_wilcox := p.adjust(wilcox_p, method = "BH")]
  ans
}

stouffer_meta <- function(stat_dt) {
  if (nrow(stat_dt) == 0) return(data.table())
  dt <- copy(stat_dt)
  dt[, aligned_z := qnorm(pmax(pmin(1 - lm_p / 2, 1 - 1e-15), 1e-15)) * sign(lm_beta_asd)]
  dt[, valid := is.finite(aligned_z)]
  dt[valid == FALSE, aligned_z := 0]
  out <- dt[, .(
    n_datasets = sum(valid),
    mean_beta = mean(lm_beta_asd, na.rm = TRUE),
    n_direction_matched = sum(sign(lm_beta_asd) == sign(mean(lm_beta_asd, na.rm = TRUE)), na.rm = TRUE),
    stouffer_z = ifelse(sum(valid) > 0, sum(aligned_z, na.rm = TRUE) / sqrt(sum(valid)), NA_real_),
    stouffer_p_two_sided = ifelse(sum(valid) > 0, 2 * pnorm(-abs(sum(aligned_z, na.rm = TRUE) / sqrt(sum(valid)))), NA_real_)
  ), by = .(broad_class, program_name, metric)]
  out[, fdr_stouffer := p.adjust(stouffer_p_two_sided, method = "BH")]
  out
}

make_primary_filter <- function(dt) {
  dt[broad_class %in% PRIMARY_CLASSES & program_name %in% c("SFARI_all", "midPrenatal_SFARI_top20")]
}

plot_box <- function(dt, value_col, outfile, title_text) {
  if (nrow(dt) == 0) return(invisible(NULL))
  d <- copy(dt)
  d[, score_value := get(value_col)]
  p <- ggplot(d, aes(x = diagnosis, y = score_value, fill = diagnosis)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.82) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.75) +
    facet_grid(broad_class ~ program_name + dataset, scales = "free_y") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = "grey80")) +
    labs(x = NULL, y = value_col, title = title_text)
  ggsave(outfile, p, width = 12, height = 6.5)
}

plot_forest <- function(dt, outfile, title_text) {
  if (nrow(dt) == 0) return(invisible(NULL))
  d <- copy(dt)
  d[, label := paste(dataset, broad_class, program_name, sep = " | ")]
  p <- ggplot(d, aes(x = lm_beta_asd, y = reorder(label, lm_beta_asd), xmin = ci_lo, xmax = ci_hi, color = dataset)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_errorbarh(height = 0.18, size = 0.5) +
    geom_point(size = 2.2) +
    theme_bw(base_size = 11) +
    labs(x = "ASD effect (linear-model beta)", y = NULL, title = title_text)
  ggsave(outfile, p, width = 11, height = max(4, 0.24 * nrow(d)))
}

log_msg("Loading principal ASD programs.")
programs <- load_programs(PROGRAMS_FILE, SFARI_FILE)
log_msg("Programs loaded: ", paste(names(programs), collapse = ", "))

psych <- read_psych_dataset()
vel <- read_velmeshev_dataset()
datasets <- list(psych, vel)

donor_score_all <- list()
pb_all <- list()
map_all <- list()
for (obj in datasets) {
  log_msg("Scoring cell-level expression-based program values for ", obj$dataset)
  cs <- cell_program_scores(obj$counts, programs)
  map_dt <- copy(cs$map)
  map_dt[, dataset := obj$dataset]
  map_all[[obj$dataset]] <- map_dt

  donor_dt <- donor_mean_scores(cs$scores, obj$meta, obj$dataset)
  donor_score_all[[obj$dataset]] <- donor_dt

  log_msg("Aggregating pseudobulk expression-based scores for ", obj$dataset)
  meta_pb <- copy(obj$meta)
  meta_pb[, dataset := obj$dataset]
  pb_dt <- pseudobulk_program_scores(obj$counts, meta_pb, programs)
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

safe_fwrite(map_dt, file.path(TABDIR, "01_program_gene_overlap.tsv"))
safe_fwrite(donor_scores, file.path(TABDIR, "02_donor_mean_expression_scores.tsv.gz"))
safe_fwrite(pb_scores, file.path(TABDIR, "03_donor_pseudobulk_expression_scores.tsv.gz"))

donor_stats <- fit_models(donor_scores, "mean_score")
pb_stats <- fit_models(pb_scores, "pseudobulk_score")
safe_fwrite(donor_stats, file.path(TABDIR, "04_donor_mean_score_dx_tests.tsv"))
safe_fwrite(pb_stats, file.path(TABDIR, "05_donor_pseudobulk_dx_tests.tsv"))

donor_meta <- stouffer_meta(donor_stats)
pb_meta <- stouffer_meta(pb_stats)
safe_fwrite(donor_meta, file.path(TABDIR, "06_donor_mean_score_meta_summary.tsv"))
safe_fwrite(pb_meta, file.path(TABDIR, "07_donor_pseudobulk_meta_summary.tsv"))

primary_donor_scores <- make_primary_filter(donor_scores)
primary_pb_scores <- make_primary_filter(pb_scores)
primary_donor_stats <- make_primary_filter(donor_stats)
primary_pb_stats <- make_primary_filter(pb_stats)
primary_donor_meta <- make_primary_filter(donor_meta)
primary_pb_meta <- make_primary_filter(pb_meta)

safe_fwrite(primary_donor_stats, file.path(TABDIR, "08_primary_panel_donor_mean_score_dx_tests.tsv"))
safe_fwrite(primary_pb_stats, file.path(TABDIR, "09_primary_panel_donor_pseudobulk_dx_tests.tsv"))
safe_fwrite(primary_donor_meta, file.path(TABDIR, "10_primary_panel_donor_mean_score_meta_summary.tsv"))
safe_fwrite(primary_pb_meta, file.path(TABDIR, "11_primary_panel_donor_pseudobulk_meta_summary.tsv"))

plot_box(primary_donor_scores, "mean_score", file.path(PLOTDIR, "01_primary_donor_mean_score_boxplots.png"),
         "Primary panel: donor-level mean expression scores")
plot_box(primary_pb_scores, "pseudobulk_score", file.path(PLOTDIR, "02_primary_donor_pseudobulk_boxplots.png"),
         "Primary panel: donor-level pseudobulk expression scores")
plot_forest(primary_donor_stats, file.path(PLOTDIR, "03_primary_donor_mean_score_forest.png"),
            "Primary panel: donor-level mean-score ASD effects")
plot_forest(primary_pb_stats, file.path(PLOTDIR, "04_primary_donor_pseudobulk_forest.png"),
            "Primary panel: donor-level pseudobulk ASD effects")

run_summary <- rbindlist(list(
  data.table(item = "datasets", value = paste(vapply(datasets, `[[`, "", "dataset"), collapse = ", ")),
  data.table(item = "programs", value = paste(names(programs), collapse = ", ")),
  data.table(item = "keep_classes", value = paste(KEEP_CLASSES, collapse = ", ")),
  data.table(item = "primary_classes", value = paste(PRIMARY_CLASSES, collapse = ", ")),
  data.table(item = "min_cells_per_donor_class", value = as.character(MIN_CELLS)),
  data.table(item = "n_donor_score_rows", value = as.character(nrow(donor_scores))),
  data.table(item = "n_pseudobulk_rows", value = as.character(nrow(pb_scores))),
  data.table(item = "best_primary_donor_mean_hit", value = ifelse(nrow(primary_donor_stats) > 0, primary_donor_stats[which.min(fdr_lm), paste(dataset, broad_class, program_name, sep = " | ")], NA_character_)),
  data.table(item = "best_primary_pseudobulk_hit", value = ifelse(nrow(primary_pb_stats) > 0, primary_pb_stats[which.min(fdr_lm), paste(dataset, broad_class, program_name, sep = " | ")], NA_character_))
), use.names = TRUE, fill = TRUE)
safe_fwrite(run_summary, file.path(TABDIR, "12_run_summary.tsv"))

log_msg("Step10C manuscript version completed successfully. Outputs written to: ", OUTDIR)
