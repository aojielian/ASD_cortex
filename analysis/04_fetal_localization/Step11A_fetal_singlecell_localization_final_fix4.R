#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

log_msg <- function(..., .file = NULL) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  if (!is.null(.file)) cat(msg, "\n", file = .file, append = TRUE)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_fwrite <- function(x, file) {
  if (grepl("\\.gz$", file, ignore.case = TRUE)) {
    fwrite(x, file = file, sep = "\t", quote = FALSE, compress = "gzip")
  } else {
    fwrite(x, file = file, sep = "\t", quote = FALSE)
  }
}

standardize_string <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null")] <- NA_character_
  x
}

standardize_gene_symbols <- function(x) {
  x <- as.character(x)
  x <- gsub("\\..*$", "", x)
  x <- trimws(x)
  toupper(x)
}

pick_col <- function(df, candidates, required = FALSE) {
  nm <- names(df)
  idx <- which(tolower(nm) %in% tolower(candidates))
  if (length(idx) > 0) return(nm[idx[1]])
  if (required) stop("Could not find required column. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

load_program_gene_sets <- function(program_table_path, requested_programs, min_genes = 5) {
  dt <- as.data.table(fread(program_table_path))
  gene_col <- pick_col(dt, c("gene_symbol", "symbol", "gene", "gene_name", "hgnc_symbol", "Gene", "SYMBOL"), required = FALSE)
  program_col <- pick_col(dt, c("program_name", "program", "gene_set", "set_name", "signature", "module"), required = FALSE)

  if (!is.na(gene_col) && !is.na(program_col)) {
    dt <- dt[!is.na(get(gene_col)) & !is.na(get(program_col))]
    dt[, gene_std := standardize_gene_symbols(get(gene_col))]
    dt[, program_std := as.character(get(program_col))]
    gs <- split(dt$gene_std, dt$program_std)
  } else {
    candidate_prog_cols <- intersect(requested_programs, names(dt))
    if (length(candidate_prog_cols) == 0) {
      stop("Could not infer gene/program columns from program table: ", program_table_path)
    }
    gs <- lapply(candidate_prog_cols, function(cc) unique(standardize_gene_symbols(dt[[cc]])))
    names(gs) <- candidate_prog_cols
  }

  gs <- gs[names(gs) %in% requested_programs]
  gs <- lapply(gs, function(v) unique(v[!is.na(v) & v != ""]))
  gs <- gs[lengths(gs) >= min_genes]
  if (length(gs) == 0) {
    stop("No usable requested programs found in ", program_table_path,
         ". Requested: ", paste(requested_programs, collapse = ", "))
  }
  gs
}

load_fetal_reference_rds <- function(rds_path, log_file) {
  log_msg("Loading fetal reference RDS: ", rds_path, .file = log_file)
  obj <- readRDS(rds_path)
  if (!is.list(obj)) stop("fetal_singlecell_reference.rds must be a list-like object.")
  if (is.null(obj$counts) || is.null(obj$meta)) {
    stop("RDS must contain at least $counts and $meta.")
  }

  counts <- obj$counts
  meta <- as.data.table(obj$meta)
  umap <- if (!is.null(obj$umap)) as.data.table(obj$umap) else NULL
  dataset_id <- if (!is.null(obj$dataset_id)) as.character(obj$dataset_id) else basename(rds_path)

  if (!inherits(counts, "Matrix") && !is.matrix(counts)) {
    stop("$counts must be a matrix or Matrix object.")
  }
  if (is.matrix(counts)) counts <- Matrix(counts, sparse = TRUE)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("$counts must have rownames (genes) and colnames (cells).")
  }

  if (!"cell_id" %in% names(meta)) {
    if (!is.null(rownames(meta))) {
      meta[, cell_id := rownames(meta)]
    } else {
      stop("$meta must contain cell_id column or rownames.")
    }
  }
  meta[, cell_id := as.character(cell_id)]
  meta <- unique(meta, by = "cell_id")

  if (!is.null(umap)) {
    if (!"cell_id" %in% names(umap)) {
      if (!is.null(rownames(umap))) {
        umap[, cell_id := rownames(umap)]
      } else {
        stop("$umap must contain cell_id column or rownames.")
      }
    }
    umap[, cell_id := as.character(cell_id)]
    xcol <- pick_col(umap, c("UMAP_1", "umap_1", "UMAP1", "umap1"), required = TRUE)
    ycol <- pick_col(umap, c("UMAP_2", "umap_2", "UMAP2", "umap2"), required = TRUE)
    setnames(umap, c(xcol, ycol), c("UMAP_1", "UMAP_2"), skip_absent = TRUE)
    umap <- unique(umap[, .(cell_id, UMAP_1, UMAP_2)], by = "cell_id")
  }

  common_cells <- intersect(colnames(counts), meta$cell_id)
  if (length(common_cells) == 0) stop("No overlapping cells between counts and meta in fetal reference RDS.")
  counts <- counts[, common_cells, drop = FALSE]
  meta <- meta[match(common_cells, cell_id)]
  stopifnot(identical(colnames(counts), meta$cell_id))

  if (!is.null(umap)) {
    common_umap <- intersect(common_cells, umap$cell_id)
    if (length(common_umap) > 0) {
      umap <- umap[match(common_cells, cell_id)]
      keep_ok <- !is.na(umap$cell_id)
      if (!all(keep_ok)) {
        umap <- umap[keep_ok]
      }
    } else {
      umap <- NULL
    }
  }

  rownames(counts) <- standardize_gene_symbols(rownames(counts))
  list(counts = counts, meta = meta, umap = umap, dataset_id = dataset_id, input_mode = "fetal_reference_rds")
}

score_ucell <- function(counts, gene_sets, log_file) {
  if (!requireNamespace("UCell", quietly = TRUE)) {
    log_msg("Package UCell not available; skipping UCell scoring.", .file = log_file)
    return(NULL)
  }
  log_msg("Running UCell scoring with matrix-first API...", .file = log_file)
  score_df <- UCell::ScoreSignatures_UCell(
    counts,
    features = gene_sets,
    maxRank = nrow(counts)
  )
  score_dt <- as.data.table(score_df, keep.rownames = "cell_id")
  setnames(score_dt, old = names(score_dt)[-1], new = paste0(names(score_dt)[-1], "__UCell"))
  score_dt
}

score_aucell <- function(counts, gene_sets, log_file, n_cores = 1L) {
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    log_msg("Package AUCell not available; skipping AUCell scoring.", .file = log_file)
    return(NULL)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    log_msg("Package SummarizedExperiment not available; skipping AUCell scoring.", .file = log_file)
    return(NULL)
  }

  log_msg("Running AUCell scoring...", .file = log_file)

  old_opt <- getOption("matrixStats.useNames.NA")
  options(matrixStats.useNames.NA = "deprecated")
  on.exit(options(matrixStats.useNames.NA = old_opt), add = TRUE)

  rankings <- tryCatch({
    AUCell::AUCell_buildRankings(counts, plotStats = FALSE, verbose = FALSE)
  }, error = function(e) {
    log_msg("AUCell_buildRankings failed; skipping AUCell. Error: ", conditionMessage(e), .file = log_file)
    return(NULL)
  })
  if (is.null(rankings)) return(NULL)

  auc <- tryCatch({
    AUCell::AUCell_calcAUC(gene_sets, rankings, nCores = n_cores, verbose = FALSE)
  }, error = function(e) {
    log_msg("AUCell_calcAUC failed; skipping AUCell. Error: ", conditionMessage(e), .file = log_file)
    return(NULL)
  })
  if (is.null(auc)) return(NULL)

  mat <- t(as.matrix(SummarizedExperiment::assay(auc)))
  if (nrow(mat) != ncol(counts)) mat <- t(mat)
  score_dt <- as.data.table(mat, keep.rownames = "cell_id")
  setnames(score_dt, old = names(score_dt)[-1], new = paste0(names(score_dt)[-1], "__AUCell"))
  score_dt
}

score_mean_expr <- function(counts, gene_sets, log_file) {
  log_msg("Running fallback mean-log-expression scoring...", .file = log_file)
  log_counts <- log1p(counts)
  out <- lapply(names(gene_sets), function(gs) {
    genes <- intersect(rownames(log_counts), gene_sets[[gs]])
    if (length(genes) == 0) {
      rep(NA_real_, ncol(log_counts))
    } else {
      Matrix::colMeans(log_counts[genes, , drop = FALSE])
    }
  })
  out <- as.data.table(out)
  names(out) <- paste0(names(gene_sets), "__MeanLogExpr")
  out[, cell_id := colnames(counts)]
  setcolorder(out, c("cell_id", setdiff(names(out), "cell_id")))
  out
}

identify_metadata_roles <- function(meta, celltype_col, stage_col, donor_col, region_col) {
  roles <- list(
    celltype = if (!is.null(celltype_col) && nzchar(celltype_col) && celltype_col %in% names(meta)) celltype_col else pick_col(meta, c("cell_type", "celltype", "CellType", "annotation", "Annotation"), required = FALSE),
    stage = if (!is.null(stage_col) && nzchar(stage_col) && stage_col %in% names(meta)) stage_col else pick_col(meta, c("stage", "Stage", "age", "Age", "developmental_stage", "pcw", "week"), required = FALSE),
    donor = if (!is.null(donor_col) && nzchar(donor_col) && donor_col %in% names(meta)) donor_col else pick_col(meta, c("donor", "donor_id", "sample", "sample_id", "subject", "individual"), required = FALSE),
    region = if (!is.null(region_col) && nzchar(region_col) && region_col %in% names(meta)) region_col else pick_col(meta, c("region", "brain_region", "area", "Area", "cortical_region"), required = FALSE)
  )
  roles
}

longify_scores <- function(scores_dt) {
  melt(
    as.data.table(scores_dt),
    id.vars = "cell_id",
    variable.name = "score_name",
    value.name = "score"
  )[, c("program", "method") := tstrsplit(score_name, "__", fixed = TRUE)][]
}

summarize_scores <- function(cell_annot, long_scores, by_cols) {
  dt <- merge(as.data.table(cell_annot), as.data.table(long_scores), by = "cell_id", all = FALSE)
  by_cols <- by_cols[!is.na(by_cols) & nzchar(by_cols) & by_cols %in% names(dt)]

  if (length(by_cols) == 0) {
    out <- dt[, .(
      n_cells = .N,
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      sd_score = sd(score, na.rm = TRUE)
    ), by = .(program, method)]
    return(out[])
  }

  for (bc in by_cols) {
    dt[[bc]] <- standardize_string(dt[[bc]])
  }
  keep <- complete.cases(as.data.frame(dt[, by_cols, with = FALSE]))
  dt <- dt[keep]
  if (nrow(dt) == 0) {
    return(data.table())
  }

  group_cols <- c(by_cols, "program", "method")
  sub <- as.data.frame(dt[, c(group_cols, "score"), with = FALSE], stringsAsFactors = FALSE)
  split_key <- interaction(sub[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE, sep = "")
  idx_list <- split(seq_len(nrow(sub)), split_key, drop = TRUE)

  out <- rbindlist(lapply(idx_list, function(idx) {
    d <- sub[idx, , drop = FALSE]
    keys <- d[1, group_cols, drop = FALSE]
    as.data.table(c(
      as.list(keys),
      list(
        n_cells = nrow(d),
        mean_score = mean(d$score, na.rm = TRUE),
        median_score = median(d$score, na.rm = TRUE),
        sd_score = sd(d$score, na.rm = TRUE)
      )
    ))
  }), fill = TRUE)

  setcolorder(out, c(group_cols, "n_cells", "mean_score", "median_score", "sd_score"))
  out[]
}

compute_method_concordance <- function(long_scores) {
  methods_by_program <- split(long_scores, long_scores$program)
  out <- rbindlist(lapply(names(methods_by_program), function(pg) {
    sub <- methods_by_program[[pg]]
    if (length(unique(sub$method)) < 2) return(NULL)
    wide <- dcast(sub, cell_id ~ method, value.var = "score")
    numeric_cols <- setdiff(names(wide), "cell_id")
    if (length(numeric_cols) < 2) return(NULL)
    combos <- combn(numeric_cols, 2, simplify = FALSE)
    rbindlist(lapply(combos, function(cc) {
      ok <- complete.cases(wide[[cc[1]]], wide[[cc[2]]])
      if (sum(ok) < 10) return(NULL)
      data.table(
        program = pg,
        method_a = cc[1],
        method_b = cc[2],
        n_cells = sum(ok),
        spearman = suppressWarnings(cor(wide[[cc[1]]][ok], wide[[cc[2]]][ok], method = "spearman")),
        pearson = suppressWarnings(cor(wide[[cc[1]]][ok], wide[[cc[2]]][ok], method = "pearson"))
      )
    }))
  }), fill = TRUE)
  if (is.null(out) || nrow(out) == 0) {
    data.table(program = character(), method_a = character(), method_b = character(), n_cells = integer(), spearman = numeric(), pearson = numeric())
  } else out
}

plot_heatmap <- function(summary_dt, group_x, group_y, out_file, title_txt) {
  if (!(group_x %in% names(summary_dt)) || !(group_y %in% names(summary_dt)) || nrow(summary_dt) == 0) return(invisible(NULL))
  p <- ggplot(summary_dt, aes(x = .data[[group_x]], y = .data[[group_y]], fill = mean_score)) +
    geom_tile() +
    facet_grid(method ~ program, scales = "free") +
    scale_fill_viridis_c(option = "C", na.value = "grey90") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title_txt, x = group_x, y = group_y, fill = "Mean score")
  ggsave(out_file, p, width = 12, height = 7)
}

plot_box_by_group <- function(cell_annot, long_scores, group_col, out_file, title_txt, max_levels = 20L) {
  if (!(group_col %in% names(cell_annot))) return(invisible(NULL))
  dt <- merge(cell_annot[, c("cell_id", group_col), with = FALSE], long_scores, by = "cell_id")
  dt[, group := standardize_string(get(group_col))]
  dt <- dt[!is.na(group)]
  if (nrow(dt) == 0) return(invisible(NULL))
  keep <- dt[, .N, by = group][order(-N)][1:min(max_levels, .N), group]
  dt <- dt[group %in% keep]
  p <- ggplot(dt, aes(x = group, y = score)) +
    geom_boxplot(outlier.size = 0.1) +
    facet_grid(method ~ program, scales = "free_y") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title_txt, x = group_col, y = "Cell score")
  ggsave(out_file, p, width = 12, height = 7)
}

plot_umap_scores <- function(umap_dt, scores_wide, score_cols, out_dir) {
  if (is.null(umap_dt) || nrow(umap_dt) == 0) return(invisible(NULL))
  dt <- merge(umap_dt, scores_wide, by = "cell_id")
  for (sc in score_cols) {
    if (!sc %in% names(dt)) next
    p <- ggplot(dt, aes(x = UMAP_1, y = UMAP_2, color = .data[[sc]])) +
      geom_point(size = 0.15, alpha = 0.8) +
      scale_color_viridis_c(option = "C", na.value = "grey85") +
      theme_bw(base_size = 11) +
      labs(title = sc, color = "Score")
    ggsave(file.path(out_dir, paste0("UMAP_", sanitize_name(sc), ".png")), p, width = 7, height = 5, dpi = 200)
  }
}

write_summary_txt <- function(path, obj, roles, overlap_dt, method_concordance, top_state_dt) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("Step11A fetal single-cell developmental-state localization summary", con)
  writeLines("", con)
  writeLines(paste0("dataset_id = ", obj$dataset_id), con)
  writeLines(paste0("input_mode = ", obj$input_mode), con)
  writeLines(paste0("n_cells = ", ncol(obj$counts)), con)
  writeLines(paste0("n_genes = ", nrow(obj$counts)), con)
  writeLines(paste0("n_metadata_columns = ", ncol(obj$meta)), con)
  writeLines(paste0("celltype_col = ", ifelse(is.na(roles$celltype), "NA", roles$celltype)), con)
  writeLines(paste0("stage_col = ", ifelse(is.na(roles$stage), "NA", roles$stage)), con)
  writeLines(paste0("donor_col = ", ifelse(is.na(roles$donor), "NA", roles$donor)), con)
  writeLines(paste0("region_col = ", ifelse(is.na(roles$region), "NA", roles$region)), con)
  writeLines("", con)
  writeLines("Program overlap with fetal dataset:", con)
  for (i in seq_len(nrow(overlap_dt))) {
    writeLines(
      paste0("  - ", overlap_dt$program[i], ": requested=", overlap_dt$n_requested[i],
             ", present=", overlap_dt$n_present[i], ", fraction=", signif(overlap_dt$fraction_present[i], 3)),
      con
    )
  }
  writeLines("", con)
  if (nrow(method_concordance) > 0) {
    writeLines("Method concordance:", con)
    for (i in seq_len(nrow(method_concordance))) {
      writeLines(
        paste0("  - ", method_concordance$program[i], " | ", method_concordance$method_a[i], " vs ", method_concordance$method_b[i],
               ": spearman=", signif(method_concordance$spearman[i], 4), ", pearson=", signif(method_concordance$pearson[i], 4),
               ", n_cells=", method_concordance$n_cells[i]),
        con
      )
    }
  }
  writeLines("", con)
  if (!is.null(top_state_dt) && nrow(top_state_dt) > 0) {
    writeLines("Top developmental states by mean score:", con)
    write.table(top_state_dt, con, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

option_list <- list(
  make_option("--base_dir", type = "character", default = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
              help = "Project base directory [default %default]"),
  make_option("--outdir", type = "character", default = NULL,
              help = "Output directory. Default: base_dir/Step11A_fetal_singlecell_localization_final"),
  make_option("--program_file", type = "character", default = NULL,
              help = "Program/gene mapping table. Default: base_dir/tables/20_midPrenatal_core_programs.tsv"),
  make_option("--fetal_rds", type = "character", default = NULL,
              help = "Path to fetal_singlecell_reference.rds. Default: base_dir/inputs/fetal_singlecell_reference.rds"),
  make_option("--programs", type = "character", default = "SFARI_all,midPrenatal_SFARI_top20,midPrenatal_SFARI_top05,midPrenatal_SFARI_top10",
              help = "Comma-separated program names to score [default %default]"),
  make_option("--min_genes", type = "integer", default = 5L,
              help = "Minimum overlapping genes required to retain a program [default %default]"),
  make_option("--n_cores", type = "integer", default = 8L,
              help = "Number of cores for AUCell if used [default %default]"),
  make_option("--celltype_col", type = "character", default = "cell_type",
              help = "Cell-type column in fetal metadata [default %default]"),
  make_option("--stage_col", type = "character", default = "stage",
              help = "Stage column in fetal metadata [default %default]"),
  make_option("--donor_col", type = "character", default = "donor",
              help = "Donor column in fetal metadata [default %default]"),
  make_option("--region_col", type = "character", default = "region",
              help = "Region column in fetal metadata [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

base_dir <- opt$base_dir
outdir <- if (is.null(opt$outdir)) file.path(base_dir, "Step11A_fetal_singlecell_localization_final") else opt$outdir
program_file <- if (is.null(opt$program_file)) file.path(base_dir, "tables", "20_midPrenatal_core_programs.tsv") else opt$program_file
fetal_rds <- if (is.null(opt$fetal_rds)) file.path(base_dir, "inputs", "fetal_singlecell_reference.rds") else opt$fetal_rds
requested_programs <- trimws(strsplit(opt$programs, ",", fixed = TRUE)[[1]])
requested_programs <- requested_programs[nzchar(requested_programs)]

safe_dir_create(outdir)
safe_dir_create(file.path(outdir, "tables"))
safe_dir_create(file.path(outdir, "plots"))
safe_dir_create(file.path(outdir, "logs"))
log_file <- file.path(outdir, "logs", "00_step11A_fetal_localization.log")
if (file.exists(log_file)) file.remove(log_file)

log_msg("Starting Step11A fetal single-cell developmental-state localization (final RDS-only version)", .file = log_file)
log_msg("base_dir=", base_dir, .file = log_file)
log_msg("outdir=", outdir, .file = log_file)
log_msg("program_file=", program_file, .file = log_file)
log_msg("fetal_rds=", fetal_rds, .file = log_file)
log_msg("requested_programs=", paste(requested_programs, collapse = ", "), .file = log_file)

if (!file.exists(program_file)) stop("Program file not found: ", program_file)
if (!file.exists(fetal_rds)) stop("Fetal RDS not found: ", fetal_rds)

obj <- load_fetal_reference_rds(fetal_rds, log_file)
counts <- obj$counts
meta <- obj$meta
umap <- obj$umap

roles <- identify_metadata_roles(
  meta = meta,
  celltype_col = opt$celltype_col,
  stage_col = opt$stage_col,
  donor_col = opt$donor_col,
  region_col = opt$region_col
)

for (nm in c("celltype", "stage", "donor", "region")) {
  if (is.na(roles[[nm]])) {
    log_msg("Metadata role unresolved: ", nm, .file = log_file)
  } else {
    log_msg("Metadata role ", nm, " -> ", roles[[nm]], .file = log_file)
  }
}

gene_sets_raw <- load_program_gene_sets(program_file, requested_programs = requested_programs, min_genes = opt$min_genes)
genes_present <- unique(rownames(counts))

gene_overlap <- rbindlist(lapply(names(gene_sets_raw), function(pg) {
  present <- intersect(gene_sets_raw[[pg]], genes_present)
  data.table(
    program = pg,
    n_requested = length(unique(gene_sets_raw[[pg]])),
    n_present = length(unique(present)),
    fraction_present = length(unique(present)) / max(1, length(unique(gene_sets_raw[[pg]])))
  )
}))
safe_fwrite(gene_overlap, file.path(outdir, "tables", "01_program_gene_overlap.tsv"))

gene_sets <- lapply(gene_sets_raw, function(gs) intersect(gs, genes_present))
gene_sets <- gene_sets[lengths(gene_sets) >= opt$min_genes]
if (length(gene_sets) == 0) stop("No program retained after overlap with fetal dataset and min_genes filter.")

program_gene_map <- rbindlist(lapply(names(gene_sets), function(pg) {
  data.table(program = pg, gene_symbol = gene_sets[[pg]])
}))
safe_fwrite(program_gene_map, file.path(outdir, "tables", "02_program_gene_map_retained.tsv.gz"))

ucell_dt <- score_ucell(counts, gene_sets, log_file)
aucell_dt <- score_aucell(counts, gene_sets, log_file, n_cores = opt$n_cores)
mean_dt <- score_mean_expr(counts, gene_sets, log_file)
score_dts <- Filter(Negate(is.null), list(ucell_dt, aucell_dt, mean_dt))
if (length(score_dts) == 0) stop("No scoring output generated.")

scores_wide <- Reduce(function(x, y) merge(x, y, by = "cell_id", all = TRUE), score_dts)
safe_fwrite(scores_wide, file.path(outdir, "tables", "03_cell_program_scores_wide.tsv.gz"))
long_scores <- longify_scores(scores_wide)
safe_fwrite(long_scores, file.path(outdir, "tables", "04_cell_program_scores_long.tsv.gz"))

input_summary <- data.table(
  metric = c("dataset_id", "input_mode", "n_cells", "n_genes", "n_metadata_columns", "celltype_col", "stage_col", "donor_col", "region_col"),
  value = c(obj$dataset_id, obj$input_mode, ncol(counts), nrow(counts), ncol(meta),
            ifelse(is.na(roles$celltype), NA_character_, roles$celltype),
            ifelse(is.na(roles$stage), NA_character_, roles$stage),
            ifelse(is.na(roles$donor), NA_character_, roles$donor),
            ifelse(is.na(roles$region), NA_character_, roles$region))
)
safe_fwrite(input_summary, file.path(outdir, "tables", "00_input_summary.tsv"))

cell_annot <- copy(meta[, .(cell_id)])
for (nm in c("celltype", "stage", "donor", "region")) {
  col_nm <- roles[[nm]]
  cell_annot[, (nm) := if (!is.na(col_nm)) as.character(meta[[col_nm]]) else NA_character_]
}

celltype_summary <- summarize_scores(cell_annot, long_scores, by_cols = c("celltype"))
stage_summary <- summarize_scores(cell_annot, long_scores, by_cols = c("stage"))
celltype_stage_summary <- summarize_scores(cell_annot, long_scores, by_cols = c("celltype", "stage"))
donor_summary <- summarize_scores(cell_annot, long_scores, by_cols = c("donor"))
region_summary <- summarize_scores(cell_annot, long_scores, by_cols = c("region"))

safe_fwrite(celltype_summary, file.path(outdir, "tables", "05_celltype_program_summary.tsv"))
safe_fwrite(stage_summary, file.path(outdir, "tables", "06_stage_program_summary.tsv"))
safe_fwrite(celltype_stage_summary, file.path(outdir, "tables", "07_celltype_stage_program_summary.tsv"))
safe_fwrite(donor_summary, file.path(outdir, "tables", "08_donor_program_summary.tsv"))
safe_fwrite(region_summary, file.path(outdir, "tables", "09_region_program_summary.tsv"))

auto_top_state <- copy(celltype_stage_summary)
if (nrow(auto_top_state) > 0) {
  setorder(auto_top_state, program, method, -mean_score)
  top_state_dt <- auto_top_state[, head(.SD, 10), by = .(program, method)]
} else {
  top_state_dt <- data.table()
}
safe_fwrite(top_state_dt, file.path(outdir, "tables", "10_top_developmental_states_by_score.tsv"))

method_concordance <- compute_method_concordance(long_scores)
safe_fwrite(method_concordance, file.path(outdir, "tables", "11_method_concordance.tsv"))

plot_dir <- file.path(outdir, "plots")
plot_heatmap(celltype_stage_summary, "stage", "celltype", file.path(plot_dir, "heatmap_stage_by_celltype_mean_scores.png"),
             "Fetal localization: mean program score by stage and cell type")
plot_box_by_group(cell_annot, long_scores, "celltype", file.path(plot_dir, "boxplot_scores_by_celltype.png"),
                  "Fetal localization by cell type")
plot_box_by_group(cell_annot, long_scores, "stage", file.path(plot_dir, "boxplot_scores_by_stage.png"),
                  "Fetal localization by developmental stage")

preferred_score_cols <- c(
  paste0(names(gene_sets), "__UCell"),
  paste0(names(gene_sets), "__AUCell"),
  paste0(names(gene_sets), "__MeanLogExpr")
)
preferred_score_cols <- preferred_score_cols[preferred_score_cols %in% names(scores_wide)]
plot_umap_scores(umap, scores_wide, preferred_score_cols, plot_dir)

write_summary_txt(
  path = file.path(outdir, "12_fetal_localization_summary.txt"),
  obj = obj,
  roles = roles,
  overlap_dt = gene_overlap,
  method_concordance = method_concordance,
  top_state_dt = top_state_dt
)

writeLines(capture.output(sessionInfo()), con = file.path(outdir, "logs", "sessionInfo.txt"))
log_msg("Step11A fetal localization completed successfully.", .file = log_file)
