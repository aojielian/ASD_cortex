#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
OUTDIR <- file.path(BASE, "step01_prepare")
METADIR <- file.path(OUTDIR, "meta")
TABDIR <- file.path(OUTDIR, "tables")
RDSDIR <- file.path(OUTDIR, "rds")
LOGDIR <- file.path(OUTDIR, "logs")

dir.create(METADIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDSDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(LOGDIR, "Step01B_BrainSpan_cortex_convergence_v2.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

normalize_nm <- function(x) {
  x <- as.character(x)
  x <- gsub("^X", "", x)
  x <- gsub("^\\.+", "", x)
  x <- gsub('"', "", x, fixed = TRUE)
  x <- trimws(x)
  x
}

age_to_stage <- function(age_vec) {
  x <- tolower(trimws(as.character(age_vec)))
  out <- rep("other", length(x))

  # prenatal pcw
  num_pcw <- suppressWarnings(as.numeric(gsub(" .*", "", sub(" pcw.*", "", x))))
  is_pcw <- grepl("pcw", x)

  out[is_pcw & !is.na(num_pcw) & num_pcw < 13] <- "early_prenatal"
  out[is_pcw & !is.na(num_pcw) & num_pcw >= 13 & num_pcw < 24] <- "mid_prenatal"
  out[is_pcw & !is.na(num_pcw) & num_pcw >= 24] <- "late_prenatal"

  out[grepl("^4 mos$|^10 mos$", x)] <- "birth"
  out[grepl("^1 yrs$|^2 yrs$|^3 yrs$", x)] <- "infancy"
  out[grepl("^4 yrs$|^8 yrs$|^11 yrs$", x)] <- "childhood"
  out[grepl("^13 yrs$|^15 yrs$|^18 yrs$|^19 yrs$", x)] <- "adolescence"
  out[grepl("^21 yrs$|^23 yrs$|^30 yrs$|^36 yrs$|^37 yrs$|^40 yrs$", x)] <- "adult"

  out
}

region_to_group <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep("other", length(y))

  out[grepl("prefrontal|frontal|orbitofrontal|dorsolateral prefrontal|ventrolateral prefrontal|medial prefrontal|inferolateral prefrontal|orbital frontal", y)] <- "frontal"
  out[grepl("primary motor|primary somatosensory|somatosensory|motor cortex|m1c|s1c|paracentral", y)] <- "sensorimotor"
  out[grepl("posterior parietal|inferior parietal|superior parietal|angular gyrus|parietal", y)] <- "parietal"
  out[grepl("superior temporal|inferolateral temporal|temporal|auditory cortex|stc", y)] <- "temporal"
  out[grepl("primary visual|visual cortex|occipital|v1c", y)] <- "occipital"

  out
}

pick_col <- function(nms, candidates = NULL, pattern = NULL, required = TRUE, label = "column") {
  hit <- character(0)
  if (!is.null(candidates)) {
    hit <- intersect(candidates, nms)
  }
  if (length(hit) == 0 && !is.null(pattern)) {
    hit <- grep(pattern, nms, ignore.case = TRUE, value = TRUE)
  }
  if (length(hit) == 0) {
    if (required) stop("Could not identify ", label, " among columns: ", paste(nms, collapse = ", "))
    return(NA_character_)
  }
  hit[1]
}

match_expression_to_metadata <- function(expr_cols, cols_meta) {
  # expr_cols excludes first feature column
  expr_orig <- as.character(expr_cols)
  expr_norm <- normalize_nm(expr_orig)

  cols_meta <- copy(cols_meta)
  cols_meta[, column_num_chr := as.character(column_num)]
  cols_meta[, column_num_norm := normalize_nm(column_num_chr)]

  match_name <- data.table(
    expr_col_index = seq_along(expr_orig),
    expr_col_original = expr_orig,
    expr_col_normalized = expr_norm,
    matched = expr_norm %in% cols_meta$column_num_norm
  )

  map_name <- cols_meta[match_name, on = .(column_num_norm = expr_col_normalized)]
  map_name[, match_method := "by_name"]
  by_name_ok <- sum(!is.na(map_name$donor_id))

  if (by_name_ok == length(expr_orig)) {
    return(list(map = map_name, method = "by_name"))
  }

  # fallback to pure positional mapping if dimensions agree
  if (nrow(cols_meta) == length(expr_orig)) {
    map_pos <- copy(cols_meta)
    map_pos[, expr_col_index := seq_len(.N)]
    map_pos[, expr_col_original := expr_orig]
    map_pos[, expr_col_normalized := expr_norm]
    map_pos[, matched := TRUE]
    map_pos[, match_method := "by_position"]
    return(list(map = map_pos, method = "by_position"))
  }

  # hybrid: keep name matches and fill unmatched by position only where unambiguous
  map_h <- copy(map_name)
  miss_idx <- which(is.na(map_h$donor_id))
  if (length(miss_idx) > 0 && nrow(cols_meta) == length(expr_orig)) {
    map_h[miss_idx, c("column_num","donor_id","donor_name","age","gender","structure_id","structure_acronym","structure_name") :=
            cols_meta[miss_idx, .(column_num, donor_id, donor_name, age, gender, structure_id, structure_acronym, structure_name)]]
    map_h[miss_idx, match_method := "hybrid_position_fill"]
  }
  if (all(!is.na(map_h$donor_id))) {
    return(list(map = map_h, method = "hybrid"))
  }

  return(list(map = map_name, method = "failed"))
}

set.seed(20260319)

sfari_file <- file.path(BASE, "sfari.genes.Sand1.txt")
brainspan_dir <- file.path(BASE, "Gencode_v3c_summarized_to_genes")
cols_file <- file.path(brainspan_dir, "columns_metadata.csv")
rows_file <- file.path(brainspan_dir, "rows_metadata.csv")
expr_file <- file.path(brainspan_dir, "expression_matrix.csv")

log_msg("Reading SFARI: ", sfari_file)
sfari <- fread(sfari_file, header = FALSE, sep = "\t", fill = TRUE)
setnames(sfari, "V1", "gene_symbol")
sfari[, gene_symbol := toupper(trimws(gene_symbol))]
sfari <- unique(sfari[gene_symbol != "" & !is.na(gene_symbol)])
safe_fwrite(sfari, file.path(TABDIR, "01_SFARI_primary_Sand1_standardized.tsv"))
log_msg("SFARI unique genes: ", nrow(sfari))

log_msg("Reading BrainSpan columns metadata: ", cols_file)
cols_meta <- fread(cols_file)
log_msg("Reading BrainSpan rows metadata: ", rows_file)
rows_meta <- fread(rows_file)

safe_fwrite(
  data.table(table = "columns_metadata", column_name = names(cols_meta)),
  file.path(METADIR, "12a_columns_metadata_dictionary.tsv")
)
safe_fwrite(
  data.table(table = "rows_metadata", column_name = names(rows_meta)),
  file.path(METADIR, "12a_rows_metadata_dictionary.tsv")
)

col_num_col   <- pick_col(names(cols_meta), candidates = c("column_num"), pattern = "^column_num$", label = "columns_metadata column_num")
donor_id_col  <- pick_col(names(cols_meta), candidates = c("donor_id"), pattern = "donor_id", label = "columns_metadata donor_id")
donor_nm_col  <- pick_col(names(cols_meta), candidates = c("donor_name"), pattern = "donor_name", required = FALSE, label = "columns_metadata donor_name")
age_col       <- pick_col(names(cols_meta), candidates = c("age"), pattern = "^age$", label = "columns_metadata age")
struct_id_col <- pick_col(names(cols_meta), candidates = c("structure_id"), pattern = "structure_id", label = "columns_metadata structure_id")
struct_acr_col<- pick_col(names(cols_meta), candidates = c("structure_acronym"), pattern = "structure_acronym", required = FALSE, label = "columns_metadata structure_acronym")
struct_nm_col <- pick_col(names(cols_meta), candidates = c("structure_name"), pattern = "structure_name", label = "columns_metadata structure_name")
gender_col    <- pick_col(names(cols_meta), candidates = c("gender"), pattern = "^gender$", required = FALSE, label = "columns_metadata gender")

gene_sym_col  <- pick_col(names(rows_meta), candidates = c("gene_symbol"), pattern = "gene_symbol", label = "rows_metadata gene_symbol")
gene_id_col   <- pick_col(names(rows_meta), candidates = c("gene_id"), pattern = "^gene_id$", required = FALSE, label = "rows_metadata gene_id")
entrez_col    <- pick_col(names(rows_meta), candidates = c("entrez_id"), pattern = "entrez", required = FALSE, label = "rows_metadata entrez_id")
ensg_col      <- pick_col(names(rows_meta), candidates = c("ensembl_gene_id"), pattern = "ensembl", required = FALSE, label = "rows_metadata ensembl_gene_id")

cols_std <- data.table(
  column_num = cols_meta[[col_num_col]],
  donor_id = cols_meta[[donor_id_col]],
  donor_name = if (!is.na(donor_nm_col)) cols_meta[[donor_nm_col]] else NA_character_,
  age = cols_meta[[age_col]],
  gender = if (!is.na(gender_col)) cols_meta[[gender_col]] else NA_character_,
  structure_id = cols_meta[[struct_id_col]],
  structure_acronym = if (!is.na(struct_acr_col)) cols_meta[[struct_acr_col]] else NA_character_,
  structure_name = cols_meta[[struct_nm_col]]
)

rows_std <- data.table(
  gene_symbol = toupper(trimws(as.character(rows_meta[[gene_sym_col]]))),
  gene_id = if (!is.na(gene_id_col)) as.character(rows_meta[[gene_id_col]]) else NA_character_,
  entrez_id = if (!is.na(entrez_col)) as.character(rows_meta[[entrez_col]]) else NA_character_,
  ensembl_gene_id = if (!is.na(ensg_col)) as.character(rows_meta[[ensg_col]]) else NA_character_
)

safe_fwrite(cols_std, file.path(METADIR, "12b_BrainSpan_columns_metadata_standardized.tsv"))
safe_fwrite(rows_std, file.path(METADIR, "12c_BrainSpan_rows_metadata_standardized.tsv"))

# Read expression header only
log_msg("Reading expression header: ", expr_file)
expr_header <- names(fread(expr_file, nrows = 0))
feature_col <- expr_header[1]
expr_cols <- expr_header[-1]

match_res <- match_expression_to_metadata(expr_cols, cols_std)
expr_map <- as.data.table(match_res$map)

safe_fwrite(expr_map, file.path(METADIR, "08b_BrainSpan_expression_column_match.tsv"))
log_msg("Expression-to-metadata match method: ", match_res$method)
log_msg("Matched columns: ", sum(!is.na(expr_map$donor_id)), " / ", nrow(expr_map))

if (match_res$method == "failed" || any(is.na(expr_map$donor_id))) {
  stop("BrainSpan expression columns still could not be mapped after fallback matching. Check 08b_BrainSpan_expression_column_match.tsv")
}

# Cortex selection
expr_map[, stage_group := age_to_stage(age)]
expr_map[, region_group := region_to_group(structure_name)]
expr_map[, is_cortex := region_group != "other"]

safe_fwrite(expr_map, file.path(TABDIR, "08_BrainSpan_samples_harmonized.tsv"))

stage_summary <- expr_map[is_cortex == TRUE, .N, by = .(stage_group)][order(stage_group)]
stage_region_summary <- expr_map[is_cortex == TRUE, .N, by = .(stage_group, region_group)][order(stage_group, region_group)]

safe_fwrite(stage_summary, file.path(TABDIR, "10_BrainSpan_cortex_stage_summary.tsv"))
safe_fwrite(stage_region_summary, file.path(TABDIR, "10b_BrainSpan_cortex_stage_region_summary.tsv"))

cortex_idx <- expr_map[is_cortex == TRUE, expr_col_index]
cortex_cols <- expr_map[is_cortex == TRUE]

log_msg("Cortex samples retained: ", length(cortex_idx))
if (length(cortex_idx) < 20) {
  stop("Too few cortex samples retained. Check region mapping in 08_BrainSpan_samples_harmonized.tsv")
}

# Read full expression matrix
log_msg("Reading full BrainSpan expression matrix. This may take a while.")
expr_dt <- fread(expr_file)
log_msg("Expression matrix dimensions: ", nrow(expr_dt), " x ", ncol(expr_dt))

# Align rows metadata
if (nrow(expr_dt) != nrow(rows_std)) {
  stop("Expression rows (", nrow(expr_dt), ") do not match rows_metadata rows (", nrow(rows_std), ").")
}

expr_gene <- cbind(rows_std, expr_dt[, -1, with = FALSE])
rm(expr_dt); gc()

# subset cortex columns
expr_cortex <- expr_gene[, c("gene_symbol", "gene_id", "entrez_id", "ensembl_gene_id", cortex_cols$expr_col_original), with = FALSE]

# collapse duplicated gene symbols by mean
sample_cols <- cortex_cols$expr_col_original
expr_long <- melt(expr_cortex,
                  id.vars = c("gene_symbol", "gene_id", "entrez_id", "ensembl_gene_id"),
                  measure.vars = sample_cols,
                  variable.name = "sample_col",
                  value.name = "expr")
expr_long[, expr := as.numeric(expr)]
expr_collapsed <- expr_long[gene_symbol != "" & !is.na(gene_symbol),
                            .(expr = mean(expr, na.rm = TRUE)),
                            by = .(gene_symbol, sample_col)]
expr_wide <- dcast(expr_collapsed, gene_symbol ~ sample_col, value.var = "expr")
expr_mat <- as.matrix(expr_wide[, -1])
rownames(expr_mat) <- expr_wide$gene_symbol

# save subset object
brainspan_obj <- list(
  expr = expr_mat,
  samples = cortex_cols,
  gene_map = rows_std
)
saveRDS(brainspan_obj, file.path(RDSDIR, "12_BrainSpan_cortex_subset.rds"))
log_msg("Saved cortex subset RDS.")

# observed convergence
risk_genes <- intersect(sfari$gene_symbol, rownames(expr_mat))
nonrisk_genes <- setdiff(rownames(expr_mat), risk_genes)

log_msg("Risk genes present in BrainSpan cortex matrix: ", length(risk_genes))

if (length(risk_genes) < 50) {
  stop("Too few SFARI genes matched in BrainSpan cortex matrix.")
}

gene_means <- rowMeans(expr_mat, na.rm = TRUE)
gene_mean_dt <- data.table(gene_symbol = names(gene_means), mean_expr = as.numeric(gene_means))
gene_mean_dt[, bin := cut(mean_expr,
                          breaks = quantile(mean_expr, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
                          include.lowest = TRUE, duplicates.ok = TRUE)]
risk_bin_dt <- gene_mean_dt[gene_symbol %in% risk_genes]
nonrisk_bin_dt <- gene_mean_dt[gene_symbol %in% nonrisk_genes]

sample_annot <- cortex_cols[, .(sample_col = expr_col_original, stage_group, region_group)]

calc_obs <- function(sample_ids) {
  col_idx <- match(sample_ids, colnames(expr_mat))
  mean(rowMeans(expr_mat[risk_genes, col_idx, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
}

perm_once <- function(sample_ids) {
  matched_pool <- character(0)
  for (b in unique(risk_bin_dt$bin)) {
    n_b <- sum(risk_bin_dt$bin == b, na.rm = TRUE)
    pool_b <- nonrisk_bin_dt[bin == b, gene_symbol]
    if (length(pool_b) < n_b) {
      pool_b <- nonrisk_bin_dt$gene_symbol
    }
    matched_pool <- c(matched_pool, sample(pool_b, n_b, replace = FALSE))
  }
  col_idx <- match(sample_ids, colnames(expr_mat))
  mean(rowMeans(expr_mat[matched_pool, col_idx, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
}

groups <- unique(sample_annot[, .(stage_group, region_group)])
groups <- groups[stage_group != "other" & region_group != "other"]

res_list <- vector("list", nrow(groups))
nperm <- 200

log_msg("Running convergence analysis with ", nperm, " permutations per stage-region group.")

for (i in seq_len(nrow(groups))) {
  stg <- groups$stage_group[i]
  reg <- groups$region_group[i]
  sids <- sample_annot[stage_group == stg & region_group == reg, sample_col]

  obs <- calc_obs(sids)
  perm_vals <- replicate(nperm, perm_once(sids))

  emp_p <- (sum(perm_vals >= obs) + 1) / (length(perm_vals) + 1)
  z <- ifelse(sd(perm_vals) > 0, (obs - mean(perm_vals)) / sd(perm_vals), NA_real_)

  res_list[[i]] <- data.table(
    stage_group = stg,
    region_group = reg,
    n_samples = length(sids),
    n_risk_genes = length(risk_genes),
    observed_mean_expr = obs,
    perm_mean_expr = mean(perm_vals),
    perm_sd_expr = sd(perm_vals),
    z_score = z,
    empirical_p = emp_p
  )
}

res <- rbindlist(res_list)
res[, fdr := p.adjust(empirical_p, method = "BH")]
setorder(res, empirical_p, -z_score)

safe_fwrite(res, file.path(TABDIR, "11_BrainSpan_stage_region_convergence.tsv"))

res_stage <- sample_annot[stage_group != "other", .(sample_col = sample_col), by = .(stage_group)]
res_stage_list <- vector("list", nrow(res_stage[, .N, by = stage_group]))

stage_only <- unique(sample_annot[, .(stage_group)])
stage_only <- stage_only[stage_group != "other"]

for (i in seq_len(nrow(stage_only))) {
  stg <- stage_only$stage_group[i]
  sids <- sample_annot[stage_group == stg, sample_col]
  obs <- calc_obs(sids)
  perm_vals <- replicate(nperm, perm_once(sids))
  emp_p <- (sum(perm_vals >= obs) + 1) / (length(perm_vals) + 1)
  z <- ifelse(sd(perm_vals) > 0, (obs - mean(perm_vals)) / sd(perm_vals), NA_real_)

  res_stage_list[[i]] <- data.table(
    stage_group = stg,
    n_samples = length(sids),
    n_risk_genes = length(risk_genes),
    observed_mean_expr = obs,
    perm_mean_expr = mean(perm_vals),
    perm_sd_expr = sd(perm_vals),
    z_score = z,
    empirical_p = emp_p
  )
}

res_stage <- rbindlist(res_stage_list)
res_stage[, fdr := p.adjust(empirical_p, method = "BH")]
setorder(res_stage, empirical_p, -z_score)

safe_fwrite(res_stage, file.path(TABDIR, "11b_BrainSpan_stage_only_convergence.tsv"))

run_summary <- rbindlist(list(
  data.table(section = "BrainSpan", metric = "expression_match_method", value = match_res$method),
  data.table(section = "BrainSpan", metric = "n_total_samples", value = nrow(cols_std)),
  data.table(section = "BrainSpan", metric = "n_cortex_samples", value = nrow(cortex_cols)),
  data.table(section = "BrainSpan", metric = "n_risk_genes_in_matrix", value = length(risk_genes)),
  data.table(section = "BrainSpan", metric = "n_stage_region_tests", value = nrow(res)),
  data.table(section = "BrainSpan", metric = "best_stage_region", value = if (nrow(res) > 0) paste(res$stage_group[1], res$region_group[1], sep = "__") else NA_character_),
  data.table(section = "BrainSpan", metric = "best_stage_region_empirical_p", value = if (nrow(res) > 0) as.character(res$empirical_p[1]) else NA_character_)
))

safe_fwrite(run_summary, file.path(METADIR, "12b_Step01B_run_summary.tsv"))

log_msg("Step01B v2 completed successfully.")
