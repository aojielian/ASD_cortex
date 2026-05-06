#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

parse_args <- function(x){
  out <- list(
    base = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    syngo_gene_file = NA_character_,
    syngo_gene_col = NA_character_,
    prepost_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/reviewer_round3/MA_round1B_prepost_CI",
    gandal_rdata = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData",
    outdir = NA_character_
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    val <- if (i < length(x)) x[[i+1L]] else NA_character_
    if (key == "--base") out$base <- val
    if (key == "--syngo_gene_file") out$syngo_gene_file <- val
    if (key == "--syngo_gene_col") out$syngo_gene_col <- val
    if (key == "--prepost_dir") out$prepost_dir <- val
    if (key == "--gandal_rdata") out$gandal_rdata <- val
    if (key == "--outdir") out$outdir <- val
    i <- i + 2L
  }
  if (is.na(out$outdir)) out$outdir <- file.path(out$base, "reviewer_round3", "MA_round1C_synaptic_nonsynaptic")
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "meta"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "MA_round1C_synaptic_nonsynaptic.log")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
safe_fwrite <- function(x, file) fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  x <- sub("\\..*$", "", x)
  x[nchar(x) == 0L] <- NA_character_
  toupper(x)
}
normalize_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub('^"+|"+$', '', x)
  x <- gsub("^'+|'+$", "", x)
  x <- sub("^(ASD|CTL)_", "", x)
  gsub("\\s+", "", x)
}
extract_attr <- function(x, key) {
  m <- regexec(paste0(key, ' "([^"]+)"'), x)
  regmatches(x, m) |> lapply(function(z) if (length(z) >= 2) z[2] else NA_character_) |> unlist(use.names = FALSE)
}
compute_logcpm <- function(count_dt) {
  mat <- as.matrix(count_dt[, -1]); rownames(mat) <- count_dt[[1]]; storage.mode(mat) <- "double"
  lib <- colSums(mat); lib[lib <= 0] <- NA_real_
  log2(t(t(mat) / lib) * 1e6 + 1)
}
read_any <- function(f) {
  if (!file.exists(f)) stop("Missing file: ", f)
  ext <- tolower(tools::file_ext(f))
  if (ext %in% c("xlsx","xls")) {
    sheets <- openxlsx::getSheetNames(f)
    if (length(sheets) < 1) stop("No sheets found in Excel file: ", f)
    return(as.data.table(openxlsx::read.xlsx(f, sheet = sheets[1], colNames = TRUE)))
  }
  if (grepl("\\.gz$", f, ignore.case = TRUE)) return(fread(cmd = paste("zcat", shQuote(f))))
  fread(f)
}
fit_ci <- function(dt, program_name, cohort, model_type, fraction_covars = character()) {
  dt <- copy(dt)
  dt <- dt[!is.na(score) & !is.na(diagnosis)]
  dt[, diagnosis := factor(as.character(diagnosis), levels = c("Control","ASD"))]
  rhs <- c("diagnosis", fraction_covars)
  rhs <- rhs[rhs %in% names(dt)]
  cc <- dt[complete.cases(dt[, ..rhs], score)]
  if (sum(cc$diagnosis == "ASD") < 3 || sum(cc$diagnosis == "Control") < 3) return(NULL)
  fml <- as.formula(paste("score ~", paste(rhs, collapse = " + ")))
  fit <- lm(fml, data = cc)
  sm <- summary(fit)$coefficients
  if (!("diagnosisASD" %in% rownames(sm))) return(NULL)
  beta <- sm["diagnosisASD","Estimate"]; se <- sm["diagnosisASD","Std. Error"]; p <- sm["diagnosisASD","Pr(>|t|)"]
  data.table(
    cohort = cohort, program_name = program_name, model_type = model_type,
    n_samples = nrow(cc), n_asd = sum(cc$diagnosis == "ASD"), n_control = sum(cc$diagnosis == "Control"),
    mean_asd = mean(cc$score[cc$diagnosis == "ASD"]), mean_control = mean(cc$score[cc$diagnosis == "Control"]),
    beta_asd = beta, se = se, ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se, p_value = p
  )
}

prog_file <- file.path(args$base, "tables", "20_midPrenatal_core_programs.tsv")
prog <- fread(prog_file)
prog[, gene_symbol := clean_symbol(gene_symbol)]
sfari_all <- unique(prog[program_name == "SFARI_all", gene_symbol])
if (length(sfari_all) == 0) stop("SFARI_all not found")

find_syngo_file <- function(base) {
  ff <- list.files(base, recursive = TRUE, full.names = TRUE)
  ff <- ff[grepl("syngo", basename(ff), ignore.case = TRUE) & grepl("\\.(tsv|txt|csv|tsv.gz|txt.gz|csv.gz|xlsx|xls)$", ff, ignore.case = TRUE)]
  if (!length(ff)) return(NA_character_)
  ff[1]
}
syn_file <- if (!is.na(args$syngo_gene_file)) args$syngo_gene_file else find_syngo_file(args$base)
if (is.na(syn_file) || !file.exists(syn_file)) stop("Could not identify SynGO gene annotation file. Provide --syngo_gene_file explicitly.")
syn <- read_any(syn_file)
if (is.na(args$syngo_gene_col)) {
  cand <- c("hgnc_symbol","gene_symbol","gene","symbol","hgnc","hgnc_id")
  hit <- cand[cand %in% names(syn)][1]
  if (is.na(hit)) stop("Could not identify gene column in SynGO file. Provide --syngo_gene_col.")
  args$syngo_gene_col <- hit
}
syn[, gene_symbol := clean_symbol(get(args$syngo_gene_col))]
syn_genes <- unique(syn[!is.na(gene_symbol), gene_symbol])

sets <- list(
  SFARI_all = sfari_all,
  SFARI_all_synaptic = intersect(sfari_all, syn_genes),
  SFARI_all_nonSynGO = setdiff(sfari_all, syn_genes)
)
set_inventory <- data.table(gene_set = names(sets), n_genes = sapply(sets, length))
safe_fwrite(set_inventory, file.path(args$outdir, "tables", "MA_round1C_gene_set_inventory.tsv"))

fraction_files <- c(
  GSE102741 = file.path(args$prepost_dir, "tables", "MA_round1B_GSE102741_nnls_fractions.tsv"),
  GSE64018 = file.path(args$prepost_dir, "tables", "MA_round1B_GSE64018_nnls_fractions.tsv"),
  Gandal2022 = file.path(args$prepost_dir, "tables", "MA_round1B_Gandal2022_nnls_fractions.tsv")
)
for (f in fraction_files) if (!file.exists(f)) stop("Missing fraction table from Round1B: ", f)
fraction_dt <- lapply(names(fraction_files), function(nm) fread(fraction_files[[nm]]))
names(fraction_dt) <- names(fraction_files)

load_gse102741 <- function(base) {
  rds_file <- file.path(base, "rds", "22_GSE102741_geneSymbol_log2cpm.rds")
  meta_file <- file.path(base, "tables", "21_GSE102741_sample_metadata.tsv")
  obj <- readRDS(rds_file)
  expr <- as.matrix(obj$log2cpm); rownames(expr) <- clean_symbol(rownames(expr))
  meta <- fread(meta_file)
  sample_col <- intersect(c("sample_id","sample","Sample.ID"), names(meta))[1]
  dx_col <- intersect(c("diagnosis","Diagnosis","dx"), names(meta))[1]
  meta <- meta[, .(sample = as.character(get(sample_col)), diagnosis = as.character(get(dx_col)))]
  meta[, diagnosis := fifelse(diagnosis %in% c("ASD","Autism"), "ASD", fifelse(diagnosis %in% c("Control","CTL","CTRL"), "Control", diagnosis))]
  keep <- intersect(colnames(expr), meta$sample)
  expr <- expr[, keep, drop = FALSE]; meta <- meta[match(keep, sample)]
  list(expr = expr, meta = meta)
}
load_gse64018 <- function(base) {
  counts_file <- file.path(base, "tables", "122_GSE64018_collapsed_counts_by_symbol.tsv.gz")
  map_file <- file.path(base, "meta", "120_GSE64018_explicit_sample_mapping.tsv")
  counts <- fread(cmd = paste("zcat", shQuote(counts_file))); expr <- compute_logcpm(counts)
  sm <- fread(map_file)
  norm_col <- intersect(c("file_sample_norm","sample","sample_norm"), names(sm))[1]
  dx_col <- intersect(c("diagnosis","group","dx"), names(sm))[1]
  title_col <- intersect(c("geo_title","title"), names(sm))[1]
  if (is.na(dx_col)) {
    if (!is.na(title_col)) { sm[, diagnosis := fifelse(startsWith(get(title_col), "ASD_"), "ASD", "Control")]; dx_col <- "diagnosis" } else stop("Could not detect diagnosis col for GSE64018")
  }
  sm[, sample := normalize_sample_id(get(norm_col))]
  sm[, diagnosis := fifelse(as.character(get(dx_col)) %in% c("ASD","Autism"), "ASD", fifelse(as.character(get(dx_col)) %in% c("Control","CTL"), "Control", as.character(get(dx_col))))]
  sm <- unique(sm[, .(sample, diagnosis)])
  meta <- data.table(sample = normalize_sample_id(colnames(expr)))
  meta <- merge(meta, sm, by = "sample", all.x = TRUE, sort = FALSE)
  colnames(expr) <- meta$sample
  keep <- !is.na(meta$diagnosis)
  list(expr = expr[, keep, drop = FALSE], meta = meta[keep])
}
load_gandal <- function(rdata_path, gtf = "/gpfs/hpc/home/lijc/lianaoj/reference/gencode.v49.basic.annotation.gtf.gz") {
  gtf_dt <- fread(cmd = paste("zcat", shQuote(gtf), "| awk '$3==\"gene\"'"), sep = "\t", header = FALSE)
  gtf_map <- unique(data.table(
    ensembl_gene_id = clean_symbol(extract_attr(gtf_dt$V9, "gene_id")),
    gene_symbol = clean_symbol(extract_attr(gtf_dt$V9, "gene_name"))
  )[ensembl_gene_id != "" & gene_symbol != ""])
  e <- new.env(parent = emptyenv()); load(rdata_path, envir = e)
  datExpr <- get("datExpr", envir = e); datMeta <- as.data.table(get("datMeta", envir = e))
  ex <- as.matrix(datExpr)
  sample_col <- intersect(c("sample_id","Sample.ID","sample","sampleID"), names(datMeta))[1]
  dx_col <- intersect(c("Diagnosis","diagnosis","dx"), names(datMeta))[1]
  sample_ids <- as.character(datMeta[[sample_col]])
  sample_cols <- intersect(colnames(ex), sample_ids)
  if (length(sample_cols) >= 10L) {
    expr <- ex[, sample_cols, drop = FALSE]; feats <- rownames(expr); meta <- datMeta[match(sample_cols, sample_ids)]; expr_samples <- sample_cols
  } else {
    sample_rows <- intersect(rownames(ex), sample_ids); expr <- t(ex[sample_rows, , drop = FALSE]); feats <- rownames(expr); meta <- datMeta[match(colnames(expr), sample_ids)]; expr_samples <- colnames(expr)
  }
  feats <- clean_symbol(feats)
  sym <- if (mean(grepl("^ENSG", feats)) > 0.5) gtf_map$gene_symbol[match(feats, gtf_map$ensembl_gene_id)] else feats
  keep <- !is.na(sym) & sym != ""
  expr <- expr[keep, , drop = FALSE]; rownames(expr) <- sym[keep]; expr <- rowsum(expr, rownames(expr))
  dx_raw <- as.character(meta[[dx_col]])
  diagnosis <- fifelse(dx_raw %in% c("ASD","Autism"), "ASD", fifelse(dx_raw %in% c("Control","CTL","CTRL"), "Control", NA_character_))
  keep2 <- !is.na(diagnosis)
  list(expr = expr[, keep2, drop = FALSE], meta = data.table(sample = expr_samples[keep2], diagnosis = diagnosis[keep2]))
}
cohorts <- list(
  GSE102741 = load_gse102741(args$base),
  GSE64018 = load_gse64018(args$base),
  Gandal2022 = load_gandal(args$gandal_rdata)
)

res_list <- list()
for (nm in names(cohorts)) {
  expr <- cohorts[[nm]]$expr
  meta <- cohorts[[nm]]$meta
  frac <- fraction_dt[[nm]]
  frac_covars <- setdiff(names(frac), c("sample","cohort","END"))
  model_dt <- merge(meta, frac, by = "sample")
  for (gs in names(sets)) {
    genes <- intersect(sets[[gs]], rownames(expr))
    if (length(genes) < 5) next
    tmp <- copy(model_dt)
    tmp[, score := colMeans(expr[genes, tmp$sample, drop = FALSE], na.rm = TRUE)]
    u <- fit_ci(tmp, gs, nm, "unadjusted", fraction_covars = character())
    a <- fit_ci(tmp, gs, nm, "nnls_adjusted", fraction_covars = frac_covars)
    if (!is.null(u)) res_list[[paste(nm, gs, "u", sep = "__")]] <- u
    if (!is.null(a)) res_list[[paste(nm, gs, "a", sep = "__")]] <- a
  }
}
res <- rbindlist(res_list, fill = TRUE)
if (nrow(res) == 0) stop("No synaptic/nonsynaptic model results generated")

wide <- dcast(
  res[, .(cohort, program_name, model_type, beta_asd, se, ci_low, ci_high, p_value, mean_asd, mean_control)],
  cohort + program_name ~ model_type,
  value.var = c("beta_asd","se","ci_low","ci_high","p_value","mean_asd","mean_control")
)
wide[, attenuation_adjusted_minus_unadjusted := beta_asd_nnls_adjusted - beta_asd_unadjusted]
wide[, direction_retained := sign(beta_asd_nnls_adjusted) == sign(beta_asd_unadjusted)]
safe_fwrite(res, file.path(args$outdir, "tables", "MA_round1C_synaptic_nonsynaptic_long.tsv"))
safe_fwrite(wide, file.path(args$outdir, "tables", "MA_round1C_synaptic_nonsynaptic_wide.tsv"))
safe_fwrite(data.table(item = c("syngo_gene_file_used", "syngo_gene_col", "n_syngo_genes_total"),
                       value = c(syn_file, args$syngo_gene_col, length(syn_genes))),
            file.path(args$outdir, "meta", "MA_round1C_run_summary.tsv"))

log_msg("Completed Round1C synaptic vs nonsynaptic sensitivity.")
