#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(openxlsx)
})

parse_args <- function(x){
  out <- list(
    base = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare",
    gandal_rdata = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData",
    syngo_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/inputs/SynGO_annotations.xlsx",
    syngo_sig_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/Step12A_SynGO_refinement/tables/11_significant_enrichment_results.tsv",
    reactome_sig_file = NA_character_,
    reactome_term_source = "msigdbr",
    top_syngo_terms = 8L,
    top_reactome_terms = 8L,
    outdir = NA_character_
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    val <- if (i < length(x)) x[[i+1L]] else NA_character_
    if (key == "--base") out$base <- val
    if (key == "--gandal_rdata") out$gandal_rdata <- val
    if (key == "--syngo_file") out$syngo_file <- val
    if (key == "--syngo_sig_file") out$syngo_sig_file <- val
    if (key == "--reactome_sig_file") out$reactome_sig_file <- val
    if (key == "--reactome_term_source") out$reactome_term_source <- val
    if (key == "--top_syngo_terms") out$top_syngo_terms <- as.integer(val)
    if (key == "--top_reactome_terms") out$top_reactome_terms <- as.integer(val)
    if (key == "--outdir") out$outdir <- val
    i <- i + 2L
  }
  if (is.na(out$outdir)) out$outdir <- file.path(out$base, "reviewer_round4", "MA_round2B_camera_roast")
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

# auto-find Reactome significant results file only if not explicitly provided
if (is.na(args$reactome_sig_file)) {
  cand <- list.files(args$base, recursive = TRUE, full.names = TRUE)
  cand <- cand[
    grepl("reactome", basename(cand), ignore.case = TRUE) &
    grepl("significant", basename(cand), ignore.case = TRUE) &
    grepl("\\.tsv(\\.gz)?$", basename(cand), ignore.case = TRUE)
  ]
  if (length(cand)) {
    args$reactome_sig_file <- cand[1]
  } else {
    stop("Could not auto-detect reactome_sig_file under base=", args$base, ". Please pass --reactome_sig_file explicitly.")
  }
}

dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "meta"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "MA_round2B_camera_roast.log")

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
norm_term <- function(x) {
  y <- toupper(as.character(x))
  y <- gsub("[^A-Z0-9]+", "_", y)
  y <- gsub("_+", "_", y)
  y
}

# program sets
prog <- fread(file.path(args$base, "tables", "20_midPrenatal_core_programs.tsv"))
prog[, gene_symbol := clean_symbol(gene_symbol)]
prog <- unique(prog[!is.na(program_name) & !is.na(gene_symbol), .(program_name, gene_symbol)])
program_sets <- split(prog$gene_symbol, prog$program_name)
program_sets <- program_sets[c("SFARI_all","midPrenatal_SFARI_top20")]
program_sets <- lapply(program_sets, unique)

# SynGO term sets
syn_annot <- read_any(args$syngo_file)
if (!all(c("hgnc_symbol","go_name") %in% names(syn_annot))) stop("SynGO file must contain hgnc_symbol and go_name")
syn_annot[, gene_symbol := clean_symbol(hgnc_symbol)]
syn_annot[, term_name := as.character(go_name)]
syn_annot <- unique(syn_annot[!is.na(gene_symbol) & !is.na(term_name), .(term_name, gene_symbol)])

syngo_sig <- read_any(args$syngo_sig_file)
if (!all(c("program","term_name") %in% names(syngo_sig))) stop("SynGO significant file must contain program and term_name")
syngo_top <- unique(syngo_sig[program == "SFARI_all"][order(q_value, -fold_enrichment), term_name])
syngo_top <- syngo_top[seq_len(min(args$top_syngo_terms, length(syngo_top)))]
syngo_sets <- lapply(syngo_top, function(tt) unique(syn_annot[term_name == tt, gene_symbol]))
names(syngo_sets) <- paste0("SYNGO__", syngo_top)

# Reactome term sets via msigdbr
if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Package msigdbr is required for Round2B but is not installed.")
reactome_raw <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
)
rn_gene <- if ("gene_symbol" %in% names(reactome_raw)) "gene_symbol" else if ("human_gene_symbol" %in% names(reactome_raw)) "human_gene_symbol" else stop("Could not identify gene symbol col in msigdbr reactome object")
rn_term <- if ("gs_name" %in% names(reactome_raw)) "gs_name" else stop("Could not identify gs_name in msigdbr reactome object")
reactome_raw <- as.data.table(reactome_raw)[, .(term_name = as.character(get(rn_term)), gene_symbol = clean_symbol(get(rn_gene)))]
reactome_raw <- unique(reactome_raw[!is.na(gene_symbol) & !is.na(term_name)])

reactome_sig <- read_any(args$reactome_sig_file)
term_col <- if ("term_name" %in% names(reactome_sig)) "term_name" else names(reactome_sig)[grepl("term", names(reactome_sig), ignore.case = TRUE)][1]
prog_col <- if ("program" %in% names(reactome_sig)) "program" else names(reactome_sig)[grepl("program", names(reactome_sig), ignore.case = TRUE)][1]
reactome_sig[, term_name_raw := as.character(get(term_col))]
reactome_sig[, program_raw := as.character(get(prog_col))]
reactome_top_raw <- unique(reactome_sig[program_raw == "midPrenatal_SFARI_top20"][order(q_value, -fold_enrichment), term_name_raw])
reactome_top_raw <- reactome_top_raw[seq_len(min(args$top_reactome_terms, length(reactome_top_raw)))]

msig_terms <- unique(reactome_raw$term_name)
msig_norm <- norm_term(msig_terms)
reactome_match <- lapply(reactome_top_raw, function(tt) {
  key <- norm_term(tt)
  key2 <- norm_term(paste0("REACTOME_", tt))
  hits <- msig_terms[msig_norm %in% c(key, key2)]
  if (!length(hits)) hits <- msig_terms[grepl(key, msig_norm, fixed = TRUE)]
  hits[1]
})
reactome_match <- unlist(reactome_match, use.names = FALSE)
names(reactome_match) <- reactome_top_raw
reactome_match <- reactome_match[!is.na(reactome_match) & reactome_match != ""]
reactome_sets <- lapply(reactome_match, function(tt) unique(reactome_raw[term_name == tt, gene_symbol]))
names(reactome_sets) <- paste0("REACTOME__", names(reactome_match))

gene_sets_all <- c(
  lapply(program_sets, unique),
  lapply(syngo_sets, unique),
  lapply(reactome_sets, unique)
)
gene_sets_all <- gene_sets_all[sapply(gene_sets_all, length) >= 5]

set_classes <- c(
  rep("program", length(program_sets)),
  rep("SynGO_term", length(syngo_sets)),
  rep("Reactome_term", length(reactome_sets))
)
set_classes <- set_classes[seq_along(gene_sets_all)]

safe_fwrite(data.table(
  set_name = names(gene_sets_all),
  n_genes = sapply(gene_sets_all, length),
  set_class = set_classes
), file.path(args$outdir, "tables", "MA_round2B_gene_set_inventory.tsv"))

# bulk loaders
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
  expr <- expr[keep, , drop = FALSE]
  rownames(expr) <- sym[keep]
  expr <- rowsum(expr, rownames(expr))
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

camera_res <- list()
roast_res <- list()

for (nm in names(cohorts)) {
  log_msg("Running camera/mroast for cohort ", nm)
  expr <- cohorts[[nm]]$expr
  meta <- cohorts[[nm]]$meta
  common_genes <- intersect(rownames(expr), unique(unlist(gene_sets_all)))
  expr2 <- expr[common_genes, , drop = FALSE]
  design <- model.matrix(~ factor(meta$diagnosis, levels = c("Control","ASD")))
  colnames(design) <- c("Intercept","ASDvsControl")
  index_list <- lapply(gene_sets_all, function(gs) which(rownames(expr2) %in% gs))
  index_list <- index_list[sapply(index_list, length) >= 5]
  cam <- camera(expr2, index = index_list, design = design, contrast = 2, sort = FALSE)
  cam_dt <- as.data.table(cam, keep.rownames = "set_name")
  cam_dt[, cohort := nm]
  cam_dt[, test_method := "camera"]
  camera_res[[nm]] <- cam_dt

  roast <- mroast(expr2, index = index_list, design = design, contrast = 2, nrot = 999)
  roast_dt <- as.data.table(roast, keep.rownames = "set_name")
  roast_dt[, cohort := nm]
  roast_dt[, test_method := "mroast"]
  roast_res[[nm]] <- roast_dt
}

camera_dt <- rbindlist(camera_res, fill = TRUE)
roast_dt <- rbindlist(roast_res, fill = TRUE)
safe_fwrite(camera_dt, file.path(args$outdir, "tables", "MA_round2B_camera_results.tsv"))
safe_fwrite(roast_dt, file.path(args$outdir, "tables", "MA_round2B_mroast_results.tsv"))

summarize_dir <- function(x) {
  xx <- as.character(x)
  fifelse(grepl("Down", xx, ignore.case = TRUE), "Down",
          fifelse(grepl("Up", xx, ignore.case = TRUE), "Up", xx))
}
cam_key <- camera_dt[, .(
  key_p = min(PValue, na.rm = TRUE),
  key_fdr = min(FDR, na.rm = TRUE),
  direction = summarize_dir(Direction[1])
), by = .(cohort, set_name)]
roast_key <- roast_dt[, .(
  key_p = min(PValue, na.rm = TRUE),
  key_fdr = min(FDR, na.rm = TRUE),
  direction = summarize_dir(Direction[1])
), by = .(cohort, set_name)]

summary_dt <- merge(cam_key, roast_key, by = c("cohort","set_name"), suffixes = c("_camera","_mroast"), all = TRUE)
safe_fwrite(summary_dt, file.path(args$outdir, "tables", "MA_round2B_camera_mroast_summary.tsv"))

run_summary <- data.table(
  item = c("syngo_sig_file", "reactome_sig_file", "n_gene_sets_tested", "n_camera_rows", "n_mroast_rows"),
  value = c(args$syngo_sig_file, args$reactome_sig_file, length(gene_sets_all), nrow(camera_dt), nrow(roast_dt))
)
safe_fwrite(run_summary, file.path(args$outdir, "meta", "MA_round2B_run_summary.tsv"))
log_msg("Completed Round2B camera/mroast sensitivity.")
