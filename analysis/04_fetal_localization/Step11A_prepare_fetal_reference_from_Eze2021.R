#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
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

`%||%` <- function(x, y) if (!is.null(x) && !is.na(x) && nzchar(x)) x else y

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  flush.console()
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) return(NA_character_)
  hit[[1L]]
}

normalize_gene <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.\\d+$", "", x)
  pipe_n <- lengths(regmatches(x, gregexpr("\\|", x, fixed = FALSE)))
  use_second <- grepl("^[A-Za-z0-9_.-]+\\|[A-Za-z0-9_.-]+$", x)
  x[use_second] <- sub("^[^|]*\\|", "", x[use_second])
  toupper(x)
}

detect_gene_col <- function(dt) {
  nms <- names(dt)
  low <- tolower(nms)
  cand <- c("gene", "gene_symbol", "symbol", "genes", "hgnc_symbol", "feature", "gene_name")
  idx <- match(cand, low)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(nms[idx[[1L]]])
  nms[[1L]]
}

detect_program_col <- function(dt) {
  nms <- names(dt)
  low <- tolower(nms)
  cand <- c("program", "set", "program_name", "gene_set", "module", "signature")
  idx <- match(cand, low)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(nms[idx[[1L]]])
  NULL
}

detect_cell_id_col <- function(dt) {
  nms <- names(dt)
  low <- tolower(nms)
  cand <- c("cell", "cell_id", "barcode", "cellid", "name")
  idx <- match(cand, low)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) return(nms[idx[[1L]]])
  nms[[1L]]
}

find_meta_like_col <- function(nms, patterns) {
  low <- tolower(nms)
  hits <- which(Reduce(`|`, lapply(patterns, function(p) grepl(p, low, perl = TRUE))))
  if (length(hits) == 0L) return(NA_character_)
  nms[[hits[[1L]]]]
}

standardize_meta <- function(meta_dt) {
  out <- copy(meta_dt)
  id_col <- detect_cell_id_col(out)
  setnames(out, id_col, "cell_id")
  recode <- list(
    cell_type = c("celltype", "cell_type", "annotation", "cluster", "class_label", "cell type", "cell_type_name"),
    stage = c("stage", "age", "pcw", "gw", "carnegie", "development"),
    donor = c("donor", "sample", "sample_id", "individual", "specimen", "donor_id"),
    region = c("region", "anatomical", "area", "tissue", "structure", "brain_region")
  )
  for (nm in names(recode)) {
    hit <- find_meta_like_col(names(out), recode[[nm]])
    if (!is.na(hit) && hit != nm && !(nm %in% names(out))) {
      setnames(out, hit, nm)
    }
  }
  out
}

read_program_genes <- function(program_file) {
  dt <- fread(program_file)
  gene_col <- detect_gene_col(dt)
  program_col <- detect_program_col(dt)
  keep_cols <- c(gene_col, program_col)
  keep_cols <- keep_cols[!is.null(keep_cols)]
  out <- unique(dt[, ..keep_cols])
  setnames(out, gene_col, "gene_raw")
  if (!is.null(program_col)) setnames(out, program_col, "program")
  if (!("program" %in% names(out))) out[, program := "program_union"]
  out[, gene_norm := normalize_gene(gene_raw)]
  out <- out[nzchar(gene_norm)]
  out
}

read_header_cells <- function(counts_file) {
  con <- gzfile(counts_file, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L)
  if (length(header) != 1L) stop("Failed to read header from counts file: ", counts_file, call. = FALSE)
  parts <- strsplit(header, "\t", fixed = TRUE)[[1L]]
  if (length(parts) < 2L) stop("Counts header has <2 columns: ", counts_file, call. = FALSE)
  list(gene_col = parts[[1L]], cell_ids = parts[-1L])
}

extract_target_sparse <- function(counts_file, target_genes, verbose_every = 2000L) {
  con <- gzfile(counts_file, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L)
  header_parts <- strsplit(header, "\t", fixed = TRUE)[[1L]]
  cell_ids <- header_parts[-1L]
  n_cells <- length(cell_ids)
  log_msg("Counts header parsed: ", format(n_cells, big.mark = ","), " cells")

  target_lookup <- unique(target_genes)
  target_lookup <- target_lookup[nzchar(target_lookup)]

  i_list <- list()
  j_list <- list()
  x_list <- list()
  kept_raw <- character()
  kept_norm <- character()
  line_no <- 0L

  repeat {
    lines <- readLines(con, n = 500L)
    if (length(lines) == 0L) break
    for (line in lines) {
      line_no <- line_no + 1L
      tab_pos <- regexpr("\t", line, fixed = TRUE)[[1L]]
      if (tab_pos < 1L) next
      gene_raw <- substr(line, 1L, tab_pos - 1L)
      gene_norm <- normalize_gene(gene_raw)
      if (!(gene_norm %in% target_lookup)) next

      vals_chr <- strsplit(substr(line, tab_pos + 1L, nchar(line)), "\t", fixed = TRUE)[[1L]]
      if (length(vals_chr) != n_cells) {
        stop(
          "Target gene line has ", length(vals_chr), " cells but header has ", n_cells,
          "; gene=", gene_raw, call. = FALSE
        )
      }
      vals_int <- suppressWarnings(as.integer(vals_chr))
      bad <- is.na(vals_int)
      if (any(bad)) {
        vals_num <- suppressWarnings(as.numeric(vals_chr[bad]))
        vals_num[is.na(vals_num)] <- 0
        vals_int[bad] <- as.integer(round(vals_num))
      }
      nz <- which(vals_int != 0L)
      row_idx <- length(kept_norm) + 1L
      kept_raw[[row_idx]] <- gene_raw
      kept_norm[[row_idx]] <- gene_norm
      if (length(nz) > 0L) {
        i_list[[row_idx]] <- rep.int(row_idx, length(nz))
        j_list[[row_idx]] <- nz
        x_list[[row_idx]] <- vals_int[nz]
      } else {
        i_list[[row_idx]] <- integer()
        j_list[[row_idx]] <- integer()
        x_list[[row_idx]] <- integer()
      }
    }
    if (line_no %% verbose_every < length(lines)) {
      log_msg("Scanned ", format(line_no, big.mark = ","), " genes; matched ", length(kept_norm), " target rows so far")
    }
  }

  if (length(kept_norm) == 0L) {
    stop("No target genes were recovered from counts file.", call. = FALSE)
  }

  sp <- sparseMatrix(
    i = unlist(i_list, use.names = FALSE),
    j = unlist(j_list, use.names = FALSE),
    x = unlist(x_list, use.names = FALSE),
    dims = c(length(kept_norm), n_cells),
    dimnames = list(kept_norm, cell_ids)
  )

  list(
    counts = sp,
    gene_raw = kept_raw,
    gene_norm = kept_norm,
    cell_ids = cell_ids,
    n_scanned_lines = line_no
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
base_dir <- args$base_dir %||% stop("--base_dir is required", call. = FALSE)
input_dir <- args$input_dir %||% file.path(base_dir, "inputs")
program_file <- args$program_file %||% stop("--program_file is required", call. = FALSE)
output_dir <- args$output_dir %||% file.path(base_dir, "Step11A_prepare_fetal_reference_from_Eze2021")
output_rds <- args$output_rds %||% file.path(input_dir, "fetal_singlecell_reference.rds")
counts_file <- args$counts_file %||% first_existing(file.path(input_dir, c("fetal_counts_exprMatrix.tsv.gz", "counts_exprMatrix.tsv.gz")))
meta_file <- args$meta_file %||% first_existing(file.path(input_dir, c("fetal_meta.tsv", "meta.tsv")))
umap_file <- args$umap_file %||% first_existing(file.path(input_dir, c("fetal_umap.coords.tsv.gz", "UMAP.coords.tsv.gz", "umap.coords.tsv.gz")))

if (is.na(counts_file) || !file.exists(counts_file)) stop("Counts file not found. Provide --counts_file", call. = FALSE)
if (is.na(meta_file) || !file.exists(meta_file)) stop("Meta file not found. Provide --meta_file", call. = FALSE)
if (is.na(umap_file) || !file.exists(umap_file)) stop("UMAP file not found. Provide --umap_file", call. = FALSE)
if (!file.exists(program_file)) stop("Program file not found: ", program_file, call. = FALSE)

ensure_dir(output_dir)
ensure_dir(dirname(output_rds))

log_msg("Reading program file: ", program_file)
program_dt <- read_program_genes(program_file)
program_gene_union <- unique(program_dt$gene_norm)
log_msg("Program gene union size: ", length(program_gene_union))

log_msg("Reading metadata: ", meta_file)
meta_dt <- fread(meta_file)
meta_dt <- standardize_meta(meta_dt)
if (anyDuplicated(meta_dt$cell_id)) {
  stop("Metadata cell_id contains duplicates.", call. = FALSE)
}
log_msg("Metadata rows: ", format(nrow(meta_dt), big.mark = ","))

log_msg("Reading UMAP coordinates: ", umap_file)
umap_dt <- fread(umap_file)
umap_id_col <- detect_cell_id_col(umap_dt)
setnames(umap_dt, umap_id_col, "cell_id")
coord_cols <- setdiff(names(umap_dt), "cell_id")
if (length(coord_cols) >= 2L) {
  setnames(umap_dt, coord_cols[1:2], c("UMAP_1", "UMAP_2"))
}
if (anyDuplicated(umap_dt$cell_id)) {
  umap_dt <- unique(umap_dt, by = "cell_id")
}
log_msg("UMAP rows: ", format(nrow(umap_dt), big.mark = ","))

header_info <- read_header_cells(counts_file)
count_cells <- header_info$cell_ids
log_msg("Counts cell columns from header: ", format(length(count_cells), big.mark = ","))

meta_match <- match(count_cells, meta_dt$cell_id)
if (anyNA(meta_match)) {
  stop(
    "Metadata is missing ", sum(is.na(meta_match)), " cell IDs present in counts header.",
    call. = FALSE
  )
}
meta_aligned <- meta_dt[meta_match]
stopifnot(identical(meta_aligned$cell_id, count_cells))
setDF(meta_aligned)
rownames(meta_aligned) <- meta_aligned$cell_id

umap_match <- match(count_cells, umap_dt$cell_id)
umap_aligned <- data.frame(cell_id = count_cells, UMAP_1 = NA_real_, UMAP_2 = NA_real_)
matched_umap <- !is.na(umap_match)
if (any(matched_umap)) {
  umap_aligned$UMAP_1[matched_umap] <- umap_dt$UMAP_1[umap_match[matched_umap]]
  umap_aligned$UMAP_2[matched_umap] <- umap_dt$UMAP_2[umap_match[matched_umap]]
}
rownames(umap_aligned) <- umap_aligned$cell_id

log_msg("Extracting only target program genes from counts matrix: ", counts_file)
extract_res <- extract_target_sparse(counts_file, target_genes = program_gene_union)
counts_sp <- extract_res$counts

kept_gene_dt <- unique(data.table(
  gene_raw_in_counts = extract_res$gene_raw,
  gene_norm = extract_res$gene_norm
))
program_map_dt <- unique(program_dt[, .(program, gene_raw, gene_norm)])
program_overlap_dt <- merge(program_map_dt, kept_gene_dt, by = "gene_norm", all.x = TRUE)
program_overlap_dt[, present_in_counts := !is.na(gene_raw_in_counts)]

missing_genes <- setdiff(program_gene_union, unique(extract_res$gene_norm))
log_msg(
  "Recovered ", nrow(kept_gene_dt), " / ", length(program_gene_union),
  " target genes from counts; missing=", length(missing_genes)
)

ref_obj <- list(
  dataset_id = "Eze2021_early_brain_UCSC_program_focused_reference",
  source = list(
    study = "Eze et al. Nat Neurosci. 2021",
    description = "Single-cell atlas of early human brain development highlights heterogeneity of human neuroepithelial cells and early radial glia",
    counts_file = normalizePath(counts_file),
    meta_file = normalizePath(meta_file),
    umap_file = normalizePath(umap_file),
    program_file = normalizePath(program_file),
    created_at = as.character(Sys.time())
  ),
  feature_space = list(
    type = "program_gene_union_only",
    n_program_genes_requested = length(program_gene_union),
    n_program_genes_recovered = nrow(kept_gene_dt),
    missing_genes = missing_genes,
    gene_map = kept_gene_dt,
    program_map = program_overlap_dt
  ),
  counts = counts_sp,
  meta = meta_aligned,
  umap = umap_aligned,
  cell_ids = count_cells
)

saveRDS(ref_obj, output_rds, compress = "xz")
log_msg("Saved fetal reference RDS: ", output_rds)

fwrite(program_overlap_dt, file.path(output_dir, "01_program_gene_overlap.tsv"), sep = "\t")
fwrite(kept_gene_dt, file.path(output_dir, "02_recovered_gene_map.tsv"), sep = "\t")
fwrite(
  data.table(
    metric = c(
      "n_cells",
      "n_meta_rows",
      "n_umap_rows",
      "n_target_genes_requested",
      "n_target_genes_recovered",
      "n_missing_genes",
      "counts_nonzero_entries"
    ),
    value = c(
      length(count_cells),
      nrow(meta_aligned),
      nrow(umap_aligned),
      length(program_gene_union),
      nrow(kept_gene_dt),
      length(missing_genes),
      length(counts_sp@x)
    )
  ),
  file.path(output_dir, "03_reference_summary.tsv"),
  sep = "\t"
)

md_lines <- c(
  "# Step11A fetal reference preparation summary",
  "",
  paste0("- Study: Eze et al. Nat Neurosci. 2021 early human brain dataset"),
  paste0("- Counts file: `", counts_file, "`"),
  paste0("- Meta file: `", meta_file, "`"),
  paste0("- UMAP file: `", umap_file, "`"),
  paste0("- Program file: `", program_file, "`"),
  paste0("- Output RDS: `", output_rds, "`"),
  "",
  "## Summary",
  paste0("- Cells in counts header: ", format(length(count_cells), big.mark = ",")),
  paste0("- Program gene union requested: ", length(program_gene_union)),
  paste0("- Program genes recovered: ", nrow(kept_gene_dt)),
  paste0("- Missing genes: ", length(missing_genes)),
  paste0("- Non-zero matrix entries stored: ", format(length(counts_sp@x), big.mark = ",")),
  "",
  "## Notes",
  "- This RDS is program-focused, not a whole-transcriptome object.",
  "- Counts are stored as a sparse matrix for the union of genes present in the supplied program table.",
  "- Metadata and UMAP are aligned to the counts header cell order."
)
writeLines(md_lines, con = file.path(output_dir, "04_prepare_fetal_reference_summary.md"))

log_msg("Done.")
