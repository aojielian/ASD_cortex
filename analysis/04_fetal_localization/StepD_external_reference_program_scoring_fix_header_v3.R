#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(argv) {
  out <- list(
    base_dir = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex",
    stepc_dir = NULL,
    outdir = NULL,
    program_file = "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex/step01_prepare/tables/20_midPrenatal_core_programs.tsv",
    sfari_all = NULL,
    top20 = NULL,
    top10 = NULL,
    top05 = NULL,
    gene_map_tsv = NULL,
    gtf = "/gpfs/hpc/home/lijc/lianaoj/reference/gencode.v49.basic.annotation.gtf.gz",
    use_rna_singlets_only = TRUE,
    use_multiome_singlets_only = TRUE,
    min_overlap = 3L,
    norm_scale_factor = 10000
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    val <- if (i < length(argv)) argv[[i + 1L]] else NA_character_
    if (key == "--base_dir") {
      out$base_dir <- val; i <- i + 2L
    } else if (key == "--stepc_dir") {
      out$stepc_dir <- val; i <- i + 2L
    } else if (key == "--outdir") {
      out$outdir <- val; i <- i + 2L
    } else if (key == "--program_file") {
      out$program_file <- val; i <- i + 2L
    } else if (key == "--sfari_all") {
      out$sfari_all <- val; i <- i + 2L
    } else if (key == "--top20") {
      out$top20 <- val; i <- i + 2L
    } else if (key == "--top10") {
      out$top10 <- val; i <- i + 2L
    } else if (key == "--top05") {
      out$top05 <- val; i <- i + 2L
    } else if (key == "--gene_map_tsv") {
      out$gene_map_tsv <- val; i <- i + 2L
    } else if (key == "--gtf") {
      out$gtf <- val; i <- i + 2L
    } else if (key == "--use_rna_singlets_only") {
      out$use_rna_singlets_only <- tolower(val) %in% c("true","t","1","yes","y"); i <- i + 2L
    } else if (key == "--use_multiome_singlets_only") {
      out$use_multiome_singlets_only <- tolower(val) %in% c("true","t","1","yes","y"); i <- i + 2L
    } else if (key == "--min_overlap") {
      out$min_overlap <- as.integer(val); i <- i + 2L
    } else if (key == "--norm_scale_factor") {
      out$norm_scale_factor <- as.numeric(val); i <- i + 2L
    } else if (key %in% c("-h", "--help")) {
      cat(
"Usage:\n",
"  Rscript StepD_external_reference_program_scoring.R [options]\n\n",
"Preferred gene-set input:\n",
"  --program_file <path>   TSV with columns program_name and gene_symbol\n",
"                         containing SFARI_all, midPrenatal_SFARI_top20,\n",
"                         midPrenatal_SFARI_top10, and midPrenatal_SFARI_top05.\n\n",
"Legacy alternative gene-set inputs:\n",
"  --sfari_all <path>\n",
"  --top20 <path>\n",
"  --top10 <path>\n",
"  --top05 <path>\n\n",
"Gene mapping options (needed when gene sets are SYMBOLS for GSE162170 RNA/multiome):\n",
"  --gene_map_tsv <path>   Two-column or named-column TSV with gene_id / gene_symbol\n",
"  --gtf <path>            GTF(.gz) from which gene_id / gene_name will be parsed\n\n",
"Other options:\n",
"  --base_dir <path>                   Base ASD_cortex directory\n",
"  --program_file <path>               Combined program file (default: step01_prepare/tables/20_midPrenatal_core_programs.tsv)\n",
"  --stepc_dir <path>                  StepC output directory (default: <base_dir>/stepC_external_reference_prep)\n",
"  --outdir <path>                     Output directory (default: <base_dir>/stepD_external_reference_program_scoring)\n",
"  --use_rna_singlets_only <T/F>       Restrict primary RNA summaries to Singlet cells (default: TRUE)\n",
"  --use_multiome_singlets_only <T/F>  Restrict primary multiome summaries to Singlet cells (default: TRUE)\n",
"  --min_overlap <int>                 Minimum mapped genes required to score a program in a dataset (default: 3)\n",
"  --norm_scale_factor <num>           Scale factor for counts-like log normalization (default: 10000)\n",
sep = "")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
  }
  if (is.null(out$stepc_dir) || is.na(out$stepc_dir) || out$stepc_dir == "") {
    out$stepc_dir <- file.path(out$base_dir, "stepC_external_reference_prep")
  }
  if (is.null(out$outdir) || is.na(out$outdir) || out$outdir == "") {
    out$outdir <- file.path(out$base_dir, "stepD_external_reference_program_scoring")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$outdir, "logs", "StepD_external_reference_program_scoring.log")
cat("", file = log_file)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_dir <- function(x) dir.create(x, recursive = TRUE, showWarnings = FALSE)

safe_fwrite <- function(x, path, sep = "\t") {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    fwrite(x, file = path, sep = sep, quote = FALSE, na = "NA", compress = "gzip")
  } else {
    fwrite(x, file = path, sep = sep, quote = FALSE, na = "NA")
  }
}

run_cmd <- function(cmd) {
  log_msg("CMD: ", cmd)
  out <- tryCatch(system(cmd, intern = TRUE), warning = function(w) system(cmd, intern = TRUE), error = function(e) e)
  if (inherits(out, "error")) stop(out$message)
  out
}

read_gz_or_plain <- function(path, header = TRUE) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    fread(cmd = sprintf("zcat %s", shQuote(path)), sep = "\t", header = header,
          data.table = FALSE, check.names = FALSE, fill = TRUE)
  } else {
    fread(path, sep = "\t", header = header, data.table = FALSE, check.names = FALSE, fill = TRUE)
  }
}

normalize_token <- function(x) {
  x <- trimws(x)
  x <- x[nzchar(x)]
  x
}

split_tab_line <- function(x) {
  strsplit(x, "	", fixed = TRUE)[[1]]
}

extract_matrix_header_cells <- function(matrix_path) {
  header_line <- run_cmd(sprintf("zcat %s | head -1", shQuote(matrix_path)))[1]
  second_line <- run_cmd(sprintf("zcat %s | head -2 | tail -1", shQuote(matrix_path)))[1]
  header_fields <- split_tab_line(header_line)
  second_fields <- split_tab_line(second_line)

  looks_like_gene_field <- nzchar(second_fields[1]) && (
    grepl("^ENSG", second_fields[1], ignore.case = TRUE) ||
    grepl("|", second_fields[1], fixed = TRUE)
  )

  has_blank_gene_placeholder <- length(header_fields) == length(second_fields) &&
    looks_like_gene_field &&
    !nzchar(header_fields[1])

  if (has_blank_gene_placeholder) {
    cell_ids <- header_fields[-1]
  } else if (length(header_fields) == (length(second_fields) - 1L)) {
    cell_ids <- header_fields
  } else {
    cell_ids <- if (length(header_fields) > 1L && grepl("_", header_fields[1])) header_fields else header_fields[-1]
  }

  list(
    cell_ids = cell_ids,
    n_cells = length(cell_ids),
    has_blank_gene_placeholder = has_blank_gene_placeholder
  )
}

load_gene_set <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- gsub("\r", "", x)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!grepl("^#", x)]
  x <- x[!tolower(x) %in% c("gene", "genes", "symbol", "gene_symbol", "gene_id")]
  unique(x)
}

load_program_table <- function(path, wanted = c("SFARI_all", "midPrenatal_SFARI_top20", "midPrenatal_SFARI_top10", "midPrenatal_SFARI_top05")) {
  dt <- read_gz_or_plain(path, header = TRUE)
  nms <- tolower(names(dt))
  pcol <- names(dt)[match("program_name", nms)]
  gcol <- if ("gene_symbol" %in% nms) {
    names(dt)[match("gene_symbol", nms)]
  } else if ("symbol" %in% nms) {
    names(dt)[match("symbol", nms)]
  } else if ("gene" %in% nms) {
    names(dt)[match("gene", nms)]
  } else if ("standardized_gene_symbol" %in% nms) {
    names(dt)[match("standardized_gene_symbol", nms)]
  } else {
    stop("Program file must contain program_name plus a gene symbol column.")
  }
  if (is.na(pcol) || is.na(gcol)) stop("Could not resolve program_name/gene symbol columns from program file: ", path)
  dt <- dt[, c(pcol, gcol), drop = FALSE]
  names(dt) <- c("program_name", "gene_symbol")
  dt$program_name <- trimws(as.character(dt$program_name))
  dt$gene_symbol <- toupper(trimws(as.character(dt$gene_symbol)))
  dt <- dt[nzchar(dt$program_name) & nzchar(dt$gene_symbol), , drop = FALSE]
  dt <- unique(dt)
  missing_programs <- setdiff(wanted, unique(dt$program_name))
  if (length(missing_programs)) {
    stop("Program file missing required programs: ", paste(missing_programs, collapse = ", "))
  }
  out <- split(dt$gene_symbol[dt$program_name %in% wanted], dt$program_name[dt$program_name %in% wanted])
  out <- out[wanted]
  out <- lapply(out, unique)
  out
}

guess_gene_type <- function(genes) {
  if (!length(genes)) return("unknown")
  frac_ensg <- mean(grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", genes, ignore.case = TRUE))
  if (frac_ensg >= 0.8) "ensg" else "symbol"
}

resolve_gene_map <- function(gene_map_tsv = NULL, gtf = NULL, outdir = NULL) {
  if (!is.null(gene_map_tsv) && !is.na(gene_map_tsv) && nzchar(gene_map_tsv)) {
    log_msg("Loading gene map TSV: ", gene_map_tsv)
    gm <- read_gz_or_plain(gene_map_tsv, header = TRUE)
    nms <- tolower(gsub("[^A-Za-z0-9]+", "_", names(gm)))
    names(gm) <- nms
    gid_col <- if ("gene_id" %in% nms) "gene_id" else names(gm)[1]
    gsym_col <- if ("gene_symbol" %in% nms) "gene_symbol" else if ("gene_name" %in% nms) "gene_name" else names(gm)[2]
    gm <- gm[, c(gid_col, gsym_col), drop = FALSE]
    names(gm) <- c("gene_id", "gene_symbol")
    gm$gene_id <- sub("\\..*$", "", trimws(gm$gene_id))
    gm$gene_symbol <- toupper(trimws(gm$gene_symbol))
    gm <- gm[nzchar(gm$gene_id) & nzchar(gm$gene_symbol), , drop = FALSE]
    gm <- unique(gm)
    if (!is.null(outdir)) safe_fwrite(gm, file.path(outdir, "resolved_gene_map.tsv.gz"))
    return(gm)
  }
  if (!is.null(gtf) && !is.na(gtf) && nzchar(gtf)) {
    log_msg("Parsing gene map from GTF: ", gtf)
    cat_cmd <- if (grepl("\\.gz$", gtf, ignore.case = TRUE)) sprintf("zcat %s", shQuote(gtf)) else sprintf("cat %s", shQuote(gtf))
    cmd <- sprintf("bash -lc %s", shQuote(sprintf("%s | awk -F '\\t' '$3==\"gene\"{print $9}'", cat_cmd)))
    attrs <- run_cmd(cmd)
    gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', attrs)
    gene_symbol <- sub('.*gene_name "([^"]+)".*', '\\1', attrs)
    gm <- data.frame(gene_id = sub("\\..*$", "", gene_id), gene_symbol = toupper(gene_symbol), stringsAsFactors = FALSE)
    gm <- gm[nzchar(gm$gene_id) & nzchar(gm$gene_symbol), , drop = FALSE]
    gm <- unique(gm)
    if (!is.null(outdir)) safe_fwrite(gm, file.path(outdir, "resolved_gene_map.tsv.gz"))
    return(gm)
  }
  NULL
}

read_dense_subset <- function(matrix_path, target_ids, header_info = NULL) {
  target_ids <- unique(target_ids[nzchar(target_ids)])
  if (!length(target_ids)) stop("No target gene IDs provided for dense subset extraction: ", matrix_path)
  if (is.null(header_info)) header_info <- extract_matrix_header_cells(matrix_path)
  header_cells <- header_info$cell_ids
  tf <- tempfile(fileext = ".genes.txt")
  writeLines(target_ids, tf)
  cmd <- sprintf(
    "bash -lc %s",
    shQuote(sprintf(
      "zcat %s | awk -F '\\t' -v GENEFILE=%s 'BEGIN{while((getline < GENEFILE) > 0) keep[$1]=1} NR>1 {gid=$1; sub(/\\..*$/, \"\", gid); if (($1 in keep) || (gid in keep)) print}'",
      shQuote(matrix_path), shQuote(tf)
    ))
  )
  on.exit(unlink(tf), add = TRUE)
  dt <- fread(cmd = cmd, sep = "\t", header = FALSE, data.table = FALSE, check.names = FALSE, fill = TRUE)
  if (!nrow(dt) || !ncol(dt)) stop("Dense subset read produced zero rows/columns: ", matrix_path)
  expected_ncol <- length(header_cells) + 1L
  if (ncol(dt) != expected_ncol) {
    stop("Subset matrix column count mismatch for ", matrix_path,
         ": observed ", ncol(dt), ", expected ", expected_ncol,
         ". This usually indicates header/data parsing drift.")
  }
  colnames(dt) <- c("gene_id", header_cells)
  dt$gene_id <- sub("\\..*$", "", dt$gene_id)
  dt
}


log_transform_counts <- function(mat, libsize, scale_factor = 10000) {
  libsize <- as.numeric(libsize)
  libsize[!is.finite(libsize) | libsize <= 0] <- NA_real_
  if (anyNA(libsize)) stop("Counts-like matrix requires positive library sizes for all cells.")
  norm <- sweep(mat, 2L, libsize / scale_factor, "/")
  log1p(norm)
}

transform_continuous <- function(mat, dataset_name = NULL) {
  if (!is.null(dataset_name) && identical(dataset_name, "Nowakowski_UCSC")) {
    return(log1p(mat))
  }
  rng <- range(mat, finite = TRUE)
  if (is.finite(rng[2]) && rng[2] > 50) {
    log1p(mat)
  } else {
    mat
  }
}

score_from_matrix <- function(expr_mat, program_to_ids) {
  out <- vector("list", length(program_to_ids))
  nm <- names(program_to_ids)
  for (i in seq_along(program_to_ids)) {
    ids <- program_to_ids[[i]]
    ids <- ids[ids %in% rownames(expr_mat)]
    if (!length(ids)) {
      out[[i]] <- rep(NA_real_, ncol(expr_mat))
    } else {
      out[[i]] <- colMeans(expr_mat[ids, , drop = FALSE])
    }
  }
  names(out) <- nm
  as.data.frame(out, check.names = FALSE)
}

summarize_long <- function(cell_df, score_cols, group_var, group_label_var = NULL, only_primary = TRUE) {
  df <- cell_df
  if (only_primary && "use_in_primary_summary" %in% names(df)) {
    df <- df[isTRUE(df$use_in_primary_summary) | (!is.na(df$use_in_primary_summary) & df$use_in_primary_summary), , drop = FALSE]
  }
  if (!nrow(df)) return(data.frame())
  if (!(group_var %in% names(df))) return(data.frame())
  groups <- unique(df[[group_var]])
  groups <- groups[!is.na(groups) & nzchar(as.character(groups))]
  if (!length(groups)) return(data.frame())
  res <- vector("list", length(groups) * length(score_cols))
  idx <- 1L
  for (g in groups) {
    sub <- df[df[[group_var]] == g, , drop = FALSE]
    glab <- if (!is.null(group_label_var) && group_label_var %in% names(sub) && any(!is.na(sub[[group_label_var]]) & nzchar(as.character(sub[[group_label_var]])))) {
      ux <- unique(as.character(sub[[group_label_var]][!is.na(sub[[group_label_var]])]))
      ux[[1]]
    } else {
      as.character(g)
    }
    for (sc in score_cols) {
      vals <- sub[[sc]]
      res[[idx]] <- data.frame(
        dataset = unique(sub$dataset)[1],
        group_var = group_var,
        group_id = as.character(g),
        group_label = glab,
        program = sc,
        n_cells = nrow(sub),
        mean_score = mean(vals, na.rm = TRUE),
        median_score = median(vals, na.rm = TRUE),
        sd_score = sd(vals, na.rm = TRUE),
        min_score = min(vals, na.rm = TRUE),
        max_score = max(vals, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  out <- rbindlist(res[seq_len(idx - 1L)], use.names = TRUE, fill = TRUE)
  out[, within_dataset_rank := frank(-mean_score, ties.method = "min"), by = .(dataset, group_var, program)]
  out[, cluster_mean_z := if (.N > 1) as.numeric(scale(mean_score)) else 0, by = .(dataset, group_var, program)]
  as.data.frame(out)
}

build_top_hits <- function(group_summary, n_top = 5L) {
  if (!nrow(group_summary)) return(data.frame())
  dt <- as.data.table(group_summary)
  dt <- dt[order(dataset, group_var, program, within_dataset_rank)]
  dt <- dt[within_dataset_rank <= n_top]
  as.data.frame(dt)
}

prepare_cell_score_output <- function(meta, score_df, dataset_name, primary_filter = NULL) {
  stopifnot(nrow(meta) == nrow(score_df))
  out <- cbind(meta, score_df)
  out$dataset <- dataset_name
  out$use_in_primary_summary <- if (is.null(primary_filter)) TRUE else primary_filter
  out
}

map_programs_to_nowakowski <- function(gene_sets, gene_type_map, now_obj, gene_map = NULL) {
  gt <- now_obj$gene_table
  if (is.null(gt) || !all(c("gene_raw", "gene_symbol") %in% names(gt))) {
    gt <- data.frame(gene_raw = rownames(now_obj$expr), gene_symbol = sub("\\|.*$", "", rownames(now_obj$expr)), stringsAsFactors = FALSE)
  }
  gt$gene_symbol_upper <- toupper(gt$gene_symbol)
  row_map <- setNames(gt$gene_raw, gt$gene_symbol_upper)
  reverse_map <- NULL
  if (!is.null(gene_map)) {
    reverse_map <- split(gene_map$gene_symbol, gene_map$gene_id)
  }
  out <- list()
  map_stats <- list()
  for (nm in names(gene_sets)) {
    genes <- gene_sets[[nm]]
    gtype <- gene_type_map[[nm]]
    if (identical(gtype, "ensg")) {
      if (is.null(reverse_map)) {
        mapped_symbols <- character()
      } else {
        mapped_symbols <- unique(toupper(unlist(reverse_map[sub("\\..*$", "", genes)], use.names = FALSE)))
      }
      mapped_rows <- unique(unname(row_map[mapped_symbols]))
      input_key_type <- "ensg"
    } else {
      mapped_symbols <- unique(toupper(genes))
      mapped_rows <- unique(unname(row_map[mapped_symbols]))
      input_key_type <- "symbol"
    }
    mapped_rows <- mapped_rows[!is.na(mapped_rows)]
    out[[nm]] <- mapped_rows
    map_stats[[nm]] <- data.frame(dataset = "Nowakowski_UCSC", program = nm,
                                  input_gene_type = input_key_type,
                                  n_input = length(unique(genes)),
                                  n_mapped = length(mapped_rows),
                                  stringsAsFactors = FALSE)
  }
  list(program_rows = out, stats = rbindlist(map_stats, use.names = TRUE, fill = TRUE))
}

map_programs_to_gene_ids <- function(gene_sets, gene_type_map, available_gene_ids, gene_map = NULL, dataset_name = NULL) {
  available_gene_ids <- unique(sub("\\..*$", "", available_gene_ids))
  avail_set <- unique(available_gene_ids)
  symbol_to_ids <- NULL
  if (!is.null(gene_map)) {
    gm <- gene_map[gene_map$gene_id %in% avail_set, , drop = FALSE]
    if (nrow(gm)) symbol_to_ids <- split(gm$gene_id, gm$gene_symbol)
  }
  out <- list()
  map_stats <- list()
  for (nm in names(gene_sets)) {
    genes <- unique(gene_sets[[nm]])
    gtype <- gene_type_map[[nm]]
    mapped_ids <- character()
    if (identical(gtype, "ensg")) {
      mapped_ids <- unique(sub("\\..*$", "", genes))
      mapped_ids <- intersect(mapped_ids, avail_set)
    } else {
      if (is.null(symbol_to_ids)) {
        mapped_ids <- character()
      } else {
        keys <- toupper(genes)
        mapped_ids <- unique(unlist(symbol_to_ids[keys], use.names = FALSE))
        mapped_ids <- intersect(mapped_ids, avail_set)
      }
    }
    out[[nm]] <- mapped_ids
    map_stats[[nm]] <- data.frame(dataset = dataset_name, program = nm,
                                  input_gene_type = gtype,
                                  n_input = length(genes),
                                  n_mapped = length(mapped_ids),
                                  stringsAsFactors = FALSE)
  }
  list(program_rows = out, stats = rbindlist(map_stats, use.names = TRUE, fill = TRUE))
}

# ---- Load gene sets ----
program_paths <- list(
  SFARI_all = args$sfari_all,
  midPrenatal_SFARI_top20 = args$top20,
  midPrenatal_SFARI_top10 = args$top10,
  midPrenatal_SFARI_top05 = args$top05
)
program_paths <- program_paths[!vapply(program_paths, function(x) is.null(x) || is.na(x) || !nzchar(x), logical(1))]

log_msg("Starting StepD external reference program scoring...")
log_msg("base_dir = ", args$base_dir)
log_msg("stepc_dir = ", args$stepc_dir)
log_msg("outdir   = ", args$outdir)

safe_dir(file.path(args$outdir, "Nowakowski"))
safe_dir(file.path(args$outdir, "GSE162170_RNA"))
safe_dir(file.path(args$outdir, "GSE162170_multiome"))

# Load gene sets
programs <- list()
program_gene_types <- list()
program_input_summary <- list()

if (!is.null(args$program_file) && !is.na(args$program_file) && nzchar(args$program_file) && file.exists(args$program_file)) {
  log_msg("Loading combined program file: ", args$program_file)
  programs <- load_program_table(args$program_file)
  for (nm in names(programs)) {
    gs <- programs[[nm]]
    gtype <- guess_gene_type(gs)
    log_msg("Loaded program ", nm, " from combined file: n=", length(gs), ", guessed_type=", gtype)
    program_gene_types[[nm]] <- gtype
    program_input_summary[[nm]] <- data.frame(
      program = nm,
      source_type = "combined_program_file",
      path = args$program_file,
      n_input = length(gs),
      guessed_type = gtype,
      stringsAsFactors = FALSE
    )
    safe_fwrite(data.frame(gene = gs, stringsAsFactors = FALSE), file.path(args$outdir, paste0(nm, "_input_genes.tsv.gz")))
  }
} else {
  if (!length(program_paths)) stop("No valid gene-set input found. Provide --program_file or legacy per-program gene-set files.")
  for (nm in names(program_paths)) {
    gs <- load_gene_set(program_paths[[nm]])
    gtype <- guess_gene_type(gs)
    log_msg("Loaded gene set ", nm, ": n=", length(gs), ", guessed_type=", gtype)
    programs[[nm]] <- gs
    program_gene_types[[nm]] <- gtype
    program_input_summary[[nm]] <- data.frame(
      program = nm,
      source_type = "legacy_gene_set_file",
      path = program_paths[[nm]],
      n_input = length(gs),
      guessed_type = gtype,
      stringsAsFactors = FALSE
    )
    safe_fwrite(data.frame(gene = gs, stringsAsFactors = FALSE), file.path(args$outdir, paste0(nm, "_input_genes.tsv.gz")))
  }
}

safe_fwrite(rbindlist(program_input_summary, use.names = TRUE, fill = TRUE), file.path(args$outdir, "StepD_gene_set_inputs.tsv"))

needs_symbol_map <- any(vapply(program_gene_types, identical, logical(1), y = "symbol"))
gene_map <- resolve_gene_map(args$gene_map_tsv, args$gtf, args$outdir)
if (needs_symbol_map && is.null(gene_map)) {
  stop("At least one gene set appears to use gene symbols, but no --gene_map_tsv or --gtf was provided.")
}

# ---- Load StepC objects ----
now_path <- file.path(args$stepc_dir, "Nowakowski", "Nowakowski_program_ready.rds")
rna_path <- file.path(args$stepc_dir, "GSE162170_RNA", "GSE162170_RNA_program_ready_light.rds")
mo_path  <- file.path(args$stepc_dir, "GSE162170_multiome", "GSE162170_multiome_program_ready_light.rds")
stopifnot(file.exists(now_path), file.exists(rna_path), file.exists(mo_path))

log_msg("Loading StepC program-ready objects...")
now_obj <- readRDS(now_path)
rna_obj <- readRDS(rna_path)
mo_obj  <- readRDS(mo_path)

mapping_summaries <- list()

# ---- Nowakowski ----
log_msg("Scoring Nowakowski...")
now_map <- map_programs_to_nowakowski(programs, program_gene_types, now_obj, gene_map)
mapping_summaries[["Nowakowski_UCSC"]] <- now_map$stats
now_valid <- names(now_map$program_rows)[vapply(now_map$program_rows, length, integer(1)) >= args$min_overlap]
if (!length(now_valid)) stop("No gene sets met min_overlap for Nowakowski.")
now_expr <- now_obj$expr
now_meta <- now_obj$meta
now_subset_genes <- unique(unlist(now_map$program_rows[now_valid], use.names = FALSE))
now_mat <- now_expr[now_subset_genes, , drop = FALSE]
mode(now_mat) <- "numeric"
now_tmat <- transform_continuous(now_mat, dataset_name = "Nowakowski_UCSC")
now_scores <- score_from_matrix(now_tmat, now_map$program_rows[now_valid])
colnames(now_scores) <- paste0("score_", colnames(now_scores))
now_primary <- rep(TRUE, nrow(now_meta))
now_cell_scores <- prepare_cell_score_output(now_meta, now_scores, "Nowakowski_UCSC", primary_filter = now_primary)
safe_fwrite(now_cell_scores, file.path(args$outdir, "Nowakowski", "Nowakowski_cell_program_scores.tsv.gz"))
now_group_summary <- summarize_long(now_cell_scores, score_cols = grep("^score_", names(now_cell_scores), value = TRUE), group_var = "cell_type", group_label_var = "cell_type", only_primary = TRUE)
now_age_summary <- summarize_long(now_cell_scores, score_cols = grep("^score_", names(now_cell_scores), value = TRUE), group_var = "age_label", group_label_var = "age_label", only_primary = TRUE)
safe_fwrite(now_group_summary, file.path(args$outdir, "Nowakowski", "Nowakowski_group_program_summary.tsv.gz"))
safe_fwrite(now_age_summary, file.path(args$outdir, "Nowakowski", "Nowakowski_age_program_summary.tsv.gz"))
saveRDS(list(dataset = "Nowakowski_UCSC", meta = now_meta, score_df = now_scores, valid_programs = now_valid),
        file.path(args$outdir, "Nowakowski", "Nowakowski_scored_object.rds"), compress = "xz")

# ---- GSE162170 RNA ----
log_msg("Scoring GSE162170 RNA...")
rna_gene_ids <- read_gz_or_plain(rna_obj$gene_ids_file, header = TRUE)[[1]]
rna_gene_ids <- sub("\\..*$", "", rna_gene_ids)
rna_map <- map_programs_to_gene_ids(programs, program_gene_types, rna_gene_ids, gene_map, dataset_name = "GSE162170_RNA")
mapping_summaries[["GSE162170_RNA"]] <- rna_map$stats
rna_valid <- names(rna_map$program_rows)[vapply(rna_map$program_rows, length, integer(1)) >= args$min_overlap]
if (!length(rna_valid)) stop("No gene sets met min_overlap for GSE162170 RNA.")
rna_union_ids <- unique(unlist(rna_map$program_rows[rna_valid], use.names = FALSE))
rna_header_info <- extract_matrix_header_cells(rna_obj$matrix_path)
rna_dt <- read_dense_subset(rna_obj$matrix_path, rna_union_ids, header_info = rna_header_info)
rna_cells_in_file <- rna_header_info$cell_ids
rna_keep_cols <- match(rna_obj$cell_ids, rna_cells_in_file)
if (anyNA(rna_keep_cols)) stop("Some StepC RNA cell IDs were not found in subset matrix header.")
rna_gene_order <- sub("\\..*$", "", rna_dt$gene_id)
rna_mat <- as.matrix(rna_dt[, rna_keep_cols + 1L, drop = FALSE])
mode(rna_mat) <- "numeric"
rownames(rna_mat) <- rna_gene_order
colnames(rna_mat) <- rna_obj$cell_ids
rna_meta <- rna_obj$meta
stopifnot(identical(rna_meta$cell_id, colnames(rna_mat)))
rna_lib <- if ("RNA_Counts" %in% names(rna_meta)) rna_meta$RNA_Counts else stop("RNA_Counts missing from GSE162170 RNA metadata.")
rna_tmat <- log_transform_counts(rna_mat, libsize = rna_lib, scale_factor = args$norm_scale_factor)
rna_scores <- score_from_matrix(rna_tmat, rna_map$program_rows[rna_valid])
colnames(rna_scores) <- paste0("score_", colnames(rna_scores))
rna_primary <- if (args$use_rna_singlets_only && "df_classification" %in% names(rna_meta)) {
  is.na(rna_meta$df_classification) | rna_meta$df_classification == "Singlet"
} else rep(TRUE, nrow(rna_meta))
rna_cell_scores <- prepare_cell_score_output(rna_meta, rna_scores, "GSE162170_RNA", primary_filter = rna_primary)
safe_fwrite(rna_cell_scores, file.path(args$outdir, "GSE162170_RNA", "GSE162170_RNA_cell_program_scores.tsv.gz"))
rna_group_summary <- summarize_long(rna_cell_scores, score_cols = grep("^score_", names(rna_cell_scores), value = TRUE), group_var = "cluster_id", group_label_var = "cluster_id", only_primary = TRUE)
rna_age_summary <- summarize_long(rna_cell_scores, score_cols = grep("^score_", names(rna_cell_scores), value = TRUE), group_var = "age_label", group_label_var = "age_label", only_primary = TRUE)
rna_sample_summary <- summarize_long(rna_cell_scores, score_cols = grep("^score_", names(rna_cell_scores), value = TRUE), group_var = "sample_id", group_label_var = "sample_id", only_primary = TRUE)
safe_fwrite(rna_group_summary, file.path(args$outdir, "GSE162170_RNA", "GSE162170_RNA_group_program_summary.tsv.gz"))
safe_fwrite(rna_age_summary, file.path(args$outdir, "GSE162170_RNA", "GSE162170_RNA_age_program_summary.tsv.gz"))
safe_fwrite(rna_sample_summary, file.path(args$outdir, "GSE162170_RNA", "GSE162170_RNA_sample_program_summary.tsv.gz"))
saveRDS(list(dataset = "GSE162170_RNA", meta = rna_meta, score_df = rna_scores, valid_programs = rna_valid),
        file.path(args$outdir, "GSE162170_RNA", "GSE162170_RNA_scored_object.rds"), compress = "xz")

# ---- GSE162170 multiome ----
log_msg("Scoring GSE162170 multiome...")
mo_gene_ids <- read_gz_or_plain(mo_obj$gene_ids_file, header = TRUE)[[1]]
mo_gene_ids <- sub("\\..*$", "", mo_gene_ids)
mo_map <- map_programs_to_gene_ids(programs, program_gene_types, mo_gene_ids, gene_map, dataset_name = "GSE162170_multiome")
mapping_summaries[["GSE162170_multiome"]] <- mo_map$stats
mo_valid <- names(mo_map$program_rows)[vapply(mo_map$program_rows, length, integer(1)) >= args$min_overlap]
if (!length(mo_valid)) stop("No gene sets met min_overlap for GSE162170 multiome.")
mo_union_ids <- unique(unlist(mo_map$program_rows[mo_valid], use.names = FALSE))
mo_header_info <- extract_matrix_header_cells(mo_obj$matrix_path)
mo_dt <- read_dense_subset(mo_obj$matrix_path, mo_union_ids, header_info = mo_header_info)
mo_cells_in_file <- mo_header_info$cell_ids
mo_keep_cols <- match(mo_obj$cell_ids, mo_cells_in_file)
if (anyNA(mo_keep_cols)) stop("Some StepC multiome cell IDs were not found in subset matrix header.")
mo_gene_order <- sub("\\..*$", "", mo_dt$gene_id)
mo_mat <- as.matrix(mo_dt[, mo_keep_cols + 1L, drop = FALSE])
mode(mo_mat) <- "numeric"
rownames(mo_mat) <- mo_gene_order
colnames(mo_mat) <- mo_obj$cell_ids
mo_meta <- mo_obj$meta
stopifnot(identical(mo_meta$cell_id, colnames(mo_mat)))
mo_tmat <- transform_continuous(mo_mat, dataset_name = "GSE162170_multiome")
mo_scores <- score_from_matrix(mo_tmat, mo_map$program_rows[mo_valid])
colnames(mo_scores) <- paste0("score_", colnames(mo_scores))
mo_primary <- if (args$use_multiome_singlets_only && "df_classification" %in% names(mo_meta)) {
  is.na(mo_meta$df_classification) | mo_meta$df_classification == "Singlet"
} else rep(TRUE, nrow(mo_meta))
mo_cell_scores <- prepare_cell_score_output(mo_meta, mo_scores, "GSE162170_multiome", primary_filter = mo_primary)
safe_fwrite(mo_cell_scores, file.path(args$outdir, "GSE162170_multiome", "GSE162170_multiome_cell_program_scores.tsv.gz"))
mo_group_var <- if ("cluster_name" %in% names(mo_cell_scores) && any(!is.na(mo_cell_scores$cluster_name))) "cluster_name" else "cluster_id"
mo_group_summary <- summarize_long(mo_cell_scores, score_cols = grep("^score_", names(mo_cell_scores), value = TRUE), group_var = mo_group_var, group_label_var = mo_group_var, only_primary = TRUE)
mo_age_summary <- summarize_long(mo_cell_scores, score_cols = grep("^score_", names(mo_cell_scores), value = TRUE), group_var = "age_label", group_label_var = "age_label", only_primary = TRUE)
mo_sample_summary <- summarize_long(mo_cell_scores, score_cols = grep("^score_", names(mo_cell_scores), value = TRUE), group_var = "sample_id", group_label_var = "sample_id", only_primary = TRUE)
safe_fwrite(mo_group_summary, file.path(args$outdir, "GSE162170_multiome", "GSE162170_multiome_group_program_summary.tsv.gz"))
safe_fwrite(mo_age_summary, file.path(args$outdir, "GSE162170_multiome", "GSE162170_multiome_age_program_summary.tsv.gz"))
safe_fwrite(mo_sample_summary, file.path(args$outdir, "GSE162170_multiome", "GSE162170_multiome_sample_program_summary.tsv.gz"))
saveRDS(list(dataset = "GSE162170_multiome", meta = mo_meta, score_df = mo_scores, valid_programs = mo_valid),
        file.path(args$outdir, "GSE162170_multiome", "GSE162170_multiome_scored_object.rds"), compress = "xz")

# ---- Combined outputs ----
map_summary <- rbindlist(mapping_summaries, use.names = TRUE, fill = TRUE)
map_summary$passed_min_overlap <- map_summary$n_mapped >= args$min_overlap
safe_fwrite(map_summary, file.path(args$outdir, "StepD_gene_set_mapping_summary.tsv"))

all_group <- rbindlist(list(now_group_summary, rna_group_summary, mo_group_summary), use.names = TRUE, fill = TRUE)
all_age   <- rbindlist(list(now_age_summary, rna_age_summary, mo_age_summary), use.names = TRUE, fill = TRUE)
all_sample <- rbindlist(list(rna_sample_summary, mo_sample_summary), use.names = TRUE, fill = TRUE)
safe_fwrite(all_group, file.path(args$outdir, "StepD_all_group_program_summary.tsv.gz"))
safe_fwrite(all_age, file.path(args$outdir, "StepD_all_age_program_summary.tsv.gz"))
safe_fwrite(all_sample, file.path(args$outdir, "StepD_all_sample_program_summary.tsv.gz"))
safe_fwrite(build_top_hits(all_group, n_top = 5L), file.path(args$outdir, "StepD_top5_groups_per_dataset_program.tsv"))

run_summary <- data.frame(
  dataset = c("Nowakowski_UCSC", "GSE162170_RNA", "GSE162170_multiome"),
  n_cells = c(nrow(now_meta), nrow(rna_meta), nrow(mo_meta)),
  n_primary_cells = c(sum(now_primary), sum(rna_primary), sum(mo_primary)),
  group_var = c("cell_type", "cluster_id", mo_group_var),
  valid_programs = c(paste(now_valid, collapse = ";"), paste(rna_valid, collapse = ";"), paste(mo_valid, collapse = ";")),
  stringsAsFactors = FALSE
)
safe_fwrite(run_summary, file.path(args$outdir, "StepD_run_summary.tsv"))

log_msg("Finished StepD external reference program scoring.")
