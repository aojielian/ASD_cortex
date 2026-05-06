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
  x <- toupper(x)
  x
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

  # wide format
  if (target_program %in% names(dt)) {
    genes <- unique(na.omit(std_gene(dt[[target_program]])))
    genes <- genes[nzchar(genes)]
    if (length(genes) > 0) {
      return(list(genes = genes, source = program_path, method = paste0("wide_column:", target_program)))
    }
  }

  # long format
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

  # fallback: any row containing target_program; take gene-like columns
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
    rn <- rownames(x); cn <- colnames(x)
    out[[length(out) + 1L]] <- list(
      path = path,
      object = x,
      nrow = NROW(x),
      ncol = NCOL(x),
      rownames_present = !is.null(rn),
      colnames_present = !is.null(cn)
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

choose_expr_from_rds <- function(obj, target_genes, outdir = NULL) {
  cands <- collect_matrix_candidates(obj, path = "rds")
  if (length(cands) == 0) stop("No numeric matrix-like candidates found in GSE102741 RDS object.")
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
  best_obj <- cands[[which(vapply(cands, function(z) identical(z$path, best_path), logical(1)))[1]]]$object
  if (!is.null(outdir)) {
    safe_write(inv, file.path(outdir, "08_GSE102741_rds_candidate_inventory.tsv"))
  }
  list(expr = safe_num_mat(best_obj), chosen_path = best_path, inventory = inv)
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

main <- function() {
  root_dir <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
  outdir <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3C_gene_harmonization")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  sfari_path <- file.path(root_dir, "step01_prepare/tables/01_SFARI_primary_Sand1_standardized.tsv")
  program_path <- file.path(root_dir, "step01_prepare/tables/20_midPrenatal_core_programs.tsv")
  syn_inventory_path <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3A_synaptic_leave1out_v6/01_synaptic_gene_inventory.tsv")
  gse102741_rds <- file.path(root_dir, "step01_prepare/rds/22_GSE102741_geneSymbol_log2cpm.rds")
  gse64018_expr <- file.path(root_dir, "step01_prepare/tables/122_GSE64018_collapsed_counts_by_symbol.tsv.gz")
  gandal_rdata <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData"

  stopifnot(file.exists(sfari_path), file.exists(program_path), file.exists(syn_inventory_path),
            file.exists(gse102741_rds), file.exists(gse64018_expr), file.exists(gandal_rdata))

  # Program definitions
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

  all_target_genes <- unique(program_inventory$gene)

  # Dataset 1: GSE102741
  gse102741_obj <- readRDS(gse102741_rds)
  chosen102 <- if (is.matrix(gse102741_obj) || is.data.frame(gse102741_obj)) {
    list(expr = safe_num_mat(gse102741_obj), chosen_path = "rds")
  } else {
    choose_expr_from_rds(gse102741_obj, target_genes = all_target_genes, outdir = outdir)
  }
  expr102 <- chosen102$expr
  genes102_raw <- rownames(expr102)
  genes102_std <- unique(na.omit(std_gene(genes102_raw)))
  genes102_std <- genes102_std[nzchar(genes102_std)]

  # Dataset 2: GSE64018
  expr64018_dt <- fread(gse64018_expr)
  gene_col64018 <- names(expr64018_dt)[1]
  genes64018_raw <- expr64018_dt[[gene_col64018]]
  genes64018_std <- unique(na.omit(std_gene(genes64018_raw)))
  genes64018_std <- genes64018_std[nzchar(genes64018_std)]

  # Dataset 3: Gandal2022
  env <- new.env(parent = emptyenv())
  load(gandal_rdata, envir = env)
  stopifnot(exists("datExpr", envir = env, inherits = FALSE))
  datExpr <- get("datExpr", envir = env)
  genesG_raw <- rownames(datExpr)
  genesG_ensg_clean <- clean_gandal_ensg(genesG_raw)
  gmap <- map_ensembl_to_symbol(genesG_ensg_clean)
  gmap[, ensembl_std := std_gene(ensembl)]
  gmap[, symbol_std := std_gene(symbol)]
  idx_dt <- data.table(raw_gene_id = genesG_raw,
                       ensembl_clean = genesG_ensg_clean,
                       ensembl_std = std_gene(genesG_ensg_clean))
  idx_dt <- merge(idx_dt, gmap[, .(ensembl_std, symbol_std)], by = "ensembl_std", all.x = TRUE)
  genesG_std <- unique(na.omit(idx_dt$symbol_std))
  genesG_std <- genesG_std[nzchar(genesG_std)]

  safe_write(idx_dt, file.path(outdir, "09_Gandal_row_mapping_table.tsv.gz"))

  # Dataset gene space summary
  gene_space_summary <- rbindlist(list(
    data.table(dataset = "GSE102741",
               raw_gene_id_type = "gene_symbol_like",
               mapping_strategy = "direct_standardized_symbol_from_chosen_RDS_object",
               chosen_object = chosen102$chosen_path,
               n_raw_gene_ids = length(genes102_raw),
               n_unique_mapped_genes = length(genes102_std),
               n_unmapped_or_dropped = sum(is.na(std_gene(genes102_raw)) | !nzchar(std_gene(genes102_raw))),
               mapped_fraction = NA_real_),
    data.table(dataset = "GSE64018",
               raw_gene_id_type = "gene_symbol_like",
               mapping_strategy = "direct_standardized_symbol_from_collapsed_counts_first_column",
               chosen_object = gene_col64018,
               n_raw_gene_ids = length(genes64018_raw),
               n_unique_mapped_genes = length(genes64018_std),
               n_unmapped_or_dropped = sum(is.na(std_gene(genes64018_raw)) | !nzchar(std_gene(genes64018_raw))),
               mapped_fraction = NA_real_),
    data.table(dataset = "Gandal2022",
               raw_gene_id_type = "ENSG.version_suffix",
               mapping_strategy = "clean_ENSG_suffix_and_version_then_ENSEMBL_to_SYMBOL_via_org.Hs.eg.db",
               chosen_object = "datExpr",
               n_raw_gene_ids = length(genesG_raw),
               n_unique_mapped_genes = length(genesG_std),
               n_unmapped_or_dropped = sum(is.na(idx_dt$symbol_std) | !nzchar(idx_dt$symbol_std)),
               mapped_fraction = mean(!is.na(idx_dt$symbol_std) & nzchar(idx_dt$symbol_std)))
  ), fill = TRUE)
  safe_write(gene_space_summary, file.path(outdir, "02_dataset_gene_space_summary.tsv"))

  # Dataset/program mapped counts
  dataset_gene_spaces <- list(
    GSE102741 = genes102_std,
    GSE64018 = genes64018_std,
    Gandal2022 = genesG_std
  )
  dataset_gene_space_raw_counts <- data.table(
    dataset = c("GSE102741","GSE64018","Gandal2022"),
    n_dataset_genes = c(length(genes102_std), length(genes64018_std), length(genesG_std))
  )

  mapped_lists <- list()
  counts_rows <- list()
  list_rows <- list()
  unmapped_rows <- list()

  for (ds in names(dataset_gene_spaces)) {
    ds_genes <- dataset_gene_spaces[[ds]]
    mapped_lists[[ds]] <- list()
    for (pg in unique(program_inventory$program)) {
      pg_genes <- unique(program_inventory[program == pg, gene])
      mapped <- sort(intersect(ds_genes, pg_genes))
      dropped <- sort(setdiff(pg_genes, ds_genes))
      mapped_lists[[ds]][[pg]] <- mapped

      counts_rows[[paste(ds, pg, sep = "::")]] <- data.table(
        dataset = ds,
        program = pg,
        n_program_input = length(pg_genes),
        n_program_mapped = length(mapped),
        n_program_dropped = length(dropped),
        mapped_fraction_of_program = ifelse(length(pg_genes) > 0, length(mapped) / length(pg_genes), NA_real_),
        n_dataset_gene_space = length(ds_genes)
      )

      if (length(mapped) > 0) {
        list_rows[[paste(ds, pg, "mapped", sep = "::")]] <- data.table(dataset = ds, program = pg, status = "mapped", gene = mapped)
      } else {
        list_rows[[paste(ds, pg, "mapped", sep = "::")]] <- data.table(dataset = character(), program = character(), status = character(), gene = character())
      }
      if (length(dropped) > 0) {
        unmapped_rows[[paste(ds, pg, "dropped", sep = "::")]] <- data.table(dataset = ds, program = pg, status = "dropped", gene = dropped)
      } else {
        unmapped_rows[[paste(ds, pg, "dropped", sep = "::")]] <- data.table(dataset = character(), program = character(), status = character(), gene = character())
      }
    }
  }

  counts_dt <- rbindlist(counts_rows, fill = TRUE)

  # Add relation to SFARI_all_synaptic within each dataset
  syn_ref <- counts_dt[program == "SFARI_all_synaptic", .(dataset, syn_mapped_n = n_program_mapped)]
  counts_dt <- merge(counts_dt, syn_ref, by = "dataset", all.x = TRUE)

  rel_rows <- list()
  for (ds in names(mapped_lists)) {
    syn_set <- mapped_lists[[ds]][["SFARI_all_synaptic"]]
    for (pg in names(mapped_lists[[ds]])) {
      set_pg <- mapped_lists[[ds]][[pg]]
      rel_rows[[paste(ds, pg, sep = "::")]] <- data.table(
        dataset = ds,
        program = pg,
        overlap_with_synaptic = length(intersect(set_pg, syn_set)),
        jaccard_with_synaptic = ifelse(length(union(set_pg, syn_set)) > 0, length(intersect(set_pg, syn_set)) / length(union(set_pg, syn_set)), NA_real_),
        mapped_set_identical_to_SFARI_all_synaptic = identical(sort(set_pg), sort(syn_set))
      )
    }
  }
  rel_dt <- rbindlist(rel_rows, fill = TRUE)

  counts_dt <- merge(counts_dt, rel_dt, by = c("dataset","program"), all.x = TRUE)

  counts_dt[, note := fifelse(program == "SFARI_all" & mapped_set_identical_to_SFARI_all_synaptic,
                              "Mapped set collapsed to SFARI_all_synaptic in this dataset",
                              fifelse(program == "SFARI_all_synaptic",
                                      "Reference synaptic mapped set",
                                      fifelse(program == "midPrenatal_SFARI_top20",
                                              "Developmental-context program",
                                              "")))]

  setorderv(counts_dt, c("dataset","program"), c(1L,1L))
  safe_write(counts_dt, file.path(outdir, "03_post_mapping_counts_table.tsv"))

  mapped_gene_list_dt <- rbindlist(list_rows, fill = TRUE)
  unmapped_gene_list_dt <- rbindlist(unmapped_rows, fill = TRUE)
  safe_write(mapped_gene_list_dt, file.path(outdir, "04_mapped_gene_lists.tsv.gz"))
  safe_write(unmapped_gene_list_dt, file.path(outdir, "05_dropped_gene_lists.tsv.gz"))

  # Pairwise overlaps per dataset
  overlap_dt <- rbindlist(lapply(names(mapped_lists), function(ds) {
    x <- make_pairwise_overlap(mapped_lists[[ds]])
    x[, dataset := ds]
    x
  }), fill = TRUE)
  setcolorder(overlap_dt, c("dataset","program1","program2","n1","n2","n_intersection","n_union","jaccard","identical_sets","p1_subset_p2","p2_subset_p1"))
  safe_write(overlap_dt, file.path(outdir, "06_pairwise_program_overlap_by_dataset.tsv"))

  # Reviewer-facing compact summary
  reviewer_dt <- copy(counts_dt)
  reviewer_dt[, mapping_comment := fifelse(dataset == "Gandal2022",
                                           "ENSG.version_suffix cleaned and mapped to SYMBOL",
                                           "Direct standardized symbol matching")]
  reviewer_dt[, collapse_flag := fifelse(mapped_set_identical_to_SFARI_all_synaptic, "YES", "NO")]
  setcolorder(reviewer_dt, c("dataset","program","n_program_input","n_program_mapped","n_program_dropped",
                             "mapped_fraction_of_program","overlap_with_synaptic","jaccard_with_synaptic",
                             "mapped_set_identical_to_SFARI_all_synaptic","collapse_flag","mapping_comment","note",
                             "n_dataset_gene_space","syn_mapped_n"))
  safe_write(reviewer_dt, file.path(outdir, "07_reviewer_facing_mapping_summary.tsv"))

  # Run info
  run_info <- data.table(
    time = timestamp(),
    root_dir = root_dir,
    outdir = outdir,
    sfari_path = sfari_path,
    program_path = program_path,
    syn_inventory_path = syn_inventory_path,
    gse102741_rds = gse102741_rds,
    gse64018_expr = gse64018_expr,
    gandal_rdata = gandal_rdata,
    n_SFARI_all_input = length(sfari_all),
    n_SFARI_all_synaptic_input = length(sfari_syn),
    n_midPrenatal_top20_input = length(top20),
    GSE102741_chosen_object = chosen102$chosen_path
  )
  safe_write(run_info, file.path(outdir, "00_run_info.txt"))

  logmsg("Done. Reviewer-facing mapping summary written to: ", file.path(outdir, "07_reviewer_facing_mapping_summary.tsv"))
}

main()
