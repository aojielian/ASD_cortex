suppressPackageStartupMessages(library(data.table))

root_dir <- "/gpfs/hpc/home/lijc/lianaoj/ASD_cortex"
expr_path <- file.path(root_dir, "step01_prepare/tables/122_GSE64018_collapsed_counts_by_symbol.tsv.gz")
sfari_path <- file.path(root_dir, "step01_prepare/tables/01_SFARI_primary_Sand1_standardized.tsv")
syn_path <- file.path(root_dir, "step01_prepare/reviewer_round5/MA_round3A_synaptic_leave1out_v6/01_synaptic_gene_inventory.tsv")

std_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NAN", "NULL")] <- NA_character_
  x
}

choose_gene_col <- function(dt) {
  cand <- intersect(names(dt), c("gene","Gene","symbol","SYMBOL","gene_symbol","std_gene_id","Gene.Symbol"))
  if (length(cand) > 0) return(cand[1])
  cand2 <- names(dt)[grepl("gene|symbol", names(dt), ignore.case = TRUE)]
  if (length(cand2) > 0) return(cand2[1])
  names(dt)[1]
}

expr_dt <- fread(expr_path)
expr_genes <- unique(na.omit(std_gene(expr_dt[[1]])))

sfari_dt <- fread(sfari_path)
sfari_all <- unique(na.omit(std_gene(sfari_dt[[choose_gene_col(sfari_dt)]])))

syn_dt <- fread(syn_path)
sfari_syn <- unique(na.omit(std_gene(syn_dt[[choose_gene_col(syn_dt)]])))

mapped_all <- sort(intersect(expr_genes, sfari_all))
mapped_syn <- sort(intersect(expr_genes, sfari_syn))

cat("n_expr_genes =", length(expr_genes), "\n")
cat("n_SFARI_all_input =", length(sfari_all), "\n")
cat("n_SFARI_syn_input =", length(sfari_syn), "\n")
cat("n_mapped_all =", length(mapped_all), "\n")
cat("n_mapped_syn =", length(mapped_syn), "\n")
cat("n_intersection =", length(intersect(mapped_all, mapped_syn)), "\n")
cat("sets_identical =", identical(mapped_all, mapped_syn), "\n")

write.table(data.frame(gene = mapped_all), "mapped_SFARI_all_GSE64018.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = mapped_syn), "mapped_SFARI_all_synaptic_GSE64018.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
