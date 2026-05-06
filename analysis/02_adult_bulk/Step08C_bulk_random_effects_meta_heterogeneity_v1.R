#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
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
BULK_TABLE <- if (!is.null(args$bulk_table)) args$bulk_table else file.path(PROJECT_ROOT, "step01_prepare", "tables", "126_bulk_threeCohort_program_effects_updated.tsv")
SELECTION_TABLE <- if (!is.null(args$selection_table)) args$selection_table else file.path(PROJECT_ROOT, "step01_prepare", "tables", "131_final_program_selection_updated.tsv")
OUTDIR <- if (!is.null(args$outdir)) args$outdir else file.path(PROJECT_ROOT, "step08C_bulk_random_effects_meta_heterogeneity")

TABDIR <- file.path(OUTDIR, "tables")
PLOTDIR <- file.path(OUTDIR, "plots")
LOGDIR <- file.path(OUTDIR, "logs")
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)
LOGFILE <- file.path(LOGDIR, "Step08C_bulk_random_effects_meta_heterogeneity_v1.log")

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

infer_se_from_beta_p <- function(beta, p) {
  if (!is.finite(beta) || !is.finite(p) || p <= 0 || p > 1) return(NA_real_)
  z <- qnorm(1 - p / 2)
  if (!is.finite(z) || z <= 0) return(NA_real_)
  if (abs(beta) == 0) return(NA_real_)
  abs(beta) / z
}

meta_one_program <- function(dt_sub) {
  x <- dt_sub$effect_beta
  se <- dt_sub$effect_se
  keep <- is.finite(x) & is.finite(se) & se > 0
  x <- x[keep]
  se <- se[keep]
  if (length(x) == 0) {
    return(data.table(
      k_valid = 0L,
      fixed_beta = NA_real_,
      fixed_se = NA_real_,
      fixed_z = NA_real_,
      fixed_p_two_sided = NA_real_,
      fixed_ci_lo = NA_real_,
      fixed_ci_hi = NA_real_,
      random_beta = NA_real_,
      random_se = NA_real_,
      random_z = NA_real_,
      random_p_two_sided = NA_real_,
      random_ci_lo = NA_real_,
      random_ci_hi = NA_real_,
      Q = NA_real_,
      Q_df = NA_integer_,
      Q_p = NA_real_,
      I2 = NA_real_,
      tau2 = NA_real_
    ))
  }

  wi <- 1 / (se ^ 2)
  fixed_beta <- sum(wi * x) / sum(wi)
  fixed_se <- sqrt(1 / sum(wi))
  fixed_z <- fixed_beta / fixed_se
  fixed_p <- 2 * pnorm(-abs(fixed_z))
  fixed_ci <- fixed_beta + c(-1.96, 1.96) * fixed_se

  Q <- sum(wi * (x - fixed_beta) ^ 2)
  df <- length(x) - 1L
  Q_p <- if (df > 0) pchisq(Q, df = df, lower.tail = FALSE) else NA_real_
  c_term <- sum(wi) - (sum(wi ^ 2) / sum(wi))
  tau2 <- if (df > 0 && is.finite(c_term) && c_term > 0) max((Q - df) / c_term, 0) else 0
  wi_re <- 1 / (se ^ 2 + tau2)
  random_beta <- sum(wi_re * x) / sum(wi_re)
  random_se <- sqrt(1 / sum(wi_re))
  random_z <- random_beta / random_se
  random_p <- 2 * pnorm(-abs(random_z))
  random_ci <- random_beta + c(-1.96, 1.96) * random_se
  I2 <- if (df > 0 && is.finite(Q) && Q > 0) max((Q - df) / Q, 0) * 100 else 0

  data.table(
    k_valid = as.integer(length(x)),
    fixed_beta = fixed_beta,
    fixed_se = fixed_se,
    fixed_z = fixed_z,
    fixed_p_two_sided = fixed_p,
    fixed_ci_lo = fixed_ci[1],
    fixed_ci_hi = fixed_ci[2],
    random_beta = random_beta,
    random_se = random_se,
    random_z = random_z,
    random_p_two_sided = random_p,
    random_ci_lo = random_ci[1],
    random_ci_hi = random_ci[2],
    Q = Q,
    Q_df = as.integer(df),
    Q_p = Q_p,
    I2 = I2,
    tau2 = tau2
  )
}

plot_forest <- function(dt, outfile, title_text) {
  d <- copy(dt)
  if (nrow(d) == 0) return(invisible(NULL))
  d[, label := ifelse(level == "meta_random", "Random-effects pooled",
               ifelse(level == "meta_fixed", "Fixed-effects pooled", cohort))]
  d[, label := factor(label, levels = rev(unique(label)))]
  p <- ggplot(d, aes(x = beta, y = label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = level), height = 0.18, size = 0.5) +
    geom_point(aes(color = level, shape = level), size = 2.2) +
    facet_wrap(~ program_name, scales = "free_y", ncol = 2) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = "grey80")) +
    labs(x = "ASD effect (delta ASD - Control)", y = NULL, title = title_text)
  ggsave(outfile, p, width = 12, height = max(5, 2.2 * ceiling(length(unique(d$program_name)) / 2)))
}

log_msg("Reading bulk table: ", BULK_TABLE)
bulk <- fread(BULK_TABLE)
log_msg("Reading selection table: ", SELECTION_TABLE)
sel <- fread(SELECTION_TABLE)

required_cols <- c("cohort", "program_name", "delta_asd_minus_control", "lm_p", "n_samples")
missing_cols <- setdiff(required_cols, names(bulk))
if (length(missing_cols) > 0) stop("Missing required columns in bulk table: ", paste(missing_cols, collapse = ", "))

bulk[, effect_beta := as.numeric(delta_asd_minus_control)]
bulk[, effect_p := as.numeric(lm_p)]
bulk[, effect_p := pmax(pmin(effect_p, 1 - 1e-15), 1e-300)]
bulk[, effect_se := mapply(infer_se_from_beta_p, effect_beta, effect_p)]
bulk[, ci_lo := effect_beta - 1.96 * effect_se]
bulk[, ci_hi := effect_beta + 1.96 * effect_se]
bulk[, se_inferred_from_p := TRUE]

safe_fwrite(bulk, file.path(TABDIR, "01_bulk_effect_inputs_with_inferred_se.tsv"))

meta_summary <- bulk[, meta_one_program(.SD), by = .(program_name)]
meta_summary <- merge(meta_summary, unique(bulk[, .(program_name, adult_expected_direction)]), by = "program_name", all.x = TRUE)
meta_summary <- merge(meta_summary, sel[, .(program_name, selection_role, rationale)], by = "program_name", all.x = TRUE)
meta_summary[, fdr_random_two_sided := p.adjust(random_p_two_sided, method = "BH")]
meta_summary[, fdr_fixed_two_sided := p.adjust(fixed_p_two_sided, method = "BH")]
setorder(meta_summary, fdr_random_two_sided, fdr_fixed_two_sided, random_p_two_sided)
safe_fwrite(meta_summary, file.path(TABDIR, "02_bulk_random_effects_meta_summary.tsv"))

hetero <- meta_summary[, .(
  program_name,
  selection_role,
  adult_expected_direction,
  k_valid,
  Q,
  Q_df,
  Q_p,
  I2,
  tau2,
  random_beta,
  random_ci_lo,
  random_ci_hi,
  random_p_two_sided,
  fdr_random_two_sided
)]
safe_fwrite(hetero, file.path(TABDIR, "03_bulk_heterogeneity_stats.tsv"))

sel_join <- merge(sel, bulk, by = "program_name", all.x = TRUE)
safe_fwrite(sel_join, file.path(TABDIR, "04_bulk_selection_joined_effects.tsv"))

forest_input <- rbindlist(lapply(unique(bulk$program_name), function(pg) {
  d <- copy(bulk[program_name == pg, .(program_name, cohort, beta = effect_beta, ci_lo, ci_hi)])
  ms <- meta_summary[program_name == pg]
  d2 <- rbindlist(list(
    d[, level := "cohort"],
    data.table(program_name = pg, cohort = "Fixed-effects pooled", beta = ms$fixed_beta, ci_lo = ms$fixed_ci_lo, ci_hi = ms$fixed_ci_hi, level = "meta_fixed"),
    data.table(program_name = pg, cohort = "Random-effects pooled", beta = ms$random_beta, ci_lo = ms$random_ci_lo, ci_hi = ms$random_ci_hi, level = "meta_random")
  ), use.names = TRUE, fill = TRUE)
  d2
}), use.names = TRUE, fill = TRUE)
safe_fwrite(forest_input, file.path(TABDIR, "05_bulk_forest_input.tsv"))

primary_programs <- sel[selection_role %in% c("primary_broad_anchor", "primary_developmental_context_program"), unique(program_name)]
if (length(primary_programs) > 0) {
  plot_forest(forest_input[program_name %in% primary_programs],
              file.path(PLOTDIR, "01_forest_primary_programs.png"),
              "Primary ASD risk-gene programs: cohort effects and pooled meta-estimates")
}
plot_forest(forest_input,
            file.path(PLOTDIR, "02_forest_all_programs.png"),
            "All retained ASD risk-gene programs: cohort effects and pooled meta-estimates")

run_summary <- rbindlist(list(
  data.table(item = "n_programs", value = uniqueN(bulk$program_name)),
  data.table(item = "n_cohorts", value = uniqueN(bulk$cohort)),
  data.table(item = "programs", value = paste(unique(bulk$program_name), collapse = ",")),
  data.table(item = "cohorts", value = paste(unique(bulk$cohort), collapse = ",")),
  data.table(item = "best_random_effects_program", value = meta_summary[which.min(fdr_random_two_sided), program_name]),
  data.table(item = "best_random_effects_fdr", value = meta_summary[which.min(fdr_random_two_sided), fdr_random_two_sided])
), use.names = TRUE, fill = TRUE)
safe_fwrite(run_summary, file.path(TABDIR, "06_run_summary.tsv"))

log_msg("Step08C completed successfully. Outputs written to: ", OUTDIR)
