library(ggplot2)
library(sf)
library(tidyr)
library(mgcv)
library(dplyr)
library(utils)
library(xtable)
library(Matrix)
library(tidyverse)
library(parallel)
library(geosphere)
library(patchwork)

source('funcs.R')
load('sim.RData')

# ---- Configuration ----
# Point this at whatever results directory you've pulled back from the cluster.
# Safe to run against a still-running cluster job: any (mechanism, option,
# method) combination without a file yet just shows up as NA/blank below,
# nothing errors.
# Usage: Rscript analysis.R [results_dir] [results_label]
# Defaults match the full run built by submit_jobs.sh.
args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) args[[1]] else "results_stability_1000reps_quadratic_alpha0.05/"
results_label <- if (length(args) >= 2L) args[[2]] else "oct2026_1000reps_alpha0.05"
n_mechanisms <- 8
mechanism_levels <- as.character(1:n_mechanisms)
option_levels <- c("linear", "nonlinear")

dir.create("images", showWarnings = FALSE)

# Ground truth (see the mechanism functions in funcs.R). In R^2 units every
# mechanism has the same truth: the confounded variation explains
# confounded_share (0.5) of Var(A). As a column count
# in the selection's search order (meaningful for the matched basis only):
# 10 for mechanisms 1, 2, 3, 7 and m - 20 for the reversed mechanisms 4, 8.
# Mechanisms 5 (Leroux CAR) and 6 (coordinate-based) have no column truth.
m_tps <- ncol(simlist$B_tps_full)
m_gl  <- ncol(simlist$B_gl_full)
stopifnot(m_tps == m_gl)
true_n_ac_by_mech <- c(`1` = 10L, `2` = 10L, `3` = 10L, `4` = m_tps - 20L,
                       `7` = 10L, `8` = m_tps - 20L)
mechanisms_with_true_nc <- as.integer(names(true_n_ac_by_mech))

# Only these methods run basis selection and produce n_uc/cor_Ac diagnostics.
# The +spatialcoord IV variants (and trueIV-spatialcoord) are dropped from
# this study -- plain spatialcoord is kept.
iv_methods <- c("IV-TPS", "IV-GraphLaplacian")

method_order <- c("oracle", "baseline", "spatialcoord", "trueIV",
                  "IV-TPS", "IV-GraphLaplacian")

# Every (mechanism, option, method) cell the full run is meant to eventually
# produce. Real data is matched onto this below rather than the other way
# around, so a still-running or not-yet-started cell shows up as an explicit
# NA/blank row instead of silently shrinking a table or (worse) crashing a
# pivot_wider() + select(all_of(...)) when an entire method has zero files
# anywhere yet.
scaffold <- expand.grid(
  confounding_mechanism = 1:n_mechanisms,
  option = option_levels,
  method = method_order,
  stringsAsFactors = FALSE
)

# ---- Discover point-estimate files that actually exist ----
# Anchored so it only matches conf{N}_{option}_{method}.csv, not the
# _ci_lower/_ci_upper/_cor_Ac/_n_uc/_n_uc_foldK auxiliary files that live in
# the same directory. Method names never contain underscores (only hyphens),
# so [^_]+ for the method segment is sufficient to exclude those suffixes.
csvs <- list.files(results_dir, pattern = "^conf[0-9]+_(linear|nonlinear)_[^_]+\\.csv$")
file_index <- data.frame(
  confounding_mechanism = as.integer(sub("^conf(\\d+)_.*$", "\\1", csvs)),
  option = sub("^conf\\d+_([^_]+)_.*$", "\\1", csvs),
  method = sub("^conf\\d+_[^_]+_([^_]+)\\.csv$", "\\1", csvs),
  csv = csvs,
  stringsAsFactors = FALSE
)
# Drop the dropped-from-the-study methods even if their files still exist on
# disk from an earlier run (trueIV-spatialcoord, IV-TPS-spatialcoord,
# IV-GraphLaplacian-spatialcoord) -- everything downstream reads from
# file_index, so filtering here is the single point of control.
file_index <- file_index %>% filter(method %in% method_order)

cat(sprintf(
  "%d of %d expected (mechanism x option x method) result files found in %s.\n",
  nrow(file_index), nrow(scaffold), results_dir
))
missing_combos <- anti_join(
  scaffold, file_index,
  by = c("confounding_mechanism", "option", "method")
) %>% arrange(confounding_mechanism, option, method)
if (nrow(missing_combos) > 0) {
  cat("Still missing (not yet run or not yet finished):\n")
  print(missing_combos, row.names = FALSE)
}

# ---- True estimand per (mechanism, option), from the oracle ----
mutrues <- expand.grid(
  confounding_mechanism = 1:n_mechanisms,
  option = option_levels,
  stringsAsFactors = FALSE
)
mutrues$theta <- NA_real_
for (i in seq_len(nrow(mutrues))) {
  f <- file.path(results_dir, sprintf(
    "conf%d_%s_oracle.csv", mutrues$confounding_mechanism[i], mutrues$option[i]
  ))
  if (file.exists(f)) {
    mutrues$theta[i] <- mean(as.vector(as.matrix(read.csv(f))), na.rm = TRUE)
  }
}

# ---- Bias / RMSE / SE of point estimates ----
# Built on the full scaffold so every method column exists in the pivoted
# tables below even if that method has no files anywhere yet.
analysisdf <- scaffold
analysisdf$bias <- NA_real_
analysisdf$RMSE <- NA_real_
analysisdf$se   <- NA_real_
for (i in seq_len(nrow(file_index))) {
  row <- file_index[i, ]
  mutrue <- mutrues$theta[mutrues$confounding_mechanism == row$confounding_mechanism &
                            mutrues$option == row$option]
  if (length(mutrue) == 0 || is.na(mutrue)) next
  muests <- as.vector(as.matrix(read.csv(file.path(results_dir, row$csv))))
  idx <- which(analysisdf$confounding_mechanism == row$confounding_mechanism &
                analysisdf$option == row$option &
                analysisdf$method == row$method)
  analysisdf$bias[idx] <- mean(muests, na.rm = TRUE) - mutrue
  analysisdf$RMSE[idx] <- sqrt(mean((muests - mutrue)^2, na.rm = TRUE))
  analysisdf$se[idx]   <- sd(muests, na.rm = TRUE)
}

analysisdf_fmt <- analysisdf %>%
  mutate(bias = round(bias * 100, 3),
        RMSE = round(RMSE * 100, 3),
        se = round(se * 100, 3)) %>%
  arrange(confounding_mechanism, option)

wide_bias <- analysisdf_fmt %>%
  select(confounding_mechanism, option, method, bias) %>%
  pivot_wider(names_from = method, values_from = bias) %>%
  select(confounding_mechanism, option, all_of(method_order)) %>%
  arrange(confounding_mechanism, option)
print(xtable(wide_bias), include.rownames = FALSE, sanitize.text.function = identity)
print(xtable(wide_bias %>% mutate(across(where(is.numeric), abs))),
      include.rownames = FALSE, sanitize.text.function = identity)

wide_rmse <- analysisdf_fmt %>%
  select(confounding_mechanism, option, method, RMSE) %>%
  pivot_wider(names_from = method, values_from = RMSE) %>%
  select(confounding_mechanism, option, all_of(method_order)) %>%
  arrange(confounding_mechanism, option)
print(xtable(wide_rmse), include.rownames = FALSE, sanitize.text.function = identity)

# ---- Boxplots of point estimates, all 8 mechanisms in plain numeric order ----
# Built only from what exists -- no scaffold needed here, a missing cell just
# means an empty/absent box in that facet panel, which ggplot already handles.
read_estimates <- function(i) {
  row <- file_index[i, ]
  dat <- read.csv(file.path(results_dir, row$csv))
  names(dat)[1] <- "estimate"
  dat$confounding_mechanism <- row$confounding_mechanism
  dat$option <- row$option
  dat$method <- row$method
  dat
}
df <- if (nrow(file_index) > 0) map_dfr(seq_len(nrow(file_index)), read_estimates) else NULL

if (!is.null(df) && nrow(df) > 0) {
  df$method[df$method == "IV-GraphLaplacian"] <- "IV-GL"
  df$method[df$method == "oracle"] <- "Oracle"
  df$method[df$method == "baseline"] <- "Baseline"
  df$method[df$method == "spatialcoord"] <- "Spatial coordinates"
  df$method[df$method == "trueIV"] <- "True IV"

  desired_order <- c("Oracle", "Baseline", "Spatial coordinates", "True IV", "IV-TPS", "IV-GL")
  df$method <- factor(df$method, levels = desired_order)
  df$confounding_mechanism <- factor(df$confounding_mechanism, levels = mechanism_levels)
  df$option <- factor(df$option, levels = option_levels)

  mutrues_f <- mutrues %>%
    filter(confounding_mechanism %in% unique(df$confounding_mechanism)) %>%
    mutate(confounding_mechanism = factor(confounding_mechanism, levels = mechanism_levels),
          option = factor(option, levels = option_levels))

  method_cols <- c(
    "Oracle"                 = "gray",
    "Baseline"               = "#D62728",
    "Spatial coordinates"    = "purple",
    "True IV"                 = "lightgreen",
    "IV-TPS"                 = "#9ECAE1",
    "IV-GL"                  = "#FDAE6B"
  )

  png(sprintf("images/boxplot_%s.png", results_label), width = 3200, height = 1700, res = 200)
  print(
    ggplot(df, aes(x = method, y = estimate, fill = method)) +
      geom_boxplot(alpha = 0.5, outliers = FALSE, staplewidth = 1) +
      #stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "blue") +
      ggh4x::facet_grid2(
        option ~ confounding_mechanism,
        scales = "free",
        independent = "all",
        labeller = labeller(confounding_mechanism = function(x) paste("Confounding \nmechanism", x))
      ) +
      geom_hline(data = mutrues_f, aes(yintercept = theta),
                 color = "red", linetype = "twodash", size = 1) +
      labs(
        x = "Method",
        y = "Truncated Exposure Effect Estimate"
      ) +
      scale_fill_manual(
        name = "Method",
        values = method_cols,
        breaks = names(method_cols)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "top"
      )
  )
  dev.off()
  
  dfsub <- df %>%
    filter(confounding_mechanism %in% 1:6) %>%
    droplevels()
  mutrues_f_sub <- mutrues_f %>%
    filter(confounding_mechanism %in% 1:6) %>%
    droplevels()
  png(sprintf("images/boxplot_%s_maintext.png", results_label), width = 2500, height = 1700, res = 200)
  print(
    ggplot(dfsub, aes(x = method, y = estimate, fill = method)) +
      geom_boxplot(alpha = 0.5, outliers = FALSE, staplewidth = 1) +
      #stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "blue") +
      ggh4x::facet_grid2(
        option ~ confounding_mechanism,
        scales = "free",
        independent = "all",
        labeller = labeller(confounding_mechanism = function(x) paste("Confounding \nmechanism", x))
      ) +
      geom_hline(data = mutrues_f_sub, aes(yintercept = theta),
                 color = "red", linetype = "twodash", size = 1) +
      labs(
        x = "Method",
        y = "Truncated Exposure Effect Estimate"
      ) +
      scale_fill_manual(
        name = "Method",
        values = method_cols,
        breaks = names(method_cols)
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "top"
      )
  )
  dev.off()
} else {
  cat("No point-estimate files found yet -- skipping boxplot.\n")
}

# ---- CI coverage and width ----
# Per-replicate truth, NOT the pooled mutrue used for bias/RMSE above.
# funcs.R's DGM (2026-09-16) redraws the confounder field fresh for every
# replicate, so each replicate has its own true population truncated-
# exposure-effect -- coverage has to be checked against each replicate's
# own truth, not one value shared across all reps (bias/RMSE against the
# pooled mutrue remain valid; only coverage needed this fix -- checking a
# narrow, correctly-centered per-replicate CI against a POOLED mutrue that
# differs from that replicate's own truth just measures between-replicate
# field variability, not miscalibration). true_psi_mech<N>.rds (one file
# per mechanism, living in results_dir itself) is produced automatically
# by submit_jobs.sh's truth_cm<N> jobs (compute_true_psi.R /
# run_job_truth.sh), submitted alongside the main cm<N> jobs every time
# submit_jobs.sh runs -- so it's always present for, and in sync with
# (same nsims/seed), whatever is in results_dir, with no separate file to
# keep up to date by hand. It replays set.seed(seed);
# compute_data_XXX(..., nsims) -- the exact call run_simfunc.R makes --
# to reconstruct each replicate's real confounder field, then runs
# estimate_truncated_effect_full() on a NOISELESS outcome (E[Y|A,U], no
# random error) using that replicate's own realized A/U. Validated
# against the real cluster output: reconstructed truth matched oracle's
# own per-replicate estimates to within its sampling noise, and in
# aggregate matched the pooled oracle mean to within ~0.1-0.5 percentage
# points across all 16 mechanism/option cells.
true_psi_files <- list.files(results_dir, pattern = "^true_psi_mech[0-9]+\\.rds$", full.names = TRUE)
if (length(true_psi_files) == 0) {
  stop("No true_psi_mech*.rds files found in ", results_dir,
       " -- did submit_jobs.sh's truth_cm<N> jobs run/finish?")
}
true_psi <- bind_rows(lapply(true_psi_files, readRDS))

ci_summary <- scaffold
ci_summary$coverage <- NA_real_
ci_summary$mean_width <- NA_real_
ci_summary$n <- 0L
for (i in seq_len(nrow(file_index))) {
  row <- file_index[i, ]
  lo_f <- file.path(results_dir, sub("\\.csv$", "_ci_lower.csv", row$csv))
  hi_f <- file.path(results_dir, sub("\\.csv$", "_ci_upper.csv", row$csv))
  if (!file.exists(lo_f) || !file.exists(hi_f)) next
  truth_sub <- true_psi[true_psi$mechanism == row$confounding_mechanism &
                          true_psi$option == row$option, ]
  if (nrow(truth_sub) == 0) next
  truth_sub <- truth_sub[order(truth_sub$replicate), ]
  lo <- as.vector(as.matrix(read.csv(lo_f)))
  hi <- as.vector(as.matrix(read.csv(hi_f)))
  if (length(lo) != nrow(truth_sub)) {
    warning(sprintf(
      "conf%d_%s: %d CI rows but %d truth rows -- using only the first %d. true_psi_mech%d.rds may be stale relative to the CSVs (re-run submit_jobs.sh's truth_cm%d job?).",
      row$confounding_mechanism, row$option, length(lo), nrow(truth_sub),
      min(length(lo), nrow(truth_sub)), row$confounding_mechanism, row$confounding_mechanism
    ))
  }
  n_common <- min(length(lo), nrow(truth_sub))
  psi_true <- truth_sub$psi_true[seq_len(n_common)]
  covers <- lo[seq_len(n_common)] <= psi_true & hi[seq_len(n_common)] >= psi_true
  idx <- which(ci_summary$confounding_mechanism == row$confounding_mechanism &
                ci_summary$option == row$option &
                ci_summary$method == row$method)
  ci_summary$coverage[idx] <- mean(covers, na.rm = TRUE)
  ci_summary$mean_width[idx] <- mean(hi[seq_len(n_common)] - lo[seq_len(n_common)], na.rm = TRUE)
  ci_summary$n[idx] <- sum(!is.na(covers))
}
ci_summary <- ci_summary %>% arrange(confounding_mechanism, option, method)
print(xtable(ci_summary, digits = 3), include.rownames = FALSE)

# ---- Error boxplots: estimate minus each replicate's own true psi ----
# The raw-estimate boxplots above compare every method with one pooled line
# (the oracle mean), so between-replicate variation in the truth shows up as
# spread, and a skewed truth distribution (e.g. mechanism 6) makes unbiased
# methods look offset (box = median, line = mean). Subtracting each
# replicate's own true psi (true_psi, as used for coverage) removes that:
# bias is the offset of the white diamond (mean error) from zero.
error_rows <- list()
for (i in seq_len(nrow(file_index))) {
  row <- file_index[i, ]
  truth_sub <- true_psi[true_psi$mechanism == row$confounding_mechanism &
                          true_psi$option == row$option, ]
  if (nrow(truth_sub) == 0) next
  truth_sub <- truth_sub[order(truth_sub$replicate), ]
  est <- as.vector(as.matrix(read.csv(file.path(results_dir, row$csv))))
  k <- min(length(est), nrow(truth_sub))
  error_rows[[length(error_rows) + 1]] <- data.frame(
    confounding_mechanism = row$confounding_mechanism, option = row$option,
    method = row$method, error = est[seq_len(k)] - truth_sub$psi_true[seq_len(k)]
  )
}
error_df <- if (length(error_rows) > 0) bind_rows(error_rows) else NULL

if (!is.null(error_df) && nrow(error_df) > 0) {
  method_labels <- c(oracle = "Oracle", baseline = "Baseline", spatialcoord = "Spatial coordinates",
                     trueIV = "True IV", `IV-TPS` = "IV-TPS", `IV-GraphLaplacian` = "IV-GL")
  error_cols <- c("Oracle" = "gray", "Baseline" = "#D62728", "Spatial coordinates" = "purple",
                  "True IV" = "lightgreen", "IV-TPS" = "#9ECAE1", "IV-GL" = "#FDAE6B")
  error_df$method <- factor(method_labels[error_df$method], levels = unname(method_labels))
  error_df$option <- factor(error_df$option, levels = option_levels)
  error_df$confounding_mechanism <- factor(error_df$confounding_mechanism, levels = mechanism_levels)
  error_df <- droplevels(error_df)

  error_plot <- function(d) {
    ggplot(d, aes(x = method, y = error, fill = method)) +
      geom_hline(yintercept = 0, color = "red", linetype = "twodash", linewidth = 1) +
      geom_boxplot(alpha = 0.5, outliers = FALSE, staplewidth = 1) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white") +
      ggh4x::facet_grid2(
        option ~ confounding_mechanism,
        scales = "free",
        independent = "all",
        labeller = labeller(confounding_mechanism = function(x) paste("Confounding \nmechanism", x))
      ) +
      labs(x = "Method", y = expression(paste("Estimate minus true ", psi, " of the replicate"))) +
      scale_fill_manual(name = "Method", values = error_cols, breaks = names(error_cols)) +
      theme_minimal(base_size = 14) +
      theme(
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "top"
      )
  }

  png(sprintf("images/error_boxplot_%s.png", results_label), width = 3200, height = 1700, res = 200)
  print(error_plot(error_df))
  dev.off()

  png(sprintf("images/error_boxplot_%s_maintext.png", results_label), width = 2500, height = 1700, res = 200)
  print(error_plot(droplevels(subset(error_df, confounding_mechanism %in% as.character(1:6)))))
  dev.off()
}

# ---- n_ac/n_uc selection recovery vs. DGP ground truth (basis-selection
# methods only) ----
# One n_ac (confounded dimension) and n_uc (= m - n_ac) selection per sim --
# no outer folds any more (simfunc() no longer cross-fits). bias_n_ac is
# mean_n_ac - true n_ac for the mechanisms with a column truth; NA for
# mechanisms 5/6. frac_no_adjust is the share of sims where selection walked
# all the way to n_ac = 1 (constant only, i.e. no spatial adjustment).
nac_summary <- expand.grid(
  confounding_mechanism = 1:n_mechanisms,
  option = option_levels,
  method = iv_methods,
  stringsAsFactors = FALSE
)
nac_summary$true_n_c <- unname(true_n_ac_by_mech[as.character(nac_summary$confounding_mechanism)])
nac_summary$mean_n_ac <- NA_real_
nac_summary$median_n_ac <- NA_real_
nac_summary$sd_n_ac <- NA_real_
nac_summary$bias_n_ac <- NA_real_
nac_summary$frac_no_adjust <- NA_real_
nac_summary$mean_n_uc <- NA_real_
nac_summary$n <- 0L
for (i in seq_len(nrow(nac_summary))) {
  mech <- nac_summary$confounding_mechanism[i]
  opt <- nac_summary$option[i]
  meth <- nac_summary$method[i]
  n_ac_file <- file.path(results_dir, sprintf("conf%d_%s_%s_n_ac.csv", mech, opt, meth))
  n_uc_file <- file.path(results_dir, sprintf("conf%d_%s_%s_n_uc.csv", mech, opt, meth))
  if (!file.exists(n_ac_file)) next
  n_ac <- unlist(read.csv(n_ac_file), use.names = FALSE)
  nac_summary$mean_n_ac[i]   <- mean(n_ac, na.rm = TRUE)
  nac_summary$median_n_ac[i] <- median(n_ac, na.rm = TRUE)
  nac_summary$sd_n_ac[i]     <- sd(n_ac, na.rm = TRUE)
  nac_summary$bias_n_ac[i]   <- nac_summary$mean_n_ac[i] - nac_summary$true_n_c[i]
  nac_summary$frac_no_adjust[i] <- mean(n_ac == 1, na.rm = TRUE)
  nac_summary$n[i]           <- sum(!is.na(n_ac))
  if (file.exists(n_uc_file)) {
    nac_summary$mean_n_uc[i] <- mean(unlist(read.csv(n_uc_file), use.names = FALSE), na.rm = TRUE)
  }
}
nac_summary <- nac_summary %>% arrange(confounding_mechanism, option, method)
print(xtable(nac_summary, digits = 2), include.rownames = FALSE)

# ---- A^c recovery (cor_Ac), basis-selection methods only ----
corac_summary <- expand.grid(
  confounding_mechanism = 1:n_mechanisms,
  option = option_levels,
  method = iv_methods,
  stringsAsFactors = FALSE
)
corac_summary$mean_cor_Ac <- NA_real_
corac_summary$sd_cor_Ac <- NA_real_
corac_summary$n <- 0L
for (i in seq_len(nrow(corac_summary))) {
  mech <- corac_summary$confounding_mechanism[i]
  opt <- corac_summary$option[i]
  meth <- corac_summary$method[i]
  f <- file.path(results_dir, sprintf("conf%d_%s_%s_cor_Ac.csv", mech, opt, meth))
  if (!file.exists(f)) next
  corac <- unlist(read.csv(f), use.names = FALSE)
  corac_summary$mean_cor_Ac[i] <- mean(corac, na.rm = TRUE)
  corac_summary$sd_cor_Ac[i] <- sd(corac, na.rm = TRUE)
  corac_summary$n[i] <- sum(!is.na(corac))
}
corac_summary <- corac_summary %>% arrange(confounding_mechanism, option, method)
print(xtable(corac_summary, digits = 3), include.rownames = FALSE)

# ---- True truncated exposure effect (oracle mean), for direct copy-paste ----
# Matches the format of the existing supplement subsection "True truncated
# exposure effects as estimated by the oracle mean".
cat("\n% ---- oracle-mean (tau*) table ----\n")
cat("\\subsection{True truncated exposure effects as estimated by the oracle mean}\n\n\n")
cat("\\renewcommand{\\arraystretch}{1}\n")
cat("\\begin{center}\\begin{tabular}{ccc}\n  \\hline\n")
cat("Confounding Mechanism & Outcome Model & $\\tau^*$ \\\\ \n  \\hline\n")
mutrues_sorted <- mutrues %>% arrange(confounding_mechanism, option)
for (i in seq_len(nrow(mutrues_sorted))) {
  cat(sprintf("%d & %s & %.4f \\\\ \n", mutrues_sorted$confounding_mechanism[i],
              mutrues_sorted$option[i], mutrues_sorted$theta[i]))
}
cat("   \\hline\n\\end{tabular}\\end{center}\n")

# ---- Coverage / CI width, wide table  ----
# Two stacked panels (coverage, then width) sharing one header row, one
# column per method. Built directly rather than through xtable, since the
# target layout (p{}-width columns, makecell headers, a multicolumn panel
# separator) isn't something xtable's automatic formatting produces.
coverage_methods <- c("oracle", "baseline", "spatialcoord", "trueIV", "IV-TPS", "IV-GraphLaplacian")
coverage_labels  <- c("Oracle", "Baseline", "Spatial\\\\coordinates", "True IV", "IV-TPS", "IV-GL")

print_wide_panel <- function(summary_df, value_col, digits = 3) {
  wide <- summary_df %>%
    select(confounding_mechanism, option, method, all_of(value_col)) %>%
    pivot_wider(names_from = method, values_from = all_of(value_col)) %>%
    select(confounding_mechanism, option, all_of(coverage_methods)) %>%
    arrange(confounding_mechanism, option)
  fmt <- paste0("%.", digits, "f")
  for (i in seq_len(nrow(wide))) {
    vals <- vapply(coverage_methods, function(m) {
      v <- wide[[m]][i]
      if (is.na(v)) "" else sprintf(fmt, v)
    }, character(1))
    cat(sprintf("%d & %s & %s \\\\ \n", wide$confounding_mechanism[i], wide$option[i],
                paste(vals, collapse = " & ")))
  }
}

cat("\n% ---- coverage / CI width table ----\n")
cat("\\begin{table}\n\\renewcommand{\\arraystretch}{1}\n")
cat("\\begin{tabular}{p{1.2cm}p{1cm}p{0.9cm}p{0.9cm}p{0.9cm}p{0.9cm}p{0.9cm}p{0.9cm}}\n")
cat(sprintf(
  "\\makecell[l]{Confounding \\\\Mechanism} & \\makecell[l]{Outcome \\\\Model} & %s & %s & \\makecell[l]{%s} & %s  & %s & %s  \\\\\n",
  coverage_labels[1], coverage_labels[2], coverage_labels[3],
  coverage_labels[4], coverage_labels[5], coverage_labels[6]
))
cat("\\hline\n\\multicolumn{8}{c}{Average coverage} \\\\ \n\\hline\n")
print_wide_panel(ci_summary, "coverage", digits = 3)
cat("  \\hline\n\\multicolumn{8}{c}{Average confidence interval width} \\\\ \n\\hline\n")
print_wide_panel(ci_summary, "mean_width", digits = 3)
cat("\\end{tabular}\n")
cat("    \\caption{Average coverage, defined as the proportion of confidence intervals containing that replicate's own true truncated exposure effect (see coverage\\_analysis/compute\\_true\\_psi.R), and average confidence interval width.}\n")
cat("\\label{tab:combined_tall_results_coverage_width}\n\\end{table}\n")

# ---- n_ac selection recovery vs. DGP ground truth (replaces the old,
# fold-based n_uc_hist_*.png) ----
# One row per sim (no outer folds any more). Pools every selected n_ac for
# IV-TPS/IV-GL, faceted by mechanism x option, with a red dashed line at the
# true confounded dimension n_c. Restricted to mechanisms 1,2,3,4,7,8 --
# mechanisms 5 (Leroux CAR) and 6 (coordinate-based) don't construct their
# confounder from a truncated basis at all, so there's no ground truth for
# these "recovery" plots to check against; they still appear in the other
# tables above (bias/RMSE/coverage/cor_Ac), just not here.
truth_mechanism_levels <- as.character(mechanisms_with_true_nc)
nac_rows <- list()
for (mech in mechanisms_with_true_nc) {
  for (opt in option_levels) {
    for (meth in iv_methods) {
      f <- file.path(results_dir, sprintf("conf%d_%s_%s_n_ac.csv", mech, opt, meth))
      if (!file.exists(f)) next
      n_ac <- unlist(read.csv(f), use.names = FALSE)
      nac_rows[[length(nac_rows) + 1]] <- data.frame(
        confounding_mechanism = mech, option = opt, method = meth, n_ac = n_ac
      )
    }
  }
}
nac_df <- if (length(nac_rows) > 0) bind_rows(nac_rows) else NULL

if (!is.null(nac_df) && nrow(nac_df) > 0) {
  nac_df$method[nac_df$method == "IV-GraphLaplacian"] <- "IV-GL"
  nac_df$confounding_mechanism <- factor(nac_df$confounding_mechanism, levels = truth_mechanism_levels)
  nac_df$option <- factor(nac_df$option, levels = option_levels)

  truth_lines <- data.frame(
    confounding_mechanism = factor(mechanisms_with_true_nc, levels = truth_mechanism_levels),
    true_n_c = unname(true_n_ac_by_mech)
  )

  # Panel 1: full distribution of selected n_ac, one histogram per
  # (mechanism, option), IV-TPS/IV-GL overlaid.
  #
  # Mechanisms 4 and 8 get their own plot: in their (reversed) search order
  # the selected n_ac sits near m - 20 (~894) rather than near 10, so a
  # shared x-range with the other mechanisms would be unreadable.
  main_mechs <- c("1", "2", "3", "7")
  nac_df_main <- nac_df[nac_df$confounding_mechanism %in% main_mechs, ] %>% droplevels()
  nac_df_reversed <- nac_df[nac_df$confounding_mechanism %in% c("4", "8"), ] %>% droplevels()
  truth_lines_main <- truth_lines[truth_lines$confounding_mechanism %in% main_mechs, ] %>% droplevels()
  truth_lines_reversed <- truth_lines[truth_lines$confounding_mechanism %in% c("4", "8"), ] %>% droplevels()

  nac_hist_theme <- theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 17),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      strip.text = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      legend.position = "top"
    )
  nac_hist_labs <- labs(
    x = expression(paste("Selected spatial basis dimension of confounding adjustment ", n[c])),
    y = "Count"
  )
  nac_hist_labeller <- labeller(confounding_mechanism = function(x) paste("Confounding \nmechanism", x))

  png(sprintf("images/n_ac_hist_%s.png", results_label), width = 2200, height = 1700, res = 200)
  print(
    ggplot(nac_df_main, aes(x = n_ac, fill = method)) +
      geom_bar(position = position_dodge(width = 3)) +
      geom_vline(data = truth_lines_main, aes(xintercept = true_n_c),
                 color = "red", linetype = "dashed", linewidth = 0.8) +
      facet_grid(option ~ confounding_mechanism, scales = "free", labeller = nac_hist_labeller) +
      scale_fill_manual(name = "Method", values = c("IV-TPS" = "#9ECAE1", "IV-GL" = "#FDAE6B")) +
      nac_hist_labs +
      nac_hist_theme
  )
  dev.off()

  # Panel 1b: mechanisms 4/8 (reversed-scale, mismatched-basis methods
  # included), own y-axis scale -- see comment above. Skipped when those
  # mechanisms were not run.
  if (nrow(nac_df_reversed) > 0) {
    png(sprintf("images/n_ac_hist_reversed_%s.png", results_label), width = 1400, height = 1700, res = 200)
    print(
      ggplot(nac_df_reversed, aes(x = n_ac, fill = method)) +
        geom_bar(position = position_dodge(width = 3)) +
        geom_vline(data = truth_lines_reversed, aes(xintercept = true_n_c),
                   color = "red", linetype = "dashed", linewidth = 0.8) +
        facet_grid(option ~ confounding_mechanism, scales = "free", labeller = nac_hist_labeller) +
        scale_fill_manual(name = "Method", values = c("IV-TPS" = "#9ECAE1", "IV-GL" = "#FDAE6B")) +
        nac_hist_labs +
        nac_hist_theme
    )
    dev.off()
  }
} else {
  cat("No n_ac files found yet -- skipping n_ac plots.\n")
}
