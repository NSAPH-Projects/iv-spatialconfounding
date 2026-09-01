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

source('../funcs.R')
load('sim.RData')

# ---- Configuration ----
# Point this at whatever results directory you've pulled back from the cluster.
# Safe to run against a still-running cluster job: any (mechanism, option,
# method) combination without a file yet just shows up as NA/blank below,
# nothing errors.
results_dir <- "results_manuscript_cluster_1000reps_quadratic_alpha1.5/"
results_label <- "cluster_1000reps_alpha1.5"
n_mechanisms <- 8
mechanism_levels <- as.character(1:n_mechanisms)
option_levels <- c("linear", "nonlinear")

dir.create("images", showWarnings = FALSE)

# Basis sizes, for computing n_max the same way make_candidate_grid() does
# (m - max(5, floor(0.02 * m))), so the n_uc diagnostics below know how close
# a selection came to the fragile endpoint for each basis.
m_tps <- ncol(simlist$B_tps_full)
m_gl  <- ncol(simlist$B_gl_full)
n_max_for <- function(method) {
  m <- if (grepl("TPS", method)) m_tps else m_gl
  m - max(5L, floor(0.02 * m))
}

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

  desired_order <- c("oracle", "baseline", "spatialcoord", "trueIV", "IV-TPS", "IV-GL")
  df$method <- factor(df$method, levels = desired_order)
  df$confounding_mechanism <- factor(df$confounding_mechanism, levels = mechanism_levels)
  df$option <- factor(df$option, levels = option_levels)

  mutrues_f <- mutrues %>%
    mutate(confounding_mechanism = factor(confounding_mechanism, levels = mechanism_levels),
          option = factor(option, levels = option_levels))

  method_cols <- c(
    "oracle"                 = "gray",
    "baseline"               = "#D62728",
    "spatialcoord"           = "purple",
    "trueIV"                 = "lightgreen",
    "IV-TPS"                 = "#9ECAE1",
    "IV-GL"                  = "#FDAE6B"
  )

  png(sprintf("images/boxplot_%s.png", results_label), width = 2500, height = 1250, res = 200)
  print(
    ggplot(df, aes(x = method, y = estimate, fill = method)) +
      geom_boxplot(alpha = 0.5, outliers = FALSE, staplewidth = 1) +
      #stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "blue") +
      ggh4x::facet_grid2(
        option ~ confounding_mechanism,
        scales = "free",
        independent = "all"
      ) +
      geom_hline(data = mutrues_f, aes(yintercept = theta),
                 color = "red", linetype = "twodash", size = 1) +
      labs(
        x = sprintf("Confounding mechanism (1-%d)", n_mechanisms),
        y = "Truncated Exposure Effect Estimate"
      ) +
      scale_fill_manual(
        name = "Method",
        values = method_cols,
        breaks = names(method_cols)
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top")
  )
  dev.off()
  
  dfsub <- df %>%
    filter(confounding_mechanism %in% 1:6) %>%
    droplevels()
  mutrues_f_sub <- mutrues_f %>%
    filter(confounding_mechanism %in% 1:6) %>%
    droplevels()
  png(sprintf("images/boxplot_%s_maintext.png", results_label), width = 1875, height = 1250, res = 200)
  print(
    ggplot(dfsub, aes(x = method, y = estimate, fill = method)) +
      geom_boxplot(alpha = 0.5, outliers = FALSE, staplewidth = 1) +
      #stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "blue") +
      ggh4x::facet_grid2(
        option ~ confounding_mechanism,
        scales = "free",
        independent = "all"
      ) +
      geom_hline(data = mutrues_f_sub, aes(yintercept = theta),
                 color = "red", linetype = "twodash", size = 1) +
      labs(
        x = sprintf("Confounding mechanism (1-%d)", 6),
        y = "Truncated Exposure Effect Estimate"
      ) +
      scale_fill_manual(
        name = "Method",
        values = method_cols,
        breaks = names(method_cols)
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top")
  )
  dev.off()
} else {
  cat("No point-estimate files found yet -- skipping boxplot.\n")
}

# ---- CI coverage and width ----
# Same scaffold pattern as the bias/RMSE table: every (mechanism, option,
# method) row is always present, NA where CI files aren't there yet.
ci_summary <- scaffold
ci_summary$coverage <- NA_real_
ci_summary$mean_width <- NA_real_
ci_summary$n <- 0L
for (i in seq_len(nrow(file_index))) {
  row <- file_index[i, ]
  lo_f <- file.path(results_dir, sub("\\.csv$", "_ci_lower.csv", row$csv))
  hi_f <- file.path(results_dir, sub("\\.csv$", "_ci_upper.csv", row$csv))
  if (!file.exists(lo_f) || !file.exists(hi_f)) next
  mutrue <- mutrues$theta[mutrues$confounding_mechanism == row$confounding_mechanism &
                            mutrues$option == row$option]
  if (length(mutrue) == 0 || is.na(mutrue)) next
  lo <- as.vector(as.matrix(read.csv(lo_f)))
  hi <- as.vector(as.matrix(read.csv(hi_f)))
  covers <- lo <= mutrue & hi >= mutrue
  idx <- which(ci_summary$confounding_mechanism == row$confounding_mechanism &
                ci_summary$option == row$option &
                ci_summary$method == row$method)
  ci_summary$coverage[idx] <- mean(covers, na.rm = TRUE)
  ci_summary$mean_width[idx] <- mean(hi - lo, na.rm = TRUE)
  ci_summary$n[idx] <- sum(!is.na(covers))
}
ci_summary <- ci_summary %>% arrange(confounding_mechanism, option, method)
print(xtable(ci_summary, digits = 3), include.rownames = FALSE)

# ---- n_uc selection stability (basis-selection methods only) ----
# frac_near_max is the share of fold-level selections landing within 10 of
# n_max for that method's basis -- the fragile endpoint region identified
# during the mechanism-2 pilot investigation. n_folds is nsims x 5 outer
# folds pooled together, not a literal count of distinct folds.
nuc_summary <- expand.grid(
  confounding_mechanism = 1:n_mechanisms,
  option = option_levels,
  method = iv_methods,
  stringsAsFactors = FALSE
)
nuc_summary$mean_n_uc <- NA_real_
nuc_summary$sd_n_uc <- NA_real_
nuc_summary$max_n_uc <- NA_integer_
nuc_summary$n_max <- vapply(nuc_summary$method, n_max_for, numeric(1))
nuc_summary$frac_near_max <- NA_real_
nuc_summary$n_folds <- 0L
for (i in seq_len(nrow(nuc_summary))) {
  mech <- nuc_summary$confounding_mechanism[i]
  opt <- nuc_summary$option[i]
  meth <- nuc_summary$method[i]
  fold_files <- file.path(results_dir, sprintf(
    "conf%d_%s_%s_n_uc_fold%d.csv", mech, opt, meth, 1:5
  ))
  if (!all(file.exists(fold_files))) next
  nuc <- unlist(lapply(fold_files, function(f) unlist(read.csv(f), use.names = FALSE)))
  nmax <- nuc_summary$n_max[i]
  nuc_summary$mean_n_uc[i] <- mean(nuc)
  nuc_summary$sd_n_uc[i] <- sd(nuc)
  nuc_summary$max_n_uc[i] <- max(nuc)
  nuc_summary$frac_near_max[i] <- mean(nuc >= nmax - 10)
  nuc_summary$n_folds[i] <- length(nuc)
}
nuc_summary <- nuc_summary %>% arrange(confounding_mechanism, option, method)
print(xtable(nuc_summary, digits = 2), include.rownames = FALSE)

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
coverage_labels  <- c("Oracle", "Baseline", "Spatial\\\\coordinates", "trueIV", "IV-TPS", "IV-GL")

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
cat("    \\caption{Average coverage, defined as the proportion of confidence intervals containing the oracle mean, and average confidence interval width.}\n")
cat("\\label{tab:combined_tall_results_coverage_width}\n\\end{table}\n")

# ---- n_uc selection histogram (replaces the old n_uc_hist_Mar27.png) ----
# Pools every fold-level n_uc selection across all reps for IV-TPS/IV-GL,
# faceted by mechanism x option
nuc_rows <- list()
for (mech in 1:n_mechanisms) {
  for (opt in option_levels) {
    for (meth in iv_methods) {
      fold_files <- file.path(results_dir, sprintf(
        "conf%d_%s_%s_n_uc_fold%d.csv", mech, opt, meth, 1:5
      ))
      if (!all(file.exists(fold_files))) next
      nuc <- unlist(lapply(fold_files, function(f) unlist(read.csv(f), use.names = FALSE)))
      nuc_rows[[length(nuc_rows) + 1]] <- data.frame(
        confounding_mechanism = mech, option = opt, method = meth, n_uc = nuc
      )
    }
  }
}
nuc_df <- if (length(nuc_rows) > 0) bind_rows(nuc_rows) else NULL

if (!is.null(nuc_df) && nrow(nuc_df) > 0) {
  nuc_df$method[nuc_df$method == "IV-GraphLaplacian"] <- "IV-GL"
  nuc_df$confounding_mechanism <- factor(nuc_df$confounding_mechanism, levels = mechanism_levels)
  nuc_df$option <- factor(nuc_df$option, levels = option_levels)

  png(sprintf("images/n_uc_hist_%s.png", results_label), width = 2500, height = 1250, res = 200)
  print(
    ggplot(nuc_df, aes(x = n_uc, fill = method)) +
      geom_bar(position = position_dodge(width = 3)) +
      facet_grid(option ~ confounding_mechanism, scales = "free") +
      scale_fill_manual(name = "Method", values = c("IV-TPS" = "#9ECAE1", "IV-GL" = "#FDAE6B")) +
      labs(
        x = expression(paste("Number of small-scale basis elements selected as instruments (", n[uc], ")")),
        y = "Count"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
  )
  dev.off()
} else {
  cat("No n_uc fold files found yet -- skipping n_uc histogram.\n")
}
