library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)
library(purrr)

source('../funcs.R')
load('sim.RData')

# ---- Configuration ----
# Sensitivity to alpha, the tuning parameter of the sequential stability criterion 
# Each directory below is a full independent rerun of the cluster simulation at a fixed alpha for the two IV methods 
alpha_dirs <- c(
  "0.5" = "results_manuscript_cluster_1000reps_quadratic_alpha0.5/",
  "1"   = "results_manuscript_cluster_1000reps_quadratic_alpha1/",
  "1.5" = "results_manuscript_cluster_1000reps_quadratic_alpha1.5/",
  "2"   = "results_manuscript_cluster_1000reps_quadratic_alpha2/"
)
n_mechanisms <- 8
mechanism_levels <- as.character(1:n_mechanisms)
option_levels <- c("linear", "nonlinear")
iv_methods <- c("IV-TPS", "IV-GraphLaplacian")

dir.create("images", showWarnings = FALSE)

m_tps <- ncol(simlist$B_tps_full)
m_gl  <- ncol(simlist$B_gl_full)
n_max_for <- function(method) {
  m <- if (grepl("TPS", method)) m_tps else m_gl
  m - max(5L, floor(0.02 * m))
}

# ---- Per-alpha summary: bias / RMSE / coverage / width / mean n_uc ----
# Ground truth (oracle mean) is computed from each alpha directory's own
# oracle files -- same convention as analysis.R 
summarize_alpha_dir <- function(results_dir) {
  scaffold <- expand.grid(
    confounding_mechanism = 1:n_mechanisms,
    option = option_levels,
    method = iv_methods,
    stringsAsFactors = FALSE
  )

  csvs <- list.files(results_dir, pattern = "^conf[0-9]+_(linear|nonlinear)_[^_]+\\.csv$")
  file_index <- data.frame(
    confounding_mechanism = as.integer(sub("^conf(\\d+)_.*$", "\\1", csvs)),
    option = sub("^conf\\d+_([^_]+)_.*$", "\\1", csvs),
    method = sub("^conf\\d+_[^_]+_([^_]+)\\.csv$", "\\1", csvs),
    csv = csvs,
    stringsAsFactors = FALSE
  )

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

  out <- scaffold
  out$bias <- NA_real_
  out$RMSE <- NA_real_
  out$se <- NA_real_
  out$coverage <- NA_real_
  out$mean_width <- NA_real_
  out$mean_n_uc <- NA_real_
  out$sd_n_uc <- NA_real_

  iv_file_index <- file_index %>% filter(method %in% iv_methods)
  for (i in seq_len(nrow(iv_file_index))) {
    row <- iv_file_index[i, ]
    mutrue <- mutrues$theta[mutrues$confounding_mechanism == row$confounding_mechanism &
                              mutrues$option == row$option]
    if (length(mutrue) == 0 || is.na(mutrue)) next
    idx <- which(out$confounding_mechanism == row$confounding_mechanism &
                  out$option == row$option &
                  out$method == row$method)

    muests <- as.vector(as.matrix(read.csv(file.path(results_dir, row$csv))))
    out$bias[idx] <- mean(muests, na.rm = TRUE) - mutrue
    out$RMSE[idx] <- sqrt(mean((muests - mutrue)^2, na.rm = TRUE))
    out$se[idx] <- sd(muests, na.rm = TRUE)

    lo_f <- file.path(results_dir, sub("\\.csv$", "_ci_lower.csv", row$csv))
    hi_f <- file.path(results_dir, sub("\\.csv$", "_ci_upper.csv", row$csv))
    if (file.exists(lo_f) && file.exists(hi_f)) {
      lo <- as.vector(as.matrix(read.csv(lo_f)))
      hi <- as.vector(as.matrix(read.csv(hi_f)))
      covers <- lo <= mutrue & hi >= mutrue
      out$coverage[idx] <- mean(covers, na.rm = TRUE)
      out$mean_width[idx] <- mean(hi - lo, na.rm = TRUE)
    }

    fold_files <- file.path(results_dir, sprintf(
      "conf%d_%s_%s_n_uc_fold%d.csv", row$confounding_mechanism, row$option, row$method, 1:5
    ))
    if (all(file.exists(fold_files))) {
      nuc <- unlist(lapply(fold_files, function(f) unlist(read.csv(f), use.names = FALSE)))
      out$mean_n_uc[idx] <- mean(nuc)
      out$sd_n_uc[idx] <- sd(nuc)
    }
  }
  out
}

alpha_summary <- imap_dfr(alpha_dirs, function(dir, alpha_label) {
  cat(sprintf("Summarizing alpha = %s (%s)...\n", alpha_label, dir))
  summarize_alpha_dir(dir) %>% mutate(alpha = as.numeric(alpha_label))
})

alpha_summary$method[alpha_summary$method == "IV-GraphLaplacian"] <- "IV-GL"
alpha_summary$confounding_mechanism <- factor(alpha_summary$confounding_mechanism, levels = mechanism_levels)
alpha_summary$option <- factor(alpha_summary$option, levels = option_levels)

missing <- alpha_summary %>% filter(is.na(bias))
if (nrow(missing) > 0) {
  cat("\nCells with no point-estimate file yet (shown as gaps in the plots):\n")
  print(missing %>% select(alpha, confounding_mechanism, option, method), row.names = FALSE)
}

method_cols <- c("IV-TPS" = "#9ECAE1", "IV-GL" = "#FDAE6B")

# ---- Figure: bias / variance tradeoff vs alpha ----
# Three stacked metrics (|bias|, empirical SE, RMSE) x alpha, one column per
# confounding mechanism, colored by method, linestyle by outcome model.
# Empirical SE -- the Monte Carlo SD of the 1000 point estimates around their
# own mean -- is the actual variance term (RMSE^2 = bias^2 + SE^2), unlike
# mean CI width, which reflects the *estimated* SE from the variance
# estimator and can fail to track the real precision cost (e.g. coverage
# collapsing while width stays flat). RMSE is included as the net effect of
# the tradeoff.
tradeoff_df <- alpha_summary %>%
  mutate(abs_bias = abs(bias) * 100, se_mc = se * 100, RMSE_pct = RMSE * 100) %>%
  select(alpha, confounding_mechanism, option, method, abs_bias, se_mc, RMSE_pct) %>%
  pivot_longer(c(abs_bias, se_mc, RMSE_pct), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("abs_bias", "se_mc", "RMSE_pct"),
                          labels = c("|Bias| (x10^2)", "Empirical SE (x10^2)", "RMSE (x10^2)")))

png("images/alpha_tradeoff_bias_variance.png", width = 3000, height = 2000, res = 200)
print(
  ggplot(tradeoff_df, aes(x = alpha, y = value, color = method, linetype = option)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE, size = 1.2) +
    ggh4x::facet_grid2(metric ~ confounding_mechanism, scales = "free_y", independent = "y") +
    scale_color_manual(name = "Method", values = method_cols) +
    scale_linetype_discrete(name = "Outcome model") +
    scale_x_continuous(breaks = as.numeric(names(alpha_dirs))) +
    labs(
      x = expression(alpha),
      y = NULL
    ) +
    theme_bw() +
    theme(legend.position = "top")
)
dev.off()

# ---- Figure: mean n_uc selected vs alpha ----
# Mechanistic link behind the tradeoff above: larger alpha permits a larger
# instrumental set n_uc, which is what drives bias (and, in principle, SE) up.
nuc_df <- alpha_summary %>%
  filter(!is.na(mean_n_uc))

png("images/alpha_tradeoff_n_uc.png", width = 3000, height = 900, res = 200)
print(
  ggplot(nuc_df, aes(x = alpha, y = mean_n_uc, color = method, linetype = option)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE, size = 1.2) +
    geom_errorbar(aes(ymin = mean_n_uc - sd_n_uc, ymax = mean_n_uc + sd_n_uc),
                  width = 0, alpha = 0.3, na.rm = TRUE) +
    facet_wrap(~ confounding_mechanism, nrow = 1) +
    scale_color_manual(name = "Method", values = method_cols) +
    scale_linetype_discrete(name = "Outcome model") +
    scale_x_continuous(breaks = as.numeric(names(alpha_dirs))) +
    labs(
      x = expression(alpha),
      y = expression(paste("Mean selected ", n[uc]))
    ) +
    theme_bw() +
    theme(legend.position = "top")
)
dev.off()

# ---- Condensed table: alpha x mechanism x option, IV-TPS/IV-GL only ----
table_df <- alpha_summary %>%
  mutate(
    bias = round(bias * 100, 2),
    se = round(se * 100, 2),
    RMSE = round(RMSE * 100, 2),
    coverage = round(coverage, 3),
    mean_width = round(mean_width, 3),
    mean_n_uc = round(mean_n_uc, 1)
  ) %>%
  arrange(confounding_mechanism, option, method, alpha) %>%
  select(confounding_mechanism, option, method, alpha, bias, se, RMSE, coverage, mean_width, mean_n_uc)

write.csv(table_df, "images/alpha_tradeoff_table.csv", row.names = FALSE)

cat("\n% ---- alpha sensitivity table (IV-TPS / IV-GL only) ----\n")
print(xtable(table_df, digits = c(0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 1)),
      include.rownames = FALSE, sanitize.text.function = identity)

cat("\nWrote images/alpha_tradeoff_bias_variance.png, images/alpha_tradeoff_n_uc.png, ",
    "images/alpha_tradeoff_table.csv\n", sep = "")
