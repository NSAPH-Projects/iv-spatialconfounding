# Bias, RMSE, coverage, and interval width for the simulation.
#
# Usage: Rscript analysis.R [results_dir] [label] [scenarios]
#   scenarios: comma-separated, in display order (default: Mechanisms 1-5,
#              LS,MS,twoU,coordnoise,stateJump). Every output, including
#              summary.csv, covers only these.
# Writes, using label (default "final"):
#   <results_dir>/summary.csv
#   tables/results_<label>.tex   (bias/RMSE and coverage/width tables, in the
#                                 format of the supplement's simulation tables)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) args[[1]] else "results_final"
label <- if (length(args) >= 2L) args[[2]] else "final"

# Supplement numbering (fixed, whatever subset is shown).
all_scenario_labels <- c(LS = "Mechanism 1", MS = "Mechanism 2", twoU = "Mechanism 3",
                         coordnoise = "Mechanism 4", stateJump = "Mechanism 5")
scenario_levels <- if (length(args) >= 3L) strsplit(args[[3]], ",")[[1]] else
  c("LS", "MS", "twoU", "coordnoise", "stateJump")
stopifnot(all(scenario_levels %in% names(all_scenario_labels)))
scenario_labels <- all_scenario_labels[scenario_levels]
option_levels <- c("linear", "tanh")
option_labels <- c(linear = "linear", tanh = "nonlinear")
method_levels <- c("oracle", "baseline", "spatialcoord", "trueIV", "IV-TPS", "IV-GL")

keep <- c("rep", "method", "psi", "ci_lower", "ci_upper", "scenario", "option", "psi_true",
          "n_ac", "n_ac0", "n_fits", "cor_ac", "cor_auc", "seconds")
# <scenario>_<option>.csv (single job), <scenario>_<option>_<group>.csv, or
# <scenario>_<option>_<group>_r<first rep>.csv (chunked runs).
files <- list.files(results_dir, pattern = "^[A-Za-z0-9]+_[a-z]+(_[a-z0-9]+)*\\.csv$", full.names = TRUE)
res <- do.call(rbind, lapply(files, function(f) {
  d <- read.csv(f)
  for (k in setdiff(keep, names(d))) d[[k]] <- NA_real_
  d[, keep]
}))
res <- res[res$scenario %in% scenario_levels & res$option %in% option_levels &
             res$method %in% method_levels, ]
res$error <- res$psi - res$psi_true
res$covers <- res$ci_lower <= res$psi_true & res$ci_upper >= res$psi_true
res$width <- res$ci_upper - res$ci_lower

# ---- Summary ----
summary_rows <- lapply(split(res, list(res$scenario, res$option, res$method), drop = TRUE), function(d) {
  d <- d[!is.na(d$psi), ]
  data.frame(
    scenario = d$scenario[1],
    option = d$option[1],
    method = d$method[1],
    n_ok = nrow(d),
    psi_true = d$psi_true[1],
    bias = mean(d$error),
    mc_se_bias = sd(d$error) / sqrt(nrow(d)),
    rmse = sqrt(mean(d$error^2)),
    error_sd_over_se = sd(d$error) / mean(d$width / (2 * qnorm(0.975))),
    coverage = mean(d$covers),
    mean_width = mean(d$width),
    median_n_ac = stats::median(d$n_ac),
    median_n_ac0 = stats::median(d$n_ac0),
    mean_fits = mean(d$n_fits),
    median_cor_ac = stats::median(d$cor_ac, na.rm = TRUE),
    median_cor_auc = stats::median(d$cor_auc, na.rm = TRUE),
    mean_seconds = mean(d$seconds)
  )
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[order(match(summary_table$scenario, scenario_levels),
                                     match(summary_table$option, option_levels),
                                     match(summary_table$method, method_levels)), ]
print(summary_table, digits = 3, row.names = FALSE)
write.csv(summary_table, file.path(results_dir, "summary.csv"), row.names = FALSE)

# ---- LaTeX tables ----
nsims <- max(summary_table$n_ok)
panel_rows <- function(value_col, fmt, scale = 1) {
  out <- character(0)
  for (sc in scenario_levels) for (op in option_levels) {
    vals <- vapply(method_levels, function(m) {
      v <- summary_table[[value_col]][summary_table$scenario == sc &
                                        summary_table$option == op &
                                        summary_table$method == m]
      if (length(v) == 0L || is.na(v)) "" else sprintf(fmt, scale * v)
    }, character(1))
    out <- c(out, sprintf("%s & %s & %s \\\\ ", sub("Mechanism ", "", scenario_labels[[sc]]), option_labels[[op]],
                          paste(vals, collapse = " & ")))
  }
  out
}
summary_table$abs_bias <- abs(summary_table$bias)

header <- c(
  "\\makecell[c]{Confounding\\\\mechanism} &",
  "\\makecell[c]{Outcome\\\\model} &",
  "\\multicolumn{1}{c}{Oracle} &",
  "\\multicolumn{1}{c}{Baseline} &",
  "\\multicolumn{1}{c}{\\makecell{Spatial\\\\coordinates}} &",
  "\\multicolumn{1}{c}{True IV} &",
  "\\multicolumn{1}{c}{IV-TPS} &",
  "\\multicolumn{1}{c}{IV-GL} \\\\"
)

bias_table <- c(
  "\\begin{table}",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\renewcommand{\\arraystretch}{1.15}",
  "",
  "\\begin{tabular}{@{}cc*{6}{S[table-format=2.2]}@{}}",
  "\\toprule",
  header,
  "\\midrule",
  "",
  "\\multicolumn{8}{l}{\\textit{Absolute bias}} \\\\",
  "\\cmidrule(lr){1-8}",
  panel_rows("abs_bias", "%.2f", 100),
  "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Root mean squared error (RMSE)}} \\\\",
  "\\cmidrule(lr){1-8}",
  panel_rows("rmse", "%.2f", 100),
  "\\bottomrule",
  "\\end{tabular}",
  sprintf("    \\caption{Absolute bias and root mean squared error of the truncated exposure effect estimates across $%d$ simulated datasets for each of the %d data-generating scenarios, relative to the true truncated exposure effect $\\psi$. All values have been multiplied by $10^2.$}", nsims, 2L * length(scenario_levels)),
  "\\label{tab:combined_tall_results}",
  "\\end{table}"
)

coverage_table <- c(
  "\\begin{table}",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\renewcommand{\\arraystretch}{1.15}",
  "",
  "\\begin{tabular}{@{}cc*{6}{S[table-format=1.3]}@{}}",
  "\\toprule",
  header,
  "\\midrule",
  "",
  "\\multicolumn{8}{l}{\\textit{Average coverage}} \\\\",
  "\\cmidrule(lr){1-8}",
  panel_rows("coverage", "%.3f"),
  "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Average confidence interval width}} \\\\",
  "\\cmidrule(lr){1-8}",
  panel_rows("mean_width", "%.3f"),
  "\\bottomrule",
  "\\end{tabular}",
  sprintf("    \\caption{Average coverage of the 95\\%% confidence intervals for the true truncated exposure effect $\\psi$, and average confidence interval width, across $%d$ simulated datasets.}", nsims),
  "\\label{tab:combined_tall_results_coverage_width}",
  "\\end{table}"
)

dir.create("tables", showWarnings = FALSE)
writeLines(c(bias_table, "", coverage_table), sprintf("tables/results_%s.tex", label))
message("Wrote ", file.path(results_dir, "summary.csv"), " and tables/results_", label, ".tex")
