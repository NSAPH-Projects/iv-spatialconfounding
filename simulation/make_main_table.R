# Main-text table: RMSE (x 10^2) of every method in every setting, from
# <results_dir>/summary.csv (written by analysis.R).
#
# Rows: Mechanisms 1-5, each with one row per outcome model (linear, then
# non-linear) in an "Outcome model" column; the mechanism number is shown on
# its first row only. Columns: oracle, true IV, baseline, spatial coordinates,
# IV-TPS, IV-GL. Settings missing from summary.csv are left out. The caption
# states 1000 datasets (the final run), whatever results_dir holds.
#
# Usage (from simulation/):
#   Rscript make_main_table.R [results_dir] [label] [options]
#   options: linear, tanh, or both (default both)
#   defaults: results_final, final, both
# Writes tables/main_rmse_<label>.tex (a table environment for \input; needs
# only booktabs, which the main text loads).

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1L) args[[1]] else "results_final"
label <- if (length(args) >= 2L) args[[2]] else "final"
options_arg <- if (length(args) >= 3L) args[[3]] else "both"
options <- if (options_arg == "both") c("linear", "tanh") else options_arg
option_names <- c(linear = "linear", tanh = "non-linear")

summ <- read.csv(file.path(results_dir, "summary.csv"))

settings <- c(LS = "1", MS = "2", twoU = "3", coordnoise = "4", stateJump = "5")
methods <- c("oracle", "trueIV", "baseline", "spatialcoord", "IV-TPS", "IV-GL")

rmse_of <- function(sc, m, o) {
  v <- summ$rmse[summ$scenario == sc & summ$method == m & summ$option == o]
  if (length(v) == 0L) NA_real_ else 100 * v
}

present <- names(settings)[names(settings) %in% summ$scenario]
row_for <- function(sc, o, first) {
  v <- vapply(methods, function(m) rmse_of(sc, m, o), numeric(1))
  sprintf("%s & %s & %s \\\\", if (first) settings[[sc]] else "", option_names[[o]],
          paste(ifelse(is.na(v), "--", sprintf("%.2f", v)), collapse = " & "))
}
body <- unlist(lapply(present, function(sc) {
  vapply(seq_along(options), function(k) row_for(sc, options[k], k == 1L), character(1))
}))

caption_scope <- if (length(options) == 1L) {
  sprintf(" under the %s outcome model", option_names[[options]])
} else ""
out <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\begin{tabular}{llcccccc}",
  "\\toprule",
  "Confounding Mechanism & Outcome model & Oracle & True IV & Baseline & Spatial coordinates & IV-TPS & IV-GL \\\\",
  "\\midrule",
  body,
  "\\bottomrule",
  "\\end{tabular}",
  sprintf(paste0(
    "\\caption{Root mean squared error ($\\times 10^2$) of the truncated exposure effect estimates ",
    "across $1000$ simulated datasets%s.}"), caption_scope),
  "\\label{tab:main-rmse}",
  "\\end{table}"
)

dir.create("tables", showWarnings = FALSE)
path <- sprintf("tables/main_rmse_%s.tex", label)
writeLines(out, path)
message("Wrote ", path, " (", length(present), " settings: ", paste(present, collapse = ", "),
        "; outcome models: ", paste(options, collapse = ", "), ")")
