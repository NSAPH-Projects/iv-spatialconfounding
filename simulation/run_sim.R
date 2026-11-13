# Run the simulation for one or more mechanisms.
#
# Usage (from simulation/, after preprocessing.R):
#   Rscript run_sim.R [key=value ...]
# Keys (defaults in parentheses):
#   nsims (100), cores (8), results (results),
#   rep_start (1; replicates rep_start, ..., rep_start + nsims - 1, so a
#              long run can be split into chunks with identical datasets),
#   scenarios (all; comma-separated, e.g. LS,twoU),
#   methods (oracle,baseline,spatialcoord,trueIV; also IV-TPS, IV-GL),
#   options (linear,tanh; outcome models, see noiseless_y()),
#   tau (0.5; exposure R^2 reference for IV-TPS / IV-GL),
#   tag (all; appended to file names so jobs running different methods on
#        the same cell write separate files).
# Example: Rscript run_sim.R nsims=2 scenarios=LS methods=oracle,trueIV
# Writes <results>/<scenario>_<option>_<tag>.csv. Datasets depend only on
# the replicate number, so files from different method groups line up.

suppressMessages({
  library(parallel)
})
source("funcs.R")
source("../R/truncated_effect.R")
source("../R/basis_selection.R")

settings <- list(
  nsims = "100", rep_start = "1", cores = "8", results = "results",
  scenarios = paste(scenario_names, collapse = ","),
  methods = "oracle,baseline,spatialcoord,trueIV",
  options = "linear,tanh",
  tag = "all", tau = "0.5"
)
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L || !kv[1] %in% names(settings)) stop("Unknown argument: ", a)
  settings[[kv[1]]] <- kv[2]
}
nsims <- as.integer(settings$nsims)
rep_start <- as.integer(settings$rep_start)
reps <- rep_start - 1L + seq_len(nsims)
n_cores <- as.integer(settings$cores)
results_dir <- settings$results
scenario_names <- strsplit(settings$scenarios, ",")[[1]]
methods <- strsplit(settings$methods, ",")[[1]]
options <- strsplit(settings$options, ",")[[1]]
tau <- as.numeric(settings$tau)
stopifnot(!is.na(nsims), !is.na(rep_start), rep_start >= 1, !is.na(n_cores),
          !is.na(tau), tau > 0, tau < 1, all(scenario_names %in% names(scenarios)))

n <- 1000L
cutoff <- 0
seed <- 20261013L

sim <- readRDS("sim.rds")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

iv_methods <- c("IV-TPS", "IV-GL")

na_iv <- c(n_ac = NA_real_, n_ac0 = NA_real_, n_fits = NA_real_,
           cor_ac = NA_real_, cor_auc = NA_real_)
na_fit <- c(psi = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, na_iv)

one_rep <- function(rep_id, scenario, option) {
  set.seed(seed + rep_id)
  points <- sample_us_points(sim$outline, n)
  d <- draw_scenario_data(sim, scenario, points)
  coords <- points / 1000  # thousands of km
  y <- noiseless_y(d, d$A, option) + rnorm(n)
  rows <- lapply(methods, function(m) {
    t0 <- proc.time()[["elapsed"]]
    # Any failure (including basis construction) gives an NA row with the
    # error message, rather than an error object in the results.
    error <- ""
    out <- tryCatch(
      if (m %in% iv_methods) {
        b <- if (m == "IV-TPS") "tps" else "gl"
        B <- if (b == "tps") build_tps_basis_full(coords) else build_gl_basis_full(coords)
        fit_iv(d, y, B, cutoff, tau = tau)
      } else {
        c(fit_method(d, y, m, cutoff), na_iv)
      },
      error = function(e) {
        error <<- conditionMessage(e)
        na_fit
      }
    )
    data.frame(rep = rep_id, method = m, t(out), seconds = proc.time()[["elapsed"]] - t0,
               error = error)
  })
  do.call(rbind, rows)
}

for (scenario in scenario_names) {
  for (option in options) {
    t0 <- Sys.time()
    psi0 <- true_psi(sim, scenario, option, cutoff)
    out <- mclapply(reps, one_rep, scenario = scenario,
                    option = option, mc.cores = n_cores)
    # A worker that crashed returns an error object instead of a data frame.
    crashed <- !vapply(out, is.data.frame, logical(1))
    if (any(crashed)) {
      stop(sprintf("%s / %s: %d replicate(s) crashed, e.g. %s", scenario, option, sum(crashed),
                   as.character(out[[which(crashed)[1]]])))
    }
    res <- do.call(rbind, out)
    if (any(res$error != "")) {
      warning(sprintf("%s / %s: %d fit(s) failed, e.g. %s", scenario, option,
                      sum(res$error != ""), res$error[res$error != ""][1]))
    }
    res$scenario <- scenario
    res$option <- option
    res$psi_true <- psi0
    res$tau <- tau
    write.csv(res, file.path(results_dir, sprintf("%s_%s_%s.csv", scenario, option, settings$tag)),
              row.names = FALSE)
    message(sprintf("%s / %s: true psi %.4f, %d reps in %.1f min", scenario, option, psi0, nsims,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
