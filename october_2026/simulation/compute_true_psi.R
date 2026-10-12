# Per-replicate coverage truth, cluster CLI version.
#
# Companion to run_simfunc.R: submit_jobs.sh submits one of these per
# mechanism (see run_job_truth.sh), alongside the main cm*/option jobs, so
# that a single `bash submit_jobs.sh` produces both the point estimates and
# the per-replicate truth needed to check CI coverage correctly.
#
# Why per-replicate truth is needed (see this session's discussion): since
# funcs.R's DGM (2026-09-16) redraws the confounder field fresh for every
# replicate, each replicate has its own true population truncated-exposure
# -effect -- coverage can't validly be checked against one pooled mutrue
# (bias/RMSE against a pooled mutrue remain fine). This script reconstructs
# each replicate's true psi without any new simulation: replay
# set.seed(seed); compute_data_XXX(..., nsims) -- the exact call
# run_simfunc.R makes, so as long as submit_jobs.sh passes the same nsims/
# seed to both (it does, by construction -- see the new loop in
# submit_jobs.sh), this reproduces the exact same confounder field the main
# job used. Then for each replicate, estimate_truncated_effect_full() is
# run on a NOISELESS outcome (E[Y|A,U], no random error term) using that
# replicate's own realized A/U -- this reuses the trusted, already-in-use
# estimator machinery rather than hand-deriving the estimand, and with zero
# outcome noise the result is accurate to many decimal places.
#
# Validated locally this session (coverage_analysis/compute_true_psi.R,
# the no-CLI-args prototype this script is adapted from): reconstructed
# truth for mechanism 1's first 5 replicates landed within 0.3-2.2
# percentage points of the ACTUAL saved oracle estimates for those same
# replicates -- exactly oracle's own sampling noise, confirming the replay
# reconstructs the correct fields.
#
# Only needs funcs.R (compute_data_*/createY) and R/truncated_effect.R
# (estimate_truncated_effect_full) -- no selection rule is involved, since
# truth conditions directly on the true confounder.
#
# One job per MECHANISM (not per mechanism x option): field generation is
# shared across both options (option only changes the outcome formula
# applied afterward), so computing both options' truth in one job avoids
# redundant field reconstruction and halves the job count.

suppressMessages({
  library(dplyr)
  library(geosphere)
  library(Matrix)
  library(parallel)
})

source("funcs.R")
source("../R/truncated_effect.R")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4L || length(args) > 5L) {
  stop(
    paste(
      "Usage: Rscript compute_true_psi.R",
      "<mechanism> <nsims> <seed> <results_dir> [n_cores]"
    )
  )
}

parse_positive_integer <- function(value, name) {
  out <- suppressWarnings(as.integer(value))
  if (length(out) != 1L || is.na(out) || out < 1L) {
    stop(name, " must be a positive integer.")
  }
  out
}

mechanism <- parse_positive_integer(args[[1]], "mechanism")
if (!mechanism %in% 1:8) {
  stop("mechanism must be between 1 and 8.")
}
nsims <- parse_positive_integer(args[[2]], "nsims")
seed <- parse_positive_integer(args[[3]], "seed")
results_dir <- args[[4]]
n_cores <- if (length(args) >= 5L) parse_positive_integer(args[[5]], "n_cores") else 8L

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(results_dir)) {
  stop("Could not create results directory: ", results_dir)
}
results_dir <- paste0(normalizePath(results_dir, mustWork = TRUE), .Platform$file.sep)

message("Truth configuration:")
print(list(mechanism = mechanism, nsims = nsims, seed = seed,
           results_dir = results_dir, n_cores = n_cores))

load("sim.RData")

cutoff <- 0
n_grid <- 25L
sl_lib <- c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction")
density_trim <- 0.01
option_levels <- c("linear", "nonlinear")

distmat <- geosphere::distm(cbind(simlist$lon, simlist$lat), fun = geosphere::distHaversine)
distmat <- distmat / 1000000

# Same call (and seed) as simfunc(), so each replicate's field matches.
gen_field <- function(mech) {
  set.seed(seed)
  generate_mechanism_data(
    mech, nsims, B_tps_full = simlist$B_tps_full, B_gl_full = simlist$B_gl_full,
    W = simlist$W, lat = simlist$lat, lon = simlist$lon, distmat = distmat
  )
}

noiseless_Y_single <- function(U, A, option) {
  if (option == "linear") -2 - 2 * U + A - A * U
  else -2 - 2 * U + A - 0.4 * A^2 - 0.5 * U * A
}
noiseless_Y_two <- function(U1, U2, A, option) {
  if (option == "linear") -2 - 2 * U1 + A - A * U1 - 1.5 * A * U2
  else -2 - 2 * U1 + A - 0.4 * A^2 - 0.5 * U1 * A - 0.5 * U2 * A
}

compute_truth_for_rep <- function(mech, dat, i, option) {
  two_conf <- mech %in% c(3, 7)
  A <- dat$Ac[, i] + dat$Auc[, i]
  if (two_conf) {
    U1 <- dat$U1[, i]; U2 <- dat$U2[, i]
    Y <- noiseless_Y_two(U1, U2, A, option)
    w <- cbind(U1 = U1, U2 = U2)
  } else {
    U <- dat$U[, i]
    Y <- noiseless_Y_single(U, A, option)
    w <- matrix(U, ncol = 1, dimnames = list(NULL, "U"))
  }
  fit <- tryCatch(
    estimate_truncated_effect_full(
      y = Y, a = A, w = w, cutoff = cutoff,
      nuisance_args = list(sl_library = sl_lib, density_trim = density_trim),
      n_grid = n_grid, constrain = TRUE, density_trim = density_trim
    ),
    error = function(e) {
      message("  truth fit failed mech=", mech, " option=", option, " rep=", i, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) NA_real_ else fit$psi
}

message("Reconstructing confounder field for mechanism ", mechanism, "...")
dat <- gen_field(mechanism)

tasks <- list()
for (option in option_levels) {
  for (i in seq_len(nsims)) {
    tasks[[length(tasks) + 1]] <- list(option = option, rep = i)
  }
}

message("Computing true psi for ", length(tasks), " (option x replicate) cells across ", n_cores, " cores...")
results <- mclapply(tasks, function(task) {
  psi_true <- compute_truth_for_rep(mechanism, dat, task$rep, task$option)
  data.frame(mechanism = mechanism, option = task$option, replicate = task$rep, psi_true = psi_true)
}, mc.cores = n_cores)

truth_df <- bind_rows(results)
out_file <- file.path(results_dir, sprintf("true_psi_mech%d.rds", mechanism))
saveRDS(truth_df, out_file)
message("\nSaved ", out_file, " (", nrow(truth_df), " rows)")

# ---- best-effort sanity check: compare mean(psi_true) per option to the
# pooled oracle mutrue, if the main job's oracle CSVs happen to exist
# already (they may not, since this runs concurrently with the main jobs
# rather than depending on them -- that's fine, this check is skipped
# gracefully if the files aren't there yet). ----
for (option in option_levels) {
  oracle_f <- file.path(results_dir, sprintf("conf%d_%s_oracle.csv", mechanism, option))
  if (!file.exists(oracle_f)) next
  oracle_mean <- mean(as.vector(as.matrix(read.csv(oracle_f))), na.rm = TRUE)
  sub <- truth_df[truth_df$option == option, ]
  message(sprintf(
    "  sanity check [%s]: mean(psi_true)=%.4f  pooled_oracle_mutrue=%.4f",
    option, mean(sub$psi_true, na.rm = TRUE), oracle_mean
  ))
}
