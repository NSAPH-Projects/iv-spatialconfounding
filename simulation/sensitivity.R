# Sensitivity of the truncated exposure effect estimate to the dimension of
# the confounding adjustment, with the sequential stability rule overlaid.
#
# For each mechanism, outcome model, basis (TPS, GL), and replicate, this runs
# select_basis_dimension(full_trajectory = TRUE) from ../R/basis_selection.R,
# which fits the estimate at every n_ac from the R^2 reference down to 1 and
# also returns the rule's selection. Datasets are those of run_sim.R (same
# seed per replicate). Each trajectory shows the reference n_uc^0 where the
# walk starts (open triangle) and the selected n_uc (diamond); the dashed red
# line is the smallest dimension whose basis captures 95% of the variance of
# A_c's spatial confounded part (U, or 0.7 U1 + 0.7 U2 in twoU), and the
# dotted grey line the true psi. The x axis is n_uc = m - n_ac, so the walk
# runs left to right. A second plot shows, along the same trajectories,
# corr(estimated A_c(n_ac), U's part of A_c), as run_sim.R's cor_ac.
#
# Usage (from simulation/):
#   Rscript sensitivity.R [key=value ...]
#   keys: reps (5), cores (SLURM_CPUS_PER_TASK or 8), tau (0.5),
#         scenarios (all), options (linear,tanh),
#         out (images/sensitivity; path prefix of the outputs)
# On the cluster: sbatch run_sensitivity.sh
# Output: <out>.rds (the trajectories), <out>.png, <out>_correlation.png.
# To redraw from the .rds without rerunning: Rscript plot_sensitivity.R <out>.rds

suppressMessages({
  library(ggplot2)
  library(parallel)
})
source("funcs.R")
source("../R/truncated_effect.R")
source("../R/basis_selection.R")

settings <- list(
  reps = "5", cores = Sys.getenv("SLURM_CPUS_PER_TASK", "8"), tau = "0.5",
  scenarios = paste(names(scenarios), collapse = ","),
  options = "linear,tanh", out = "images/sensitivity"
)
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L || !kv[1] %in% names(settings)) stop("Unknown argument: ", a)
  settings[[kv[1]]] <- kv[2]
}
n_reps <- as.integer(settings$reps)
n_cores <- as.integer(settings$cores)
tau <- as.numeric(settings$tau)
scenario_names <- strsplit(settings$scenarios, ",")[[1]]
options_ <- strsplit(settings$options, ",")[[1]]

n <- 1000L
cutoff <- 0
seed <- 20261013L  # as in run_sim.R, so replicate r is the same dataset
stopifnot(all(scenario_names %in% names(scenarios)))
sim <- readRDS("sim.rds")

jobs <- expand.grid(rep = seq_len(n_reps), basis = c("TPS", "GL"), option = options_,
                    scenario = scenario_names, stringsAsFactors = FALSE)

run_job <- function(j) {
  job <- jobs[j, ]
  set.seed(seed + job$rep)
  b <- tolower(job$basis)
  points <- sample_us_points(sim$outline, n)
  d <- draw_scenario_data(sim, job$scenario, points)
  y <- noiseless_y(d, d$A, job$option) + rnorm(n)
  coords <- points / 1000
  B <- if (b == "tps") build_tps_basis_full(coords) else build_gl_basis_full(coords)
  sel <- select_basis_dimension(
    y = y, a = d$A, x = d$X, B = B, cutoff = cutoff,
    effect_estimator = estimate_truncated_effect_full, effect_args = effect_settings,
    r2_reference = tau, family_alpha = 0.05, seed = 12001L, full_trajectory = TRUE
  )
  target <- confounded_part(d)
  r2_target <- exposure_r2_curve(B, target)
  # Split of A at every n_ac on the trajectory: A^c(q) = projection of A on
  # the first q columns, via one QR (the first q columns of Q span them).
  ns <- sel$diagnostics$n_ac
  Q <- qr.Q(qr(B))[, seq_len(max(ns)), drop = FALSE]
  coef <- drop(crossprod(Q, d$A))
  Ac_hat <- t(apply(sweep(Q, 2, coef, "*"), 1, cumsum))[, ns, drop = FALSE]
  cor_ac <- vapply(seq_along(ns), function(k) safe_cor(Ac_hat[, k], target), numeric(1))
  cor_auc <- vapply(seq_along(ns), function(k) safe_cor(d$A - Ac_hat[, k], d$Auc), numeric(1))
  data.frame(
    scenario = job$scenario, option = job$option, basis = job$basis, rep = job$rep,
    m = ncol(B), n_ac = sel$diagnostics$n_ac, psi = sel$diagnostics$psi,
    selected = sel$selected_n_ac, reference = sel$settings$n_ac0,
    n_ac_95 = which(r2_target >= 0.95)[1],
    cor_ac = cor_ac, cor_auc = cor_auc
  )
}

res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(j) {
  tryCatch(run_job(j), error = function(e) {
    message(sprintf("job %d (%s %s %s rep %d) failed: %s", j, jobs$scenario[j], jobs$option[j],
                    jobs$basis[j], jobs$rep[j], conditionMessage(e)))
    NULL
  })
}, mc.cores = n_cores, mc.preschedule = FALSE))

truth <- do.call(rbind, lapply(scenario_names, function(sc) do.call(rbind, lapply(options_, function(o) {
  psi <- true_psi(sim, sc, o, cutoff)
  data.frame(scenario = sc, option = o, psi_true = psi)
}))))

dir.create(dirname(settings$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(trajectories = res, truth = truth, tau = tau), paste0(settings$out, ".rds"))

# ---- Summary and plot ----
summ <- unique(res[, c("scenario", "option", "basis", "rep", "m", "reference", "selected", "n_ac_95")])
summ$reference_n_uc <- summ$m - summ$reference
summ$selected_n_uc <- summ$m - summ$selected
summ$n_uc_95 <- summ$m - summ$n_ac_95
print(summ[, c("scenario", "option", "basis", "rep", "reference_n_uc", "selected_n_uc",
               "n_uc_95")], row.names = FALSE)

rds_path <- paste0(settings$out, ".rds")
png_path <- paste0(settings$out, ".png")
plot_scenarios <- scenario_names
# A plotting error must not fail the job: the .rds is already saved and can
# be plotted locally with plot_sensitivity.R.
tryCatch(source("plot_sensitivity.R"), error = function(e) {
  message("Plotting failed (", conditionMessage(e), "); the trajectories are in ", rds_path,
          ". Run: Rscript plot_sensitivity.R ", rds_path)
})
