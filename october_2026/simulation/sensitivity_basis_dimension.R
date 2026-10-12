# Sensitivity of the truncated exposure effect estimate to the dimension n_c
# of the confounding adjustment, with the sequential stability rule overlaid.
#
# For mechanisms 1, 2, 3, 7 (matched basis: TPS for 1 and 3, GL for 2 and 7),
# both outcome models, and n_reps replicates, this fits the estimate at every
# n_c = 1, ..., max(50, reference) and plots the trajectories. Each replicate
# shows:
#   - the selected n_c (diamond), obtained by applying the rule in
#     ../R/basis_selection.R to these fits: same reference (smallest B^c with
#     exposure R^2 >= r2_reference), same one-column walk, same Bonferroni
#     test, same seed per fit, so it matches select_basis_dimension();
#   - the R^2 reference where the walk starts (open triangle);
#   - the true n_c = 10 (dashed red line).
#
# Run from october_2026/simulation/:  Rscript sensitivity_basis_dimension.R
# Output: images/sensitivity_basis_dimension.png (+ .rds with the fits).

suppressMessages({
  library(ggplot2)
  library(Matrix)
  library(parallel)
})
source("funcs.R")
source("../R/basis_selection.R")
source("../R/truncated_effect.R")
load("sim.RData")

n_reps       <- 5L
max_nc       <- 50L     # fit every n_c from 1 to max(max_nc, reference)
true_nc      <- 10L
r2_reference <- 0.7
family_alpha <- 0.05
fit_seed     <- 12001L   # same per-fit seed as select_basis_dimension()
n_cores      <- 8L
options_     <- c("linear", "nonlinear")

mechanisms <- list(
  list(id = 1, method = "IV-TPS", B = simlist$B_tps_full),
  list(id = 2, method = "IV-GL",  B = simlist$B_gl_full),
  list(id = 3, method = "IV-TPS", B = simlist$B_tps_full),
  list(id = 7, method = "IV-GL",  B = simlist$B_gl_full)
)

# Estimator settings: simfunc()'s defaults.
effect_args <- list(
  nuisance_args = list(sl_library = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction"),
                       density_trim = 0.01),
  n_grid = 25L, density_trim = 0.01, constrain = TRUE
)
intercept <- matrix(1, nrow(simlist$B_tps_full), 1, dimnames = list(NULL, "Intercept"))

# Outcome, as in simfunc(): two confounders for mechanisms 3 and 7.
make_y <- function(mech, dat, i, option) {
  A <- dat$Ac[, i] + dat$Auc[, i]
  n <- length(A)
  if (mech %in% c(3, 7)) {
    U1 <- dat$U1[, i]; U2 <- dat$U2[, i]
    eta <- if (option == "linear") -2 - 2*U1 + A - A*U1 - 1.5*A*U2
           else -2 - 2*U1 + A - 0.4*A^2 - 0.5*U1*A - 0.5*U2*A
    rnorm(n, eta, 1)
  } else {
    createY(matrix(dat$U[, i]), matrix(A), option)[, 1]
  }
}

# ---- Simulate data (fixed seed per mechanism) and list the jobs ----
jobs <- list()
for (mc in mechanisms) {
  set.seed(20261012 + mc$id)
  dat <- generate_mechanism_data(
    mc$id, n_reps, B_tps_full = simlist$B_tps_full, B_gl_full = simlist$B_gl_full,
    W = simlist$W, lat = simlist$lat, lon = simlist$lon, distmat = NULL
  )
  for (option in options_) for (i in seq_len(n_reps)) {
    jobs[[length(jobs) + 1]] <- list(
      mech = mc$id, method = mc$method, B = mc$B, option = option, rep = i,
      a = dat$Ac[, i] + dat$Auc[, i], y = make_y(mc$id, dat, i, option)
    )
  }
}

# ---- Fit every n_c and apply the selection rule ----
run_job <- function(job) {
  r2 <- exposure_r2_curve(job$B, job$a)
  n0 <- which(r2 >= r2_reference)[1]            # reference (walk starts here)
  fits <- lapply(seq_len(max(max_nc, n0)), function(q) {
    set.seed(fit_seed)
    w <- cbind(Ac = project_Ac_full(job$B, job$a, q)$Ac, intercept)
    do.call(estimate_truncated_effect_full,
            c(list(y = job$y, a = job$a, w = w, cutoff = 0), effect_args))
  })

  # Same rule as select_basis_dimension(): walk down from n0, stop at the
  # first Bonferroni departure, keep the previous n_c; n_c = 1 if none.
  J <- n0 - 1L
  z <- qnorm(1 - family_alpha / (2 * J))
  selected <- 1L
  prev <- n0
  for (q in rev(seq_len(J))) {
    d  <- fits[[q]]$psi - fits[[n0]]$psi
    se <- paired_difference_se(fits[[q]]$influence_curve, fits[[n0]]$influence_curve)
    if (abs(d) > z * se) { selected <- prev; break }
    prev <- q
  }

  data.frame(
    mechanism = job$mech, method = job$method, option = job$option, rep = job$rep,
    n_c = seq_along(fits), psi = vapply(fits, `[[`, numeric(1), "psi"),
    exposure_R2 = r2[seq_along(fits)], reference = n0, selected = selected
  )
}

res <- do.call(rbind, mclapply(jobs, function(j) {
  tryCatch(suppressMessages(run_job(j)), error = function(e) {
    message(sprintf("mechanism %d %s rep %d failed: %s", j$mech, j$option, j$rep, conditionMessage(e)))
    NULL
  })
}, mc.cores = n_cores, mc.preschedule = FALSE))

dir.create("images", showWarnings = FALSE)
saveRDS(res, "images/sensitivity_basis_dimension.rds")

# ---- Plot ----
res$rep <- factor(res$rep)
res$option <- factor(res$option, levels = options_)
res$panel <- factor(
  sprintf("Confounding\nmechanism %d (%s)", res$mechanism, res$method),
  levels = sprintf("Confounding\nmechanism %d (%s)",
                   sapply(mechanisms, `[[`, "id"), sapply(mechanisms, `[[`, "method"))
)
selected_pts  <- subset(res, n_c == selected)
reference_pts <- subset(res, n_c == reference)

# Lines run from n_c = 1 to max(50, reference) for each replicate.
p <- ggplot(res, aes(n_c, psi, colour = rep, group = rep)) +
  geom_vline(xintercept = true_nc, colour = "red", linetype = "dashed") +
  geom_line(alpha = 0.8) +
  geom_point(size = 0.6) +
  geom_point(data = reference_pts, shape = 2, size = 2.5, stroke = 0.9, show.legend = FALSE) +
  geom_point(data = selected_pts, shape = 18, size = 4, show.legend = FALSE) +
  facet_grid(option ~ panel, scales = "free_y") +
  scale_colour_brewer(palette = "Set1", name = "Replicate") +
  labs(x = expression(paste("Spatial basis dimension of confounding adjustment ", n[c])),
       y = "Truncated exposure effect estimate",
       caption = paste("Diamond: selected n_c. Triangle: R\u00b2 reference where the walk starts.",
                       "Dashed line: true n_c = 10.")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"))

ggsave("images/sensitivity_basis_dimension.png", p, width = 15, height = 8, dpi = 200)

cat("\nSelected n_c (truth 10) and reference by mechanism / option:\n")
print(unique(res[, c("mechanism", "option", "rep", "reference", "selected")]), row.names = FALSE)
