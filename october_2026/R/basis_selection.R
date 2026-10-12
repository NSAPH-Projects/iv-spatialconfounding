# Sequential stability selection of the spatial-adjustment dimension
#
# Shared by october_2026/simulation and october_2026/data_application.
#
# B is an ordered basis (column 1 = largest scale / constant, column m =
# smallest scale; reverse it beforehand for scale-reversed settings).
# B^c(n_ac) is the first n_ac columns; A^c = P_{B^c} A is added to the
# adjustment set. Equivalently, n_uc = m - n_ac.
#
# Rule (manuscript Section 4):
#   1. Reference: the smallest n_ac0 whose B^c explains at least r2_reference
#      of Var(A). This is the heaviest adjustment considered; it is assumed to
#      contain all confounded exposure variation.
#   2. Candidates: n_ac = n_ac0 - 1, n_ac0 - 2, ..., 1, i.e. one basis
#      function at a time is moved from B^c to B^uc. J = n_ac0 - 1.
#   3. For each candidate, test H0: psi(candidate) = psi(reference) using the
#      paired influence-function standard error, rejecting when
#      |d| > z_{1 - alpha / (2J)} * se(d) (Bonferroni over all J candidates).
#   4. Stop at the first rejection and select the candidate before it (the
#      reference itself if the first candidate is rejected). If no candidate
#      is rejected, select the last candidate, n_ac = 1.
# The rule only uses A to define the candidates, never Y.

project_Ac_full <- function(B, A, n_ac) {
  B <- as.matrix(B)
  A <- as.numeric(A)
  n_ac <- as.integer(n_ac)

  if (length(A) != nrow(B)) stop("A must have one value per row of B.")
  if (n_ac < 1L || n_ac > ncol(B)) {
    stop("n_ac must satisfy 1 <= n_ac <= ncol(B).")
  }

  Bc <- B[, seq_len(n_ac), drop = FALSE]
  Ac <- drop(Bc %*% qr.coef(qr(Bc), A))
  Auc <- A - Ac

  exposure_R2 <- 1 - sum(Auc^2) / sum((A - mean(A))^2)

  list(
    Ac = Ac,
    Auc = Auc,
    exposure_R2 = min(max(exposure_R2, 0), 1),
    n_ac = n_ac,
    n_uc = ncol(B) - n_ac
  )
}


# Exposure R^2 of B^c(q) for every q = 1, ..., m, via one QR of B (valid for
# any full-rank B, orthonormal or not).
exposure_r2_curve <- function(B, A) {
  B <- as.matrix(B)
  A <- as.numeric(A)
  qty <- qr.qty(qr(B), A)[seq_len(ncol(B))]
  1 - (sum(A^2) - cumsum(qty^2)) / sum((A - mean(A))^2)
}


paired_difference_se <- function(phi_candidate, phi_reference) {
  if (length(phi_candidate) != length(phi_reference)) {
    stop("The candidate and reference influence curves must have equal length.")
  }
  stats::sd(phi_candidate - phi_reference) / sqrt(length(phi_reference))
}


select_basis_dimension <- function(
    y,
    a,
    x,
    B,
    cutoff,
    effect_estimator,
    effect_args = list(),
    r2_reference = 0.7,
    family_alpha = 0.05,
    seed = 12001L) {

  B <- as.matrix(B)
  x <- as.matrix(x)
  n <- length(y)
  m <- ncol(B)

  if (length(a) != n || nrow(x) != n || nrow(B) != n) {
    stop("y, a, x, and B must contain the same observations.")
  }
  if (r2_reference <= 0 || r2_reference >= 1) {
    stop("r2_reference must be strictly between 0 and 1.")
  }

  fit_at <- function(n_ac) {
    projection <- project_Ac_full(B, a, n_ac)
    set.seed(seed)  # same Super Learner CV splits for every fit
    fit <- do.call(
      effect_estimator,
      c(list(y = y, a = a, w = cbind(Ac = projection$Ac, x), cutoff = cutoff),
        effect_args)
    )
    list(projection = projection, fit = fit)
  }

  # 1. Reference.
  r2 <- exposure_r2_curve(B, a)
  reached <- which(r2 >= r2_reference)
  reached_reference <- length(reached) > 0L
  n_ac0 <- if (reached_reference) reached[1L] else m
  if (!reached_reference) {
    warning(sprintf(
      "The full basis explains only %.3f of Var(A) < r2_reference = %.2f; using all %d columns as the reference.",
      r2[m], r2_reference, m
    ))
  }
  reference <- fit_at(n_ac0)

  # 2-3. Candidates and Bonferroni threshold (fixed upfront from A alone).
  candidates <- if (n_ac0 > 1L) seq.int(n_ac0 - 1L, 1L) else integer(0)
  J <- length(candidates)
  critical_value <- if (J > 0L) stats::qnorm(1 - family_alpha / (2 * J)) else NA_real_

  row_for <- function(n_ac, current, difference, paired_se, departs) {
    data.frame(
      n_ac = n_ac,
      n_uc = m - n_ac,
      exposure_R2 = current$projection$exposure_R2,
      psi = current$fit$psi,
      ci_lower = current$fit$confidence_interval[1L],
      ci_upper = current$fit$confidence_interval[2L],
      standard_error = current$fit$standard_error,
      difference_from_reference = difference,
      paired_se = paired_se,
      departure_threshold = critical_value * paired_se,
      departs_from_reference = departs
    )
  }

  rows <- list(row_for(n_ac0, reference, 0, NA_real_, FALSE))
  selected <- reference
  selected_n_ac <- n_ac0
  rejected <- FALSE

  # 4. Walk down; stop at the first rejection.
  for (q in candidates) {
    current <- tryCatch(fit_at(q), error = function(e) {
      message("  candidate fit failed for n_ac = ", q, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(current)) next

    difference <- current$fit$psi - reference$fit$psi
    paired_se <- paired_difference_se(
      current$fit$influence_curve, reference$fit$influence_curve
    )
    departs <- abs(difference) > critical_value * paired_se
    rows[[length(rows) + 1L]] <- row_for(q, current, difference, paired_se, departs)

    if (departs) {
      rejected <- TRUE
      break
    }
    selected <- current
    selected_n_ac <- q
  }

  diagnostics <- do.call(rbind, rows)
  diagnostics$selected <- diagnostics$n_ac == selected_n_ac

  list(
    selected_n_ac = as.integer(selected_n_ac),
    selected_n_uc = as.integer(m - selected_n_ac),
    selected_fit = selected$fit,
    selected_projection = selected$projection,
    reference_fit = reference$fit,
    diagnostics = diagnostics,
    settings = list(
      cutoff = cutoff,
      r2_reference = r2_reference,
      reached_reference = reached_reference,
      n_ac0 = n_ac0,
      J = J,
      family_alpha = family_alpha,
      bonferroni_critical = critical_value,
      rejected = rejected
    )
  )
}
