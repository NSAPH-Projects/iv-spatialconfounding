# Helper functions for spatial basis selection.

.basis_stop <- function(message) {
  stop(message, call. = FALSE)
}

.is_integerish_scalar <- function(x) {
  length(x) == 1L && is.numeric(x) && is.finite(x) &&
    abs(x - round(x)) < sqrt(.Machine$double.eps)
}

.validate_index <- function(index, n, name) {
  if (!is.numeric(index) || length(index) == 0L || any(!is.finite(index)) ||
      any(index != as.integer(index)) || any(index < 1L) || any(index > n) ||
      anyDuplicated(index)) {
    .basis_stop(sprintf(
      "%s must contain unique integer row indices between 1 and %d.",
      name, n
    ))
  }
  as.integer(index)
}

.validate_basis <- function(B) {
  if (!is.matrix(B) || !is.numeric(B) || nrow(B) < 2L || ncol(B) < 2L ||
      any(!is.finite(B))) {
    .basis_stop("B must be a finite numeric matrix with at least two rows and columns.")
  }
  invisible(B)
}

# Construct the candidate grid in Supplement Section 4.
#
# n_uc_star = floor(core_fraction * m) unless supplied explicitly.
# n_max = m - max(5, floor(0.02 * m)) unless min_confounded is supplied.
# The sequence follows the displayed formula exactly and therefore does not
# append n_max when it is not reached by an integer number of `step`s.
make_candidate_grid <- function(m,
                                core_fraction = 0.8,
                                step = 5L,
                                min_confounded = NULL,
                                n_uc_star = NULL) {
  if (!.is_integerish_scalar(m) || m < 2L) {
    .basis_stop("m must be an integer greater than or equal to 2.")
  }
  m <- as.integer(m)

  if (!is.numeric(core_fraction) || length(core_fraction) != 1L ||
      !is.finite(core_fraction) || core_fraction <= 0 || core_fraction >= 1) {
    .basis_stop("core_fraction must be a finite scalar strictly between 0 and 1.")
  }
  if (!.is_integerish_scalar(step) || step < 1L) {
    .basis_stop("step must be a positive integer.")
  }
  step <- as.integer(step)

  if (is.null(min_confounded)) {
    min_confounded <- max(5L, floor(0.02 * m))
  }
  if (!.is_integerish_scalar(min_confounded) ||
      min_confounded < 1L || min_confounded >= m) {
    .basis_stop("min_confounded must be a positive integer smaller than m.")
  }
  min_confounded <- as.integer(min_confounded)

  if (is.null(n_uc_star)) {
    n_uc_star <- floor(core_fraction * m)
  }
  if (!.is_integerish_scalar(n_uc_star) || n_uc_star < 1L) {
    .basis_stop("n_uc_star must be a positive integer.")
  }
  n_uc_star <- as.integer(n_uc_star)

  n_max <- m - min_confounded
  if (n_uc_star > n_max) {
    .basis_stop(sprintf(
      paste0(
        "The core set (%d) exceeds n_max (%d). Increase m, reduce the core ",
        "fraction, or reduce min_confounded."
      ),
      n_uc_star, n_max
    ))
  }

  as.integer(seq.int(n_uc_star, n_max, by = step))
}

# Return column indices for B^c and B^uc.
#
# `instrument_scale = "small"` is the primary: B is
# ordered from largest to smallest spatial scale and its last n_uc columns are
# candidate instruments. `"large"` supports the reversed-scale simulation.
partition_basis_columns <- function(m, n_uc,
                                    instrument_scale = c("small", "large")) {
  instrument_scale <- match.arg(instrument_scale)
  if (!.is_integerish_scalar(m) || m < 2L) {
    .basis_stop("m must be an integer greater than or equal to 2.")
  }
  if (!.is_integerish_scalar(n_uc) || n_uc < 1L || n_uc >= m) {
    .basis_stop("n_uc must be a positive integer smaller than m.")
  }
  m <- as.integer(m)
  n_uc <- as.integer(n_uc)

  if (instrument_scale == "small") {
    confounded <- seq_len(m - n_uc)
    instruments <- seq.int(m - n_uc + 1L, m)
  } else {
    instruments <- seq_len(n_uc)
    confounded <- seq.int(n_uc + 1L, m)
  }

  list(confounded = confounded, instruments = instruments)
}

# Estimate gamma^c on training rows and apply the same coefficients elsewhere.
#
# This deliberately uses a QR least-squares fit. Although the complete basis B
# may be orthonormal, B[train_idx, ] generally is not, so t(B_train) %*% A_train
# is not the fold-specific coefficient estimator in Supplement equation (4.7).
project_Ac <- function(B, A, n_uc, train_idx,
                       predict_idx = seq_len(nrow(B)),
                       instrument_scale = c("small", "large"),
                       tol = 1e-10) {
  .validate_basis(B)
  instrument_scale <- match.arg(instrument_scale)
  n <- nrow(B)
  m <- ncol(B)

  if (!is.numeric(A) || length(A) != n || any(!is.finite(A))) {
    .basis_stop("A must be a finite numeric vector with one value per row of B.")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    .basis_stop("tol must be a positive finite scalar.")
  }

  train_idx <- .validate_index(train_idx, n, "train_idx")
  predict_idx <- .validate_index(predict_idx, n, "predict_idx")
  columns <- partition_basis_columns(m, n_uc, instrument_scale)
  B_train <- B[train_idx, columns$confounded, drop = FALSE]

  if (nrow(B_train) <= ncol(B_train)) {
    .basis_stop(sprintf(
      "The training set has %d rows but B^c has %d columns; the projection is not identified.",
      nrow(B_train), ncol(B_train)
    ))
  }

  B_qr <- qr(B_train, tol = tol)
  if (B_qr$rank < ncol(B_train)) {
    .basis_stop("B^c is rank deficient on the training rows.")
  }
  gamma <- drop(qr.coef(B_qr, A[train_idx]))
  Ac <- drop(B[predict_idx, columns$confounded, drop = FALSE] %*% gamma)

  structure(
    list(
      Ac = Ac,
      gamma = gamma,
      n_uc = as.integer(n_uc),
      train_idx = train_idx,
      predict_idx = predict_idx,
      confounded_cols = columns$confounded,
      instrument_cols = columns$instruments,
      instrument_scale = instrument_scale
    ),
    class = "Ac_projection"
  )
}

# Paired standard error from Supplement Section 4.3.
paired_difference_se <- function(phi_candidate, phi_core) {
  if (!is.numeric(phi_candidate) || !is.numeric(phi_core) ||
      length(phi_candidate) != length(phi_core) ||
      length(phi_candidate) < 2L ||
      any(!is.finite(phi_candidate)) || any(!is.finite(phi_core))) {
    .basis_stop(
      "phi_candidate and phi_core must be finite numeric vectors of equal length (at least 2)."
    )
  }

  phi_difference <- phi_candidate - phi_core
  n <- length(phi_difference)
  centered <- phi_difference - mean(phi_difference)
  sqrt(sum(centered^2) / (n * (n - 1L)))
}

# Apply the core-referenced, contiguous selection rule.
select_contiguous_candidate <- function(n_uc_cands, psi, influence_curves,
                                        alpha = 1) {
  if (!is.numeric(n_uc_cands) || length(n_uc_cands) < 1L ||
      any(!is.finite(n_uc_cands)) || any(n_uc_cands != as.integer(n_uc_cands)) ||
      is.unsorted(n_uc_cands, strictly = TRUE)) {
    .basis_stop("n_uc_cands must be a strictly increasing integer vector.")
  }
  n_uc_cands <- as.integer(n_uc_cands)
  n_candidates <- length(n_uc_cands)

  if (!is.numeric(psi) || length(psi) != n_candidates || any(!is.finite(psi))) {
    .basis_stop("psi must contain one finite estimate per candidate.")
  }
  if (is.list(influence_curves)) {
    influence_curves <- do.call(cbind, influence_curves)
  }
  if (!is.matrix(influence_curves) || !is.numeric(influence_curves) ||
      ncol(influence_curves) != n_candidates || nrow(influence_curves) < 2L ||
      any(!is.finite(influence_curves))) {
    .basis_stop(
      "influence_curves must be a finite matrix with observations in rows and candidates in columns."
    )
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0) {
    .basis_stop("alpha must be a positive finite scalar.")
  }

  difference <- psi - psi[1L]
  paired_se <- threshold <- rep(NA_real_, n_candidates)
  passes <- rep(FALSE, n_candidates)
  passes[1L] <- TRUE

  if (n_candidates > 1L) {
    for (j in 2:n_candidates) {
      paired_se[j] <- paired_difference_se(
        influence_curves[, j], influence_curves[, 1L]
      )
      threshold[j] <- alpha * paired_se[j]
      passes[j] <- abs(difference[j]) <= threshold[j]
    }
  }

  in_contiguous_sequence <- as.logical(cumprod(as.integer(passes)))
  selected_index <- max(which(in_contiguous_sequence))
  diagnostics <- data.frame(
    n_uc = n_uc_cands,
    psi = psi,
    difference_from_core = difference,
    paired_se = paired_se,
    threshold = threshold,
    passes = passes,
    in_contiguous_sequence = in_contiguous_sequence
  )

  structure(
    list(
      selected_n_uc = n_uc_cands[selected_index],
      selected_index = selected_index,
      alpha = alpha,
      diagnostics = diagnostics
    ),
    class = "basis_selection"
  )
}

# Uniform random folds, as specified by Supplement Algorithm 1.
# Call set.seed() before this function when a reproducible split is required.
make_uniform_folds <- function(n, K = 5L) {
  if (!.is_integerish_scalar(n) || n < 2L) {
    .basis_stop("n must be an integer greater than or equal to 2.")
  }
  if (!.is_integerish_scalar(K) || K < 2L || K > n) {
    .basis_stop("K must be an integer between 2 and n.")
  }
  n <- as.integer(n)
  K <- as.integer(K)
  sample(rep_len(seq_len(K), n), size = n, replace = FALSE)
}

# Select n_uc for one outer fold using only that fold's training observations.
#
# `estimate_candidate` is called once per candidate with these named arguments:
#   y, a, x, Ac, cutoff, folds, n_uc
# It must use the supplied two-fold split and return:
#   list(psi = <scalar>, influence_curve = <length(training rows) vector>)
# The influence curve must be for the final truncated effect estimate, after the
# delta-method combination of its four components.
select_basis_outer_fold <- function(y, a, x = NULL, B, outer_train_idx,
                                    n_uc_cands, cutoff, estimate_candidate,
                                    inner_folds = NULL, alpha = 1,
                                    instrument_scale = c("small", "large"), ...) {
  .validate_basis(B)
  instrument_scale <- match.arg(instrument_scale)
  n <- nrow(B)

  if (!is.numeric(y) || length(y) != n || any(!is.finite(y))) {
    .basis_stop("y must be a finite numeric vector with one value per row of B.")
  }
  if (!is.numeric(a) || length(a) != n || any(!is.finite(a))) {
    .basis_stop("a must be a finite numeric vector with one value per row of B.")
  }
  if (is.null(x)) {
    x <- matrix(numeric(0), nrow = n, ncol = 0L)
  } else {
    x <- as.matrix(x)
    if (!is.numeric(x) || nrow(x) != n || any(!is.finite(x))) {
      .basis_stop("x must be NULL or a finite numeric matrix with one row per row of B.")
    }
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff)) {
    .basis_stop("cutoff must be a finite numeric scalar.")
  }
  if (!is.function(estimate_candidate)) {
    .basis_stop("estimate_candidate must be a function.")
  }

  outer_train_idx <- .validate_index(outer_train_idx, n, "outer_train_idx")
  n_train <- length(outer_train_idx)
  if (n_train < 4L) {
    .basis_stop("At least four outer-training observations are required.")
  }

  if (is.null(inner_folds)) {
    inner_folds <- make_uniform_folds(n_train, K = 2L)
  }
  if (!is.numeric(inner_folds) || length(inner_folds) != n_train ||
      any(!is.finite(inner_folds)) || length(unique(inner_folds)) != 2L) {
    .basis_stop("inner_folds must assign every outer-training row to one of exactly two folds.")
  }

  if (!is.numeric(n_uc_cands) || length(n_uc_cands) < 1L ||
      any(!is.finite(n_uc_cands)) ||
      any(n_uc_cands != as.integer(n_uc_cands)) ||
      any(n_uc_cands < 1L) || any(n_uc_cands >= ncol(B)) ||
      is.unsorted(n_uc_cands, strictly = TRUE)) {
    .basis_stop(
      "n_uc_cands must be a strictly increasing integer vector between 1 and ncol(B) - 1."
    )
  }
  n_uc_cands <- as.integer(n_uc_cands)

  candidate_results <- vector("list", length(n_uc_cands))
  influence_curves <- matrix(NA_real_, nrow = n_train, ncol = length(n_uc_cands))
  psi <- rep(NA_real_, length(n_uc_cands))
  extra_args <- list(...)

  for (j in seq_along(n_uc_cands)) {
    projection <- project_Ac(
      B = B,
      A = a,
      n_uc = n_uc_cands[j],
      train_idx = outer_train_idx,
      predict_idx = outer_train_idx,
      instrument_scale = instrument_scale
    )

    result <- do.call(
      estimate_candidate,
      c(
        list(
          y = y[outer_train_idx],
          a = a[outer_train_idx],
          x = x[outer_train_idx, , drop = FALSE],
          Ac = projection$Ac,
          cutoff = cutoff,
          folds = inner_folds,
          n_uc = n_uc_cands[j]
        ),
        extra_args
      )
    )

    if (!is.list(result) || !is.numeric(result$psi) || length(result$psi) != 1L ||
        !is.finite(result$psi) || !is.numeric(result$influence_curve) ||
        length(result$influence_curve) != n_train ||
        any(!is.finite(result$influence_curve))) {
      .basis_stop(
        "estimate_candidate must return a finite scalar psi and a finite influence_curve for every training row."
      )
    }

    psi[j] <- result$psi
    influence_curves[, j] <- result$influence_curve
    candidate_results[[j]] <- result
  }

  selection <- select_contiguous_candidate(
    n_uc_cands = n_uc_cands,
    psi = psi,
    influence_curves = influence_curves,
    alpha = alpha
  )

  selected_projection <- project_Ac(
    B = B,
    A = a,
    n_uc = selection$selected_n_uc,
    train_idx = outer_train_idx,
    predict_idx = seq_len(n),
    instrument_scale = instrument_scale
  )

  structure(
    list(
      selected_n_uc = selection$selected_n_uc,
      selection = selection,
      selected_projection = selected_projection,
      candidate_results = candidate_results,
      influence_curves = influence_curves,
      outer_train_idx = outer_train_idx,
      inner_folds = inner_folds
    ),
    class = "outer_fold_basis_selection"
  )
}
