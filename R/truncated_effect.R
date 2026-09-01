# Doubly robust estimation of the truncated exposure effect.
#
# This file implements the candidate estimator used inside each outer training
# fold during basis selection.
# nuisance estimation is restricted to A >= cutoff, pseudo-outcomes are
# cross-fitted, nu(c) is estimated by a one-sided local quadratic regression, (quadratic better than linear for boundary)
# and the final influence curve combines four components through the delta method.

.effect_stop <- function(message) {
  stop(message, call. = FALSE)
}

.effect_integerish <- function(x) {
  length(x) == 1L && is.numeric(x) && is.finite(x) &&
    abs(x - round(x)) < sqrt(.Machine$double.eps)
}

.validate_nuisance_result <- function(result, n, n_grid) {
  required <- c("pihat", "pihat_grid", "muhat", "muhat_grid")
  if (!is.list(result) || !all(required %in% names(result))) {
    .effect_stop(sprintf(
      "The nuisance estimator must return: %s.",
      paste(required, collapse = ", ")
    ))
  }
  if (!is.numeric(result$pihat) || length(result$pihat) != n ||
      any(!is.finite(result$pihat)) || any(result$pihat < 0) ||
      !is.numeric(result$muhat) || length(result$muhat) != n ||
      any(!is.finite(result$muhat))) {
    .effect_stop("pihat and muhat must be finite numeric vectors of the expected length.")
  }
  if (!is.matrix(result$pihat_grid) || !is.numeric(result$pihat_grid) ||
      !identical(dim(result$pihat_grid), c(n, n_grid)) ||
      any(!is.finite(result$pihat_grid)) || any(result$pihat_grid < 0) ||
      !is.matrix(result$muhat_grid) || !is.numeric(result$muhat_grid) ||
      !identical(dim(result$muhat_grid), c(n, n_grid)) ||
      any(!is.finite(result$muhat_grid))) {
    .effect_stop(paste0(
      "pihat_grid and muhat_grid must be finite n-by-length(a_grid) numeric ",
      "matrices, and density estimates cannot be negative."
    ))
  }
  invisible(result)
}

# Cap the inverse-density weight ratio (varpi / pihat) at its own upper
# empirical quantile.
#
# The pihat floor below is a quantile of pihat's own distribution, so when a
# candidate's conditional-density estimate is degenerate across the board
# (e.g. very few B^c columns near n_max), the floor is itself tiny and
# varpi / pihat can still explode. Truncating the ratio directly bounds the
# pseudo-outcome regardless of why the density estimate degenerated.
.truncate_ipw_ratio <- function(ratio, ipw_ratio_trim) {
  if (!is.numeric(ratio) || length(ratio) < 1L || any(!is.finite(ratio)) ||
      any(ratio < 0)) {
    .effect_stop("ratio must be a finite, non-negative numeric vector.")
  }
  if (!is.numeric(ipw_ratio_trim) || length(ipw_ratio_trim) != 1L ||
      !is.finite(ipw_ratio_trim) || ipw_ratio_trim < 0 ||
      ipw_ratio_trim >= 0.5) {
    .effect_stop("ipw_ratio_trim must be a finite scalar in [0, 0.5).")
  }
  if (ipw_ratio_trim == 0 || length(ratio) < 2L) {
    return(ratio)
  }
  cap <- as.numeric(stats::quantile(
    ratio,
    probs = 1 - ipw_ratio_trim,
    names = FALSE
  ))
  pmin(ratio, cap)
}

.trapezoid_rows <- function(values, grid) {
  if (!is.matrix(values) || ncol(values) != length(grid) || length(grid) < 2L) {
    .effect_stop("values must be a matrix with one column per grid point.")
  }
  increments <- diff(grid)
  rowSums(
    sweep(
      (values[, -1L, drop = FALSE] +
         values[, -ncol(values), drop = FALSE]) / 2,
      2L,
      increments,
      `*`
    )
  )
}

.local_linear_at <- function(a, outcome, target, bandwidth,
                             kernel = stats::dnorm, degree = 1L) {
  if (!degree %in% c(1L, 2L)) {
    .effect_stop("degree must be 1 (local linear) or 2 (local quadratic).")
  }
  u <- (a - target) / bandwidth
  weights <- kernel(u) / bandwidth
  design <- if (degree == 1L) cbind(1, u) else cbind(1, u, u^2)
  D <- crossprod(design, design * weights) / length(a)

  if (any(!is.finite(D)) || rcond(D) < 1e-10) {
    .effect_stop("The local-linear design is singular at the requested bandwidth.")
  }
  moment <- colMeans(design * (weights * outcome))
  beta <- drop(solve(D, moment))

  list(beta = beta, D = D, u = u, weights = weights, design = design, degree = degree)
}

.select_local_linear_bandwidth <- function(a, outcome, bw_seq,
                                           kernel = stats::dnorm) {
  if (!requireNamespace("KernSmooth", quietly = TRUE)) {
    .effect_stop("KernSmooth is required for automatic bandwidth selection.")
  }
  if (!is.numeric(bw_seq) || length(bw_seq) < 1L ||
      any(!is.finite(bw_seq)) || any(bw_seq <= 0)) {
    .effect_stop("bw_seq must contain positive finite bandwidths.")
  }

  evaluation_grid <- seq(min(a), max(a), length.out = min(100L, length(a)))

  hat_values <- function(bandwidth) {
    values <- vapply(evaluation_grid, function(target) {
      u <- (a - target) / bandwidth
      weights <- kernel(u) / bandwidth
      s0 <- mean(weights)
      s1 <- mean(u * weights)
      s2 <- mean(u^2 * weights)
      denominator <- s0 * s2 - s1^2
      if (!is.finite(denominator) || abs(denominator) < 1e-10) {
        return(NA_real_)
      }
      s2 * (kernel(0) / bandwidth) / (denominator * length(a))
    }, numeric(1))
    stats::approx(evaluation_grid, values, xout = a, rule = 2)$y
  }

  fitted_values <- function(bandwidth) {
    fit <- tryCatch(
      KernSmooth::locpoly(
        x = a,
        y = outcome,
        bandwidth = bandwidth,
        degree = 1L
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(rep(NA_real_, length(a)))
    }
    stats::approx(fit$x, fit$y, xout = a, rule = 2)$y
  }

  risk <- vapply(bw_seq, function(bandwidth) {
    hats <- hat_values(bandwidth)
    fitted <- fitted_values(bandwidth)
    denominator <- 1 - hats
    valid <- is.finite(hats) & is.finite(fitted) & abs(denominator) > 1e-8
    if (mean(valid) < 0.5) {
      return(NA_real_)
    }
    mean(((outcome[valid] - fitted[valid]) / denominator[valid])^2)
  }, numeric(1))

  if (!any(is.finite(risk))) {
    .effect_stop("Bandwidth selection failed: no candidate had finite risk.")
  }

  selected <- which.min(replace(risk, !is.finite(risk), Inf))
  list(
    bandwidth = bw_seq[selected],
    risk = data.frame(bandwidth = bw_seq, risk = risk)
  )
}

# Cross-fitted nuisance functions using the Super Learner library reported in
# the supplement. This function is called only on the A >= cutoff subpopulation.
estimate_nuisance_superlearner <- function(
    y, a, w, folds, a_grid,
    sl_library = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction"),
    density_trim = 0,
    cv_control = list(V = 2L),
    variance_floor = 1e-4) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    .effect_stop("SuperLearner is required for the nuisance estimator.")
  }

  w <- as.data.frame(w)
  n <- length(a)
  n_grid <- length(a_grid)
  fold_values <- sort(unique(folds))
  if (length(fold_values) != 2L) {
    .effect_stop("The nuisance estimator requires exactly two cross-fitting folds.")
  }
  if (!is.numeric(density_trim) || length(density_trim) != 1L ||
      !is.finite(density_trim) || density_trim < 0 || density_trim >= 0.5) {
    .effect_stop("density_trim must be a finite scalar in [0, 0.5).")
  }
  if (!is.numeric(variance_floor) || length(variance_floor) != 1L ||
      !is.finite(variance_floor) || variance_floor <= 0) {
    .effect_stop("variance_floor must be a positive finite scalar.")
  }

  pihat <- muhat <- rep(NA_real_, n)
  pihat_grid <- muhat_grid <- matrix(NA_real_, nrow = n, ncol = n_grid)

  for (fold in fold_values) {
    train <- which(folds != fold)
    holdout <- which(folds == fold)
    n_holdout <- length(holdout)
    if (length(train) < 4L || n_holdout < 1L) {
      .effect_stop("Each nuisance fold needs at least four training rows and one holdout row.")
    }

    w_grid <- w[holdout[rep(seq_len(n_holdout), times = n_grid)], , drop = FALSE]
    w_new <- rbind(w, w_grid)

    mean_fit <- SuperLearner::SuperLearner(
      Y = a[train],
      X = w[train, , drop = FALSE],
      newX = w_new,
      SL.library = sl_library,
      cvControl = cv_control,
      env = asNamespace("SuperLearner")
    )
    mean_prediction <- drop(mean_fit$SL.predict)

    training_residual <- a[train] - mean_prediction[train]
    variance_fit <- SuperLearner::SuperLearner(
      Y = log(pmax(training_residual^2, .Machine$double.eps)),
      X = w[train, , drop = FALSE],
      newX = w_new,
      SL.library = sl_library,
      cvControl = cv_control,
      env = asNamespace("SuperLearner")
    )
    variance_prediction <- pmax(
      exp(drop(variance_fit$SL.predict)),
      variance_floor
    )

    wa_observed <- data.frame(w, a = a, check.names = FALSE)
    wa_grid <- data.frame(
      w_grid,
      a = rep(a_grid, each = n_holdout),
      check.names = FALSE
    )
    wa_new <- rbind(wa_observed, wa_grid)
    outcome_fit <- SuperLearner::SuperLearner(
      Y = y[train],
      X = wa_observed[train, , drop = FALSE],
      newX = wa_new,
      SL.library = sl_library,
      cvControl = cv_control,
      env = asNamespace("SuperLearner")
    )
    outcome_prediction <- drop(outcome_fit$SL.predict)

    grid_positions <- n + seq_len(n_holdout * n_grid)
    standardized_train <- training_residual /
      sqrt(variance_prediction[train])
    standardized_holdout <- (a[holdout] - mean_prediction[holdout]) /
      sqrt(variance_prediction[holdout])
    standardized_grid <- (
      rep(a_grid, each = n_holdout) - mean_prediction[grid_positions]
    ) / sqrt(variance_prediction[grid_positions])

    density_range <- range(
      c(standardized_train, standardized_holdout, standardized_grid),
      finite = TRUE
    )
    residual_density <- stats::density(
      standardized_train,
      from = density_range[1L],
      to = density_range[2L]
    )

    pihat[holdout] <- stats::approx(
      residual_density$x,
      residual_density$y,
      xout = standardized_holdout,
      rule = 2
    )$y / sqrt(variance_prediction[holdout])
    pihat_grid[holdout, ] <- matrix(
      stats::approx(
        residual_density$x,
        residual_density$y,
        xout = standardized_grid,
        rule = 2
      )$y / sqrt(variance_prediction[grid_positions]),
      nrow = n_holdout,
      ncol = n_grid
    )
    muhat[holdout] <- outcome_prediction[holdout]
    muhat_grid[holdout, ] <- matrix(
      outcome_prediction[grid_positions],
      nrow = n_holdout,
      ncol = n_grid
    )
  }

  density_floor <- max(
    as.numeric(stats::quantile(pihat, probs = density_trim, names = FALSE)),
    .Machine$double.eps
  )

  list(
    pihat = pmax(pihat, density_floor),
    pihat_grid = pihat_grid,
    muhat = muhat,
    muhat_grid = muhat_grid,
    density_floor = density_floor
  )
}

# Estimate the truncated exposure effect on one outer-training sample.
#
# This signature matches the callback expected by select_basis_outer_fold().
# `n_uc` is accepted for diagnostics but does not otherwise alter estimation;
# the candidate-specific A^c has already been supplied in `Ac`.
estimate_truncated_effect_candidate <- function(
    y, a, x, Ac, cutoff, folds, n_uc = NULL,
    nuisance_estimator = estimate_nuisance_superlearner,
    nuisance_args = list(),
    a_grid = NULL,
    n_grid = 100L,
    bandwidth = NULL,
    bw_seq = NULL,
    kernel = stats::dnorm,
    density_trim = 0,
    ipw_ratio_trim = 0,
    constrain = FALSE,
    local_degree = 1L) {
  n <- length(y)
  if (!is.numeric(y) || !is.numeric(a) || length(a) != n || n < 8L ||
      any(!is.finite(y)) || any(!is.finite(a))) {
    .effect_stop("y and a must be finite numeric vectors of equal length (at least 8).")
  }
  if (missing(x) || is.null(x)) {
    x <- matrix(numeric(0), nrow = n, ncol = 0L)
  } else {
    x <- as.matrix(x)
    if (!is.numeric(x) || nrow(x) != n || any(!is.finite(x))) {
      .effect_stop("x must be NULL or a finite numeric matrix with one row per observation.")
    }
  }
  if (!is.numeric(Ac) || length(Ac) != n || any(!is.finite(Ac))) {
    .effect_stop("Ac must be a finite numeric vector with one value per observation.")
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff)) {
    .effect_stop("cutoff must be a finite numeric scalar.")
  }
  if (!is.numeric(folds) || length(folds) != n || any(!is.finite(folds)) ||
      length(unique(folds)) != 2L) {
    .effect_stop("folds must assign every observation to one of exactly two folds.")
  }
  if (!is.function(nuisance_estimator)) {
    .effect_stop("nuisance_estimator must be a function.")
  }
  if (!is.list(nuisance_args)) {
    .effect_stop("nuisance_args must be a list.")
  }
  if (!is.function(kernel)) {
    .effect_stop("kernel must be a function.")
  }
  if (!is.logical(constrain) || length(constrain) != 1L || is.na(constrain)) {
    .effect_stop("constrain must be TRUE or FALSE.")
  }
  if (!is.numeric(density_trim) || length(density_trim) != 1L ||
      !is.finite(density_trim) || density_trim < 0 || density_trim >= 0.5) {
    .effect_stop("density_trim must be a finite scalar in [0, 0.5).")
  }
  if (!is.numeric(ipw_ratio_trim) || length(ipw_ratio_trim) != 1L ||
      !is.finite(ipw_ratio_trim) || ipw_ratio_trim < 0 ||
      ipw_ratio_trim >= 0.5) {
    .effect_stop("ipw_ratio_trim must be a finite scalar in [0, 0.5).")
  }

  above <- a >= cutoff
  below <- a < cutoff
  n_above <- sum(above)
  n_below <- sum(below)
  if (n_above < 6L || n_below < 2L || length(unique(folds[above])) != 2L) {
    .effect_stop(
      "The cutoff must leave at least six observations above, two below, and both nuisance folds above."
    )
  }
  if (diff(range(a[above])) <= 0) {
    .effect_stop("Exposure must vary among observations with A >= cutoff.")
  }

  if (is.null(a_grid)) {
    if (!.effect_integerish(n_grid) || n_grid < 5L) {
      .effect_stop("n_grid must be an integer greater than or equal to 5.")
    }
    a_grid <- seq(cutoff, max(a[above]), length.out = as.integer(n_grid))
  } else {
    if (!is.numeric(a_grid) || length(a_grid) < 5L || any(!is.finite(a_grid)) ||
        is.unsorted(a_grid, strictly = TRUE) ||
        abs(a_grid[1L] - cutoff) > sqrt(.Machine$double.eps) ||
        max(a_grid) < max(a[above])) {
      .effect_stop(
        paste0(
          "a_grid must be a strictly increasing finite grid beginning at cutoff ",
          "and spanning max(A | A >= cutoff)."
        )
      )
    }
  }

  w <- cbind(Ac = Ac, x)
  w_names <- colnames(w)
  unnamed <- is.na(w_names) | !nzchar(w_names)
  w_names[unnamed] <- paste0("x", which(unnamed))
  colnames(w) <- make.unique(make.names(w_names))
  w_above <- w[above, , drop = FALSE]
  nuisance <- do.call(
    nuisance_estimator,
    c(
      list(
        y = y[above],
        a = a[above],
        w = w_above,
        folds = folds[above],
        a_grid = a_grid
      ),
      nuisance_args
    )
  )
  .validate_nuisance_result(nuisance, n_above, length(a_grid))

  density_floor <- max(
    as.numeric(stats::quantile(
      nuisance$pihat,
      probs = density_trim,
      names = FALSE
    )),
    .Machine$double.eps
  )
  pihat <- pmax(nuisance$pihat, density_floor)
  varpi_grid <- pmax(colMeans(nuisance$pihat_grid), .Machine$double.eps)
  m_grid <- colMeans(nuisance$muhat_grid)
  varpi_observed <- pmax(
    stats::approx(a_grid, varpi_grid, xout = a[above], rule = 2)$y,
    .Machine$double.eps
  )
  m_observed <- stats::approx(a_grid, m_grid, xout = a[above], rule = 2)$y

  ipw_ratio_raw <- varpi_observed / pihat
  ipw_ratio <- .truncate_ipw_ratio(ipw_ratio_raw, ipw_ratio_trim)
  ipw_ratio_cap <- if (ipw_ratio_trim > 0) {
    as.numeric(stats::quantile(
      ipw_ratio_raw,
      probs = 1 - ipw_ratio_trim,
      names = FALSE
    ))
  } else {
    Inf
  }
  pseudo_raw <- (y[above] - nuisance$muhat) * ipw_ratio +
    m_observed
  pseudo_for_fit <- pseudo_raw
  if (constrain) {
    outcome_range <- range(y)
    pseudo_for_fit <- pmin(pmax(pseudo_for_fit, outcome_range[1L]), outcome_range[2L])
  }

  if (is.null(bandwidth)) {
    if (is.null(bw_seq)) {
      exposure_sd <- stats::sd(a[above])
      bw_seq <- seq(exposure_sd / 10, exposure_sd, length.out = 100L)
    }
    bandwidth_fit <- .select_local_linear_bandwidth(
      a = a[above],
      outcome = pseudo_for_fit,
      bw_seq = bw_seq,
      kernel = kernel
    )
    bandwidth <- bandwidth_fit$bandwidth
    bandwidth_risk <- bandwidth_fit$risk
  } else {
    if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
        !is.finite(bandwidth) || bandwidth <= 0) {
      .effect_stop("bandwidth must be a positive finite scalar.")
    }
    bandwidth_risk <- NULL
  }

  local_fit <- .local_linear_at(
    a = a[above],
    outcome = pseudo_for_fit,
    target = cutoff,
    bandwidth = bandwidth,
    kernel = kernel,
    degree = local_degree
  )
  theta1 <- local_fit$beta[1L]

  grid_u <- (a_grid - cutoff) / bandwidth
  grid_kernel <- kernel(grid_u) / bandwidth
  centered_mu <- sweep(nuisance$muhat_grid, 2L, m_grid, `-`)
  correction_base <- sweep(centered_mu, 2L, varpi_grid * grid_kernel, `*`)
  correction_intercept <- .trapezoid_rows(correction_base, a_grid)
  correction_slope <- .trapezoid_rows(
    sweep(correction_base, 2L, grid_u, `*`),
    a_grid
  )

  if (local_degree == 1L) {
    local_residual <- pseudo_for_fit -
      local_fit$beta[1L] - local_fit$beta[2L] * local_fit$u
    phi_moments <- rbind(
      local_fit$weights * local_residual + correction_intercept,
      local_fit$u * local_fit$weights * local_residual + correction_slope
    )
  } else {
    correction_quadratic <- .trapezoid_rows(
      sweep(correction_base, 2L, grid_u^2, `*`),
      a_grid
    )
    local_residual <- pseudo_for_fit -
      local_fit$beta[1L] - local_fit$beta[2L] * local_fit$u -
      local_fit$beta[3L] * local_fit$u^2
    phi_moments <- rbind(
      local_fit$weights * local_residual + correction_intercept,
      local_fit$u * local_fit$weights * local_residual + correction_slope,
      local_fit$u^2 * local_fit$weights * local_residual + correction_quadratic
    )
  }
  phi1_subpopulation <- drop(t(solve(local_fit$D, phi_moments))[, 1L])

  theta2 <- mean(below)
  theta3 <- mean(y[below])
  theta4 <- mean(y)
  if (!is.finite(theta4) || abs(theta4) < .Machine$double.eps) {
    .effect_stop("The truncated-effect ratio is undefined because mean(y) is zero.")
  }

  psi <- unname((theta1 * (1 - theta2) + theta2 * theta3) / theta4)
  phi1 <- rep(0, n)
  phi1[above] <- (n / n_above) * phi1_subpopulation
  phi2 <- as.numeric(below) - theta2
  phi3 <- rep(0, n)
  phi3[below] <- (y[below] - theta3) / theta2
  phi4 <- y - theta4
  phi_components <- cbind(phi1 = phi1, phi2 = phi2, phi3 = phi3, phi4 = phi4)

  gradient <- c(
    (1 - theta2) / theta4,
    (-theta1 + theta3) / theta4,
    theta2 / theta4,
    -(theta1 * (1 - theta2) + theta2 * theta3) / theta4^2
  )
  influence_curve <- drop(phi_components %*% gradient)
  variance <- stats::var(influence_curve) / n
  local_effective_n <- if (sum(local_fit$weights^2) > 0) {
    sum(local_fit$weights)^2 / sum(local_fit$weights^2)
  } else {
    0
  }
  diagnostics <- list(
    density_trim = density_trim,
    density_floor = density_floor,
    min_pihat = min(pihat),
    max_pihat = max(pihat),
    ipw_ratio_trim = ipw_ratio_trim,
    ipw_ratio_cap = ipw_ratio_cap,
    max_ipw_ratio_raw = max(ipw_ratio_raw),
    max_ipw_ratio_used = max(ipw_ratio),
    n_ipw_ratio_truncated = sum(ipw_ratio < ipw_ratio_raw),
    max_abs_pseudo_outcome_raw = max(abs(pseudo_raw)),
    bandwidth = bandwidth,
    local_design_condition = kappa(local_fit$D, exact = TRUE),
    local_effective_n = local_effective_n
  )

  structure(
    list(
      psi = psi,
      influence_curve = influence_curve,
      standard_error = sqrt(variance),
      variance = variance,
      theta = c(theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4),
      gradient = gradient,
      phi_components = phi_components,
      pseudo_outcome = pseudo_for_fit,
      pseudo_outcome_raw = pseudo_raw,
      above_index = which(above),
      bandwidth = bandwidth,
      bandwidth_risk = bandwidth_risk,
      nuisance = nuisance,
      diagnostics = diagnostics,
      n_uc = n_uc
    ),
    class = "truncated_effect_candidate"
  )
}
