# Doubly robust estimator of the truncated exposure effect
#
# Nuisance functions are estimated on the full supplied sample. 
# The estimator returns an observation-level influence
# curve, which is later used to compare nested spatial adjustments.

.name_adjustment_matrix <- function(w) {
  w <- as.matrix(w)
  if (!is.numeric(w) || nrow(w) < 1L || ncol(w) < 1L ||
      any(!is.finite(w))) {
    stop("w must be a finite numeric matrix with at least one column.")
  }
  
  if (is.null(colnames(w))) {
    colnames(w) <- paste0("x", seq_len(ncol(w)))
  }
  colnames(w) <- make.unique(make.names(colnames(w)))
  w
}


.trapezoid_rows <- function(values, grid) {
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
                             kernel = stats::dnorm) {
  u <- (a - target) / bandwidth
  weights <- kernel(u) / bandwidth
  design <- cbind(intercept = 1, slope = u)
  
  D <- crossprod(design, design * weights) / length(a)
  if (rcond(D) < 1e-10) {
    stop("The local-linear design is singular at this bandwidth.")
  }
  
  moment <- colMeans(design * (weights * outcome))
  beta <- drop(solve(D, moment))
  
  list(beta = beta, D = D, u = u, weights = weights)
}


.select_local_linear_bandwidth <- function(a, outcome, bw_seq,
                                           kernel = stats::dnorm) {
  if (!requireNamespace("KernSmooth", quietly = TRUE)) {
    stop("Install KernSmooth for automatic bandwidth selection.")
  }
  
  evaluation_grid <- seq(
    min(a), max(a),
    length.out = min(100L, length(a))
  )
  
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
      
      s2 * (kernel(0) / bandwidth) /
        (denominator * length(a))
    }, numeric(1))
    
    stats::approx(
      evaluation_grid,
      values,
      xout = a,
      rule = 2
    )$y
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
    
    if (is.null(fit)) return(rep(NA_real_, length(a)))
    stats::approx(fit$x, fit$y, xout = a, rule = 2)$y
  }
  
  risk <- vapply(bw_seq, function(bandwidth) {
    hats <- hat_values(bandwidth)
    fitted <- fitted_values(bandwidth)
    denominator <- 1 - hats
    usable <- is.finite(hats) & is.finite(fitted) &
      abs(denominator) > 1e-8
    
    if (mean(usable) < 0.5) return(NA_real_)
    mean(((outcome[usable] - fitted[usable]) / denominator[usable])^2)
  }, numeric(1))
  
  if (!any(is.finite(risk))) {
    stop("Bandwidth selection failed for every candidate bandwidth.")
  }
  
  selected_index <- which.min(replace(risk, !is.finite(risk), Inf))
  
  list(
    bandwidth = bw_seq[selected_index],
    risk = data.frame(bandwidth = bw_seq, risk = risk),
    at_lower_boundary = selected_index == 1L,
    at_upper_boundary = selected_index == length(bw_seq)
  )
}


# Full-sample Super Learner nuisance estimation.
estimate_nuisance_superlearner_full <- function(
    y,
    a,
    w,
    a_grid,
    sl_library = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction"),
    density_trim = 0.01,
    cv_control = list(V = 2L),
    variance_floor = 1e-4,
    variance_ceiling_multiplier = 100) {
  
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Install SuperLearner for nuisance estimation.")
  }
  
  w <- as.data.frame(.name_adjustment_matrix(w))
  n <- length(a)
  n_grid <- length(a_grid)
  
  # Prediction rows: observed exposure first, then every exposure-grid value.
  grid_rows <- rep(seq_len(n), times = n_grid)
  w_grid <- w[grid_rows, , drop = FALSE]
  w_new <- rbind(w, w_grid)
  observed_positions <- seq_len(n)
  grid_positions <- n + seq_len(n * n_grid)
  
  # Conditional exposure mean.
  mean_fit <- SuperLearner::SuperLearner(
    Y = a,
    X = w,
    newX = w_new,
    SL.library = sl_library,
    cvControl = cv_control,
    env = asNamespace("SuperLearner")
  )
  exposure_mean <- drop(mean_fit$SL.predict)
  exposure_residual <- a - exposure_mean[observed_positions]
  
  # Conditional exposure variance.
  variance_fit <- SuperLearner::SuperLearner(
    Y = log(pmax(exposure_residual^2, .Machine$double.eps)),
    X = w,
    newX = w_new,
    SL.library = sl_library,
    cvControl = cv_control,
    env = asNamespace("SuperLearner")
  )
  exposure_variance <- exp(drop(variance_fit$SL.predict))
  variance_ceiling <- max(
    variance_floor,
    variance_ceiling_multiplier * stats::var(exposure_residual)
  )
  exposure_variance <- pmin(
    pmax(exposure_variance, variance_floor),
    variance_ceiling
  )
  
  # Conditional outcome mean at observed and grid exposure values.
  wa_observed <- data.frame(w, a = a, check.names = FALSE)
  wa_grid <- data.frame(
    w_grid,
    a = rep(a_grid, each = n),
    check.names = FALSE
  )
  wa_new <- rbind(wa_observed, wa_grid)
  
  outcome_fit <- SuperLearner::SuperLearner(
    Y = y,
    X = wa_observed,
    newX = wa_new,
    SL.library = sl_library,
    cvControl = cv_control,
    env = asNamespace("SuperLearner")
  )
  outcome_mean <- drop(outcome_fit$SL.predict)
  
  # Location-scale estimate of the conditional exposure density.
  standardized_observed <- exposure_residual /
    sqrt(exposure_variance[observed_positions])
  standardized_grid <- (
    rep(a_grid, each = n) - exposure_mean[grid_positions]
  ) / sqrt(exposure_variance[grid_positions])
  
  density_range <- range(
    c(standardized_observed, standardized_grid),
    finite = TRUE
  )
  residual_density <- stats::density(
    standardized_observed,
    from = density_range[1L],
    to = density_range[2L]
  )
  evaluate_density <- function(values) {
    stats::approx(
      residual_density$x,
      residual_density$y,
      xout = values,
      rule = 2
    )$y
  }
  
  density_observed_raw <- evaluate_density(standardized_observed) /
    sqrt(exposure_variance[observed_positions])
  density_grid <- matrix(
    evaluate_density(standardized_grid) /
      sqrt(exposure_variance[grid_positions]),
    nrow = n,
    ncol = n_grid
  )
  
  density_floor <- max(
    as.numeric(stats::quantile(
      density_observed_raw,
      probs = density_trim,
      names = FALSE
    )),
    .Machine$double.eps
  )
  
  list(
    pihat = pmax(density_observed_raw, density_floor),
    pihat_raw = density_observed_raw,
    pihat_grid = density_grid,
    muhat = outcome_mean[observed_positions],
    muhat_grid = matrix(
      outcome_mean[grid_positions],
      nrow = n,
      ncol = n_grid
    ),
    density_floor = density_floor
  )
}


estimate_truncated_effect_full <- function(
    y,
    a,
    w,
    cutoff,
    nuisance_estimator = estimate_nuisance_superlearner_full,
    nuisance_args = list(),
    a_grid = NULL,
    n_grid = 25L,
    bandwidth = NULL,
    bw_seq = NULL,
    kernel = stats::dnorm,
    density_trim = 0.01,
    constrain = TRUE) {
  
  y <- as.numeric(y)
  a <- as.numeric(a)
  w <- .name_adjustment_matrix(w)
  n <- length(y)
  
  if (length(a) != n || nrow(w) != n) {
    stop("y, a, and w must contain the same number of observations.")
  }
  
  # These conventions match the primary analysis.
  above <- a >= cutoff - 1
  below <- a <= cutoff
  n_above <- sum(above)
  
  if (n_above < 6L || sum(below) < 2L) {
    stop("Too few observations remain around this cutoff.")
  }
  
  if (is.null(a_grid)) {
    a_grid <- seq(cutoff - 1, cutoff + 1, length.out = n_grid)
  }
  
  nuisance <- do.call(
    nuisance_estimator,
    utils::modifyList(
      list(
        y = y[above],
        a = a[above],
        w = w[above, , drop = FALSE],
        a_grid = a_grid,
        density_trim = density_trim
      ),
      nuisance_args
    )
  )
  
  pihat <- pmax(nuisance$pihat, nuisance$density_floor)
  marginal_density_grid <- pmax(
    colMeans(nuisance$pihat_grid),
    .Machine$double.eps
  )
  marginal_outcome_grid <- colMeans(nuisance$muhat_grid)
  
  marginal_density_observed <- pmax(
    stats::approx(
      a_grid,
      marginal_density_grid,
      xout = a[above],
      rule = 2
    )$y,
    nuisance$density_floor
  )
  marginal_outcome_observed <- stats::approx(
    a_grid,
    marginal_outcome_grid,
    xout = a[above],
    rule = 2
  )$y
  
  ipw_ratio <- marginal_density_observed / pihat
  pseudo_raw <-
    (y[above] - nuisance$muhat) * ipw_ratio +
    marginal_outcome_observed
  
  pseudo <- pseudo_raw
  if (constrain) {
    outcome_range <- range(y)
    pseudo <- pmin(pmax(pseudo, outcome_range[1L]), outcome_range[2L])
  }
  
  if (is.null(bandwidth)) {
    if (is.null(bw_seq)) {
      exposure_sd <- stats::sd(a)
      # Same grid in ../simulation and ../data_application.
      bw_seq <- seq(exposure_sd / 10, exposure_sd, length.out = 100L)
    }
    bandwidth_fit <- .select_local_linear_bandwidth(
      a = a[above],
      outcome = pseudo,
      bw_seq = bw_seq,
      kernel = kernel
    )
    bandwidth <- bandwidth_fit$bandwidth
  } else {
    bandwidth_fit <- NULL
  }
  
  local_fit <- .local_linear_at(
    a = a[above],
    outcome = pseudo,
    target = cutoff,
    bandwidth = bandwidth,
    kernel = kernel
  )
  theta1 <- local_fit$beta[1L]
  
  # Influence-function contribution for the local-linear intercept.
  grid_u <- (a_grid - cutoff) / bandwidth
  grid_kernel <- kernel(grid_u) / bandwidth
  centered_mu <- sweep(
    nuisance$muhat_grid,
    2L,
    marginal_outcome_grid,
    `-`
  )
  correction_base <- sweep(
    centered_mu,
    2L,
    marginal_density_grid * grid_kernel,
    `*`
  )
  correction_intercept <- .trapezoid_rows(correction_base, a_grid)
  correction_slope <- .trapezoid_rows(
    sweep(correction_base, 2L, grid_u, `*`),
    a_grid
  )
  
  local_residual <- pseudo -
    local_fit$beta[1L] -
    local_fit$beta[2L] * local_fit$u
  phi_moments <- rbind(
    local_fit$weights * local_residual + correction_intercept,
    local_fit$u * local_fit$weights * local_residual + correction_slope
  )
  phi1_above <- drop(
    t(solve(local_fit$D, phi_moments))[, 1L]
  )
  
  theta2 <- mean(below)
  theta3 <- mean(y[below])
  theta4 <- mean(y)
  
  psi <- unname(
    (theta1 * (1 - theta2) + theta2 * theta3) / theta4
  )
  
  phi1 <- rep(0, n)
  phi1[above] <- (n / n_above) * phi1_above
  phi2 <- as.numeric(below) - theta2
  phi3 <- rep(0, n)
  phi3[below] <- (y[below] - theta3) / theta2
  phi4 <- y - theta4
  
  gradient <- c(
    (1 - theta2) / theta4,
    (-theta1 + theta3) / theta4,
    theta2 / theta4,
    -(theta1 * (1 - theta2) + theta2 * theta3) / theta4^2
  )
  
  influence_curve <- drop(
    cbind(phi1, phi2, phi3, phi4) %*% gradient
  )
  standard_error <- stats::sd(influence_curve) / sqrt(n)
  confidence_interval <- psi +
    c(-1, 1) * stats::qnorm(0.975) * standard_error
  
  local_effective_n <-
    sum(local_fit$weights)^2 / sum(local_fit$weights^2)
  
  list(
    psi = psi,
    standard_error = standard_error,
    confidence_interval = confidence_interval,
    influence_curve = influence_curve,
    bandwidth = bandwidth,
    diagnostics = list(
      n = n,
      n_above = n_above,
      density_floor = nuisance$density_floor,
      max_ipw_ratio = max(ipw_ratio),
      proportion_constrained = mean(pseudo != pseudo_raw),
      local_effective_n = local_effective_n,
      bandwidth_at_lower_boundary = if (is.null(bandwidth_fit)) {
        NA
      } else {
        bandwidth_fit$at_lower_boundary
      },
      bandwidth_at_upper_boundary = if (is.null(bandwidth_fit)) {
        NA
      } else {
        bandwidth_fit$at_upper_boundary
      }
    )
  )
}
