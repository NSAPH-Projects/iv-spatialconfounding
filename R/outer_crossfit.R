# Five-fold CV for estimation of the truncated exposure effect.
#
# Dependencies:
#   R/basis_selection.R
#   R/truncated_effect.R

.validate_external_nuisance_result <- function(result, n_target, n_grid) {
  required <- c("pihat", "muhat", "muhat_grid", "varpi_grid", "m_grid")
  if (!is.list(result) || !all(required %in% names(result))) {
    .effect_stop(sprintf(
      "The outer-fold nuisance estimator must return: %s.",
      paste(required, collapse = ", ")
    ))
  }
  if (!is.numeric(result$pihat) || length(result$pihat) != n_target ||
      any(!is.finite(result$pihat)) || any(result$pihat < 0) ||
      !is.numeric(result$muhat) || length(result$muhat) != n_target ||
      any(!is.finite(result$muhat))) {
    .effect_stop(
      "Outer-fold pihat and muhat must be finite target-length vectors and pihat cannot be negative."
    )
  }
  if (!is.matrix(result$muhat_grid) || !is.numeric(result$muhat_grid) ||
      !identical(dim(result$muhat_grid), c(n_target, n_grid)) ||
      any(!is.finite(result$muhat_grid))) {
    .effect_stop("Outer-fold muhat_grid must be a finite target-by-grid matrix.")
  }
  if (!is.numeric(result$varpi_grid) || length(result$varpi_grid) != n_grid ||
      any(!is.finite(result$varpi_grid)) || any(result$varpi_grid < 0) ||
      !is.numeric(result$m_grid) || length(result$m_grid) != n_grid ||
      any(!is.finite(result$m_grid))) {
    .effect_stop(
      "Outer-fold varpi_grid and m_grid must be finite grid-length vectors and varpi_grid cannot be negative."
    )
  }
  invisible(result)
}

.name_adjustment_matrix <- function(w) {
  w <- as.matrix(w)
  w_names <- colnames(w)
  if (is.null(w_names)) {
    w_names <- rep("", ncol(w))
  }
  unnamed <- is.na(w_names) | !nzchar(w_names)
  w_names[unnamed] <- paste0("x", which(unnamed))
  colnames(w) <- make.unique(make.names(w_names))
  w
}

# Fit the nuisance functions on one outer training set and predict the held-out
# observations. Marginal density and outcome regressions are averaged over the
# training covariate distribution, matching the sums over T_k in the supplement.
estimate_nuisance_superlearner_external <- function(
    y_train, a_train, w_train, a_target, w_target, a_grid,
    sl_library = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction"),
    density_trim = 0,
    cv_control = list(V = 2L),
    variance_floor = 1e-4) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    .effect_stop("SuperLearner is required for the nuisance estimator.")
  }

  w_train <- as.data.frame(.name_adjustment_matrix(w_train))
  w_target <- as.data.frame(.name_adjustment_matrix(w_target))
  if (!identical(names(w_train), names(w_target))) {
    .effect_stop("w_train and w_target must have identical columns.")
  }
  n_train <- length(a_train)
  n_target <- length(a_target)
  n_grid <- length(a_grid)
  if (n_train < 6L || n_target < 1L || length(y_train) != n_train ||
      nrow(w_train) != n_train || nrow(w_target) != n_target ||
      any(!is.finite(y_train)) || any(!is.finite(a_train)) ||
      any(!is.finite(a_target))) {
    .effect_stop(
      "The outer nuisance fit requires at least six training rows and one finite target row."
    )
  }
  if (!is.numeric(density_trim) || length(density_trim) != 1L ||
      !is.finite(density_trim) || density_trim < 0 || density_trim >= 0.5) {
    .effect_stop("density_trim must be a finite scalar in [0, 0.5).")
  }
  if (!is.numeric(variance_floor) || length(variance_floor) != 1L ||
      !is.finite(variance_floor) || variance_floor <= 0) {
    .effect_stop("variance_floor must be a positive finite scalar.")
  }

  target_grid <- w_target[
    rep(seq_len(n_target), times = n_grid),
    ,
    drop = FALSE
  ]
  reference_grid <- w_train[
    rep(seq_len(n_train), times = n_grid),
    ,
    drop = FALSE
  ]
  w_new <- rbind(w_train, w_target, target_grid, reference_grid)

  train_positions <- seq_len(n_train)
  target_positions <- n_train + seq_len(n_target)
  target_grid_positions <- n_train + n_target + seq_len(n_target * n_grid)
  reference_grid_positions <- n_train + n_target + n_target * n_grid +
    seq_len(n_train * n_grid)

  mean_fit <- SuperLearner::SuperLearner(
    Y = a_train,
    X = w_train,
    newX = w_new,
    SL.library = sl_library,
    cvControl = cv_control,
    env = asNamespace("SuperLearner")
  )
  mean_prediction <- drop(mean_fit$SL.predict)
  training_residual <- a_train - mean_prediction[train_positions]

  variance_fit <- SuperLearner::SuperLearner(
    Y = log(pmax(training_residual^2, .Machine$double.eps)),
    X = w_train,
    newX = w_new,
    SL.library = sl_library,
    cvControl = cv_control,
    env = asNamespace("SuperLearner")
  )
  variance_prediction <- pmax(
    exp(drop(variance_fit$SL.predict)),
    variance_floor
  )

  new_exposure <- c(
    a_train,
    a_target,
    rep(a_grid, each = n_target),
    rep(a_grid, each = n_train)
  )
  wa_train <- data.frame(w_train, a = a_train, check.names = FALSE)
  wa_new <- data.frame(w_new, a = new_exposure, check.names = FALSE)
  outcome_fit <- SuperLearner::SuperLearner(
    Y = y_train,
    X = wa_train,
    newX = wa_new,
    SL.library = sl_library,
    cvControl = cv_control,
    env = asNamespace("SuperLearner")
  )
  outcome_prediction <- drop(outcome_fit$SL.predict)

  standardized_train <- training_residual /
    sqrt(variance_prediction[train_positions])
  standardized_target <- (
    a_target - mean_prediction[target_positions]
  ) / sqrt(variance_prediction[target_positions])
  standardized_target_grid <- (
    rep(a_grid, each = n_target) - mean_prediction[target_grid_positions]
  ) / sqrt(variance_prediction[target_grid_positions])
  standardized_reference_grid <- (
    rep(a_grid, each = n_train) - mean_prediction[reference_grid_positions]
  ) / sqrt(variance_prediction[reference_grid_positions])

  density_range <- range(
    c(
      standardized_train,
      standardized_target,
      standardized_target_grid,
      standardized_reference_grid
    ),
    finite = TRUE
  )
  residual_density <- stats::density(
    standardized_train,
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

  pihat_train <- evaluate_density(standardized_train) /
    sqrt(variance_prediction[train_positions])
  density_floor <- max(
    as.numeric(stats::quantile(
      pihat_train,
      probs = density_trim,
      names = FALSE
    )),
    .Machine$double.eps
  )
  pihat <- pmax(
    evaluate_density(standardized_target) /
      sqrt(variance_prediction[target_positions]),
    density_floor
  )
  pihat_reference_grid <- matrix(
    evaluate_density(standardized_reference_grid) /
      sqrt(variance_prediction[reference_grid_positions]),
    nrow = n_train,
    ncol = n_grid
  )
  muhat_grid <- matrix(
    outcome_prediction[target_grid_positions],
    nrow = n_target,
    ncol = n_grid
  )
  muhat_reference_grid <- matrix(
    outcome_prediction[reference_grid_positions],
    nrow = n_train,
    ncol = n_grid
  )

  list(
    pihat = pihat,
    muhat = outcome_prediction[target_positions],
    muhat_grid = muhat_grid,
    varpi_grid = colMeans(pihat_reference_grid),
    m_grid = colMeans(muhat_reference_grid),
    density_floor = density_floor
  )
}

# Run the complete selection-and-estimation algorithm across outer folds.
estimate_selected_truncated_effect <- function(
    y, a, x = NULL, B = NULL, cutoff,
    n_uc_cands = if (!is.null(B)) make_candidate_grid(ncol(B)) else NULL,
    outer_folds = NULL,
    K = 5L,
    alpha = 1,
    density_trim = 0,
    ipw_ratio_trim = 0,
    instrument_scale = c("small", "large"),
    candidate_estimator = estimate_truncated_effect_candidate,
    candidate_args = list(),
    final_nuisance_estimator = estimate_nuisance_superlearner_external,
    final_nuisance_args = list(),
    a_grid = NULL,
    n_grid = 100L,
    bandwidth = NULL,
    bw_seq = NULL,
    kernel = stats::dnorm,
    constrain = FALSE,
    local_degree = 1L,
    keep_fold_details = FALSE,
    fixed_adjustment = NULL) {
  use_fixed_adjustment <- !is.null(fixed_adjustment)
  if (!use_fixed_adjustment) {
    .validate_basis(B)
  }
  instrument_scale <- match.arg(instrument_scale)
  if (use_fixed_adjustment) {
    fixed_adjustment <- .name_adjustment_matrix(fixed_adjustment)
    n <- nrow(fixed_adjustment)
    if (!is.numeric(fixed_adjustment) || ncol(fixed_adjustment) < 1L ||
        any(!is.finite(fixed_adjustment))) {
      .effect_stop(
        "fixed_adjustment must be a finite numeric matrix with at least one column."
      )
    }
  } else {
    n <- nrow(B)
  }

  if (!is.numeric(y) || !is.numeric(a) || length(y) != n || length(a) != n ||
      any(!is.finite(y)) || any(!is.finite(a))) {
    .effect_stop("y and a must be finite numeric vectors with one value per row of B.")
  }
  if (is.null(x)) {
    x <- matrix(numeric(0), nrow = n, ncol = 0L)
  } else {
    x <- as.matrix(x)
    if (!is.numeric(x) || nrow(x) != n || any(!is.finite(x))) {
      .effect_stop("x must be NULL or a finite numeric matrix with one row per row of B.")
    }
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff)) {
    .effect_stop("cutoff must be a finite numeric scalar.")
  }
  if ((!use_fixed_adjustment && !is.function(candidate_estimator)) ||
      !is.function(final_nuisance_estimator)) {
    .effect_stop(
      "candidate_estimator and final_nuisance_estimator must be functions when used."
    )
  }
  if (!is.list(candidate_args) || !is.list(final_nuisance_args)) {
    .effect_stop("candidate_args and final_nuisance_args must be lists.")
  }
  if (!is.function(kernel)) {
    .effect_stop("kernel must be a function.")
  }
  if (!is.logical(constrain) || length(constrain) != 1L || is.na(constrain) ||
      !is.logical(keep_fold_details) || length(keep_fold_details) != 1L ||
      is.na(keep_fold_details)) {
    .effect_stop("constrain and keep_fold_details must each be TRUE or FALSE.")
  }
  if (!use_fixed_adjustment &&
      (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
       alpha <= 0)) {
    .effect_stop("alpha must be a positive finite scalar.")
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

  if (is.null(outer_folds)) {
    outer_folds <- make_uniform_folds(n, K = K)
  }
  if (!is.numeric(outer_folds) || length(outer_folds) != n ||
      any(!is.finite(outer_folds)) || length(unique(outer_folds)) < 2L) {
    .effect_stop("outer_folds must assign every observation to at least two folds.")
  }
  fold_values <- sort(unique(outer_folds))

  if (!use_fixed_adjustment) {
    if (!is.numeric(n_uc_cands) || length(n_uc_cands) < 1L ||
        any(!is.finite(n_uc_cands)) ||
        any(n_uc_cands != as.integer(n_uc_cands)) ||
        any(n_uc_cands < 1L) || any(n_uc_cands >= ncol(B)) ||
        is.unsorted(n_uc_cands, strictly = TRUE)) {
      .effect_stop(
        "n_uc_cands must be a strictly increasing integer vector between 1 and ncol(B) - 1."
      )
    }
    n_uc_cands <- as.integer(n_uc_cands)
  }

  above <- a >= cutoff
  below <- a < cutoff
  n_above <- sum(above)
  if (n_above < 6L || sum(below) < 2L || diff(range(a[above])) <= 0) {
    .effect_stop(
      "The cutoff must leave varying exposure for at least six observations above and two below."
    )
  }
  if (any(vapply(fold_values, function(fold) {
    sum(outer_folds == fold & above) < 1L ||
      sum(outer_folds != fold & above) < 6L
  }, logical(1)))) {
    .effect_stop(
      "Every outer fold needs at least one held-out and six training observations with A >= cutoff."
    )
  }

  if (is.null(a_grid)) {
    if (!.effect_integerish(n_grid) || n_grid < 5L) {
      .effect_stop("n_grid must be an integer greater than or equal to 5.")
    }
    a_grid <- seq(cutoff, max(a[above]), length.out = as.integer(n_grid))
  } else if (!is.numeric(a_grid) || length(a_grid) < 5L ||
             any(!is.finite(a_grid)) || is.unsorted(a_grid, strictly = TRUE) ||
             abs(a_grid[1L] - cutoff) > sqrt(.Machine$double.eps) ||
             max(a_grid) < max(a[above])) {
    .effect_stop(
      "a_grid must be strictly increasing, begin at cutoff, and span max(A | A >= cutoff)."
    )
  }

  selected_n_uc <- if (use_fixed_adjustment) {
    NULL
  } else {
    stats::setNames(integer(length(fold_values)), as.character(fold_values))
  }
  selection_diagnostics <- if (use_fixed_adjustment) {
    NULL
  } else {
    stats::setNames(vector("list", length(fold_values)), names(selected_n_uc))
  }
  fold_records <- vector("list", length(fold_values))
  names(fold_records) <- as.character(fold_values)
  Ac_crossfit <- if (use_fixed_adjustment) NULL else rep(NA_real_, n)
  pseudo_raw_full <- pseudo_fit_full <- rep(NA_real_, n)

  for (fold_index in seq_along(fold_values)) {
    fold <- fold_values[fold_index]
    train <- which(outer_folds != fold)
    holdout <- which(outer_folds == fold)
    holdout_above <- holdout[above[holdout]]
    train_above <- train[above[train]]
    if (use_fixed_adjustment) {
      w_fold <- fixed_adjustment
    } else {
      inner_folds <- make_uniform_folds(length(train), K = 2L)
      selection <- do.call(
        select_basis_outer_fold,
        c(
          list(
            y = y,
            a = a,
            x = x,
            B = B,
            outer_train_idx = train,
            n_uc_cands = n_uc_cands,
            cutoff = cutoff,
            estimate_candidate = candidate_estimator,
            inner_folds = inner_folds,
            alpha = alpha,
            instrument_scale = instrument_scale
          ),
          candidate_args
        )
      )
      selected_n_uc[fold_index] <- selection$selected_n_uc
      candidate_diagnostics <- do.call(rbind, lapply(
        selection$candidate_results,
        function(candidate) {
          values <- candidate$diagnostics
          value_or_na <- function(name) {
            value <- values[[name]]
            if (is.numeric(value) && length(value) == 1L) value else NA_real_
          }
          data.frame(
            density_floor = value_or_na("density_floor"),
            min_pihat = value_or_na("min_pihat"),
            max_ipw_ratio_raw = value_or_na("max_ipw_ratio_raw"),
            max_ipw_ratio_used = value_or_na("max_ipw_ratio_used"),
            n_ipw_ratio_truncated = value_or_na("n_ipw_ratio_truncated"),
            max_abs_pseudo_outcome_raw =
              value_or_na("max_abs_pseudo_outcome_raw"),
            bandwidth = value_or_na("bandwidth"),
            local_design_condition = value_or_na("local_design_condition"),
            local_effective_n = value_or_na("local_effective_n")
          )
        }
      ))
      selection_diagnostics[[fold_index]] <- cbind(
        selection$selection$diagnostics,
        candidate_diagnostics
      )
      Ac_fold <- selection$selected_projection$Ac
      Ac_crossfit[holdout] <- Ac_fold[holdout]
      w_fold <- .name_adjustment_matrix(cbind(Ac = Ac_fold, x))
    }
    nuisance <- do.call(
      final_nuisance_estimator,
      c(
        list(
          y_train = y[train_above],
          a_train = a[train_above],
          w_train = w_fold[train_above, , drop = FALSE],
          a_target = a[holdout_above],
          w_target = w_fold[holdout_above, , drop = FALSE],
          a_grid = a_grid
        ),
        final_nuisance_args
      )
    )
    .validate_external_nuisance_result(
      nuisance,
      n_target = length(holdout_above),
      n_grid = length(a_grid)
    )

    pihat <- pmax(nuisance$pihat, .Machine$double.eps)
    varpi_grid <- pmax(nuisance$varpi_grid, .Machine$double.eps)
    varpi_observed <- pmax(
      stats::approx(
        a_grid,
        varpi_grid,
        xout = a[holdout_above],
        rule = 2
      )$y,
      .Machine$double.eps
    )
    m_observed <- stats::approx(
      a_grid,
      nuisance$m_grid,
      xout = a[holdout_above],
      rule = 2
    )$y
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
    pseudo_raw <- (
      y[holdout_above] - nuisance$muhat
    ) * ipw_ratio + m_observed
    pseudo_fit <- pseudo_raw
    if (constrain) {
      outcome_range <- range(y)
      pseudo_fit <- pmin(pmax(pseudo_fit, outcome_range[1L]), outcome_range[2L])
    }
    pseudo_raw_full[holdout_above] <- pseudo_raw
    pseudo_fit_full[holdout_above] <- pseudo_fit

    fold_records[[fold_index]] <- list(
      fold = fold,
      train_index = train,
      holdout_index = holdout,
      holdout_above_index = holdout_above,
      muhat_grid = nuisance$muhat_grid,
      m_grid = nuisance$m_grid,
      varpi_grid = varpi_grid,
      nuisance_summary = list(
        density_floor = nuisance$density_floor,
        min_pihat = min(pihat),
        max_ipw_ratio_raw = max(ipw_ratio_raw),
        max_ipw_ratio_used = max(ipw_ratio),
        ipw_ratio_cap = ipw_ratio_cap,
        n_ipw_ratio_truncated = sum(ipw_ratio < ipw_ratio_raw),
        max_abs_pseudo_outcome_raw = max(abs(pseudo_raw)),
        n_train_above = length(train_above),
        n_holdout_above = length(holdout_above)
      )
    )
  }

  if ((!use_fixed_adjustment && any(!is.finite(Ac_crossfit))) ||
      any(!is.finite(pseudo_raw_full[above])) ||
      any(!is.finite(pseudo_fit_full[above]))) {
    .effect_stop("Outer-fold pooling failed to produce complete cross-fitted values.")
  }

  if (is.null(bandwidth)) {
    if (is.null(bw_seq)) {
      exposure_sd <- stats::sd(a[above])
      bw_seq <- seq(exposure_sd / 10, exposure_sd, length.out = 100L)
    }
    bandwidth_fit <- .select_local_linear_bandwidth(
      a = a[above],
      outcome = pseudo_fit_full[above],
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
    outcome = pseudo_fit_full[above],
    target = cutoff,
    bandwidth = bandwidth,
    kernel = kernel,
    degree = local_degree
  )
  theta1 <- local_fit$beta[1L]
  above_index <- which(above)
  correction_intercept <- correction_slope <- rep(NA_real_, n_above)
  correction_quadratic <- if (local_degree == 2L) rep(NA_real_, n_above) else NULL
  grid_u <- (a_grid - cutoff) / bandwidth
  grid_kernel <- kernel(grid_u) / bandwidth

  for (record in fold_records) {
    positions <- match(record$holdout_above_index, above_index)
    centered_mu <- sweep(record$muhat_grid, 2L, record$m_grid, `-`)
    correction_base <- sweep(
      centered_mu,
      2L,
      record$varpi_grid * grid_kernel,
      `*`
    )
    correction_intercept[positions] <- .trapezoid_rows(correction_base, a_grid)
    correction_slope[positions] <- .trapezoid_rows(
      sweep(correction_base, 2L, grid_u, `*`),
      a_grid
    )
    if (local_degree == 2L) {
      correction_quadratic[positions] <- .trapezoid_rows(
        sweep(correction_base, 2L, grid_u^2, `*`),
        a_grid
      )
    }
  }

  if (local_degree == 1L) {
    local_residual <- pseudo_fit_full[above] -
      local_fit$beta[1L] - local_fit$beta[2L] * local_fit$u
    phi_moments <- rbind(
      local_fit$weights * local_residual + correction_intercept,
      local_fit$u * local_fit$weights * local_residual + correction_slope
    )
  } else {
    local_residual <- pseudo_fit_full[above] -
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
  standard_error <- sqrt(variance)
  local_effective_n <- if (sum(local_fit$weights^2) > 0) {
    sum(local_fit$weights)^2 / sum(local_fit$weights^2)
  } else {
    0
  }

  fold_summaries <- lapply(fold_records, function(record) {
    list(
      fold = record$fold,
      n_train = length(record$train_index),
      n_holdout = length(record$holdout_index),
      n_train_above = record$nuisance_summary$n_train_above,
      n_holdout_above = record$nuisance_summary$n_holdout_above,
      density_floor = record$nuisance_summary$density_floor,
      min_pihat = record$nuisance_summary$min_pihat,
      max_ipw_ratio_raw = record$nuisance_summary$max_ipw_ratio_raw,
      max_ipw_ratio_used = record$nuisance_summary$max_ipw_ratio_used,
      ipw_ratio_cap = record$nuisance_summary$ipw_ratio_cap,
      n_ipw_ratio_truncated = record$nuisance_summary$n_ipw_ratio_truncated,
      max_abs_pseudo_outcome_raw =
        record$nuisance_summary$max_abs_pseudo_outcome_raw
    )
  })
  diagnostics <- list(
    density_trim = density_trim,
    ipw_ratio_trim = ipw_ratio_trim,
    max_abs_pseudo_outcome_raw = max(abs(pseudo_raw_full[above])),
    bandwidth = bandwidth,
    local_design_condition = kappa(local_fit$D, exact = TRUE),
    local_effective_n = local_effective_n
  )

  result <- structure(
    list(
      psi = psi,
      standard_error = standard_error,
      variance = variance,
      confidence_interval = psi + c(-1, 1) * stats::qnorm(0.975) * standard_error,
      theta = c(theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4),
      gradient = gradient,
      influence_curve = influence_curve,
      phi_components = phi_components,
      selected_n_uc = selected_n_uc,
      selection_diagnostics = selection_diagnostics,
      Ac_crossfit = Ac_crossfit,
      pseudo_outcome = pseudo_fit_full[above],
      pseudo_outcome_raw = pseudo_raw_full[above],
      above_index = above_index,
      outer_folds = outer_folds,
      bandwidth = bandwidth,
      bandwidth_risk = bandwidth_risk,
      diagnostics = diagnostics,
      fold_summaries = fold_summaries,
      fold_details = if (keep_fold_details) fold_records else NULL
    ),
    class = if (use_fixed_adjustment) {
      "crossfit_truncated_effect"
    } else {
      "selected_truncated_effect"
    }
  )
  result
}

# Apply the same five-fold nuisance, pseudo-outcome, local-linear, and
# influence-function estimator without basis selection. This is used by the
# oracle and other non-IV simulation methods so that methods differ only in
# their adjustment variables.
estimate_crossfit_truncated_effect <- function(
    y, a, w, cutoff,
    outer_folds = NULL,
    K = 5L,
    density_trim = 0,
    ipw_ratio_trim = 0,
    final_nuisance_estimator = estimate_nuisance_superlearner_external,
    final_nuisance_args = list(),
    a_grid = NULL,
    n_grid = 100L,
    bandwidth = NULL,
    bw_seq = NULL,
    kernel = stats::dnorm,
    constrain = FALSE,
    local_degree = 1L,
    keep_fold_details = FALSE) {
  estimate_selected_truncated_effect(
    y = y,
    a = a,
    x = NULL,
    B = NULL,
    cutoff = cutoff,
    n_uc_cands = NULL,
    outer_folds = outer_folds,
    K = K,
    density_trim = density_trim,
    ipw_ratio_trim = ipw_ratio_trim,
    final_nuisance_estimator = final_nuisance_estimator,
    final_nuisance_args = final_nuisance_args,
    a_grid = a_grid,
    n_grid = n_grid,
    bandwidth = bandwidth,
    bw_seq = bw_seq,
    kernel = kernel,
    constrain = constrain,
    local_degree = local_degree,
    keep_fold_details = keep_fold_details,
    fixed_adjustment = w
  )
}
