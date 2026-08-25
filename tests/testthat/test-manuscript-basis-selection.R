testthat::test_that("candidate grid follows the supplement formula", {
  testthat::expect_identical(
    make_candidate_grid(m = 100L),
    c(80L, 85L, 90L, 95L)
  )

  # TPS data-application example: 50 radial columns plus 3 affine columns.
  testthat::expect_identical(
    make_candidate_grid(m = 53L),
    c(42L, 47L)
  )

  grid <- make_candidate_grid(m = 250L)
  testthat::expect_true(all(250L - grid >= max(5L, floor(0.02 * 250L))))
})

testthat::test_that("candidate grid rejects an incompatible core and minimum", {
  testthat::expect_error(
    make_candidate_grid(m = 20L),
    "core set"
  )
  testthat::expect_error(
    make_candidate_grid(m = 50L, n_uc_star = 49L),
    "exceeds n_max"
  )
})

testthat::test_that("basis columns are partitioned in both scale directions", {
  small <- partition_basis_columns(m = 10L, n_uc = 3L, "small")
  testthat::expect_identical(small$confounded, 1:7)
  testthat::expect_identical(small$instruments, 8:10)

  large <- partition_basis_columns(m = 10L, n_uc = 3L, "large")
  testthat::expect_identical(large$instruments, 1:3)
  testthat::expect_identical(large$confounded, 4:10)
})

testthat::test_that("A-c projection is fit on training rows and has no holdout leakage", {
  set.seed(11)
  n <- 30L
  m <- 8L
  B <- qr.Q(qr(matrix(rnorm(n * m), nrow = n, ncol = m)))
  A <- drop(B %*% seq_len(m)) + rnorm(n, sd = 0.1)
  train <- 1:22

  fit <- project_Ac(B, A, n_uc = 3L, train_idx = train)
  Bc_train <- B[train, 1:5, drop = FALSE]
  gamma_manual <- drop(qr.solve(Bc_train, A[train]))

  testthat::expect_equal(fit$gamma, gamma_manual, tolerance = 1e-12)
  testthat::expect_equal(
    fit$Ac,
    drop(B[, 1:5, drop = FALSE] %*% gamma_manual),
    tolerance = 1e-12
  )

  A_perturbed <- A
  A_perturbed[-train] <- A_perturbed[-train] + 1e6
  perturbed_fit <- project_Ac(B, A_perturbed, n_uc = 3L, train_idx = train)

  testthat::expect_equal(perturbed_fit$gamma, fit$gamma, tolerance = 0)
  testthat::expect_equal(perturbed_fit$Ac, fit$Ac, tolerance = 0)
})

testthat::test_that("paired SE matches the centered influence-function formula", {
  phi_core <- c(-1, 0, 1, 0)
  phi_candidate <- c(-0.5, 0.5, 1.5, 0.5)
  difference <- phi_candidate - phi_core

  testthat::expect_equal(
    paired_difference_se(phi_candidate, phi_core),
    stats::sd(difference) / sqrt(length(difference)),
    tolerance = 1e-14
  )
})

testthat::test_that("selection is core-referenced and stops at the first failure", {
  n <- 20L
  signs <- rep(c(-1, 1), length.out = n)
  phi <- cbind(
    core = rep(0, n),
    first = 0.2 * signs,
    second = 0.4 * signs,
    failure = 0.1 * signs
  )
  selection <- select_contiguous_candidate(
    n_uc_cands = c(80L, 85L, 90L, 95L),
    psi = c(1, 1.01, 1.02, 1.4),
    influence_curves = phi,
    alpha = 1
  )

  testthat::expect_identical(selection$selected_n_uc, 90L)
  testthat::expect_identical(
    selection$diagnostics$in_contiguous_sequence,
    c(TRUE, TRUE, TRUE, FALSE)
  )

  noncontiguous <- select_contiguous_candidate(
    n_uc_cands = c(80L, 85L, 90L),
    psi = c(1, 1.2, 1.001),
    influence_curves = cbind(rep(0, n), 0.1 * signs, signs),
    alpha = 1
  )
  testthat::expect_identical(noncontiguous$selected_n_uc, 80L)
  testthat::expect_identical(
    noncontiguous$diagnostics$in_contiguous_sequence,
    c(TRUE, FALSE, FALSE)
  )
})

testthat::test_that("outer-fold selector exposes only training data to candidates", {
  set.seed(19)
  n <- 24L
  m <- 10L
  B <- qr.Q(qr(matrix(rnorm(n * m), nrow = n, ncol = m)))
  a <- drop(B %*% seq_len(m)) + rnorm(n, sd = 0.1)
  y <- 2 + a + rnorm(n)
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  train <- 1:18
  inner <- rep(1:2, length.out = length(train))
  observed_sizes <- integer(0)

  candidate_estimator <- function(y, a, x, Ac, cutoff, folds, n_uc) {
    observed_sizes <<- c(observed_sizes, length(y))
    testthat::expect_length(a, length(train))
    testthat::expect_identical(nrow(x), length(train))
    testthat::expect_length(Ac, length(train))
    testthat::expect_identical(folds, inner)

    signs <- rep(c(-1, 1), length.out = length(y))
    psi <- switch(
      as.character(n_uc),
      `5` = 1,
      `6` = 1.01,
      `7` = 1.5
    )
    scale <- switch(
      as.character(n_uc),
      `5` = 0,
      `6` = 0.2,
      `7` = 0.1
    )
    list(psi = psi, influence_curve = scale * signs)
  }

  result <- select_basis_outer_fold(
    y = y,
    a = a,
    x = x,
    B = B,
    outer_train_idx = train,
    n_uc_cands = 5:7,
    cutoff = 0,
    estimate_candidate = candidate_estimator,
    inner_folds = inner,
    alpha = 1
  )

  testthat::expect_identical(observed_sizes, rep(length(train), 3L))
  testthat::expect_identical(result$selected_n_uc, 6L)
  testthat::expect_length(result$selected_projection$Ac, n)

  a_perturbed <- a
  a_perturbed[-train] <- a_perturbed[-train] + 1e6
  result_perturbed <- select_basis_outer_fold(
    y = y,
    a = a_perturbed,
    x = x,
    B = B,
    outer_train_idx = train,
    n_uc_cands = 5:7,
    cutoff = 0,
    estimate_candidate = candidate_estimator,
    inner_folds = inner,
    alpha = 1
  )

  testthat::expect_equal(
    result_perturbed$selected_projection$gamma,
    result$selected_projection$gamma,
    tolerance = 0
  )
})

testthat::test_that("uniform folds are balanced", {
  set.seed(23)
  folds <- make_uniform_folds(103L, K = 5L)
  testthat::expect_identical(sort(unique(folds)), 1:5)
  testthat::expect_lte(max(tabulate(folds)) - min(tabulate(folds)), 1L)
})

testthat::test_that("outer-fold candidates cannot be silently truncated", {
  set.seed(29)
  B <- qr.Q(qr(matrix(rnorm(80), nrow = 10, ncol = 8)))
  estimator <- function(y, a, x, Ac, cutoff, folds, n_uc) {
    list(psi = 1, influence_curve = rep(0, length(y)))
  }

  testthat::expect_error(
    select_basis_outer_fold(
      y = rnorm(10),
      a = rnorm(10),
      B = B,
      outer_train_idx = 1:8,
      n_uc_cands = c(5.5, 6.5),
      cutoff = 0,
      estimate_candidate = estimator,
      inner_folds = rep(1:2, 4)
    ),
    "integer vector"
  )
})
