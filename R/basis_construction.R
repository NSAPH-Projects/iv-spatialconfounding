# Ordered spatial bases
#
# Every returned basis is ordered from the largest/smoothest spatial scale
# to smaller/more localized scales. The first column is the constant vector.

build_gl_basis_truncated <- function(L, k_max) {
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("Install RSpectra to construct the graph-Laplacian basis.")
  }
  
  n <- nrow(L)
  if (ncol(L) != n) stop("L must be square.")
  if (k_max < 2L || k_max >= n) stop("k_max must satisfy 2 <= k_max < n.")
  if (!isTRUE(isSymmetric(L))) stop("L must be symmetric.")
  
  # RSpectra works most reliably with a general compressed sparse matrix.
  if (methods::is(L, "symmetricMatrix")) {
    L <- methods::as(L, "generalMatrix")
  }
  if (methods::is(L, "sparseMatrix") &&
      !methods::is(L, "CsparseMatrix")) {
    L <- methods::as(L, "CsparseMatrix")
  }
  
  eig <- RSpectra::eigs_sym(
    A = L,
    k = as.integer(k_max),
    which = "SM",
    opts = list(tol = 1e-8, maxitr = 10000)
  )
  
  order_index <- order(eig$values)
  eigenvalues <- eig$values[order_index]
  B <- eig$vectors[, order_index, drop = FALSE]
  
  zero_tolerance <- 1e-8 * max(1, max(abs(eigenvalues)))
  eigenvalues[abs(eigenvalues) < zero_tolerance] <- 0
  
  if (sum(eigenvalues == 0) != 1L) {
    stop("The graph must be connected and have exactly one zero eigenvalue.")
  }
  
  # Fix arbitrary eigenvector signs so repeated runs are easier to compare.
  for (j in seq_len(ncol(B))) {
    pivot <- which.max(abs(B[, j]))
    if (B[pivot, j] < 0) B[, j] <- -B[, j]
  }
  
  colnames(B) <- NULL
  attr(B, "roughness") <- eigenvalues
  attr(B, "intercept_col") <- 1L
  B
}


build_tps_basis_truncated <- function(lat, lon, k_max) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install mgcv to construct the thin-plate-spline basis.")
  }
  
  lat <- as.numeric(lat)
  lon <- as.numeric(lon)
  n <- length(lat)
  
  if (length(lon) != n) stop("lat and lon must have the same length.")
  if (k_max < 4L || k_max >= n) stop("k_max must satisfy 4 <= k_max < n.")
  
  smooth <- mgcv::smoothCon(
    mgcv::s(lat, lon, k = as.integer(k_max), bs = "tp", fx = FALSE),
    data = data.frame(lat = lat, lon = lon),
    absorb.cons = FALSE
  )[[1L]]
  
  X <- as.matrix(smooth$X)
  S <- as.matrix(smooth$S[[1L]])
  null_dimension <- as.integer(smooth$null.space.dim)
  
  if (null_dimension != 3L) {
    stop("Expected the TPS affine null space to have dimension three.")
  }
  
  # Demmler-Reinsch construction: eigenvalues measure basis roughness.
  gram_chol <- chol(crossprod(X))
  gram_chol_inverse <- backsolve(gram_chol, diag(k_max))
  S_orthogonal <- crossprod(
    gram_chol_inverse,
    S %*% gram_chol_inverse
  )
  S_orthogonal <- (S_orthogonal + t(S_orthogonal)) / 2
  
  eig <- eigen(S_orthogonal, symmetric = TRUE)
  order_index <- order(eig$values)
  roughness_all <- eig$values[order_index]
  B_dr <- X %*% gram_chol_inverse %*%
    eig$vectors[, order_index, drop = FALSE]
  
  # Use a transparent affine basis for the three unpenalized directions.
  affine <- cbind(constant = 1, latitude = lat, longitude = lon)
  Q_affine <- qr.Q(qr(affine))
  
  B_penalized <- B_dr[, seq.int(null_dimension + 1L, k_max), drop = FALSE]
  B_penalized <- B_penalized -
    Q_affine %*% crossprod(Q_affine, B_penalized)
  B_penalized <- qr.Q(qr(B_penalized))
  
  B <- cbind(Q_affine, B_penalized)
  roughness <- c(
    rep(0, null_dimension),
    pmax(roughness_all[seq.int(null_dimension + 1L, k_max)], 0)
  )
  
  for (j in seq_len(ncol(B))) {
    pivot <- which.max(abs(B[, j]))
    if (B[pivot, j] < 0) B[, j] <- -B[, j]
  }
  
  colnames(B) <- NULL
  attr(B, "roughness") <- roughness
  attr(B, "null_space_dim") <- null_dimension
  attr(B, "intercept_col") <- 1L
  attr(B, "initial_spatial_block") <- 2L:3L
  B
}


validate_ordered_basis <- function(B, roughness = attr(B, "roughness")) {
  if (!is.matrix(B) || !is.numeric(B) || any(!is.finite(B))) {
    stop("B must be a finite numeric matrix.")
  }
  if (length(roughness) != ncol(B) ||
      any(diff(roughness) < -1e-8)) {
    stop("Basis roughness must be finite and nondecreasing.")
  }
  
  gram_error <- max(abs(crossprod(B) - diag(ncol(B))))
  if (gram_error > 1e-5) {
    stop(sprintf("The basis is not orthonormal; Gram error = %.3g.", gram_error))
  }
  
  invisible(TRUE)
}
