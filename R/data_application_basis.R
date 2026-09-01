# Scalable, truncated basis construction for the data application.
#
# build_tps_basis()/eigen(K_proj) in funcs.R compute the FULL n x n TPS
# eigenbasis and are only tractable at the scale of the simulation (n = 915). 
# At the zip code scale (n ~ 33,000) the computation of the full basis is infeasible
# The functions here compute only the leading k_max columns of the
# same ordered (largest to smallest spatial scale) basis, using solvers whose
# cost scales with k_max rather than with n^3.
#
# Both return an n x k_max matrix in the same "largest to smallest spatial
# scale" column order that B expects, so
# partition_basis_columns()/project_Ac()/select_basis_outer_fold() can still operate as usual

# Truncated thin-plate-regression-spline basis via mgcv::smoothCon(bs = "tp").
# This is Wood (2003)'s low-rank approximation to the same r^2*log(r) radial
# kernel basis build_tps_basis() computes exactly, restricted up front to a
# k_max-dimensional subspace -- it never forms the full n x n kernel matrix.
#
# Returns an n x k_max matrix. Columns are NOT guaranteed to be individually
# orthonormal the way build_tps_basis()'s eigenvectors are (mgcv's TPRS basis
# is a smoothness-ordered basis, not literally an eigenbasis), so callers
# that need orthonormal columns (e.g. anything doing principal-angle
# comparisons across bases) should orthonormalize via qr.Q(qr(B)) first.
#
# Column order: mgcv::smoothCon() returns X as
#   [ penalised range space: smoothest -> most localized | null space: 1, lat, lon ]
# So move the affine terms first
build_tps_basis_truncated <- function(lat, lon, k_max) {
  # lat, lon : coordinate vectors, length n
  # k_max    : number of basis columns to keep (including the affine part)
  # returns: basis B
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("mgcv is required for build_tps_basis_truncated().", call. = FALSE)
  }
  n <- length(lat)
  if (length(lon) != n) {
    stop("lat and lon must have the same length.", call. = FALSE)
  }
  if (!is.numeric(k_max) || length(k_max) != 1L || k_max < 4L || k_max >= n) {
    stop("k_max must be a scalar with 4 <= k_max < n (mgcv requires k >= 4 for bs = 'tp').",
         call. = FALSE)
  }
  k_max <- as.integer(k_max)

  data <- data.frame(lat = lat, lon = lon)
  sm <- mgcv::smoothCon(
    mgcv::s(lat, lon, k = k_max, bs = "tp", fx = TRUE),
    data = data,
    knots = NULL,
    absorb.cons = FALSE
  )[[1]]

  B <- sm$X
  if (!identical(dim(B), c(n, k_max))) {
    stop(sprintf(
      "mgcv::smoothCon returned a %d x %d basis; expected %d x %d.",
      nrow(B), ncol(B), n, k_max
    ), call. = FALSE)
  }

  # Reorder: affine / largest-scale block first,
  # then the penalised range space from smoothest to most localized. mgcv
  # returns the M-dimensional null space as the trailing M columns.
  m_null <- sm$null.space.dim
  if (!is.numeric(m_null) || length(m_null) != 1L || is.na(m_null) ||
      m_null < 1L || m_null >= k_max) {
    stop(sprintf(
      "Unexpected mgcv null.space.dim (%s) for a %d-column tp basis.",
      format(m_null), k_max
    ), call. = FALSE)
  }
  m_null <- as.integer(m_null)
  B <- B[, c(seq.int(k_max - m_null + 1L, k_max), seq_len(k_max - m_null)),
         drop = FALSE]

  colnames(B) <- NULL
  B
}

# Truncated graph Laplacian eigenbasis via RSpectra::eigs_sym(..., which = "SM").
# Computes only the k_max smallest-eigenvalue (smoothest, largest-scale)
# eigenpairs of the sparse graph Laplacian L = D - W, using an iterative
# Lanczos solver that only needs sparse matrix-vector products
#
# Returns an n x k_max matrix with orthonormal columns (a true truncated
# eigenbasis, unlike the TPS constructor above).
build_gl_basis_truncated <- function(L, k_max) {
  # L     : sparse n x n graph Laplacian (e.g. Diagonal(rowSums(W)) - W)
  # k_max : number of eigenvectors to keep, ordered smallest eigenvalue first
  #         (largest spatial scale first, matching the TPS ordering convention)
  # returns: basis B
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("RSpectra is required for build_gl_basis_truncated().", call. = FALSE)
  }
  n <- nrow(L)
  if (ncol(L) != n) {
    stop("L must be a square matrix.", call. = FALSE)
  }
  if (!is.numeric(k_max) || length(k_max) != 1L || k_max < 1L || k_max >= n) {
    stop("k_max must be a scalar with 1 <= k_max < n.", call. = FALSE)
  }
  k_max <- as.integer(k_max)

  # RSpectra::eigs_sym() only dispatches on general sparse (dgCMatrix) or
  # dense matrices -- it errors on Matrix's symmetric-tagged dsCMatrix class,
  # which Matrix() / forceSymmetric() / D - W can silently produce for a
  # matrix that happens to be numerically symmetric. Coerce defensively so
  # this works regardless of how the caller built L.
  if (methods::is(L, "sparseMatrix") && !methods::is(L, "CsparseMatrix")) {
    L <- methods::as(L, "CsparseMatrix")
  }
  if (methods::is(L, "symmetricMatrix")) {
    L <- methods::as(L, "generalMatrix")
  }

  eig <- RSpectra::eigs_sym(
    L, k = k_max, which = "SM",
    opts = list(tol = 1e-3, maxitr = 10000)
  )
  # eigs_sym does not guarantee ascending eigenvalue order; enforce it so
  # column 1 is the smoothest/largest-scale direction, matching the TPS
  # constructor's ordering convention.
  ord <- order(eig$values)
  B <- eig$vectors[, ord, drop = FALSE]
  colnames(B) <- NULL
  B
}
