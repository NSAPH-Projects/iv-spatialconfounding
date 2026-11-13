# Superpopulation simulation with a simple causal graph.
#
# Each replicate draws n locations S_1, ..., S_n uniformly over the
# contiguous US, together with independent non-spatial noise, so the
# observations (S, X, U, A, Y) are iid draws from one population with a
# single true psi. The graph is
#
#   S -> X,  S -> U,  X -> A_uc,  (X, U) -> A_c,  A = A_c + A_uc,
#   (X, A, U) -> Y,
#
# with every variable also having its own independent noise, except U, which
# is a fixed function of location (a spatial confounder). In every mechanism
#
#   X    = 0.6 F(S) + 0.8 eta                  F medium-scale, Var(X) = 1
#   U    = G(S)                                G mechanism-specific
#   A_c  = U + 0.3 X + 0.5 nu
#   A_uc = 0.3 X + eps
#   A    = A^c + A^uc
#   Y    = mu(A, X, U) + N(0, 1)               (noiseless_y())
#
# with eta, nu, eps iid N(0, 1). No two variables are equal: A_c carries its
# own noise nu, so the oracle ({U, X}) and true IV ({A_c, X}) adjustment sets
# differ, and both are valid: given (X, U) or (X, A_c), A varies only through
# noise (nu + eps, or eps) that is independent of U. The confounded exposure
# variation U has variance 1, about 38% of Var(A), so the IV methods' tau = 0.5
# exposure R^2 reference lies beyond the basis functions capturing 95% of it
# (the selection rule's assumption).
# The surfaces are nearly uncorrelated over the US (|r| < 0.1, preprocessing.R),
# so X and U are independent apart from sharing the parent S.
#
# Mechanisms differ only in G (two confounders in twoU):
#   LS          G large-scale (1200-2500 km wavelengths; see preprocessing.R)
#   MS          G medium-scale (800-1500 km)
#   twoU        U1 = G1(S), U2 = G2(S), independent large-scale surfaces
#               (G1 = LS's G), A_c = 0.7 U1 + 0.7 U2 + 0.3 X + 0.5 nu
#   coordnoise  G = standardized{sin(pi x) + cos(pi y) + x y}, (x, y) the
#               coordinates rescaled to [0, 1] over the US bounding box
#   stateJump   G = standardized{G_LS(S) + b_state(S)}: LS's surface plus a
#               fixed offset b_r ~ N(0, 1) per state (preprocessing.R), so U
#               jumps at state borders (state policies, health systems)

# CONUS Albers equal-area, in km.
us_crs_km <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +datum=NAD83 +units=km"

n_truth_points <- 2e6

# Wavelength bands (km) of the cosine-wave surfaces (preprocessing.R).
wavelength_bands <- list(
  large = c(1200, 2500),   # G in LS and stateJump, U1 and U2 in twoU
  medium = c(800, 1500),   # G in MS
  x = c(500, 1500),        # F, the spatial part of X
  small = c(200, 500)      # used only in the seed rule for F (preprocessing.R)
)

# Which fixed surface plays G in each mechanism (names in sim$surfaces).
scenarios <- list(
  LS = list(u = "large"),
  MS = list(u = "medium_conf"),
  twoU = list(u = c("large", "large2")),
  coordnoise = list(u = "coord"),
  stateJump = list(u = "jump")
)
scenario_names <- names(scenarios)

# ---- Locations ----

# Uniform points inside the outline, by rejection from its bounding box.
sample_us_points <- function(outline, n) {
  bb <- sf::st_bbox(outline)
  out <- matrix(numeric(0), 0, 2)
  while (nrow(out) < n) {
    m <- ceiling(1.8 * (n - nrow(out))) + 100L
    cand <- cbind(runif(m, bb["xmin"], bb["xmax"]), runif(m, bb["ymin"], bb["ymax"]))
    pts <- sf::st_as_sf(as.data.frame(cand), coords = c(1, 2), crs = sf::st_crs(outline))
    keep <- lengths(sf::st_intersects(pts, outline)) > 0
    out <- rbind(out, cand[keep, , drop = FALSE])
  }
  out <- out[seq_len(n), , drop = FALSE]
  colnames(out) <- c("x", "y")
  out
}

# ---- Surfaces ----

# Sum of n_waves cosines with random directions and phases; wavelengths are
# log-uniform within band (km). Standardized over the reference sample.
make_band_surface <- function(band, reference_points, n_waves = 30L, seed) {
  set.seed(seed)
  wavelength <- exp(runif(n_waves, log(band[1]), log(band[2])))
  angle <- runif(n_waves, 0, 2 * pi)
  surface <- list(
    omega = cbind(cos(angle), sin(angle)) / wavelength,
    phase = runif(n_waves, 0, 2 * pi),
    center = 0,
    scale = 1
  )
  ref <- eval_surface(surface, reference_points)
  surface$center <- mean(ref)
  surface$scale <- sd(ref)
  surface
}

eval_surface <- function(surface, points, chunk = 2e5) {
  n <- nrow(points)
  out <- numeric(n)
  for (start in seq(1L, n, by = chunk)) {
    idx <- start:min(n, start + chunk - 1L)
    arg <- 2 * pi * points[idx, , drop = FALSE] %*% t(surface$omega)
    arg <- sweep(arg, 2L, surface$phase, `+`)
    out[idx] <- rowSums(cos(arg))
  }
  (out - surface$center) / surface$scale
}


# Unstandardized coordinate function of mechanism coordnoise; (x, y) rescaled
# to [0, 1] over the US bounding box. Center and scale from preprocessing.R.
coord_surface_raw <- function(points, bbox) {
  x <- (points[, 1] - bbox[["xmin"]]) / (bbox[["xmax"]] - bbox[["xmin"]])
  y <- (points[, 2] - bbox[["ymin"]]) / (bbox[["ymax"]] - bbox[["ymin"]])
  sin(pi * x) + cos(pi * y) + x * y
}

# ---- Data ----

coord_surface <- function(sim, points) {
  (coord_surface_raw(points, sim$coord_u$bbox) - sim$coord_u$center) / sim$coord_u$scale
}

# Index (into sim$states) of the state containing each point. Points on a
# border take the first match; the few in no polygon (slivers between the
# state polygons and the dissolved outline) take the nearest state. The
# reference points' states are cached in sim$jump$ref_state.
state_of <- function(sim, points) {
  if (nrow(points) == nrow(sim$reference_points) && identical(points, sim$reference_points)) {
    return(sim$jump$ref_state)
  }
  pts <- sf::st_as_sf(as.data.frame(points), coords = c(1, 2), crs = sf::st_crs(sim$states))
  hits <- sf::st_intersects(pts, sim$states)
  idx <- vapply(hits, function(h) if (length(h)) h[1] else NA_integer_, integer(1))
  if (anyNA(idx)) idx[is.na(idx)] <- sf::st_nearest_feature(pts[is.na(idx), ], sim$states)
  idx
}

# Mechanism stateJump: LS's surface plus a fixed offset per state.
jump_surface <- function(sim, points) {
  raw <- eval_surface(sim$surfaces$large, points) + sim$jump$offsets[state_of(sim, points)]
  (raw - sim$jump$center) / sim$jump$scale
}

eval_named_surface <- function(sim, name, points) {
  switch(name,
         coord = coord_surface(sim, points),
         jump = jump_surface(sim, points),
         eval_surface(sim$surfaces[[name]], points))
}

draw_scenario_data <- function(sim, scenario, points) {
  sc <- scenarios[[scenario]]
  n <- nrow(points)
  d <- list(s = points)
  d$X <- cbind(X = 0.6 * eval_surface(sim$surfaces$x, points) + 0.8 * rnorm(n))
  if (length(sc$u) == 2L) {
    d$U1 <- eval_named_surface(sim, sc$u[1], points)
    d$U2 <- eval_named_surface(sim, sc$u[2], points)
    u_part <- 0.7 * d$U1 + 0.7 * d$U2
  } else {
    d$U <- eval_named_surface(sim, sc$u, points)
    u_part <- d$U
  }
  d$Ac <- u_part + 0.3 * d$X[, 1] + 0.5 * rnorm(n)
  d$Auc <- 0.3 * d$X[, 1] + rnorm(n)
  d$A <- d$Ac + d$Auc
  d
}

# E(Y | A = a, X, U). Intercept 8 keeps E(Y) > 0, so psi reads as a relative
# outcome rate; tanh is nonlinear in A but has no curvature at the cutoff 0.
# twoU: -2 U1 - U2 and interaction -(U1 + 0.5 U2) a (halved in tanh).
noiseless_y <- function(d, a, option) {
  if (is.null(d$U1)) {
    base <- 8 - 2 * d$U
    slope_u <- -d$U
  } else {
    base <- 8 - 2 * d$U1 - d$U2
    slope_u <- -(d$U1 + 0.5 * d$U2)
  }
  base <- base + 0.5 * d$X[, 1]
  switch(
    option,
    linear = base + a + slope_u * a,
    tanh = base + 2 * tanh(a / 2) + 0.5 * slope_u * a,
    stop("Unknown outcome option: ", option)
  )
}

# ---- Truth ----

# psi = E{Y(min(A, c))} / E(Y) (Corollary 2):
#   psi = (theta1 (1 - theta2) + theta2 theta3) / theta4,
#   theta1 = E[ E(Y | A = c, X, U) | A > c ],
#   theta2 = P(A <= c), theta3 = E(Y | A <= c), theta4 = E(Y),
# by Monte Carlo over the 2e6 reference points and fresh noise.
true_psi <- function(sim, scenario, option, cutoff, seed = 999L) {
  set.seed(seed)
  d <- draw_scenario_data(sim, scenario, sim$reference_points)
  m <- noiseless_y(d, d$A, option)
  m_at_cutoff <- noiseless_y(d, cutoff, option)
  below <- d$A <= cutoff
  theta1 <- mean(m_at_cutoff[!below])
  theta2 <- mean(below)
  theta3 <- mean(m[below])
  theta4 <- mean(m)
  (theta1 * (1 - theta2) + theta2 * theta3) / theta4
}

# ---- Spatial bases (rebuilt for every dataset's sampled locations) ----

# Full-rank thin plate spline basis, ordered from smoothest to roughest:
# the affine null space (constant, x, y), then the Demmler-Reinsch basis of
# the TPS kernel E_ij = r_ij^2 log r_ij projected off the affine space. For
# M = Z' E Z = V diag(lambda) V', the column Z v_j has roughness 1 / lambda_j,
# so columns are ordered by decreasing lambda. Equivalent to
# build_tps_basis_truncated(k_max = n - 1) in ../R/basis_construction.R but
# one eigendecomposition instead of mgcv's construction (~5 s vs ~13 s at
# n = 1000); the leading columns span the same subspaces. Drops the roughest
# column, giving m = n - 1 columns. coords are in thousands of km.
build_tps_basis_full <- function(coords) {
  n <- nrow(coords)
  r <- as.matrix(stats::dist(coords))
  E <- ifelse(r > 0, r^2 * log(r), 0)
  Q <- qr.Q(qr(cbind(1, coords)), complete = TRUE)
  affine <- Q[, 1:3]
  Z <- Q[, -(1:3)]
  eig <- eigen(crossprod(Z, E %*% Z), symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  B <- cbind(affine, Z %*% eig$vectors[, ord])[, seq_len(n - 1L)]
  for (j in seq_len(ncol(B))) {
    pivot <- which.max(abs(B[, j]))
    if (B[pivot, j] < 0) B[, j] <- -B[, j]
  }
  attr(B, "roughness") <- c(0, 0, 0, 1 / eig$values[ord])[seq_len(n - 1L)]
  B
}

# Graph Laplacian basis of the symmetric k-nearest-neighbour graph, with the
# smallest k (up to k_max) that makes the graph connected, as in
# ../data_application/0_preprocessing.R. Ordered by eigenvalue; drops the
# roughest column (m = n - 1). Base R only (no RANN / igraph), so cluster
# jobs do not depend on extra packages.
build_gl_basis_full <- function(coords, k_max = 30L) {
  n <- nrow(coords)
  dist_matrix <- as.matrix(stats::dist(coords))
  diag(dist_matrix) <- Inf
  nn <- t(apply(dist_matrix, 1L, function(row) order(row)[seq_len(k_max)]))
  is_connected <- function(W) {
    seen <- rep(FALSE, n)
    frontier <- 1L
    seen[1L] <- TRUE
    while (length(frontier) > 0L) {
      nxt <- which(colSums(W[frontier, , drop = FALSE]) > 0 & !seen)
      seen[nxt] <- TRUE
      frontier <- nxt
    }
    all(seen)
  }
  for (k in seq_len(k_max)) {
    W <- matrix(0, n, n)
    W[cbind(rep(seq_len(n), each = k), as.vector(t(nn[, seq_len(k), drop = FALSE])))] <- 1
    W <- pmax(W, t(W))
    if (is_connected(W)) break
  }
  if (!is_connected(W)) stop("kNN graph not connected for k <= ", k_max)
  eig <- eigen(diag(rowSums(W)) - W, symmetric = TRUE)
  ord <- order(eig$values)
  B <- eig$vectors[, ord][, seq_len(n - 1L)]
  for (j in seq_len(ncol(B))) {
    pivot <- which.max(abs(B[, j]))
    if (B[pivot, j] < 0) B[, j] <- -B[, j]
  }
  attr(B, "roughness") <- pmax(eig$values[ord], 0)[seq_len(n - 1L)]
  attr(B, "knn_k") <- k
  B
}

# ---- Fits ----

# Each method's own covariates plus X (baseline = X only).
method_covariates <- function(d, method) {
  own <- switch(
    method,
    oracle = if (is.null(d$U1)) cbind(U = d$U) else cbind(U1 = d$U1, U2 = d$U2),
    baseline = cbind(Intercept = rep(1, length(d$A))),
    spatialcoord = {
      # Coordinates scaled to [0, 1].
      x <- d$s[, 1]
      y <- d$s[, 2]
      cbind(x = (x - min(x)) / diff(range(x)), y = (y - min(y)) / diff(range(y)))
    },
    trueIV = cbind(Ac = d$Ac),
    stop("Unknown method: ", method)
  )
  if (method == "baseline") d$X else cbind(own, d$X)
}

# Super Learner library includes SL.earth (MARS with pairwise interactions),
# which can represent interactions such as A x (A^c)^2 that the other
# candidates cannot.
effect_settings <- list(
  nuisance_args = list(
    sl_library = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction", "SL.earth"),
    density_trim = 0.01
  ),
  n_grid = 25L,
  density_trim = 0.01,
  constrain = TRUE
)

fit_method <- function(d, y, method, cutoff) {
  # Same Super Learner CV splits whatever ran before (as in select_basis_dimension()).
  set.seed(12001L)
  # Errors propagate to run_sim.R, which records the message.
  fit <- do.call(
    estimate_truncated_effect_full,
    c(list(y = y, a = d$A, w = method_covariates(d, method), cutoff = cutoff), effect_settings)
  )
  c(psi = fit$psi, ci_lower = fit$confidence_interval[1L], ci_upper = fit$confidence_interval[2L])
}

# IV-TPS / IV-GL: sequential stability selection of n_ac on the dataset's
# basis (../R/basis_selection.R), then the estimate at the selected n_ac.
# Also returns how well the selected split recovers the confounding:
# cor_ac = corr(selected A^c, U's part of A^c), cor_auc = corr(selected
# A^uc, true A^uc). The basis cannot capture X or the noise, so A^c's spatial
# confounded part (U, or 0.7 U1 + 0.7 U2) is the target of cor_ac. NA when one
# side is constant (e.g. n_ac = 1, where A^c is the mean of A).
fit_iv <- function(d, y, B, cutoff, tau = 0.5) {
  x <- d$X
  sel <- select_basis_dimension(
    y = y, a = d$A, x = x, B = B, cutoff = cutoff,
    effect_estimator = estimate_truncated_effect_full,
    effect_args = effect_settings,
    r2_reference = tau, family_alpha = 0.05, seed = 12001L
  )
  c(psi = sel$selected_fit$psi,
    ci_lower = sel$selected_fit$confidence_interval[1L],
    ci_upper = sel$selected_fit$confidence_interval[2L],
    n_ac = sel$selected_n_ac,
    n_ac0 = sel$settings$n_ac0,
    n_fits = nrow(sel$diagnostics),
    cor_ac = safe_cor(sel$selected_projection$Ac, confounded_part(d)),
    cor_auc = safe_cor(sel$selected_projection$Auc, d$Auc))
}

# Spatial confounded part of A^c.
confounded_part <- function(d) if (is.null(d$U1)) d$U else 0.7 * d$U1 + 0.7 * d$U2

safe_cor <- function(u, v) {
  if (stats::sd(u) < 1e-12 || stats::sd(v) < 1e-12) NA_real_ else stats::cor(u, v)
}
