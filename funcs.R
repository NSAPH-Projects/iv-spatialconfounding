# Function used to estimate the exposure-response curve
ctseff <- function(y, a, x, bw.seq, n.pts = 100, a.rng = c(min(a), max(a)),
                   sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"),
                   constrain = T, trim = 0.01,
                   folds = NULL) {
  # y is outcome
  # a is exposure
  # x is covariate matrix
  # bw.seq is a sequence of bandwidth values
  # a.rng is the range of exposure values to evaluate the ERF
  # n.pts is the number of points within a.rng at which to evaluate the ERF
  # sl.lib is the library of SuperLearner algorithms to use
  # constrain is a boolean indicating whether pseudo-outcome is restricted to (min(Y), max(Y))
  # trim is the quantile at which to trim the density estimate to prevent extreme IPW weights
  # returns a list of two dataframes and a list
  
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require(KernSmooth)
  kern <- function(t) {
    dnorm(t)
  }
  
  n <- nrow(x)

  # set up evaluation points & matrices for predictions
  a.min <- a.rng[1]
  a.max <- a.rng[2]
  a.vals <- seq(a.min, a.max, length.out = n.pts)
  
  xa.new <- rbind(cbind(x, a), cbind(x[rep(1:n, length(a.vals)), ], 
                                     a = rep(a.vals, rep(n, length(a.vals)))))
  colnames(xa.new) <- c(colnames(x), "a") # sophie's change.
  x.new <- xa.new[, -dim(xa.new)[2]]
  x <- data.frame(x)
  x.new <- data.frame(x.new)
  colnames(x.new) <- colnames(x) # sophie's change
  xa.new <- data.frame(xa.new)
  #print('created evaluation points and matrices for predictions')
  
  # estimate nuisance functions via super learner
  # note: other methods could be used here instead
  if (is.null(folds)) {
    # In-sample estimation (used inside select_t and when folds not provided).
    pimod      <- SuperLearner(Y = a, X = data.frame(x), SL.library = sl.lib, newX = x.new)
    pimod.vals <- pimod$SL.predict
    pi2mod     <- SuperLearner(Y = log((a - pimod.vals[1:n])^2), X = x,
                               SL.library = sl.lib, newX = x.new)
    pi2mod.vals <- exp(pi2mod$SL.predict)
    mumod      <- SuperLearner(Y = y, X = cbind(x, a), SL.library = sl.lib, newX = xa.new)
    muhat.vals <- mumod$SL.predict

    a.std      <- (xa.new$a - pimod.vals) / sqrt(pi2mod.vals)
    pihat.vals <- approx(density(a.std[1:n], from = min(a.std), to = max(a.std))$x,
                         density(a.std[1:n], from = min(a.std), to = max(a.std))$y,
                         xout = a.std)$y / sqrt(pi2mod.vals)
    pihat      <- pihat.vals[1:n]
    pihat      <- pmax(pihat, quantile(pihat, trim))
    pihat.mat  <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
    muhat      <- muhat.vals[1:n]
    muhat.mat  <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))

  } else {
    # K-fold cross-fitting of nuisance functions (supplement Section 3, step 3d-e).
    # For each fold k, pimod/mumod are trained on the other K-1 folds and
    # predicted on the held-out fold, giving out-of-fold nuisance estimates.
    K         <- max(folds)
    pihat     <- numeric(n)
    pihat.mat <- matrix(NA_real_, n, length(a.vals))
    muhat     <- numeric(n)
    muhat.mat <- matrix(NA_real_, n, length(a.vals))

    for (k in seq_len(K)) {
      tr   <- which(folds != k)
      ho   <- which(folds == k)
      n_ho <- length(ho)

      x_tr  <- data.frame(x[tr, , drop = FALSE])
      xa_tr <- data.frame(cbind(x[tr, , drop = FALSE], a = a[tr]))
      colnames(xa_tr) <- c(colnames(x), "a")

      # newX: all n obs (needed for density support from training residuals)
      #       + held-out fold x a.vals grid (for pihat.mat, muhat.mat)
      x_ho_grid  <- data.frame(x[ho[rep(seq_len(n_ho), length(a.vals))], , drop = FALSE])
      colnames(x_ho_grid) <- colnames(x)
      x_new_k    <- rbind(data.frame(x), x_ho_grid)

      xa_ho_grid <- data.frame(cbind(
        x[ho[rep(seq_len(n_ho), length(a.vals))], , drop = FALSE],
        a = rep(a.vals, each = n_ho)
      ))
      colnames(xa_ho_grid) <- c(colnames(x), "a")
      xa_all   <- data.frame(cbind(x, a = a))
      colnames(xa_all) <- c(colnames(x), "a")
      xa_new_k <- rbind(xa_all, xa_ho_grid)

      pimod_k  <- SuperLearner(Y = a[tr], X = x_tr, SL.library = sl.lib, newX = x_new_k)
      pi_k     <- pimod_k$SL.predict          # length n + n_ho*n.pts
      pi2mod_k <- SuperLearner(Y = log((a[tr] - pi_k[tr])^2),
                               X = x_tr, SL.library = sl.lib, newX = x_new_k)
      pi2_k    <- exp(pi2mod_k$SL.predict)    # length n + n_ho*n.pts
      mumod_k  <- SuperLearner(Y = y[tr], X = xa_tr, SL.library = sl.lib, newX = xa_new_k)
      mu_k     <- mumod_k$SL.predict          # length n + n_ho*n.pts

      # Density of standardised residuals, estimated from training fold
      a_std_k  <- (a - pi_k[1:n]) / sqrt(pi2_k[1:n])
      dens_k   <- density(a_std_k[tr], from = min(a_std_k), to = max(a_std_k))

      # pihat at held-out observed values (out-of-fold)
      pihat[ho] <- approx(dens_k$x, dens_k$y, xout = a_std_k[ho])$y / sqrt(pi2_k[ho])

      # pihat.mat and muhat.mat at held-out fold x a.vals grid (out-of-fold)
      grid_idx   <- (n + 1):(n + n_ho * length(a.vals))
      a_std_grid <- (rep(a.vals, each = n_ho) - pi_k[grid_idx]) / sqrt(pi2_k[grid_idx])
      pihat.mat[ho, ] <- matrix(
        approx(dens_k$x, dens_k$y, xout = a_std_grid)$y / sqrt(pi2_k[grid_idx]),
        nrow = n_ho, ncol = length(a.vals)
      )
      muhat.mat[ho, ] <- matrix(mu_k[grid_idx], nrow = n_ho, ncol = length(a.vals))
      muhat[ho]       <- mu_k[ho]
    }

    pihat <- pmax(pihat, quantile(pihat, trim))
  }

  # construct varpi/m from pihat.mat and muhat.mat (common to both paths)
  varpihat     <- predict(smooth.spline(a.vals, apply(pihat.mat, 2, mean)), x = a)$y
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = T, nrow = n)
  mhat         <- predict(smooth.spline(a.vals, apply(muhat.mat, 2, mean)), x = a)$y
  mhat.mat     <- matrix(rep(apply(muhat.mat, 2, mean), n), byrow = T, nrow = n)
  
  
  # form adjusted/pseudo outcome xi
  pseudo.out <- (y - muhat) / (pihat / varpihat) + mhat
  # Keep unconstrained copy for influence function computation: clipping pseudo.out
  # before computing phis shrinks the residuals and biases the variance downward.
  pseudo.out.if <- pseudo.out
  if (constrain){
    pseudo.out[pseudo.out > max(y)] <- max(y)
    pseudo.out[pseudo.out < min(y)] <- min(y)
  }

  #print('calculated pseudo.out')
  
  # leave-one-out cross-validation to select bandwidth
  w.fn <- function(bw, a.vals) { # sophie's change
    w.avals <- NULL
    for (a.val in a.vals) {
      a.std <- (a - a.val) / bw
      kern.std <- kern(a.std) / bw
      denom <- mean(kern.std) * mean(a.std^2 * kern.std) - mean(a.std * kern.std)^2
      w.avals <- c(w.avals, if (abs(denom) < 1e-10) NA_real_
                            else mean(a.std^2 * kern.std) * (kern(0) / bw) / denom)
    }
    return(w.avals / n)
  }
  
  hatvals <- function(bw) {
    asubset = seq(min(a), max(a), length.out = 100)
    approx(asubset, w.fn(bw, a.vals = asubset), xout = a)$y # sophie's change
  }
  cts.eff.fn <- function(out, bw) {
    approx(locpoly(a, out, bandwidth = bw), xout = a)$y 
  }
  # note: choice of bandwidth range depends on specific problem,
  # make sure to inspect plot of risk as function of bandwidth
  risk.fn <- function(h) {
    hats <- hatvals(h)
    mean(((pseudo.out - cts.eff.fn(pseudo.out, bw = h)) / (1 - hats))^2)
  } 
  risk.est <- sapply(bw.seq, risk.fn)
  h.opt <- bw.seq[which.min(risk.est)]
  bw.risk <- data.frame(bw = bw.seq, risk = risk.est)
  #print('calculated h.opt')
  
  # alternative approach:
  # h.opt <- optimize(function(h){ hats <- hatvals(h); mean( ((pseudo.out[a > a.min & a < a.max]-cts.eff.fn(pseudo.out,bw=h))/(1-hats))^2) } ,
  #  bw.seq, tol=0.01)$minimum
  
  # estimate effect curve with optimal bandwidth
  est <- approx(locpoly(a, pseudo.out, bandwidth = h.opt), xout = a.vals)$y
  #print('calculated est')
  
  phis <- list()
  ix <- 1
  for (a.val in a.vals) {
    a.std <- (a - a.val) / h.opt
    kern.std <- kern(a.std) / h.opt
    beta <- coef(lm(pseudo.out ~ a.std, weights = kern.std))
    Dh <- matrix(c(
      mean(kern.std), mean(kern.std * a.std),
      mean(kern.std * a.std), mean(kern.std * a.std^2)
    ), nrow = 2)
    kern.mat <- matrix(rep(kern((a.vals - a.val) / h.opt) / h.opt, n), byrow = T, nrow = n)
    g2 <- matrix(rep((a.vals - a.val) / h.opt, n), byrow = T, nrow = n)
    intfn1.mat <- kern.mat * (muhat.mat - mhat.mat) * varpihat.mat
    intfn2.mat <- g2 * kern.mat * (muhat.mat - mhat.mat) * varpihat.mat
    int1 <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]),n),
                         byrow=T,nrow=n)*(intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)]) / 2, 1,sum)
    int2 <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]),n),
                         byrow=T,nrow=n)* ( intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)]) /2, 1,sum)
    # Use unconstrained pseudo.out for influence functions to avoid downward bias in variance
    if (rcond(Dh) < 1e-10) {
      phis[[ix]] <- rep(NA_real_, n)
    } else {
      phi_both <- t(solve(Dh) %*%
                      rbind(
                        kern.std * (pseudo.out.if - beta[1] - beta[2] * a.std) + int1,
                        a.std * kern.std * (pseudo.out.if - beta[1] - beta[2] * a.std) + int2
                      ))
      phis[[ix]] <- phi_both[,1]
    }
    ix <- ix + 1
  }
  
  res <- data.frame(a.vals, est)
  
  return(invisible(list(res = res, bw.risk = bw.risk, phi = phis)))
}

# Function to create outcome
createY <- function(Us, As, option = c('linear', 'nonlinear')){
  # Us is a n x nreps matrix of simulated unmeasured confounder
  # As is a n x nreps matrix of simulated exposure
  # option is a string indicating the form of the outcome model
  # returns a vector of outcome
  option <- match.arg(option)
  n <- nrow(Us)
  nreps <- ncol(Us)
  Ys <- matrix(NA, n, nreps)
  # linear outcome model
  if (option == 'linear'){
    for (i in 1:nreps){
      Ys[,i] <- rnorm(n, -0.5 + (-1)*Us[,i] + As[,i] - 0.5*As[,i]*Us[,i], 1) 
    }
  }
  # nonlinear outcome model
  if (option == 'nonlinear'){
    for (i in 1:nreps){
      eta <- -0.5 - 0.5*Us[,i] +
        tanh(1.5*As[,i]) - 0.2*Us[,i]*tanh(As[,i]) + 
        0.1*tanh(As[,i])^2
      Ys[,i] <- rnorm(n, eta, 1)
    }
  } 
  return(Ys)
}

# Function to plot the variables titled "names" in the dataframe "df"
plotfunc <- function(df, names, labels=names,
                    xlimits = c(-125, -65), ylimits = c(25, 50)){
  # df is a sf dataframe including columns names
  # names is a vector of strings with the names of the columns to be plotted
  # labels is a vector of strings to title the ggplots
  # returns a list of ggplots

  K <- length(names)
  gs <- list()
  for (k in 1:K){
    # extract the column with names[k] from df
    var <- df[[names[k]]]
    qs <- round(quantile(var, probs = c(0.1, 0.3, 0.5, 0.7, 0.9)),2)
    # Turn qs into a vector of strings
    qschar <- as.character(qs)
    gs[[k]] <- ggplot(df) +
      xlim(xlimits[1],xlimits[2]) + 
      ylim(ylimits[1], ylimits[2]) +
      geom_sf(aes_string(fill = names[k]), color=NA, size = 0.005) +
      scale_fill_gradient2(low = "#1e90ff", 
                           mid = "white", 
                           high = "#8b0000", 
                           midpoint = 0,
                           breaks = qs,
                           labels = qschar,
                           limits = c(min(var), max(var)),
                           na.value = "white") +
      theme_minimal() +
      theme(plot.title = element_text(size = 24 * 2,hjust = 0.5),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal", 
            legend.text.align = 0.75,
            legend.key.width = unit(100, "points"),
            panel.grid.major = element_line(colour = "transparent"),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 25)
            ) + 
      ggtitle(labels[k])
  }
  return(gs)
}

# Function to calculate average absolute bias, avg RMSE, avg coverage
metrics <- function(a.vals, muests, mutrue){
  # a.vals is a vector of exposure values for which we want to estimate ERF
  # muests is a matrix of estimated ERFs, columns correspond to diff sims
  # mutrue is a vector of true ERF
  # returns a list with avgabsbias, avgRMSE, avgse
  
  avgabsbias <- mean(abs(rowMeans(muests - mutrue, na.rm = T)), na.rm = T)
  avgRMSE <- mean(sqrt(rowMeans((muests - mutrue)^2, na.rm = T)), na.rm = T)
  avgse <- mean(apply(muests,1,sd,na.rm = T), na.rm = T)
  return(list(avgabsbias = avgabsbias, 
              avgRMSE = avgRMSE,
              avgse = avgse))
}

# Build an orthonormal TPS basis ordered from largest to smallest spatial scale.
# Columns 1-3 are the orthonormal affine basis {1, lat, lon} (largest scale).
# Columns 4..n are eigenvectors of the projected r^2*log(r) kernel matrix,
# ordered from largest eigenvalue (smoothest) to smallest eigenvalue (most
# localized). Larger eigenvalue of the TPS kernel <=> smoother pattern <=>
# larger spatial scale. Returns an n x n matrix.
build_tps_basis <- function(lat, lon) {
  n     <- length(lat)
  T_mat <- cbind(1, lat, lon)
  Q_T   <- qr.Q(qr(T_mat))                 # n x 3 orthonormal affine basis
  P_perp <- diag(n) - tcrossprod(Q_T)      # projection off affine subspace

  # Radial kernel: K_ij = r_ij^2 * log(r_ij), Euclidean distance on (lat, lon)
  r_mat  <- as.matrix(dist(cbind(lat, lon)))
  K      <- matrix(0.0, n, n)
  nz     <- r_mat > 0
  K[nz]  <- r_mat[nz]^2 * log(r_mat[nz])

  # Project K onto the null-space complement; result has rank n-3.
  K_proj <- P_perp %*% K %*% P_perp

  # eigen() returns eigenvalues in decreasing order.
  # Larger eigenvalue <=> smoother <=> larger spatial scale.
  # Skip the 3 near-zero eigenvalues from the affine null space.
  E      <- eigen(K_proj, symmetric = TRUE)
  B_rad  <- E$vectors[, seq_len(n - 3L)]   # n x (n-3): largest eigenvalue first

  cbind(Q_T, B_rad)   # n x n: col 1-3 affine, col 4..n radial large->small scale
}

# Build an orthonormal GL basis ordered from largest to smallest spatial scale.
# Column 1 is the eigenvector with the smallest eigenvalue of L (constant,
# largest scale); column n is the eigenvector with the largest eigenvalue
# (most localized, smallest scale). Returns an n x n matrix.
build_gl_basis <- function(W) {
  n  <- nrow(W)
  D  <- diag(rowSums(as.matrix(W)))
  L  <- D - as.matrix(W)
  E  <- eigen(L, symmetric = TRUE)   # decreasing eigenvalues: largest lambda first
  # Largest eigenvalue of L = most localized = smallest scale -> currently in col 1.
  # We want largest scale first, so reverse the column order.
  E$vectors[, n:1]                   # n x n: col 1 = smoothest, col n = most localized
}

# Partition n observations into K spatial blocks via k-means on approximately
# equal-distance coordinates. Longitude is scaled by cos(mean_lat) so that
# one unit in each direction corresponds to roughly the same arc length.
# Returns an integer vector of fold assignments (values 1..K).
make_spatial_folds <- function(lat, lon, K = 5L) {
  cos_lat <- cos(mean(lat) * pi / 180)
  coords  <- cbind(lat, lon * cos_lat)
  km <- kmeans(coords, centers = K, nstart = 25L)
  print("computed spatial folds")
  km$cluster
}

# Sequential equivalence-based instrument selection (Supplement Algorithm 1).
#
# B_full  : n x m orthonormal basis, col 1 = largest scale, col m = smallest scale.
#           B^c(n_uc) = B_full[, 1:(m - n_uc)]   (large-scale, confounded)
#           B^uc(n_uc) = B_full[, (m - n_uc + 1):m] (small-scale, instruments)
# n_uc_star  : core number of instruments (small; exogeneity most plausible).
# n_uc_cands : increasing integer sequence of candidate instrument counts,
#              starting at n_uc_star.
# alpha      : tolerance multiplier; delta = alpha * se(psi_hat(n_uc_star)).
# folds      : integer vector (length n) of spatial fold assignments (1..K).
# cutoff     : truncation cutoff c.
# x_extra    : optional matrix of additional covariates (n x p); NULL if none.
# delta_bw   : half-bandwidth for the local ERF window around cutoff.
#
# Returns an integer vector of length K: fold-specific selected instrument counts.
select_t <- function(y, a, B_full, n_uc_star, n_uc_cands, alpha = 1.0,
                     folds, cutoff, x_extra = NULL, delta_bw = 0.05) {
  K       <- max(folds)
  m       <- ncol(B_full)
  n_cand  <- length(n_uc_cands)
  psi_by_fold <- matrix(NA_real_, nrow = K, ncol = n_cand)

  for (k in seq_len(K)) {
    tr   <- which(folds != k)
    y_tr <- y[tr];  a_tr <- a[tr]
    B_tr <- B_full[tr, , drop = FALSE]
    x_tr <- if (!is.null(x_extra)) x_extra[tr, , drop = FALSE] else NULL

    # Precompute projection scores once per fold: avoids recomputing t(B_c) %*% a_tr
    # for every candidate. Each candidate just slices the first (m - n_uc) elements.
    scores <- drop(t(B_tr) %*% a_tr)

    # Initialise Ac for the first (smallest) n_uc candidate.
    n_cols_c <- m - n_uc_cands[1L]
    Ac <- drop(B_tr[, seq_len(n_cols_c), drop = FALSE] %*% scores[seq_len(n_cols_c)])

    for (j in seq_len(n_cand)) {
      # Incremental update: moving 'by' columns from confounded → unconfounded
      # costs O(n × step) instead of O(n × (m - n_uc)) per candidate.
      if (j > 1L) {
        drop_cols <- seq(m - n_uc_cands[j] + 1L, m - n_uc_cands[j - 1L])
        Ac <- Ac - drop(B_tr[, drop_cols, drop = FALSE] %*% scores[drop_cols])
      }

      xmat <- if (!is.null(x_tr)) cbind(Ac, x_tr) else matrix(Ac, ncol = 1L)
      colnames(xmat)[1L] <- 'Ac'

      psi_by_fold[k, j] <- tryCatch({
        sub    <- a_tr > cutoff - delta_bw
        a_sub  <- a_tr[sub]
        xsub   <- matrix(xmat[sub, , drop = FALSE], ncol = ncol(xmat))
        colnames(xsub) <- colnames(xmat)
        erf  <- ctseff(
          y      = y_tr[sub],
          a      = a_sub,
          x      = xsub,
          n.pts  = 5L,
          a.rng  = c(cutoff - delta_bw, cutoff + delta_bw),
          bw.seq = seq(sd(a_sub) / 10, sd(a_sub), length.out = 50L)
        )
        ix <- which.min(abs(erf$res$a.vals - cutoff))
        (erf$res$est[ix] * mean(a_tr > cutoff) +
           mean(y_tr[a_tr <= cutoff]) * mean(a_tr <= cutoff)) / mean(y_tr)
      }, error = function(e) NA_real_)
    }
  }

  # For each fold k: select t_k using psi_by_fold[k, ] (estimated on training folds != k),
  # with SE calibrated from the other K-1 folds' core estimates.
  t_by_fold <- integer(K)
  for (k in seq_len(K)) {
    other_core <- psi_by_fold[-k, 1L]
    n_other    <- sum(!is.na(other_core))
    se_k       <- if (n_other > 1L)
                    sd(other_core, na.rm = TRUE) / sqrt(n_other)
                  else
                    sd(psi_by_fold[, 1L], na.rm = TRUE) / sqrt(K)
    delta_k <- alpha * se_k

    t_k <- n_uc_cands[1L]
    for (j in seq_len(n_cand - 1L)) {
      d <- psi_by_fold[k, j + 1L] - psi_by_fold[k, j]
      if (!is.na(d) && abs(d) < delta_k) {
        t_k <- n_uc_cands[j + 1L]
      } else {
        break
      }
    }
    t_by_fold[k] <- t_k
    print(paste0("fold ", k, ": selected n_uc = ", t_k))
  }
  t_by_fold
}

# Function that simulates data, estimates truncated exposure effect using different methods, and saves results to csvs
simfunc <- function(nsims,
                   lat,
                   lon,
                   confounding_mechanism,
                   option = c('linear', 'nonlinear'),
                   methods = c(
                     'baseline',
                     'oracle',
                     'spatialcoord',
                     'IV-TPS',
                     'IV-GraphLaplacian',
                     'IV-TPS-spatialcoord',
                     'IV-GraphLaplacian-spatialcoord',
                     'trueIV',
                     'trueIV-spatialcoord'
                   ),
                   B_tps_full,
                   B_gl_full,
                   statemat,
                   W = NULL,
                   cutoff = 0.5,
                   select_basis = TRUE,
                   n_cores = 1L,
                   results_dir = "results_Mar27/")
{
  # nsims is the number of simulations
  # lat is a vector of latitudes
  # lon is a vector of longitudes
  # option is the form of the outcome model
  # methods are the methods used to estimate truncated exposure effect
  # statemat is the matrix of state-level indicators
  # cutoff is c
  # select_basis: if TRUE, use select_t() to choose n_uc; otherwise use floor(0.07*n)
  # results_dir: directory to write output CSVs

  confounding_mechanism <- as.integer(confounding_mechanism)
  option <- match.arg(option)
  
  ################# GENERATE DATA #################
  
  # Compute distance matrix
  distmat <- geosphere::distm(cbind(lon, lat), 
                              fun = distHaversine)
  distmat <- distmat/1000000 # scale so range (0,2)
  n <- length(lat)
  # Mechanism 1: baseline GP
  if (confounding_mechanism == 1){
    dat <- compute_data_GP(nsims = nsims, distmat = distmat)
  }
  # Mechanism 3: within-state GP
  if (confounding_mechanism == 3){
    dat <- compute_data_GP_state(nsims = nsims, distmat = distmat, statemat = statemat)
  }
  # Mechanism 6: coordinate-based nonlinear
  if (confounding_mechanism == 6){
    dat <- compute_data_spatialcoord(lat = lat, long = lon, nsims = nsims, distmat = distmat)
  }
  # Mechanism 4: two-confounder GP
  if (confounding_mechanism == 4){
    dat <- compute_data_GP_2U(nsims = nsims, distmat = distmat)
    Ac  <- dat$Ac
    Auc <- dat$Auc
    U1  <- dat$U1
    U2  <- dat$U2
    A   <- Ac + Auc
    Y   <- matrix(NA, n, nsims)
    if (option == 'linear'){
      for (i in 1:nsims){
        Y[, i] <- rnorm(n, -0.5 + (-1)*U1[,i] + A[,i] - 0.5*A[,i]*U1[,i] - 0.75*A[,i]*U2[,i], 1)
      }
    }
    if (option == 'nonlinear'){
      for (i in 1:nsims){
        eta <- -0.5 - 0.5*U1[,i] +
          tanh(1.5*A[,i]) - 0.2*U2[,i]*tanh(A[,i]) +
          0.1*tanh(A[,i])^2
        Y[, i] <- rnorm(n, eta, 1)
      }
    }
  }
  # Mechanism 5: reversed-scale GP (B_uc is large-scale, B_c is small-scale)
  if (confounding_mechanism == 5){
    dat <- compute_data_GP(nsims = nsims, distmat = distmat,
                           theta_B1uc = 0.80, theta_B2uc = 0.50,
                           theta_B1c  = 0.10, theta_B2c  = 0.05)
  }

  if (confounding_mechanism == 2){
    stopifnot(!is.null(W))
    dat <- compute_data_leroux(W = W, nsims = nsims)
  }

  if (confounding_mechanism != 4){
    Ac <- dat$Ac
    Auc <- dat$Auc
    U <- dat$U
    A <- Ac + Auc # all have dimension n x nsims
    Y <- createY(Us=U, As=A, option = option)
  }
  
  
  ################# FIT MODELS #################

  # Pre-compute spatial folds once; used by select_t() for basis selection
  # and by ctseff() for cross-fitting nuisance functions.
  folds <- make_spatial_folds(lat, lon)

  t_simfunc_start <- proc.time()["elapsed"]

  for (method in methods){
    # Create filename for csvs containing estimates
    filename <- paste0(results_dir, 'conf', confounding_mechanism, '_', option, '_', method, '.csv')
  
    t_method_start <- proc.time()["elapsed"]

    # Run sims in parallel (n_cores = 1 → sequential, identical behaviour).
    results_list <- parallel::mclapply(seq_len(nsims), function(sim) {
      t_sim_start <- proc.time()["elapsed"]
      message(paste(method, sim))

      # ---- build xmat ----
      n_uc_sel <- NA_integer_
      if (method == 'baseline'){
        xmat <- matrix(rep(1, n), ncol = 1)
        colnames(xmat) <- 'Intercept'
      }
      if (method == 'oracle'){
        if (confounding_mechanism != 4){
          xmat <- matrix(U[, sim], ncol = 1)
          colnames(xmat) <- 'U'
        } else {
          xmat <- cbind(U1[, sim], U2[, sim])
          colnames(xmat) <- c('U1', 'U2')
        }
      }
      if (method == 'spatialcoord'){
        xmat <- cbind(lat, lon)
        colnames(xmat) <- c('Latitude', 'Longitude')
      }
      if (method == 'IV-TPS'){
        m      <- ncol(B_tps_full)
        Ac_vec <- numeric(n)
        if (select_basis) {
          t_by_fold_k <- select_t(y = Y[, sim], a = A[, sim], B_full = B_tps_full,
                                  n_uc_star  = floor(0.9 * n),
                                  n_uc_cands = seq(floor(0.9 * n), m - 3L, by = 3L),
                                  alpha = 1.0, folds = folds, cutoff = cutoff)
          n_uc_sel <- round(mean(t_by_fold_k))
          for (k in seq_len(max(folds))) {
            ho     <- which(folds == k)
            n_uc_k <- t_by_fold_k[k]
            if (confounding_mechanism != 5) {
              B_c_k  <- B_tps_full[, seq_len(m - n_uc_k), drop = FALSE]
              Ac_vec[ho] <- drop(B_c_k[ho, , drop = FALSE] %*% (t(B_c_k) %*% A[, sim]))
            } else {
              B_uc_k <- B_tps_full[, seq(m - n_uc_k + 1L, m), drop = FALSE]
              Ac_vec[ho] <- drop(B_uc_k[ho, , drop = FALSE] %*% (t(B_uc_k) %*% A[, sim]))
            }
          }
        } else {
          n_uc     <- floor(0.9 * n)
          n_uc_sel <- n_uc
          if (confounding_mechanism != 5) {
            B_c    <- B_tps_full[, seq_len(m - n_uc), drop = FALSE]
            Ac_vec <- drop(B_c %*% (t(B_c) %*% A[, sim]))
          } else {
            B_uc   <- B_tps_full[, seq(m - n_uc + 1L, m), drop = FALSE]
            Ac_vec <- drop(B_uc %*% (t(B_uc) %*% A[, sim]))
          }
        }
        xmat <- matrix(Ac_vec, ncol = 1)
        colnames(xmat) <- if (confounding_mechanism != 5) 'Ac-TPS' else 'Ac-TPS-reverse'
      }
      if (method == 'IV-GraphLaplacian'){
        m      <- ncol(B_gl_full)
        Ac_vec <- numeric(n)
        if (select_basis) {
          t_by_fold_k <- select_t(y = Y[, sim], a = A[, sim], B_full = B_gl_full,
                                  n_uc_star  = floor(0.9 * n),
                                  n_uc_cands = seq(floor(0.9 * n), m - 3L, by = 3L),
                                  alpha = 1.0, folds = folds, cutoff = cutoff)
          n_uc_sel <- round(mean(t_by_fold_k))
          for (k in seq_len(max(folds))) {
            ho     <- which(folds == k)
            n_uc_k <- t_by_fold_k[k]
            if (confounding_mechanism != 5) {
              B_c_k  <- B_gl_full[, seq_len(m - n_uc_k), drop = FALSE]
              Ac_vec[ho] <- drop(B_c_k[ho, , drop = FALSE] %*% (t(B_c_k) %*% A[, sim]))
            } else {
              B_uc_k <- B_gl_full[, seq(m - n_uc_k + 1L, m), drop = FALSE]
              Ac_vec[ho] <- drop(B_uc_k[ho, , drop = FALSE] %*% (t(B_uc_k) %*% A[, sim]))
            }
          }
        } else {
          n_uc     <- floor(0.9 * n)
          n_uc_sel <- n_uc
          if (confounding_mechanism != 5) {
            B_c    <- B_gl_full[, seq_len(m - n_uc), drop = FALSE]
            Ac_vec <- drop(B_c %*% (t(B_c) %*% A[, sim]))
          } else {
            B_uc   <- B_gl_full[, seq(m - n_uc + 1L, m), drop = FALSE]
            Ac_vec <- drop(B_uc %*% (t(B_uc) %*% A[, sim]))
          }
        }
        xmat <- matrix(Ac_vec, ncol = 1)
        colnames(xmat) <- if (confounding_mechanism != 5) 'Ac-GraphLaplacian' else 'Ac-GraphLaplacian-reverse'
      }
      if (method == 'IV-TPS-spatialcoord'){
        m      <- ncol(B_tps_full)
        Ac_vec <- numeric(n)
        if (select_basis) {
          t_by_fold_k <- select_t(y = Y[, sim], a = A[, sim], B_full = B_tps_full,
                                  n_uc_star  = floor(0.9 * n),
                                  n_uc_cands = seq(floor(0.9 * n), m - 3L, by = 3L),
                                  alpha = 1.0, folds = folds, cutoff = cutoff,
                                  x_extra = cbind(lat, lon))
          n_uc_sel <- round(mean(t_by_fold_k))
          for (k in seq_len(max(folds))) {
            ho     <- which(folds == k)
            n_uc_k <- t_by_fold_k[k]
            if (confounding_mechanism != 5) {
              B_c_k  <- B_tps_full[, seq_len(m - n_uc_k), drop = FALSE]
              Ac_vec[ho] <- drop(B_c_k[ho, , drop = FALSE] %*% (t(B_c_k) %*% A[, sim]))
            } else {
              B_uc_k <- B_tps_full[, seq(m - n_uc_k + 1L, m), drop = FALSE]
              Ac_vec[ho] <- drop(B_uc_k[ho, , drop = FALSE] %*% (t(B_uc_k) %*% A[, sim]))
            }
          }
        } else {
          n_uc     <- floor(0.9 * n)
          n_uc_sel <- n_uc
          if (confounding_mechanism != 5) {
            B_c    <- B_tps_full[, seq_len(m - n_uc), drop = FALSE]
            Ac_vec <- drop(B_c %*% (t(B_c) %*% A[, sim]))
          } else {
            B_uc   <- B_tps_full[, seq(m - n_uc + 1L, m), drop = FALSE]
            Ac_vec <- drop(B_uc %*% (t(B_uc) %*% A[, sim]))
          }
        }
        Achat <- matrix(Ac_vec, ncol = 1)
        xmat  <- cbind(Achat, lat, lon)
        colnames(xmat) <- c(if (confounding_mechanism != 5) 'Ac-TPS' else 'Ac-TPS-reverse',
                            'Latitude', 'Longitude')
      }
      if (method == 'IV-GraphLaplacian-spatialcoord'){
        m      <- ncol(B_gl_full)
        Ac_vec <- numeric(n)
        if (select_basis) {
          t_by_fold_k <- select_t(y = Y[, sim], a = A[, sim], B_full = B_gl_full,
                                  n_uc_star  = floor(0.9 * n),
                                  n_uc_cands = seq(floor(0.9 * n), m - 3L, by = 3L),
                                  alpha = 1.0, folds = folds, cutoff = cutoff,
                                  x_extra = cbind(lat, lon))
          n_uc_sel <- round(mean(t_by_fold_k))
          for (k in seq_len(max(folds))) {
            ho     <- which(folds == k)
            n_uc_k <- t_by_fold_k[k]
            if (confounding_mechanism != 5) {
              B_c_k  <- B_gl_full[, seq_len(m - n_uc_k), drop = FALSE]
              Ac_vec[ho] <- drop(B_c_k[ho, , drop = FALSE] %*% (t(B_c_k) %*% A[, sim]))
            } else {
              B_uc_k <- B_gl_full[, seq(m - n_uc_k + 1L, m), drop = FALSE]
              Ac_vec[ho] <- drop(B_uc_k[ho, , drop = FALSE] %*% (t(B_uc_k) %*% A[, sim]))
            }
          }
        } else {
          n_uc     <- floor(0.9 * n)
          n_uc_sel <- n_uc
          if (confounding_mechanism != 5) {
            B_c    <- B_gl_full[, seq_len(m - n_uc), drop = FALSE]
            Ac_vec <- drop(B_c %*% (t(B_c) %*% A[, sim]))
          } else {
            B_uc   <- B_gl_full[, seq(m - n_uc + 1L, m), drop = FALSE]
            Ac_vec <- drop(B_uc %*% (t(B_uc) %*% A[, sim]))
          }
        }
        Achat <- matrix(Ac_vec, ncol = 1)
        xmat  <- cbind(Achat, lat, lon)
        colnames(xmat) <- c(if (confounding_mechanism != 5) 'Ac-GraphLaplacian' else 'Ac-GraphLaplacian-reverse',
                            'Latitude', 'Longitude')
      }
      if (method == 'trueIV-spatialcoord'){
        xmat <- cbind(matrix(Ac[, sim], ncol = 1), lat, lon)
        colnames(xmat) <- c('Ac-TPS', 'Latitude', 'Longitude')
      }
      if (method == 'trueIV'){
        xmat <- matrix(Ac[, sim], ncol = 1)
        colnames(xmat) <- 'Ac-TPS'
      }

      # ---- fit ERF ----
      delta <- 0.05
      y <- Y[, sim]
      a <- A[, sim]
      a_sub <- a[a > cutoff - delta]
      out <- tryCatch({
        sub_idx  <- which(a > cutoff - delta)
        xsub     <- matrix(xmat[sub_idx, , drop = FALSE], ncol = ncol(xmat))
        colnames(xsub) <- colnames(xmat)
        folds_sub <- folds[sub_idx]
        folds_sub <- match(folds_sub, sort(unique(folds_sub)))  # relabel to 1..K
        ctseff(
          y     = y[sub_idx],
          a     = a_sub,
          x     = xsub,
          n.pts = 5,
          a.rng = c(cutoff - delta, cutoff + delta),
          bw.seq = seq(sd(a_sub) / 10, sd(a_sub), length.out = 100),
          folds  = folds_sub
        )
      }, error = function(e) {
        message("Error encountered: ", e$message)
        NA
      })

      # ---- estimate + CI ----
      ix_cut <- which.min(abs(out$res$a.vals - cutoff))
      muest <- (out$res$est[ix_cut] * mean(a > cutoff) +
                  mean(y[a <= cutoff]) * mean(a <= cutoff)) / mean(y)
      ci <- c(NA_real_, NA_real_)
      if (is.list(out)) {
        asym_var <- tryCatch(
          asymptotic_variance_delta(y = y, a = a, erfest = out,
                                    cutoff = cutoff, delta = delta),
          error = function(e) NA_real_
        )
        se_est <- sqrt(as.numeric(asym_var) / n)
        ci <- muest + c(-1, 1) * 1.96 * se_est
      }

      list(muest  = muest,
           ci     = ci,
           n_uc   = n_uc_sel,
           time_s = proc.time()["elapsed"] - t_sim_start)

    }, mc.cores = n_cores)

    # ---- aggregate results ----
    safe <- function(r, field, default)
      if (is.list(r)) r[[field]] else default

    muests    <- sapply(results_list, safe, "muest",  NA_real_)
    cis       <- do.call(rbind, lapply(results_list, safe, "ci",
                                       c(NA_real_, NA_real_)))
    n_uc_sels <- as.integer(sapply(results_list, safe, "n_uc", NA_integer_))
    sim_times <- sapply(results_list, safe, "time_s", NA_real_)

    t_method_elapsed <- proc.time()["elapsed"] - t_method_start
    t_total_elapsed  <- proc.time()["elapsed"] - t_simfunc_start
    message(sprintf("[%s] method=%s  sims=%d  method_time=%.1fs  total_time=%.1fs  mean_sim_time=%.2fs",
                    format(Sys.time(), "%H:%M:%S"), method, nsims,
                    t_method_elapsed, t_total_elapsed, mean(sim_times, na.rm = TRUE)))

    # Create dataframe whose first column is a.vals and the rest of cols are muests
    df <- muests
    
    # write results to file
    # Check if file exists
    if (file.exists(filename)){
      # write new sims to file as new columns
      olddf <- read.csv(filename)
      newdf <- cbind(olddf, muests)
      write.csv(newdf, filename, row.names = FALSE)
    }
    # if file for estimates does not exist create it and write results
    else{
      write.csv(df, filename, row.names = FALSE)
    }

    # Save confidence intervals as two separate CSVs (lower / upper bounds),
    # matching the estimates format: nsims rows, one column per batch.
    if (any(!is.na(cis))) {
      filename_ci_lower <- paste0(results_dir, 'conf', confounding_mechanism,
                                  '_', option, '_', method, '_ci_lower.csv')
      filename_ci_upper <- paste0(results_dir, 'conf', confounding_mechanism,
                                  '_', option, '_', method, '_ci_upper.csv')
      df_lower <- data.frame(lower = cis[, 1])
      df_upper <- data.frame(upper = cis[, 2])
      if (file.exists(filename_ci_lower)) {
        write.csv(cbind(read.csv(filename_ci_lower), df_lower),
                  filename_ci_lower, row.names = FALSE)
        write.csv(cbind(read.csv(filename_ci_upper), df_upper),
                  filename_ci_upper, row.names = FALSE)
      } else {
        write.csv(df_lower, filename_ci_lower, row.names = FALSE)
        write.csv(df_upper, filename_ci_upper, row.names = FALSE)
      }
    }

    # Save per-sim wall times (seconds)
  #   filename_time <- paste0(results_dir, 'conf', confounding_mechanism,
  #                           '_', option, '_', method, '_time.csv')
  #   df_time <- data.frame(time_s = sim_times)
  #   if (file.exists(filename_time)) {
  #     write.csv(cbind(read.csv(filename_time), df_time), filename_time, row.names = FALSE)
  #   } else {
  #     write.csv(df_time, filename_time, row.names = FALSE)
  #   }
  # 
  #   # Save n_uc selections for methods that use basis selection
    if (select_basis && any(!is.na(n_uc_sels))) {
      filename_n_uc <- paste0(results_dir, 'conf', confounding_mechanism,
                              '_', option, '_', method, '_n_uc.csv')
      df_n_uc <- data.frame(n_uc = n_uc_sels)
      if (file.exists(filename_n_uc)) {
        olddf_n_uc <- read.csv(filename_n_uc)
        write.csv(cbind(olddf_n_uc, df_n_uc), filename_n_uc, row.names = FALSE)
      } else {
        write.csv(df_n_uc, filename_n_uc, row.names = FALSE)
      }
    }
  }

  invisible(filename)
}

# Function that computes the covariance matrix of the GP
compute_Sigma_GP <- function(distmat, 
                        kappa=2, 
                        rangeu, 
                        rangec,
                        rho = 0.95,
                        sigu = 1, 
                        sigc = 1, 
                        sigz = 1){
  # distmat is the distance matrix
  # kappa is the smoothness parameter
  # rangeu is the range of the GP for the unconfounded part of exposure
  # rangec is the range of the GP for the confounded part of exposure
  # rho is the correlation between the exposure and unmeasured confounder
  # sigu, sigc, sigz are the standard deviations of the Auc, Ac, and U
  # returns the covariance matrix of the GP
  
  n <- nrow(distmat)
  phiu <- rangeu/(2*sqrt(kappa)) # to match Paciorek implementation
  phic <- rangec/(2*sqrt(kappa)) 
  Sigmau <- geoR::matern(u=distmat, phi=phiu, kappa=kappa)
  Sigmac <- geoR::matern(u=distmat, phi=phic, kappa=kappa)
  Sigma <- matrix(0, nrow = 3*n, ncol = 3*n)
  Sigma[1:n, 1:n] <- sigu^2*Sigmau
  Sigma[(n+1):(2*n), (n+1):(2*n)] <- sigc^2*Sigmac
  Sigma[(2*n+1):(3*n), (2*n+1):(3*n)] <- sigz^2*Sigmac
  # # Auc is uncorrelated + indep of Ac and U
  # Sigma[1:n, (n+1):(3*n)] <- 0
  # Sigma[(n+1):(3*n), 1:n] <- 0
  # Ac and U are highly dependent
  Sigma[(n+1):(2*n), (2*n+1):(3*n)] <- rho*sigc*sigz*Sigmac
  Sigma[(2*n+1):(3*n), (n+1):(2*n)] <- rho*sigc*sigz*Sigmac
  return(Sigma)
}

# Function that computes the covariance matrix for two confounders
compute_Sigma_GP_2U <- function(distmat,
                                kappa = 2,
                                rangeu, rangec, rangez1, rangez2,
                                rho1 = 0.9, rho2 = 0.7,
                                sigu = 1, sigc = 1, sigz1 = 1, sigz2 = 1) {
  stopifnot(abs(rho1) <= 1, abs(rho2) <= 1)
  n <- nrow(distmat)
  phi <- function(r) r / (2 * sqrt(kappa))
  
  Ku  <- geoR::matern(u = distmat, phi = phi(rangeu),  kappa = kappa)
  Kc  <- geoR::matern(u = distmat, phi = phi(rangec),  kappa = kappa)
  Kz1 <- geoR::matern(u = distmat, phi = phi(rangez1), kappa = kappa)
  Kz2 <- geoR::matern(u = distmat, phi = phi(rangez2), kappa = kappa)
  
  Sigma <- matrix(0, nrow = 4*n, ncol = 4*n)
  iAuc <- 1:n; iAc <- (n+1):(2*n); iU1 <- (2*n+1):(3*n); iU2 <- (3*n+1):(4*n)
  
  # Auc (independent)
  Sigma[iAuc, iAuc] <- sigu^2 * Ku
  
  # Coregionalized part on Kc (Ac, U1, U2 share it)
  Sigma[iAc, iAc]   <- sigc^2 * Kc
  Sigma[iU1, iU1]   <- (rho1^2) * sigz1^2 * Kc
  Sigma[iU2, iU2]   <- (rho2^2) * sigz2^2 * Kc
  
  Sigma[iAc, iU1] <- rho1 * sigc * sigz1 * Kc
  Sigma[iU1, iAc] <- t(Sigma[iAc, iU1])
  
  Sigma[iAc, iU2] <- rho2 * sigc * sigz2 * Kc
  Sigma[iU2, iAc] <- t(Sigma[iAc, iU2])
  
  # U1–U2 correlation induced by sharing Kc
  Sigma[iU1, iU2] <- (rho1 * rho2 * sigz1 * sigz2) * Kc
  Sigma[iU2, iU1] <- t(Sigma[iU1, iU2])
  
  # Idiosyncratic scales for U1, U2 (their own kernels)
  Sigma[iU1, iU1] <- Sigma[iU1, iU1] + (1 - rho1^2) * sigz1^2 * Kz1
  Sigma[iU2, iU2] <- Sigma[iU2, iU2] + (1 - rho2^2) * sigz2^2 * Kz2
  
  return(Sigma)
}


# Baseline GP data-generating process (supplement Section 4.2, mechanism 1).
# Also used for the reversed-scale mechanism (5) via different theta values.
# Generates 4 independent mean-zero GP fields, then assembles:
#   Auc = 0.6*B1_uc + 0.4*B2_uc
#   Ac  = 1.0*B1_c  + 0.8*B2_c
#   U   = rho1*B1_c + rho2*B2_c + Zu
compute_data_GP <- function(nsims, distmat,
                            theta_B1uc = 0.10, theta_B2uc = 0.05,
                            theta_B1c  = 0.80, theta_B2c  = 0.50,
                            theta_u    = 1.00,
                            rho1 = 0.8, rho2 = 0.6,
                            kappa = 2) {
  n   <- nrow(distmat)
  phi <- function(r) r / (2 * sqrt(kappa))
  # Sample nsims independent draws from N(0, K); always return n x nsims matrix.
  gp  <- function(K) {
    s <- MASS::mvrnorm(nsims, rep(0, n), K)
    if (nsims == 1L) matrix(s, n, 1L) else t(s)
  }
  B1uc <- gp(geoR::matern(u = distmat, phi = phi(theta_B1uc), kappa = kappa))
  B2uc <- gp(geoR::matern(u = distmat, phi = phi(theta_B2uc), kappa = kappa))
  B1c  <- gp(geoR::matern(u = distmat, phi = phi(theta_B1c),  kappa = kappa))
  B2c  <- gp(geoR::matern(u = distmat, phi = phi(theta_B2c),  kappa = kappa))
  Zu   <- gp(geoR::matern(u = distmat, phi = phi(theta_u),    kappa = kappa))
  list(
    Auc = 0.6 * B1uc + 0.4 * B2uc,
    Ac  = 1.0 * B1c  + 0.8 * B2c,
    U   = rho1 * B1c + rho2 * B2c + Zu
  )
}

# Two-confounder GP data-generating process (supplement Section 4.2, mechanism 4).
# Same 4 latent fields as mechanism 1; two confounders each linked to B1_c and B2_c.
#   U1 = rho1*B1_c + rho2*B2_c + Zu1  (Zu1 ~ GP(R(theta_u1)))
#   U2 = rho3*B1_c + rho4*B2_c + Zu2  (Zu2 ~ GP(R(theta_u2)))
compute_data_GP_2U <- function(nsims, distmat,
                               theta_B1uc = 0.10, theta_B2uc = 0.05,
                               theta_B1c  = 0.80, theta_B2c  = 0.50,
                               theta_u1   = 0.50, theta_u2   = 0.30,
                               rho1 = 0.8, rho2 = 0.6,
                               rho3 = 0.5, rho4 = 0.4,
                               kappa = 2) {
  n   <- nrow(distmat)
  phi <- function(r) r / (2 * sqrt(kappa))
  gp  <- function(K) {
    s <- MASS::mvrnorm(nsims, rep(0, n), K)
    if (nsims == 1L) matrix(s, n, 1L) else t(s)
  }
  B1uc <- gp(geoR::matern(u = distmat, phi = phi(theta_B1uc), kappa = kappa))
  B2uc <- gp(geoR::matern(u = distmat, phi = phi(theta_B2uc), kappa = kappa))
  B1c  <- gp(geoR::matern(u = distmat, phi = phi(theta_B1c),  kappa = kappa))
  B2c  <- gp(geoR::matern(u = distmat, phi = phi(theta_B2c),  kappa = kappa))
  Zu1  <- gp(geoR::matern(u = distmat, phi = phi(theta_u1),   kappa = kappa))
  Zu2  <- gp(geoR::matern(u = distmat, phi = phi(theta_u2),   kappa = kappa))
  list(
    Auc = 0.6 * B1uc + 0.4 * B2uc,
    Ac  = 1.0 * B1c  + 0.8 * B2c,
    U1  = rho1 * B1c + rho2 * B2c + Zu1,
    U2  = rho3 * B1c + rho4 * B2c + Zu2
  )
}


# Within-state GP mechanism (supplement Section 4.2, mechanism 3).
# Applies mechanism 1 independently within each state so fields are smooth
# within states but discontinuous across administrative boundaries.
compute_data_GP_state <- function(nsims, distmat, statemat) {
  nobs <- nrow(distmat)
  nst  <- ncol(statemat)
  out  <- list(Auc = matrix(NA_real_, nobs, nsims),
               Ac  = matrix(NA_real_, nobs, nsims),
               U   = matrix(NA_real_, nobs, nsims))
  for (i in seq_len(nst)) {
    ixs    <- which(statemat[, i] == 1)
    dat_st <- compute_data_GP(nsims = nsims, distmat = distmat[ixs, ixs])
    out$Auc[ixs, ] <- dat_st$Auc
    out$Ac[ixs, ]  <- dat_st$Ac
    out$U[ixs, ]   <- dat_st$U
  }
  return(out)
}

# Coordinate-based nonlinear mechanism (supplement Section 4.2, mechanism 6).
# U = sin(2*pi*lat*long) + lat + long  (deterministic, normalized coordinates)
# B1_c = U, B2_c = 0;  B1_uc, B2_uc ~ independent GP fields
# Auc = 0.6*B1_uc + 0.4*B2_uc
# Ac  = 1.0*U  (= 1.0*B1_c + 0.8*0)
compute_data_spatialcoord <- function(lat, long, nsims, distmat) {
  n      <- length(lat)
  lat_n  <- (lat  - min(lat))  / (max(lat)  - min(lat))
  long_n <- (long - min(long)) / (max(long) - min(long))
  U      <- sin(2 * pi * lat_n * long_n) + lat_n + long_n
  phi    <- function(r, kappa = 2) r / (2 * sqrt(kappa))
  gp     <- function(K) {
    s <- MASS::mvrnorm(nsims, rep(0, n), K)
    if (nsims == 1L) matrix(s, n, 1L) else t(s)
  }
  B1uc  <- gp(geoR::matern(u = distmat, phi = phi(0.10), kappa = 2))
  B2uc  <- gp(geoR::matern(u = distmat, phi = phi(0.05), kappa = 2))
  list(
    Auc = 0.6 * B1uc + 0.4 * B2uc,
    Ac  = matrix(U, n, nsims),   # 1.0*B1_c + 0.8*0, B1_c = U
    U   = matrix(U, n, nsims)
  )
}

asymptotic_variance_delta <- function(y, a, erfest, cutoff, delta){
  # Calculate parameters
  ix_cut <- which.min(abs(erfest$res$a.vals - cutoff))
  theta1 <- erfest$res$est[ix_cut]
  theta2 <- mean(a <= cutoff)
  theta3 <- mean(y[a <= cutoff])
  theta4 <- mean(y)
  
  # Calculate estimated influence functions
  # phi1 comes from ctseff run on the n_sub-observation subset (a > cutoff - delta).
  # Padding with zeros dilutes its variance by n_sub/n; rescale by n/n_sub so that
  # cov(phi1,...)/n correctly estimates Var(hat_theta1) = Var(phi1_Kennedy)/n_sub.
  n <- length(a)
  n_sub <- sum(a > cutoff - delta)
  phi1 <- rep(0, n)
  phi1[a > cutoff - delta] <- erfest$phi[[ix_cut]] * (n / n_sub)
  phi2 <- 1*(a < cutoff) - mean(a < cutoff)
  phi3 <- rep(NA, length(a))
  phi3[a <= cutoff] <- (y[a <= cutoff] - mean(y[a <= cutoff])) / theta2
  phi3[a > cutoff] <- 0
  phi4 <- y - mean(y)
  
  # Calculate the 4 x 4 covariance matrix of the IFs
  Sigma <- cov(cbind(phi1, phi2, phi3, phi4))
  
  # Calculate partial derivatives of (theta1(1-theta2) + theta3*theta2)/theta4
  grad <- c((1-theta2)/theta4,
            (-theta1 + theta3)/theta4,
            theta2/theta4,
            -(theta1*(1-theta2) + theta2*theta3)/theta4^2)
  # Return asymptotic variance via delta method
  return(t(grad) %*% Sigma %*% grad)
}

hausdorff_distance <- function(interval1, interval2){
  dist1 <- abs(interval1[1] - interval2[1])
  dist2 <- abs(interval1[2] - interval2[2])
  return(max(dist1,dist2))
}

rleroux_univariate <- function(W, rho_sp, tau_sp = 1) {
  n  <- nrow(W)
  D  <- Diagonal(n, rowSums(W))
  Qs <- tau_sp * ((1 - rho_sp) * Diagonal(n) + rho_sp * (D - W))
  Qs <- forceSymmetric(Qs)
  cf <- Cholesky(Qs, LDL = FALSE, perm = TRUE, super = TRUE)
  Z  <- rnorm(n)
  Yp <- solve(cf, Z, system = "L")
  Xp <- solve(cf, Yp, system = "Lt")
  Xp[order(cf@perm)]
}

rleroux_bivariate <- function(W, rho_sp, tau_sp, R){
  n  <- nrow(W)
  D  <- Diagonal(n, rowSums(W))
  Qs <- tau_sp * ((1 - rho_sp) * Diagonal(n) + rho_sp * (D - W))
  Qs <- forceSymmetric(Qs)
  
  Q  <- kronecker(R, Qs)           # joint precision (2n x 2n), sparse
  cf <- Cholesky(Q, LDL=FALSE, perm=TRUE, super=TRUE)
  
  Z  <- rnorm(2*n)
  Yp <- solve(cf, Z, system="L")
  Xp <- solve(cf, Yp, system="Lt")
  X  <- Xp[order(cf@perm)]
  list(phi1 = X[1:n], phi2 = X[(n+1):(2*n)])
}

# Bivariate Leroux CAR mechanism (supplement Section 4.2, mechanism 2).
#   (U, B1_c) ~ bivariate Leroux with rho_c=0.8, cross-corr 0.7
#   B2_c ~ univariate Leroux(rho_c) independently
#   B1_uc, B2_uc ~ univariate Leroux(rho_uc=0.1) independently
#   Auc = 0.6*B1_uc + 0.4*B2_uc,  Ac = 1.0*B1_c + 0.8*B2_c
compute_data_leroux <- function(W, nsims) {
  n      <- nrow(W)
  stopifnot(ncol(W) == n)
  rho_c  <- 0.8
  rho_uc <- 0.1
  R_c    <- solve(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  out <- list(Auc = matrix(NA_real_, n, nsims),
              Ac  = matrix(NA_real_, n, nsims),
              U   = matrix(NA_real_, n, nsims))
  for (i in seq_len(nsims)) {
    smp    <- rleroux_bivariate(W = W, rho_sp = rho_c, tau_sp = 1, R = R_c)
    U_i    <- smp$phi1
    B1c_i  <- smp$phi2
    B2c_i  <- rleroux_univariate(W = W, rho_sp = rho_c)
    B1uc_i <- rleroux_univariate(W = W, rho_sp = rho_uc)
    B2uc_i <- rleroux_univariate(W = W, rho_sp = rho_uc)
    out$U[, i]   <- U_i
    out$Ac[, i]  <- 1.0 * B1c_i  + 0.8 * B2c_i
    out$Auc[, i] <- 0.6 * B1uc_i + 0.4 * B2uc_i
  }
  return(out)
}

