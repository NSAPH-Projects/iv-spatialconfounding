# Function used to estimate the exposure-response curve
ctseff <- function(y, a, x, bw.seq, n.pts = 100, a.rng = c(min(a), max(a)),
                   sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean", "SL.earth"),
                   constrain = T, trim = 0.01,
                   folds = NULL,
                   y_full = NULL, a_full = NULL, x_full = NULL,
                   folds_full = NULL, sub_rows = NULL) {
  # y is outcome
  # a is exposure
  # x is covariate matrix
  # bw.seq is a sequence of bandwidth values
  # a.rng is the range of exposure values to evaluate the ERF
  # n.pts is the number of points within a.rng at which to evaluate the ERF
  # sl.lib is the library of SuperLearner algorithms to use
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

  # Helper: run one cross-fitting pass on (y_cf, a_cf, x_cf) with given folds.
  # Returns list(pihat, pihat.mat, muhat, muhat.mat) of length n_cf and n_cf x n.pts.
  .crossfit <- function(y_cf, a_cf, x_cf, folds_cf) {
    n_cf      <- length(a_cf)
    K_cf      <- max(folds_cf)
    pihat_cf     <- numeric(n_cf)
    pihat_mat_cf <- matrix(NA_real_, n_cf, length(a.vals))
    muhat_cf     <- numeric(n_cf)
    muhat_mat_cf <- matrix(NA_real_, n_cf, length(a.vals))
    cv2 <- list(V = 2L)

    for (k in seq_len(K_cf)) {
      tr   <- which(folds_cf != k)
      ho   <- which(folds_cf == k)
      n_ho <- length(ho)

      x_tr  <- data.frame(x_cf[tr, , drop = FALSE])
      xa_tr <- data.frame(cbind(x_cf[tr, , drop = FALSE], a = a_cf[tr]))
      colnames(xa_tr) <- c(colnames(x_cf), "a")

      # newX: all n_cf obs + held-out fold x a.vals grid
      x_ho_grid  <- data.frame(x_cf[ho[rep(seq_len(n_ho), length(a.vals))], , drop = FALSE])
      colnames(x_ho_grid) <- colnames(x_cf)
      x_new_k    <- rbind(data.frame(x_cf), x_ho_grid)

      xa_ho_grid <- data.frame(cbind(
        x_cf[ho[rep(seq_len(n_ho), length(a.vals))], , drop = FALSE],
        a = rep(a.vals, each = n_ho)
      ))
      colnames(xa_ho_grid) <- c(colnames(x_cf), "a")
      xa_all_cf   <- data.frame(cbind(x_cf, a = a_cf))
      colnames(xa_all_cf) <- c(colnames(x_cf), "a")
      xa_new_k <- rbind(xa_all_cf, xa_ho_grid)

      # cvControl = V=2 reduces SuperLearner's internal CV fold count, preventing
      # GAM/glm.interaction from running out of residual df on the smaller fold samples.
      pimod_k  <- SuperLearner(Y = a_cf[tr], X = x_tr, SL.library = sl.lib,
                               newX = x_new_k, cvControl = cv2)
      pi_k     <- pimod_k$SL.predict
      pi2mod_k <- SuperLearner(Y = log(pmax((a_cf[tr] - pi_k[tr])^2, .Machine$double.eps)),
                               X = x_tr, SL.library = sl.lib,
                               newX = x_new_k, cvControl = cv2)
      pi2_k    <- pmax(pmin(exp(pi2mod_k$SL.predict),
                           100 * var(a_cf[tr] - pi_k[tr])),
                      1e-4)
      mumod_k  <- SuperLearner(Y = y_cf[tr], X = xa_tr, SL.library = sl.lib,
                               newX = xa_new_k, cvControl = cv2)
      mu_k     <- mumod_k$SL.predict

      # Density of standardised residuals, estimated from training fold.
      # Compute a_std_grid first so the density range covers both observed
      # and grid values, preventing approx() from returning NA out-of-range.
      a_std_k    <- (a_cf - pi_k[1:n_cf]) / sqrt(pi2_k[1:n_cf])
      grid_idx   <- (n_cf + 1):(n_cf + n_ho * length(a.vals))
      a_std_grid <- (rep(a.vals, each = n_ho) - pi_k[grid_idx]) / sqrt(pi2_k[grid_idx])
      a_std_all  <- c(a_std_k, a_std_grid)
      dens_k     <- density(a_std_k[tr], from = min(a_std_all), to = max(a_std_all))

      pihat_cf[ho] <- approx(dens_k$x, dens_k$y, xout = a_std_k[ho], rule = 2)$y /
                        sqrt(pi2_k[ho])
      pihat_mat_cf[ho, ] <- matrix(
        approx(dens_k$x, dens_k$y, xout = a_std_grid, rule = 2)$y / sqrt(pi2_k[grid_idx]),
        nrow = n_ho, ncol = length(a.vals)
      )
      muhat_mat_cf[ho, ] <- matrix(mu_k[grid_idx], nrow = n_ho, ncol = length(a.vals))
      muhat_cf[ho]       <- mu_k[ho]
    }
    pihat_cf <- pmax(pihat_cf, quantile(pihat_cf, trim))
    list(pihat = pihat_cf, pihat.mat = pihat_mat_cf,
         muhat = muhat_cf, muhat.mat = muhat_mat_cf)
  }

  if (!is.null(a_full) && !is.null(folds_full) && !is.null(sub_rows)) {
    # Preferred path: cross-fit nuisances on the full dataset (n_full obs), then
    # extract the subset rows for ERF fitting.  Estimating on all data rather than
    # the truncated subset (a > cutoff - delta) gives larger training folds and
    # avoids selection bias in the nuisance models.
    cf <- .crossfit(y_full, a_full, data.frame(x_full), folds_full)
    pihat     <- cf$pihat[sub_rows]
    pihat.mat <- cf$pihat.mat[sub_rows, , drop = FALSE]
    muhat     <- cf$muhat[sub_rows]
    muhat.mat <- cf$muhat.mat[sub_rows, , drop = FALSE]

  } else if (!is.null(folds)) {
    # Fallback: cross-fit on the subset passed to ctseff.
    # K-fold cross-fitting of nuisance functions (supplement Section 3, step 3d-e).
    cf <- .crossfit(y, a, data.frame(x), folds)
    pihat     <- cf$pihat
    pihat.mat <- cf$pihat.mat
    muhat     <- cf$muhat
    muhat.mat <- cf$muhat.mat

  } else {
    # In-sample estimation (used inside select_t and when folds not provided).
    pimod      <- SuperLearner(Y = a, X = data.frame(x), SL.library = sl.lib, newX = x.new)
    pimod.vals <- pimod$SL.predict
    pi2mod     <- SuperLearner(Y = log(pmax((a - pimod.vals[1:n])^2, .Machine$double.eps)),
                               X = x, SL.library = sl.lib, newX = x.new)
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
  }

  # construct varpi/m from pihat.mat and muhat.mat (common to both paths)
  varpihat     <- predict(smooth.spline(a.vals, apply(pihat.mat, 2, mean)), x = a)$y
  # Smooth-spline extrapolates beyond a.vals and can go negative for observations
  # with a >> max(a.vals). Clip from below at the same floor used for pihat so
  # the ratio pihat/varpihat retains the correct sign and stays bounded.
  varpihat     <- pmax(varpihat, quantile(pihat, trim))
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = T, nrow = n)
  mhat         <- predict(smooth.spline(a.vals, apply(muhat.mat, 2, mean)), x = a)$y
  mhat.mat     <- matrix(rep(apply(muhat.mat, 2, mean), n), byrow = T, nrow = n)
  
  
  # form adjusted/pseudo outcome xi
  pseudo.out <- (y - muhat) / (pihat / varpihat) + mhat
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
    wvals <- w.fn(bw, a.vals = asubset)
    tryCatch(
      approx(asubset, wvals, xout = a)$y, # sophie's change
      error = function(e) rep(NA_real_, length(a))
    )
  }
  cts.eff.fn <- function(out, bw) {
    tryCatch(
      approx(locpoly(a, out, bandwidth = bw), xout = a)$y,
      error = function(e) rep(NA_real_, length(a))
    )
  }
  # note: choice of bandwidth range depends on specific problem,
  # make sure to inspect plot of risk as function of bandwidth
  risk.fn <- function(h) {
    hats <- hatvals(h)
    mean(((pseudo.out - cts.eff.fn(pseudo.out, bw = h)) / (1 - hats))^2)
  } 
  risk.est <- sapply(bw.seq, risk.fn)
  if (mean(is.finite(risk.est)) < 0.5)
    stop("bandwidth selection failed: fewer than half of bandwidths gave finite risk")
  h.opt <- bw.seq[which.min(risk.est)]
  bw.risk <- data.frame(bw = bw.seq, risk = risk.est)
  #print('calculated h.opt')
  
  # alternative approach:
  # h.opt <- optimize(function(h){ hats <- hatvals(h); mean( ((pseudo.out[a > a.min & a < a.max]-cts.eff.fn(pseudo.out,bw=h))/(1-hats))^2) } ,
  #  bw.seq, tol=0.01)$minimum
  
  # estimate effect curve with optimal bandwidth
  est <- tryCatch(
    approx(locpoly(a, pseudo.out, bandwidth = h.opt), xout = a.vals)$y,
    error = function(e) rep(NA_real_, length(a.vals))
  )
  est <- pmin(pmax(est, min(pseudo.out, na.rm = TRUE)), max(pseudo.out, na.rm = TRUE))
  # print(summary(est))
  
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
      Ys[,i] <- rnorm(n, -2 + (-1)*Us[,i] + As[,i] - 0.5*As[,i]*Us[,i], 1) 
    }
  }
  # nonlinear outcome model
  if (option == 'nonlinear'){
    for (i in 1:nreps){
      eta <- -2 - Us[,i] + As[,i] - 0.4*As[,i]^2 - 0.25*Us[,i]*As[,i]
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
            legend.key.width = unit(80, "points"),
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

# Partition n observations into K random (non-spatial) folds.
# Returns an integer vector of fold assignments (values 1..K).
make_random_folds <- function(n, K = 5L) {
  sample(rep_len(1:K, n))
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
# Returns an integer vector of length K (all equal): globally selected instrument count
# chosen by the one-SE rule on cross-fold variance of the psi path.
select_t <- function(y, a, B_full, n_uc_star, n_uc_cands, alpha = 5.0,
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
          bw.seq = seq(sd(a) / 10, sd(a), length.out = 50L)
        )
        ix <- which.min(abs(erf$res$a.vals - cutoff))
        (erf$res$est[ix] * mean(a_tr > cutoff) +
           mean(y_tr[a_tr <= cutoff]) * mean(a_tr <= cutoff)) / mean(y_tr)
      }, error = function(e) NA_real_)
    }
  }

  # Fold-specific one-SE rule (reference-based).
  # For fold k, use only the OTHER K-1 rows of psi_by_fold to select n_uc_k.
  # This ensures the selection for fold k is independent of fold k's data entirely
  # (neither its holdout nor its training set), supporting valid post-selection inference.

  # Guard: exclude candidates where Ac has too few columns (near-zero IV).
  # n_uc_cands[n_cand] = m gives Ac with 0 columns and a degenerate (likely NA) psi.
  min_cols_c <- max(5L, floor(0.02 * m))

  t_by_fold <- integer(K)
  for (k in seq_len(K)) {
    other       <- psi_by_fold[-k, , drop = FALSE]          # (K-1) x n_cand
    psi_bar_k   <- colMeans(other, na.rm = TRUE)
    psi_se_k    <- apply(other, 2L,
                         function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

    valid_k <- which((m - n_uc_cands) >= min_cols_c & !is.na(psi_bar_k))
    if (length(valid_k) == 0L) {
      t_by_fold[k] <- n_uc_cands[1L]
      next
    }

    ref_j_k     <- max(valid_k)
    dist_k      <- abs(psi_bar_k - psi_bar_k[ref_j_k])
    thr_k       <- alpha * psi_se_k[ref_j_k]

    eligible_k  <- valid_k[dist_k[valid_k] <= thr_k]
    j_star_k    <- if (length(eligible_k) > 0L) min(eligible_k) else ref_j_k
    t_by_fold[k] <- n_uc_cands[j_star_k]

    message(sprintf("select_t fold %d: ref_j=%d  n_uc_ref=%d  thr=%.4f  n_uc_sel=%d",
                    k, ref_j_k, n_uc_cands[ref_j_k], thr_k, t_by_fold[k]))
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
                     # 'baseline',
                      'oracle',
                     # 'spatialcoord',
                     # 'trueIV',
                     # 'trueIV-spatialcoord',
                     'IV-TPS',
                     'IV-GraphLaplacian'#,
                     # 'IV-TPS-spatialcoord',
                     # 'IV-GraphLaplacian-spatialcoord'
                   ),
                   B_tps_full,
                   B_gl_full,
                   statemat,
                   W = NULL,
                   cutoff = 1,
                   select_basis = TRUE,
                   n_cores = 1L,
                   spatial_folds = FALSE,
                   results_dir = "results_Apr9/",
                   iv_control = list())
{
  # nsims is the number of simulations
  # lat is a vector of latitudes
  # lon is a vector of longitudes
  # option is the form of the outcome model
  # methods are the methods used to estimate truncated exposure effect
  # statemat is the matrix of state-level indicators
  # cutoff is c
  # select_basis: if TRUE, use the manuscript selection rule; otherwise use a
  # single fixed candidate containing 90% of the ordered basis elements.
  # results_dir: directory to write output CSVs

  iv_defaults <- list(
    alpha = 1,
    density_trim = 0,
    ipw_ratio_trim = 0,
    core_fraction = 0.9,
    step = 5L,
    min_confounded = NULL,
    fixed_n_uc_fraction = 0.9,
    n_grid = 25L,
    sl_library = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction"),
    candidate_bandwidth = NULL,
    candidate_bw_seq = NULL,
    final_bandwidth = NULL,
    final_bw_seq = NULL,
    constrain = TRUE,
    save_diagnostics = TRUE
  )
  if (!is.list(iv_control)) {
    stop("iv_control must be a list.")
  }
  iv_control <- utils::modifyList(iv_defaults, iv_control)
  iv_methods <- c(
    "IV-TPS",
    "IV-GraphLaplacian",
    "IV-TPS-spatialcoord",
    "IV-GraphLaplacian-spatialcoord"
  )
  required_estimator_functions <- "estimate_crossfit_truncated_effect"
  if (any(methods %in% iv_methods)) {
    required_estimator_functions <- c(
      required_estimator_functions,
      "make_candidate_grid",
      "estimate_truncated_effect_candidate",
      "estimate_selected_truncated_effect"
    )
  }
  missing_estimator_functions <- required_estimator_functions[
    !vapply(required_estimator_functions, exists, logical(1), mode = "function")
  ]
  if (length(missing_estimator_functions) > 0L) {
      stop(
        paste0(
          "Source the manuscript-aligned R modules before simfunc(): ",
          paste(missing_estimator_functions, collapse = ", ")
        )
      )
  }

  confounding_mechanism <- as.integer(confounding_mechanism)
  option <- match.arg(option)
  
  ################# GENERATE DATA #################
  
  # Compute distance matrix
  distmat <- geosphere::distm(cbind(lon, lat), 
                              fun = geosphere::distHaversine)
  distmat <- distmat/1000000 # scale so range (0,2)
  n <- length(lat)
  # Mechanism 1: TPS-basis confounding (IV-TPS recovers Ac exactly)
  if (confounding_mechanism == 1){
    dat <- compute_data_TPS_basis(B_tps_full = B_tps_full, nsims = nsims)
  }
  # Mechanism 2: GL-basis confounding (IV-GraphLaplacian recovers Ac exactly)
  if (confounding_mechanism == 2){
    dat <- compute_data_GL_basis(B_gl_full = B_gl_full, nsims = nsims)
  }
  # Mechanism 3: two-confounder TPS-basis
  if (confounding_mechanism == 3){
    dat <- compute_data_TPS_2U(B_tps_full = B_tps_full, nsims = nsims)
    Ac  <- dat$Ac;  Auc <- dat$Auc;  U1 <- dat$U1;  U2 <- dat$U2
    A   <- Ac + Auc
    Y   <- matrix(NA, n, nsims)
    if (option == 'linear'){
      for (i in 1:nsims)
        Y[, i] <- rnorm(n, -2 + (-1)*U1[,i] + A[,i] - 0.5*A[,i]*U1[,i] - 0.75*A[,i]*U2[,i], 1)
    }
    if (option == 'nonlinear'){
      for (i in 1:nsims){
        eta <- -2 - U1[,i] + A[,i] - 0.4*A[,i]^2 - 0.25*U1[,i]*A[,i] - 0.25*U2[,i]*A[,i]
        Y[, i] <- rnorm(n, eta, 1)
      }
    }
  }
  # Mechanism 4: reversed-scale TPS (confounding in high-frequency TPS components)
  if (confounding_mechanism == 4){
    dat <- compute_data_TPS_reversed(B_tps_full = B_tps_full, nsims = nsims)
  }
  # Mechanism 5: bivariate Leroux CAR
  if (confounding_mechanism == 5){
    stopifnot(!is.null(W))
    dat <- compute_data_leroux(W = W, nsims = nsims)
  }
  # Mechanism 6: two-confounder GL-basis
  if (confounding_mechanism == 6){
    dat <- compute_data_GL_2U(B_gl_full = B_gl_full, nsims = nsims)
    Ac  <- dat$Ac;  Auc <- dat$Auc;  U1 <- dat$U1;  U2 <- dat$U2
    A   <- Ac + Auc
    Y   <- matrix(NA, n, nsims)
    if (option == 'linear'){
      for (i in 1:nsims)
        Y[, i] <- rnorm(n, -2 + (-1)*U1[,i] + A[,i] - 0.5*A[,i]*U1[,i] - 0.75*A[,i]*U2[,i], 1)
    }
    if (option == 'nonlinear'){
      for (i in 1:nsims){
        eta <- -2 - U1[,i] + A[,i] - 0.4*A[,i]^2 - 0.25*U1[,i]*A[,i] - 0.25*U2[,i]*A[,i]
        Y[, i] <- rnorm(n, eta, 1)
      }
    }
  }
  # Mechanism 7: coordinate-based nonlinear
  if (confounding_mechanism == 7){
    dat <- compute_data_spatialcoord(lat = lat, long = lon, nsims = nsims, distmat = distmat)
  }
  # Mechanism 8: reversed-scale GL (confounding in high-frequency GL components)
  if (confounding_mechanism == 8){
    dat <- compute_data_GL_reversed(B_gl_full = B_gl_full, nsims = nsims)
  }

  if (!confounding_mechanism %in% c(3, 6)){
    Ac <- dat$Ac
    Auc <- dat$Auc
    U <- dat$U
    A <- Ac + Auc # all have dimension n x nsims
    Y <- createY(Us=U, As=A, option = option)
  }
  
  
  ################# FIT MODELS #################

  # Normalize coordinates to [0,1] once; used by all spatialcoord methods so
  # that SuperLearner splines operate on the same scale as the DGP basis functions.
  lat_n_m <- (lat  - min(lat))  / (max(lat)  - min(lat))
  lon_n_m <- (lon - min(lon)) / (max(lon) - min(lon))

  # Pre-compute outer folds once. The default is the uniform random split in
  # Supplement Algorithm 1; spatial folds remain available as a sensitivity run.
  folds <- if (spatial_folds) make_spatial_folds(lat, lon) else make_random_folds(n)
  n_outer_folds <- length(unique(folds))

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
      n_uc_folds <- rep(NA_integer_, n_outer_folds)
      selection_diagnostics <- NULL
      cor_Ac   <- NA_real_
      y <- Y[, sim]
      a <- A[, sim]

      # The IV methods use the complete manuscript algorithm: fold-specific
      # selection, training-only A^c projection, held-out pseudo-outcomes, and
      # pooled doubly robust estimation. Non-IV methods continue through the
      # legacy ctseff() path below.
      if (method %in% iv_methods) {
        B_method <- if (method %in% c("IV-TPS", "IV-TPS-spatialcoord")) {
          B_tps_full
        } else {
          B_gl_full
        }
        x_extra <- if (method %in% c(
          "IV-TPS-spatialcoord",
          "IV-GraphLaplacian-spatialcoord"
        )) {
          out <- cbind(lat_n_m, lon_n_m)
          colnames(out) <- c("Latitude", "Longitude")
          out
        } else {
          NULL
        }
        instrument_scale <- if (confounding_mechanism %in% c(4L, 8L)) {
          "large"
        } else {
          "small"
        }
        m <- ncol(B_method)
        n_uc_cands <- if (select_basis) {
          make_candidate_grid(
            m = m,
            core_fraction = iv_control$core_fraction,
            step = iv_control$step,
            min_confounded = iv_control$min_confounded
          )
        } else {
          as.integer(floor(iv_control$fixed_n_uc_fraction * m))
        }

        iv_fit <- tryCatch(
          estimate_selected_truncated_effect(
            y = y,
            a = a,
            x = x_extra,
            B = B_method,
            cutoff = cutoff,
            n_uc_cands = n_uc_cands,
            outer_folds = folds,
            alpha = iv_control$alpha,
            density_trim = iv_control$density_trim,
            ipw_ratio_trim = iv_control$ipw_ratio_trim,
            instrument_scale = instrument_scale,
            candidate_args = list(
              nuisance_args = list(
                sl_library = iv_control$sl_library,
                density_trim = iv_control$density_trim
              ),
              n_grid = iv_control$n_grid,
              bandwidth = iv_control$candidate_bandwidth,
              bw_seq = iv_control$candidate_bw_seq,
              constrain = iv_control$constrain,
              density_trim = iv_control$density_trim,
              ipw_ratio_trim = iv_control$ipw_ratio_trim
            ),
            final_nuisance_args = list(
              sl_library = iv_control$sl_library,
              density_trim = iv_control$density_trim
            ),
            n_grid = iv_control$n_grid,
            bandwidth = iv_control$final_bandwidth,
            bw_seq = iv_control$final_bw_seq,
            constrain = iv_control$constrain
          ),
          error = function(e) {
            message("Manuscript IV estimator error: ", e$message)
            NULL
          }
        )

        if (is.null(iv_fit)) {
          return(list(
            muest = NA_real_,
            ci = c(NA_real_, NA_real_),
            n_uc = NA_integer_,
            n_uc_folds = n_uc_folds,
            selection_diagnostics = NULL,
            cor_Ac = NA_real_,
            time_s = proc.time()["elapsed"] - t_sim_start
          ))
        }

        n_uc_folds <- as.integer(iv_fit$selected_n_uc)
        n_uc_sel <- as.integer(round(mean(n_uc_folds)))
        selection_diagnostics <- iv_fit$selection_diagnostics
        attr(selection_diagnostics, "final_diagnostics") <- iv_fit$diagnostics
        attr(selection_diagnostics, "fold_summaries") <- iv_fit$fold_summaries
        cor_Ac <- stats::cor(iv_fit$Ac_crossfit, Ac[, sim])
        return(list(
          muest = iv_fit$psi,
          ci = iv_fit$confidence_interval,
          n_uc = n_uc_sel,
          n_uc_folds = n_uc_folds,
          selection_diagnostics = selection_diagnostics,
          cor_Ac = cor_Ac,
          time_s = proc.time()["elapsed"] - t_sim_start
        ))
      }

      if (method == 'baseline'){
        xmat <- matrix(rep(1, n), ncol = 1)
        colnames(xmat) <- 'Intercept'
      }
      if (method == 'oracle'){
        if (!confounding_mechanism %in% c(3, 6)){
          xmat <- matrix(U[, sim], ncol = 1)
          colnames(xmat) <- 'U'
        } else {
          xmat <- cbind(U1[, sim], U2[, sim])
          colnames(xmat) <- c('U1', 'U2')
        }
      }
      if (method == 'spatialcoord'){
        xmat <- cbind(lat_n_m, lon_n_m)
        colnames(xmat) <- c('Latitude', 'Longitude')
      }
      if (method == 'trueIV-spatialcoord'){
        xmat <- cbind(matrix(Ac[, sim], ncol = 1), lat_n_m, lon_n_m)
        colnames(xmat) <- c('Ac_TPS', 'Latitude', 'Longitude')
      }
      if (method == 'trueIV'){
        xmat <- matrix(Ac[, sim], ncol = 1)
        colnames(xmat) <- 'Ac_TPS'
      }

      non_iv_fit <- tryCatch(
        estimate_crossfit_truncated_effect(
          y = y,
          a = a,
          w = xmat,
          cutoff = cutoff,
          outer_folds = folds,
          density_trim = iv_control$density_trim,
          ipw_ratio_trim = iv_control$ipw_ratio_trim,
          final_nuisance_args = list(
            sl_library = iv_control$sl_library,
            density_trim = iv_control$density_trim
          ),
          n_grid = iv_control$n_grid,
          bandwidth = iv_control$final_bandwidth,
          bw_seq = iv_control$final_bw_seq,
          constrain = iv_control$constrain
        ),
        error = function(e) {
          message("Manuscript non-IV estimator error: ", e$message)
          NULL
        }
      )

      muest <- if (is.null(non_iv_fit)) NA_real_ else non_iv_fit$psi
      ci <- if (is.null(non_iv_fit)) {
        c(NA_real_, NA_real_)
      } else {
        non_iv_fit$confidence_interval
      }

      list(muest  = muest,
           ci     = ci,
           n_uc   = n_uc_sel,
           n_uc_folds = n_uc_folds,
           selection_diagnostics = selection_diagnostics,
           cor_Ac = cor_Ac,
           time_s = proc.time()["elapsed"] - t_sim_start)

    }, mc.cores = n_cores)

    # ---- aggregate results ----
    safe <- function(r, field, default)
      if (is.list(r)) r[[field]] else default

    muests    <- sapply(results_list, safe, "muest",  NA_real_)
    cis       <- do.call(rbind, lapply(results_list, safe, "ci",
                                       c(NA_real_, NA_real_)))
    n_uc_sels <- as.integer(sapply(results_list, safe, "n_uc", NA_integer_))
    n_uc_by_fold <- do.call(rbind, lapply(
      results_list,
      safe,
      field = "n_uc_folds",
      default = rep(NA_integer_, n_outer_folds)
    ))
    selection_diagnostics <- lapply(
      results_list,
      safe,
      field = "selection_diagnostics",
      default = NULL
    )
    cor_Ac_sels <- sapply(results_list, safe, "cor_Ac", NA_real_)
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

      # Preserve the legacy mean-selection file above and additionally save the
      # actual fold-specific selections required by the manuscript algorithm.
      for (fold_index in seq_len(ncol(n_uc_by_fold))) {
        filename_n_uc_fold <- paste0(
          results_dir, 'conf', confounding_mechanism, '_', option, '_', method,
          '_n_uc_fold', fold_index, '.csv'
        )
        df_n_uc_fold <- data.frame(n_uc = n_uc_by_fold[, fold_index])
        if (file.exists(filename_n_uc_fold)) {
          write.csv(
            cbind(read.csv(filename_n_uc_fold), df_n_uc_fold),
            filename_n_uc_fold,
            row.names = FALSE
          )
        } else {
          write.csv(df_n_uc_fold, filename_n_uc_fold, row.names = FALSE)
        }
      }

      if (isTRUE(iv_control$save_diagnostics) &&
          any(lengths(selection_diagnostics) > 0L)) {
        diagnostic_tag <- paste0(
          format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()
        )
        saveRDS(
          selection_diagnostics,
          file = paste0(
            results_dir, 'conf', confounding_mechanism, '_', option, '_', method,
            '_selection_diagnostics_', diagnostic_tag, '.rds'
          )
        )
      }
    }

    # Save per-sim correlation between estimated Ac and true Ac (IV-TPS and IV-GL only)
    if (any(!is.na(cor_Ac_sels))) {
      filename_cor_Ac <- paste0(results_dir, 'conf', confounding_mechanism,
                                '_', option, '_', method, '_cor_Ac.csv')
      df_cor_Ac <- data.frame(cor_Ac = cor_Ac_sels)
      if (file.exists(filename_cor_Ac)) {
        write.csv(cbind(read.csv(filename_cor_Ac), df_cor_Ac),
                  filename_cor_Ac, row.names = FALSE)
      } else {
        write.csv(df_cor_Ac, filename_cor_Ac, row.names = FALSE)
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
                            theta_B1uc = 0.025, theta_B2uc = 0.01,
                            theta_B1c  = 0.20, theta_B2c  = 0.10,
                            theta_u    = 0.30,
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
    Auc = 1.1 * B1uc + 0.7 * B2uc,
    Ac  = 1.0 * B1c  + 0.8 * B2c,
    U   = rho1 * B1c + rho2 * B2c + Zu
  )
}

# Two-confounder GP data-generating process (supplement Section 4.2, mechanism 4).
# Same 4 latent fields as mechanism 1; two confounders each linked to B1_c and B2_c.
#   U1 = rho1*B1_c + rho2*B2_c + Zu1  (Zu1 ~ GP(R(theta_u1)))
#   U2 = rho3*B1_c + rho4*B2_c + Zu2  (Zu2 ~ GP(R(theta_u2)))
compute_data_GP_2U <- function(nsims, distmat,
                               theta_B1uc = 0.025, theta_B2uc = 0.01,
                               theta_B1c  = 0.20, theta_B2c  = 0.10,
                               theta_u1   = 0.30, theta_u2   = 0.20,
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
    Auc = 1.1 * B1uc + 0.7 * B2uc,
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

# Coordinate-based nonlinear mechanism (supplement Section 4.2, mechanism 7).
# U = sin(2*pi*lat*long) + lat + long  (deterministic, normalized coordinates)
# B1_c = U, B2_c = 0;  B1_uc, B2_uc ~ independent GP fields
# Auc = 0.6*B1_uc + 0.4*B2_uc
# Ac  = 1.0*U  (= 1.0*B1_c + 0.8*0)
compute_data_spatialcoord <- function(lat, long, nsims, distmat) {
  n      <- length(lat)
  lat_n  <- (lat  - min(lat))  / (max(lat)  - min(lat))
  long_n <- (long - min(long)) / (max(long) - min(long))
  # Coordinate basis functions (fixed in space)
  f1 <- sin(pi * lat_n)
  f2 <- cos(pi * long_n)
  f3 <- lat_n * long_n
  # Random coefficients per simulation -> U varies across sims
  a1 <- rnorm(nsims); a2 <- rnorm(nsims); a3 <- rnorm(nsims)
  U_raw <- outer(f1, a1) + outer(f2, a2) + outer(f3, a3)  # n x nsims
  U <- sweep(U_raw, 2, colMeans(U_raw), "-")               # center each column
  phi_fn <- function(r, kappa = 2) r / (2 * sqrt(kappa))
  gp     <- function(K) {
    s <- MASS::mvrnorm(nsims, rep(0, n), K)
    if (nsims == 1L) matrix(s, n, 1L) else t(s)
  }
  B1uc  <- gp(geoR::matern(u = distmat, phi = phi_fn(0.10), kappa = 2))
  B2uc  <- gp(geoR::matern(u = distmat, phi = phi_fn(0.05), kappa = 2))
  list(
    Auc = 1.1 * B1uc + 0.7 * B2uc,
    Ac  = U,
    U   = U
  )
}

# Estimate the effective spatial range of the phi1 influence function for
# adaptive HAC bandwidth selection. Uses the unpadded phi values for the
# a > cutoff - delta subset to avoid the artificial autocorrelation created
# by zero-padding (which would bias the range estimate upward).
# Returns the first binned distance at which the empirical autocorrelation of
# phi1 drops below `target` * var(phi1), or `dmax` if still correlated throughout.
# Estimate the HAC bandwidth from the unpadded phi1 influence function values
# (for the a > cutoff - delta subsample only, to avoid zero-padding artifacts).
# Uses `target = 0`: finds the first distance bin where the empirical
# spatial autocorrelation of phi1 turns NEGATIVE. This distinguishes:
#   - Long-range positive processes (GP, mech 1/4): autocorrelation never goes
#     negative -> returns max pairwise distance -> large h -> captures positive
#     long-range correlations -> SE increases toward 0.95 from under-coverage.
#   - CAR / within-state / deterministic processes (mech 2/3/6): negative
#     autocorrelation appears at some intermediate distance -> h = that distance
#     -> includes the negative long-range contributions -> HAC decreases -> SE
#     decreases toward 0.95 from over-coverage.
hac_adaptive_cutoff <- function(phi1_raw, sub_idx, distmat,
                                 target = 0, n_bins = 20L, fallback = 0.05) {
  ok <- !is.na(phi1_raw)
  if (sum(ok) < 10L) return(fallback)

  phi1_use <- phi1_raw[ok]
  phi1_c   <- phi1_use - mean(phi1_use)
  var0     <- var(phi1_c)
  if (var0 < 1e-12) return(fallback)

  sub_ok <- sub_idx[ok]
  dsub   <- distmat[sub_ok, sub_ok]
  ut     <- upper.tri(dsub)
  dvec   <- dsub[ut]
  cprod  <- outer(phi1_c, phi1_c)[ut]

  pos_d  <- dvec[dvec > 0]
  if (length(pos_d) == 0L) return(fallback)
  dmax   <- max(pos_d)              # use full range so negative bins aren't missed
  breaks <- seq(0, dmax, length.out = n_bins + 1L)
  grp    <- findInterval(dvec, breaks, rightmost.closed = TRUE)

  acov  <- tapply(cprod, grp, mean)
  mid_d <- (breaks[-1] + breaks[-length(breaks)]) / 2

  below <- which(!is.na(acov) & acov < target)   # first bin with negative autocov
  if (length(below) == 0L) return(dmax)           # never negative -> use full range
  mid_d[min(below)]
}

asymptotic_variance_delta <- function(y, a, erfest, cutoff, delta, distmat,
                                      hac_cutoff = NULL){
  # Calculate parameters
  ix_cut <- which.min(abs(erfest$res$a.vals - cutoff))
  theta1 <- erfest$res$est[ix_cut]
  theta2 <- mean(a <= cutoff)
  theta3 <- mean(y[a <= cutoff])
  theta4 <- mean(y)

  # Calculate estimated influence functions
  # phi1 comes from ctseff run on the n_sub-observation subset (a > cutoff - delta).
  # Padding with zeros dilutes its variance by n_sub/n; rescale by n/n_sub so that
  # Sigma_HAC/n correctly estimates Var(hat_theta1) = Var(phi1_Kennedy)/n_sub.
  n <- length(a)
  n_sub <- sum(a > cutoff - delta)
  phi1 <- rep(0, n)
  phi1[a > cutoff - delta] <- erfest$phi[[ix_cut]] * (n / n_sub)
  phi2 <- 1*(a < cutoff) - mean(a < cutoff)
  phi3 <- rep(NA, length(a))
  phi3[a <= cutoff] <- (y[a <= cutoff] - mean(y[a <= cutoff])) / theta2
  phi3[a > cutoff] <- 0
  phi4 <- y - mean(y)

  # Adaptive HAC bandwidth: estimate the effective spatial range of phi1 from
  # data. Mechanisms with long-range confounding produce phi1 values correlated
  # over long distances -> larger cutoff -> larger SE. Short-range mechanisms
  # produce rapidly-decaying phi1 autocorrelation -> smaller cutoff -> smaller
  # SE. This prevents systematic over/under-coverage across mechanisms with
  # different spatial scales. Pass an explicit hac_cutoff to override.
  # if (is.null(hac_cutoff)) {
  #   sub_idx    <- which(a > cutoff - delta)
  #   phi1_raw   <- erfest$phi[[ix_cut]]   # unpadded, unscaled
  #   hac_cutoff <- hac_adaptive_cutoff(phi1_raw, sub_idx, distmat)
  # }

  # Conley (1999) spatial HAC covariance matrix.
  # Bartlett kernel: K(d/h) = max(0, 1 - d/h).
  # Sigma_HAC = (1/n) * t(phi_c) %*% W %*% phi_c  -->  se = sqrt(asym_var / n)
  phi_mat <- cbind(phi1, phi2, phi3, phi4)
  phi_c   <- sweep(phi_mat, 2, colMeans(phi_mat))
  # W       <- matrix(pmax(0, 1 - distmat / hac_cutoff), nrow(distmat), ncol(distmat))
  # Sigma   <- t(phi_c) %*% W %*% phi_c / n
  Sigma   <- cov(phi_c)

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

# Bivariate Leroux CAR mechanism (supplement Section 4.2, mechanism 5).
#   (U, B1_c) ~ bivariate Leroux with rho_c=0.8, cross-corr 0.7
#   B2_c ~ univariate Leroux(rho_c) independently
#   B1_uc, B2_uc ~ univariate Leroux(rho_uc=0.1) independently
#   Auc = 0.6*B1_uc + 0.4*B2_uc,  Ac = 1.0*B1_c + 0.8*B2_c
compute_data_leroux <- function(W, nsims) {
  n      <- nrow(W)
  stopifnot(ncol(W) == n)
  rho_c  <- 0.6
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
    out$Auc[, i] <- 1.1 * B1uc_i + 0.7 * B2uc_i
  }
  return(out)
}

# TPS-basis confounding mechanism (mechanism 1).
# Ac and U are generated from the n_c lowest-frequency TPS eigenvectors;
# Auc is generated from the remaining high-frequency TPS eigenvectors.
# Because B_tps_full is orthonormal, the IV-TPS projection B_c %*% t(B_c) %*% A
# recovers Ac exactly (Auc is orthogonal to Col(B_c) by construction).
# Coefficients are scaled by sqrt(n/n_c) so each component has average
# marginal variance ≈ 1, matching the GP mechanisms.
#   Ac  = B_c  %*% (sc_c  * a1)
#   Auc = B_uc %*% (sc_uc * au)
#   U   = B_c  %*% (sc_c  * (rho1*a1 + rho2*a2))
compute_data_TPS_basis <- function(B_tps_full, nsims,
                                   n_c  = NULL,
                                   rho1 = 0.8, rho2 = 0.6) {
  n   <- nrow(B_tps_full)
  m   <- ncol(B_tps_full)
  if (is.null(n_c)) n_c <- m - floor(0.9 * n)
  n_uc <- m - n_c
  B_c  <- B_tps_full[, seq_len(n_c),      drop = FALSE]
  B_uc <- B_tps_full[, seq(n_c + 1L, m), drop = FALSE]
  sc_c  <- sqrt(n / n_c)
  sc_uc <- sqrt(n / n_uc)
  Ac  <- matrix(NA_real_, n, nsims)
  Auc <- matrix(NA_real_, n, nsims)
  U   <- matrix(NA_real_, n, nsims)
  for (i in seq_len(nsims)) {
    a1 <- rnorm(n_c);  a2 <- rnorm(n_c);  au <- rnorm(n_uc)
    Ac[, i]  <- B_c  %*% (sc_c  * a1)
    Auc[, i] <- B_uc %*% (sc_uc * au)
    U[, i]   <- B_c  %*% (sc_c  * (rho1 * a1 + rho2 * a2))
  }
  list(Ac = Ac, Auc = Auc, U = U)
}

# GL-basis confounding mechanism (mechanism 2).
# Identical structure to mechanism 1 but using the Graph Laplacian eigenvectors
# instead of TPS eigenvectors.  IV-GraphLaplacian recovers Ac exactly.
compute_data_GL_basis <- function(B_gl_full, nsims,
                                  n_c  = NULL,
                                  rho1 = 0.8, rho2 = 0.6) {
  n   <- nrow(B_gl_full)
  m   <- ncol(B_gl_full)
  if (is.null(n_c)) n_c <- m - floor(0.9 * n)
  n_uc <- m - n_c
  B_c  <- B_gl_full[, seq_len(n_c),      drop = FALSE]
  B_uc <- B_gl_full[, seq(n_c + 1L, m), drop = FALSE]
  sc_c  <- sqrt(n / n_c)
  sc_uc <- sqrt(n / n_uc)
  Ac  <- matrix(NA_real_, n, nsims)
  Auc <- matrix(NA_real_, n, nsims)
  U   <- matrix(NA_real_, n, nsims)
  for (i in seq_len(nsims)) {
    a1 <- rnorm(n_c);  a2 <- rnorm(n_c);  au <- rnorm(n_uc)
    Ac[, i]  <- B_c  %*% (sc_c  * a1)
    Auc[, i] <- B_uc %*% (sc_uc * au)
    U[, i]   <- B_c  %*% (sc_c  * (rho1 * a1 + rho2 * a2))
  }
  list(Ac = Ac, Auc = Auc, U = U)
}

# Two-confounder TPS-basis mechanism (mechanism 3).
# Two independent confounded "fields" B1_c, B2_c are drawn from the low-frequency
# TPS subspace; Auc from the high-frequency TPS subspace.
# U1 and U2 are each correlated with both B1_c and B2_c via different weights,
# mirroring compute_data_GP_2U but using TPS projections instead of GP fields.
# IV-TPS recovers Ac exactly (Ac in Col(B_c), Auc orthogonal to Col(B_c)).
compute_data_TPS_2U <- function(B_tps_full, nsims,
                                n_c  = NULL,
                                rho1 = 0.8, rho2 = 0.6,
                                rho3 = 0.5, rho4 = 0.4) {
  n   <- nrow(B_tps_full)
  m   <- ncol(B_tps_full)
  if (is.null(n_c)) n_c <- m - floor(0.9 * n)
  n_uc <- m - n_c
  B_c  <- B_tps_full[, seq_len(n_c),      drop = FALSE]
  B_uc <- B_tps_full[, seq(n_c + 1L, m), drop = FALSE]
  sc_c  <- sqrt(n / n_c)
  sc_uc <- sqrt(n / n_uc)
  Ac  <- matrix(NA_real_, n, nsims)
  Auc <- matrix(NA_real_, n, nsims)
  U1  <- matrix(NA_real_, n, nsims)
  U2  <- matrix(NA_real_, n, nsims)
  for (i in seq_len(nsims)) {
    a1  <- rnorm(n_c);  a2  <- rnorm(n_c)
    au1 <- rnorm(n_c);  au2 <- rnorm(n_c)
    auc <- rnorm(n_uc)
    B1c <- B_c %*% (sc_c * a1)
    B2c <- B_c %*% (sc_c * a2)
    Ac[, i]  <- 1.0 * B1c + 0.8 * B2c
    Auc[, i] <- B_uc %*% (sc_uc * auc)
    U1[, i]  <- rho1 * B1c + rho2 * B2c + B_c %*% (sc_c * au1)
    U2[, i]  <- rho3 * B1c + rho4 * B2c + B_c %*% (sc_c * au2)
  }
  list(Ac = Ac, Auc = Auc, U1 = U1, U2 = U2)
}

# Two-confounder GL-basis mechanism (mechanism 6).
# Identical structure to mechanism 3 but using Graph Laplacian eigenvectors.
# IV-GraphLaplacian recovers Ac exactly.
compute_data_GL_2U <- function(B_gl_full, nsims,
                               n_c  = NULL,
                               rho1 = 0.8, rho2 = 0.6,
                               rho3 = 0.5, rho4 = 0.4) {
  n   <- nrow(B_gl_full)
  m   <- ncol(B_gl_full)
  if (is.null(n_c)) n_c <- m - floor(0.9 * n)
  n_uc <- m - n_c
  B_c  <- B_gl_full[, seq_len(n_c),      drop = FALSE]
  B_uc <- B_gl_full[, seq(n_c + 1L, m), drop = FALSE]
  sc_c  <- sqrt(n / n_c)
  sc_uc <- sqrt(n / n_uc)
  Ac  <- matrix(NA_real_, n, nsims)
  Auc <- matrix(NA_real_, n, nsims)
  U1  <- matrix(NA_real_, n, nsims)
  U2  <- matrix(NA_real_, n, nsims)
  for (i in seq_len(nsims)) {
    a1  <- rnorm(n_c);  a2  <- rnorm(n_c)
    au1 <- rnorm(n_c);  au2 <- rnorm(n_c)
    auc <- rnorm(n_uc)
    B1c <- B_c %*% (sc_c * a1)
    B2c <- B_c %*% (sc_c * a2)
    Ac[, i]  <- 1.0 * B1c + 0.8 * B2c
    Auc[, i] <- B_uc %*% (sc_uc * auc)
    U1[, i]  <- rho1 * B1c + rho2 * B2c + B_c %*% (sc_c * au1)
    U2[, i]  <- rho3 * B1c + rho4 * B2c + B_c %*% (sc_c * au2)
  }
  list(Ac = Ac, Auc = Auc, U1 = U1, U2 = U2)
}

# Reversed-scale TPS mechanism (mechanism 4).
# Same confounded/unconfounded proportions as compute_data_TPS_basis (a small
# n_c-sized confounded block, ~10% of the basis), but Ac and U now live in the
# HIGH-frequency (small-scale, last n_c columns) subspace instead of the
# low-frequency (first n_c columns) subspace, and Auc lives in the remaining
# low-frequency columns. This matches the instrument_scale = "large" partition
# used at estimation time (partition_basis_columns(): confounded = last n_uc
# columns' complement, i.e. the high-frequency tail), so IV-TPS with the
# reversed-scale flip applied should recover Ac about as well as it recovers
# Ac under compute_data_TPS_basis for the standard (non-reversed) mechanism.
compute_data_TPS_reversed <- function(B_tps_full, nsims,
                                      n_c  = NULL,
                                      rho1 = 0.8, rho2 = 0.6) {
  n    <- nrow(B_tps_full)
  m    <- ncol(B_tps_full)
  if (is.null(n_c)) n_c <- m - floor(0.9 * n)  # small confounded count
  n_uc <- m - n_c                             # large instrument count
  B_c  <- B_tps_full[, seq(n_uc + 1L, m), drop = FALSE]  # last n_c cols (high-freq): confounded
  B_uc <- B_tps_full[, seq_len(n_uc),     drop = FALSE]  # first n_uc cols (low-freq): instruments
  sc_c  <- sqrt(n / n_c)
  sc_uc <- sqrt(n / n_uc)
  Ac  <- matrix(NA_real_, n, nsims)
  Auc <- matrix(NA_real_, n, nsims)
  U   <- matrix(NA_real_, n, nsims)
  for (i in seq_len(nsims)) {
    a1 <- rnorm(n_c);  a2 <- rnorm(n_c);  au <- rnorm(n_uc)
    Ac[, i]  <- B_c  %*% (sc_c  * a1)
    Auc[, i] <- B_uc %*% (sc_uc * au)
    U[, i]   <- B_c  %*% (sc_c  * (rho1 * a1 + rho2 * a2))
  }
  list(Ac = Ac, Auc = Auc, U = U)
}

# Reversed-scale GL mechanism (mechanism 8).
# Identical structure to compute_data_TPS_reversed but using Graph Laplacian
# eigenvectors: a small, high-frequency confounded block (last n_c columns)
# matching the instrument_scale = "large" partition at estimation time.
compute_data_GL_reversed <- function(B_gl_full, nsims,
                                     n_c  = NULL,
                                     rho1 = 0.8, rho2 = 0.6) {
  n    <- nrow(B_gl_full)
  m    <- ncol(B_gl_full)
  if (is.null(n_c)) n_c <- m - floor(0.9 * n)
  n_uc <- m - n_c
  B_c  <- B_gl_full[, seq(n_uc + 1L, m), drop = FALSE]
  B_uc <- B_gl_full[, seq_len(n_uc),     drop = FALSE]
  sc_c  <- sqrt(n / n_c)
  sc_uc <- sqrt(n / n_uc)
  Ac  <- matrix(NA_real_, n, nsims)
  Auc <- matrix(NA_real_, n, nsims)
  U   <- matrix(NA_real_, n, nsims)
  for (i in seq_len(nsims)) {
    a1 <- rnorm(n_c);  a2 <- rnorm(n_c);  au <- rnorm(n_uc)
    Ac[, i]  <- B_c  %*% (sc_c  * a1)
    Auc[, i] <- B_uc %*% (sc_uc * au)
    U[, i]   <- B_c  %*% (sc_c  * (rho1 * a1 + rho2 * a2))
  }
  list(Ac = Ac, Auc = Auc, U = U)
}
