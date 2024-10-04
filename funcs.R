# npcausal ctseff (edited so can deal with 1 covariate)
ctseff <- function(y, a, x, bw.seq, n.pts = 100, a.rng = c(min(a), max(a)),
                   sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean", "SL.ranger")) {
  require("SuperLearner")
  require("earth")
  require("gam")
  require("ranger")
  require(KernSmooth)
  kern <- function(t) {
    dnorm(t)
  }
  
  n <- dim(x)[1]
  
  # set up evaluation points & matrices for predictions
  a.min <- a.rng[1]
  a.max <- a.rng[2]
  a.vals <- seq(a.min, a.max, length.out = n.pts)
  xa.new <- rbind(cbind(x, a), cbind(x[rep(1:n, length(a.vals)), ], a = rep(a.vals, rep(n, length(a.vals)))))
  colnames(xa.new) <- c(colnames(x), "a") # sophie's change.
  x.new <- xa.new[, -dim(xa.new)[2]]
  x <- data.frame(x)
  x.new <- data.frame(x.new)
  colnames(x.new) <- colnames(x) # sophie's change
  xa.new <- data.frame(xa.new)
  
  # estimate nuisance functions via super learner
  # note: other methods could be used here instead
  pimod <- SuperLearner(Y = a, X = data.frame(x), SL.library = sl.lib, newX = x.new)
  pimod.vals <- pimod$SL.predict
  pi2mod <- SuperLearner(Y = log((a - pimod.vals[1:n])^2), X = x, SL.library = sl.lib, newX = x.new)
  pi2mod.vals <- exp(pi2mod$SL.predict) # sophie's change
  mumod <- SuperLearner(Y = y, X = cbind(x, a), SL.library = sl.lib, newX = xa.new)
  muhat.vals <- mumod$SL.predict
  
  # construct estimated pi/varpi and mu/m values
  a.std <- (xa.new$a - pimod.vals) / sqrt(pi2mod.vals)
  pihat.vals <- approx(density(a.std[1:n],from = min(a.std), to =max(a.std))$x, 
                       density(a.std[1:n],from = min(a.std), to =max(a.std))$y, 
                       xout = a.std)$y / sqrt(pi2mod.vals) # sophie's change
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  varpihat <- predict(smooth.spline(a.vals, apply(pihat.mat, 2, mean)), x = a)$y
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = T, nrow = n)
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  mhat <- predict(smooth.spline(a.vals, apply(muhat.mat, 2, mean)), x = a)$y
  mhat.mat <- matrix(rep(apply(muhat.mat, 2, mean), n), byrow = T, nrow = n)
  
  
  # form adjusted/pseudo outcome xi
  pseudo.out <- (y - muhat) / (pihat / varpihat) + mhat
  
  # leave-one-out cross-validation to select bandwidth
  w.fn <- function(bw, a.vals) { # sophie's change
    w.avals <- NULL
    for (a.val in a.vals) {
      a.std <- (a - a.val) / bw
      kern.std <- kern(a.std) / bw
      w.avals <- c(w.avals, mean(a.std^2 * kern.std) * (kern(0) / bw) /
                     (mean(kern.std) * mean(a.std^2 * kern.std) - mean(a.std * kern.std)^2))
    }
    return(w.avals / n)
  }
  hatvals <- function(bw) {
    #approx(a.vals, w.fn(bw), xout = a[a > a.min & a < a.max])$y # sophie's change
    w.fn(bw, a.vals = a) # sophie's change
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
  
  # alternative approach:
  # h.opt <- optimize(function(h){ hats <- hatvals(h); mean( ((pseudo.out[a > a.min & a < a.max]-cts.eff.fn(pseudo.out,bw=h))/(1-hats))^2) } ,
  #  bw.seq, tol=0.01)$minimum
  
  # estimate effect curve with optimal bandwidth
  est <- approx(locpoly(a, pseudo.out, bandwidth = h.opt), xout = a.vals)$y
  
  # estimate pointwise confidence band
  # note: other methods could also be used
  se <- NULL
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
    sigma <- cov(t(solve(Dh) %*%
                     rbind(
                       kern.std * (pseudo.out - beta[1] - beta[2] * a.std) + int1,
                       a.std * kern.std * (pseudo.out - beta[1] - beta[2] * a.std) + int2
                     )))
    se <- c(se, sqrt(sigma[1, 1]))
  }
  
  ci.ll <- est - 1.96 * se / sqrt(n)
  ci.ul <- est + 1.96 * se / sqrt(n)
  res <- data.frame(a.vals, est, se, ci.ll, ci.ul)
  
  return(invisible(list(res = res, bw.risk = bw.risk)))
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
      Ys[,i] <- rnorm(n, 2*Us[,i] + As[,i] - 0.2*As[,i]*Us[,i], 1)
    }
  }
  # nonlinear outcome model
  if (option == 'nonlinear'){
    for (i in 1:nreps){
      Ys[,i] <- rnorm(n, 2*Us[,i] + As[,i] - 0.2*As[,i]*Us[,i] - 0.1*As[,i]^2 + 0.05*As[,i]^2*Us[,i]
                     , 1)
    }
  } 
  return(Ys)
}

# Function to compute true ERF
computemutrue <- function(a.vals, 
                         option = c('linear', 'nonlinear')){
  # a.vals is a vector of exposure values
  # option is a string indicating the form of the outcome model (see createY)
  # returns a vector of true ERF
  
  n <- length(a.vals)
  option <- match.arg(option)

  mutrue <- rep(NA, n)
  
  for (i in 1:n){
    aval <- a.vals[i]
    if (option == 'linear'){
      mutrue[i] <- 2*0.3 + aval - 0.2*aval*0.3
    }
    if (option == 'nonlinear'){
      mutrue[i] <- 2*0.3 + aval - 0.2*aval*0.3 - 0.1*aval^2 + 0.05*aval^2*0.3
    }
  }
  return(mutrue)
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

# Function that simulates data, estimates ERF using different methods, and saves results to csvs
simfunc <- function(nsims,
                   lat,
                   lon,
                   rangeu = c('tinyscale', 'smallscale'),
                   option = c('linear', 'nonlinear'),
                   methods = c(
                     'baseline',
                     'spatialcoord',
                     'IV-TPS',
                     'IV-GraphLaplacian'
                   ),
                   GFT_conf,
                   statemat,
                   within_state_GP = F) {
  # nsims is the number of simulations
  # lat is a vector of latitudes
  # lon is a vector of longitudes
  # rangeu is the scale of the unconfounded component of exposure
  # option is the form of the outcome model
  # methods are the methods used to predict ERF
  # GFT_conf are the matrix of eigenvectors of the Graph Fourier to adjust for
  # statemat is the matrix of state-level indicators
  # within_state_GP is a boolean indicating whether we're in confounding mech 3 or not

  # writes estimated ERFs to a csv file named filename
  
  rangeu <- match.arg(rangeu)
  option <- match.arg(option)
  n <- length(lat)
  if (rangeu == 'tinyscale'){
    rangeu <- 0.05
  }
  if (rangeu == 'smallscale') {
    rangeu <- 0.1
  }
  
  ################# GENERATE DATA #################
  
  # Compute distance matrix
  distmat <- geosphere::distm(cbind(lon, lat), 
                  fun = distHaversine)
  distmat <- distmat/1000000 # scale so range (0,2)
  # confounding mech 1-2
  if (!within_state_GP){ 
    # Compute variance of GP
    Sigma_GP <- compute_Sigma_GP(distmat = distmat,
                                rangeu = rangeu, 
                                rangec = 0.5)
    # Simulate nsims of data according to GP
    dat <- compute_data_GP(n = nsims, Sigma_GP = Sigma_GP)
  }
  # confounding mech 3
  else{  
    # Simulate data as GPs within each state
    dat <- compute_data_GP_state(distmat = distmat,
                                rangeu = rangeu, 
                                rangec = 0.5,
                                n = nsims,
                                statemat = statemat)
  }
  
  Ac <- dat$Ac 
  Auc <- dat$Auc
  U <- dat$U
  A <- Ac + Auc # all have dimension n x nsims
  avals <- seq(-2, 2, length.out = 100)
  Y <- createY(Us=U, As=A, option = option)
  
  ################# FIT MODELS #################
  
  for (method in methods){
    # Create filename for csvs containing estimated erfs
    if (!within_state_GP){
      filename <- paste0('results_Oct1/', rangeu, '_', option, '_', method, '.csv')
    }
    else{
      filename <- paste0('results_Oct1/within_state/', rangeu, '_', option, '_', method, '.csv')
    }
    
    # Create storage for estimated erfs
    muests <- matrix(NA, nrow = length(avals), ncol = nsims)
    
    for (sim in 1:nsims){
      print(c(method, sim))
      
      # Create xmat, the confounding adjustment.
      if (method == 'baseline'){
        xmat <- matrix(rep(1,n), ncol = 1)
        colnames(xmat) <- 'Intercept'
      }
      
      if (method == 'spatialcoord'){
        xmat <- cbind(lat, lon)
        colnames(xmat) <- c('Latitude', 'Longitude')
      }
      if (method == 'IV-TPS'){
        mod <- mgcv::gam(A[,sim] ~ s(lat,lon,k=floor(0.2*n),fx=T)) # unpenalized
        xmat <- matrix(predict(mod), ncol = 1)
        colnames(xmat) <- 'Ac-TPS'
      }
      if (method == 'IV-GraphLaplacian'){
        mod <- lm(A[,sim] ~ GFT_conf)
        xmat <- matrix(predict(mod), ncol = 1)
        colnames(xmat) <- 'Ac-GraphLaplacian'
      }
      
      # Fit the ERF adjusting for xmat.
      muests[, sim] <- tryCatch({
        erfest <- ctseff(
          y = Y[, sim],
          a = A[, sim],
          x = xmat,
          a.rng = c(min(avals), max(avals)),
          n.pts = length(avals),
          bw.seq = seq(0.2, 2, length.out = 100)
        )
        erfest$res$est
      }, error = function(e) {
        message("Error encountered: ", e$message)
        NA  # Set muests[,sim] to NA if an error occurs
      })
    } # (There shouldn't be errors but in case)
    
    # Create dataframe whose first column is a.vals and the rest of cols are muests
    df <- cbind(avals, muests)
    
    # write results to file
    # Check if file exists
    if (file.exists(filename)){
      # write new sims to file as new columns
      olddf <- read.csv(filename)
      newdf <- cbind(olddf, muests)
      write.csv(newdf, filename, row.names = FALSE)
    }
    # if file for ERF ests does not exist create it and write results
    else{
      write.csv(df, filename, row.names = FALSE)
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
  Sigma <- matrix(NA, nrow = 3*n, ncol = 3*n)
  Sigma[1:n, 1:n] <- sigu^2*Sigmau
  Sigma[(n+1):(2*n), (n+1):(2*n)] <- sigc^2*Sigmac
  Sigma[(2*n+1):(3*n), (2*n+1):(3*n)] <- sigz^2*Sigmac
  # Auc is uncorrelated + indep of Ac and U
  Sigma[1:n, (n+1):(3*n)] <- 0
  Sigma[(n+1):(3*n), 1:n] <- 0
  # Ac and U are highly dependent
  Sigma[(n+1):(2*n), (2*n+1):(3*n)] <- rho*sigc*sigz*Sigmac
  Sigma[(2*n+1):(3*n), (n+1):(2*n)] <- rho*sigc*sigz*Sigmac
  return(Sigma)
}

# Function that computes the data from the GP given the covariance matrix
compute_data_GP <- function(n, 
                           Sigma_GP,
                           mu = c(rep(0.1, nrow(Sigma_GP)/3), 
                                  rep(-0.2, nrow(Sigma_GP)/3), 
                                  rep(0.3, nrow(Sigma_GP)/3))
                           ){
  # n is the number of observations
  # Sigma_GP is the covariance matrix
  # mu is the mean vector
  # returns a list with the data Auc,Ac,U
  
  stopifnot(nrow(Sigma_GP) %% 3 == 0)
  dat <- matrix(MASS::mvrnorm(n=n, mu = mu, Sigma=Sigma_GP), 
               nrow = nrow(Sigma_GP), ncol = n, 
               byrow = TRUE)
  k <- nrow(Sigma_GP)/3
  return(list('Auc' = dat[1:k,],
              'Ac' = dat[(k+1):(2*k),],
              'U' = dat[(2*k+1):(3*k),]))
}

# Function that computes the data for confounding mechanism 3
compute_data_GP_state <- function(distmat,
                                  rangeu,
                                  rangec,
                                  n,
                                  statemat){

  # distmat is the distance matrix
  # kappa is the smoothness parameter
  # rangeu is the range of the GP for the unconfounded part of exposure
  # rangec is the range of the GP for the confounded part of exposure
  # n is the number of simulations NOT sample size
  # statemat is the indicator matrix for states (sample size x number of states)
  # returns a list with the data Auc,Ac,U
  
  c <- ncol(statemat)
  out <- list('Auc' = matrix(NA, nrow = nrow(distmat), ncol = n),
             'Ac' = matrix(NA, nrow = nrow(distmat), ncol = n),
             'U' = matrix(NA, nrow = nrow(distmat), ncol = n))
  
  # Loop through states
  for (i in 1:c){
    ixs <- which(statemat[,i] == 1)
    distmat_st <- distmat[ixs, ixs]
    # Compute GP covariance and data within each state
    Sigma <- compute_Sigma_GP(distmat = distmat_st,
                             rangeu=rangeu,
                             rangec=rangec)
    date_state <- compute_data_GP(n = n, Sigma_GP = Sigma)
    # assign to dat
    out$Auc[ixs,] <- date_state$Auc
    out$Ac[ixs,] <- date_state$Ac
    out$U[ixs,] <- date_state$U
  }
  return(out)
}

