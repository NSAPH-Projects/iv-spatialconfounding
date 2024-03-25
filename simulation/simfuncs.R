# Function used to calculate nested decomposition
nested_decomp_mats = function(groups, append_identity = TRUE) {
  if (!is.matrix(groups)) {
    groups = as.matrix(groups, ncol = 1)
  }
  if (append_identity) {
    groups = cbind(groups, 1:nrow(groups))
  }
  n = nrow(groups)
  L = ncol(groups)
  
  proj_mats = list()
  decomp_mats = list()
  for (l in 1:L) {
    P_l = matrix(0, n, n)
    for (group in unique(groups[, l])) {
      ix = which(groups[, l] == group)
      P_l[ix, ix] = 1 / length(ix)
    }
    proj_mats[[l]] = P_l
    if (l > 1) {
      decomp_mats[[l]] = P_l - proj_mats[[l - 1]]
    } else {
      decomp_mats[[l]] = P_l
    }
  }
  return(list(
    proj_mats = proj_mats,
    decomp_mats = decomp_mats
  ))
}

# Function to create unmeasured confounder
createU = function(latnorm, longnorm){
  # latnorm and longnorm are vectors of normalized latitude and longitude
  # returns a vector of unmeasured confounder
  n = length(latnorm)
  U = cos(2*pi*latnorm*longnorm) + 2*latnorm - longnorm + rnorm(n, mean = 0, sd = 0.5)
  return(as.numeric(U))
}

# Function to create measured confounder
createX = function(U, latnorm, longnorm){
  # U is a vector of unmeasured confounder
  # latnorm and longnorm are vectors of normalized latitude and longitude
  # returns a vector of measured confounder
  n = length(U)
  X = rnorm(n, 1 + 0.1*U + sin(pi*latnorm), 1) # correlation with U? 
  return(as.numeric(X))
}

# Function to create exposure
createA = function(U, X, projmat){
  # U is a vector of unmeasured confounder
  # X is a vector of measured confounder
  # projmat is the projection matrix onto the subspace of confounding
  # returns a list of Ac, Auc, and A
  n = length(U)
  Ac = projmat %*% rnorm(n, X*0.25*U^2, 1)
  Auc = (diag(n) - projmat) %*% rnorm(n, X, 1)
  A = Ac + Auc
  return(list(Ac = as.numeric(Ac), Auc = as.numeric(Auc), A = as.numeric(A)))
}

# Function to create outcome
createY = function(U, X, A, option = c('linear', 'interaction', 'nonlinear')){
  # U is a vector of unmeasured confounder
  # X is a vector of measured confounder
  # A is a vector of exposure
  # option is a string indicating the form of the outcome model
  # returns a vector of outcome
  option = match.arg(option)
  n = length(U)
  if (option == 'linear'){
    Y = rnorm(n, 2*U + A + X , 1)
  }
  if (option == 'interaction'){
    Y = rnorm(n, 2 + A + U - X - 3*A*U + A*X - 0.5*U*X, 1)
  }
  if (option == 'nonlinear'){
    Y = rnorm(n, X + log(1+A^2+U^2) + 0.5*X*A*cos(A/(1+U^2)), 1)
  }
  return(as.numeric(Y))
}

computemutrue = function(a.vals, U, latnorm, longnorm, option = c('linear', 'interaction', 'nonlinear')){
  # a.vals is a vector of exposure values
  # U is a vector of unmeasured confounder
  # latnorm and longnorm are vectors of normalized latitude and longitude
  # option is a string indicating the form of the outcome model (see createY)
  # returns a vector of true ERF
  n = length(a.vals)
  option = match.arg(option)
  # To calculate true curve, write all expectations in terms of U and then approximate using observed dist of U
  mutrue = rep(NA, n)
  for (i in 1:n){
    aval = a.vals[i]
    if (option == 'linear'){
      mutrue[i] = 2*mean(U) + aval + (1 + 0.1*mean(U) + mean(sin(pi*latnorm)))
    }
    if (option == 'interaction'){
      mutrue[i] = 2 + aval + mean(U) - (1 + 0.1*mean(U) + mean(sin(pi*latnorm))) - 3*aval*mean(U) + 
        aval*(1 + 0.1*mean(U) + mean(sin(pi*latnorm))) - 0.5*mean(U*(1 + 0.1*U + sin(pi*latnorm)))
    }
    if (option == 'nonlinear'){
      mutrue[i] = (1 + 0.1*mean(U) + mean(sin(pi*latnorm))) + mean(log(1 + aval^2 +U^2)) + 
        0.5*aval*mean((1 + 0.1*U + sin(pi*latnorm))*cos(aval/(1+U^2)))
    }
  }
  return(as.numeric(mutrue))
}

plotfunc = function(df, boundaries, names, labels=NULL){
  # df is a sf dataframe including columns names
  # boundaries is a sf dataframe with the boundaries of the region
  # names is a vector of strings with the names of the columns to be plotted
  # labels is a vector of strings to title the ggplots
  
  # returns a list of ggplots
  
  if (is.null(labels)){
    labels = names
  }
  K = length(names)
  gs = list()
  for (k in 1:K){
    # extract the column with names[k] from df
    var = df[[names[k]]]
    qs = round(quantile(var, probs = c(0.1, 0.3, 0.5, 0.7, 0.9)),2)
    # Turn qs into a vector of strings
    qschar = as.character(qs)
    gs[[k]] = ggplot(df) +
      xlim(-125, -65) + 
      ylim(25, 50) +
      geom_sf(aes_string(fill = names[k]), color=NA, size = 0.005) +
      geom_sf(data = boundaries, fill = NA, color = "black", size = 3) +
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
            panel.grid.major = element_line(colour = "transparent")) + 
      ggtitle(labels[k])
  }
  return(gs)
}

# Estimate density of A using kernel density estimation, but could also do this analytically probably
densityA = function(A, a.vals){
  # A is a vector of observed exposure
  # a.vals is a vector of exposure values for which we want to estimate density
  # returns a vector of density estimates
  out = density(A)
  dens = approx(out$x, out$y, xout = a.vals)
  return(dens)
}

# Function to calculate average absolute bias, avg RMSE, avg coverage
metrics = function(a.vals, muests, mutrue, A, cils=NULL, cius=NULL){
  # a.vals is a vector of exposure values
  # muests is a matrix of estimated ERFs, columns correspond to diff sims
  # mutrue is a vector of true ERF
  stopifnot(length(a.vals) == nrow(muests))
  stopifnot(length(a.vals) == nrow(mutrue))
  adens = densityA(A, a.vals)
  avgabsbias = sum(rowMeans(abs(muests - mutrue), na.rm = T) * adens$y, na.rm = T)/sum(adens$y, na.rm = T)
  avgRMSE = sum(sqrt(rowMeans((muests - mutrue)^2, na.rm = T))*adens$y, na.rm = T)/sum(adens$y, na.rm = T)
  if (is.null(cils) | is.null(cius)){
    return(list(avgabsbias = avgabsbias, avgRMSE = avgRMSE))
  }
  avgcoverage = sum(rowMeans((cils < mutrue) & (mutrue < cius), na.rm = T)*adens$y, na.rm = T)/sum(adens$y, na.rm = T)
  return(list(avgabsbias = avgabsbias, avgRMSE = avgRMSE, avgcoverage = avgcoverage))
}

spatialplus = function(y,
                       a,
                       x,
                       latnorm,
                       longnorm,
                       a.vals,
                       k_sp=100,
                       method = 'GCV.Cp',
                       model_fx = T) {
  xmod = mgcv::gam(a ~ x + s(latnorm,longnorm,k=k_sp,fx=model_fx),method=method)
  f_a_hat = xmod$fitted.values 
  r_a = a - f_a_hat
  mod=mgcv::gam(y ~ r_a + x + s(latnorm,longnorm,k=k_sp, fx=model_fx),method=method)
  erf = rep(NA, length(a.vals))
  # lconf = rep(NA, length(a.vals))
  # uconf = rep(NA, length(a.vals))
  
  # Predict f_a_hat for each a.val in a.vals
  for (i in 1:length(a.vals)) {
    df = data.frame(
      a = a.vals[i],
      x = x,
      latnorm = latnorm,
      longnorm = longnorm
    )
    fahat = mean(predict(
      xmod,
      newdata = df
    ))
    df$r_a = a.vals[i] - fahat
    preds = predict(
      mod,
      newdata = df
    )
    erf[i] = mean(preds)
    # lconf[i] = erf[i] - 1.96 * sd(preds)/sqrt(length(preds)) # first 
    # uconf[i] = erf[i] + 1.96 * sd(preds)/sqrt(length(preds))
  }

  # Predict ERF for each ra.val in ra.vals
  return(list('erf' = erf))#, 'lconf'= lconf, 'uconf' = uconf))
}

bobb_exposurepenalized = function(y,
                                  a,
                                  x,
                                  latnorm,
                                  longnorm,
                                  a.vals,
                                  k_sp=100) {
  mod0 = mgcv::gam(a ~ x + s(latnorm, longnorm, k = k_sp))
  mod = mgcv::gam(y ~ a + x + s(latnorm, longnorm, k = k_sp), sp = mod0$sp)
  erf = rep(NA, length(a.vals))
  # lconf = rep(NA, length(a.vals))
  # uconf = rep(NA, length(a.vals))
  for (i in 1:length(a.vals)) {
    df = data.frame(
      a = a.vals[i],
      x = x,
      latnorm = latnorm,
      longnorm = longnorm
    )
    preds = predict(
      mod,
      newdata = df
    )
    erf[i] = mean(preds)
    # lconf[i] = erf[i] - 1.96 * sd(preds)/sqrt(length(preds)) # first 
    # uconf[i] = erf[i] + 1.96 * sd(preds)/sqrt(length(preds))
  }
  return(list('erf' = erf))#, 'lconf'= lconf, 'uconf' = uconf))
}

spectral_discrete = function(y,
                             a,
                             x,
                             a.vals,
                             adj,
                             Ls = c(10)) {
  fits = list()
  dics = list()
  for (i in 1:length(Ls)) {
    fit = semipar.eCAR.Leroux(y, x=a, W=adj,
                              E=NULL,
                              C=as.matrix(x),
                              names.covariates='x',
                              L=Ls[i],
                              verbose = F)
    fits[[i]] = fit
    dics[[i]] = fit$DIC
  }
  # Return the fit with the lowest DIC
  fit = fits[[which.min(dics)]]
  omegas = fit$beta_omega[,"omega"]
  betas = fit$beta_omega[,"beta.mn"]
  #cils = fit$beta_omega[,"beta.q025"]
  #cius = fit$beta_omega[,"beta.q975"]
  
  # Extract beta of largest omega
  beta = betas[length(betas)]
  coeffs = fit$regrcoef[1:2,1] # coeffs of covariate
  
  n = length(y)
  erf = rep(NA, length(a.vals))
  # lconf = rep(NA, length(a.vals))
  # uconf = rep(NA, length(a.vals))
  for (i in 1:length(a.vals)) {
    preds = cbind(1,x) %*% coeffs + beta*rep(a.vals[i],n)
    erf[i] = mean(preds)
    # lconf[i] = erf[i] - 1.96 * sd(preds)/sqrt(length(preds))
    # uconf[i] = erf[i] + 1.96 * sd(preds)/sqrt(length(preds))
  }
  return(list('erf' = erf))#, 'lconf'= lconf, 'uconf' = uconf))
}

keller_szpiro_selectingscale = function(y,
                                        a,
                                        x,
                                        a.vals,
                                        latnorm,
                                        longnorm,
                                        fitmethod = c('outcome', 'preadjustment'),
                                        selection = c('AIC-NE', 'BIC-NE'),
                                        ks = seq(2,100,by=5)) {
  fitmethod = match.arg(fitmethod)
  selection = match.arg(selection)
  # Select the degrees of freedom in thin plate regression spline of no-exposure model
  
  criteria = rep(NA, length(ks))

  for (i in 1:length(ks)){
    k = ks[i]
    fit = mgcv::gam(y ~ x + s(latnorm, longnorm, bs = 'tp', k = k))
    if (selection == 'AIC-NE'){
      criteria[i] = AIC(fit)
    }
    if (selection == 'BIC-NE'){
      criteria[i] = BIC(fit)
    }
  }
  mink = ks[which.min(criteria)]
  erf = rep(NA, length(a.vals))
  # lconf = rep(NA, length(a.vals))
  # uconf = rep(NA, length(a.vals))
  if (fitmethod == 'outcome'){
    fit = mgcv::gam(y ~ a + x + s(latnorm, longnorm, bs = 'tp', k = ks[which.min(criteria)]))
    for (i in 1:length(a.vals)) {
      df = data.frame(
        a = a.vals[i],
        x = x,
        latnorm = latnorm,
        longnorm = longnorm
      )
      preds = predict(
        fit,
        newdata = df
      )
      erf[i] = mean(preds)
      # lconf[i] = erf[i] - 1.96 * sd(preds)/sqrt(length(preds)) # first 
      # uconf[i] = erf[i] + 1.96 * sd(preds)/sqrt(length(preds))
    }
  }
  if (fitmethod == 'preadjustment'){
    spl = mgcv::smoothCon(mgcv::s(latnorm, longnorm, k = mink), 
                          data = cbind.data.frame('latnorm' = latnorm, 'longnorm' = longnorm), knots = NULL)[[1]]
    mat = mgcv::PredictMat(spl,cbind.data.frame('latnorm' = latnorm, 'longnorm' = longnorm))
    Hmmat = mat %*% solve(t(mat) %*% mat) %*% t(mat)
    a1 = Hmmat %*% a
    a2 = a - a1
    fit = lm(y ~ x + a1 + a2)
    n = length(y)
    for (i in 1:length(a.vals)) {
      preds = cbind(1,x) %*% coefficients(fit)[1:2] + beta*rep(a.vals[i],n)
      erf[i] = mean(preds)
      # uconf[i] = erf[i] + 1.96 * sd(preds)/sqrt(length(preds))
      # lconf[i] = erf[i] - 1.96 * sd(preds)/sqrt(length(preds))
    }
  }
  return(list('erf' = erf))#, 'lconf' = lconf, 'uconf' = uconf))
}

# Non oracle function
# nonoracle = function(nsims,
#                      a.vals,
#                      latnorm,
#                      longnorm,
#                      projmats,
#                      binsize){
#   # For each matrix in projmats,
#   # for each binsize in binsizes,
#   # we'll fit ERF using one of the methods (KennedyERC or GPCERF)
#   
# }

# Function that runs a simulation for a given method and outcome model
simfunc = function(nsims,
                   a.vals,
                   latnorm,
                   longnorm,
                   projmat,
                   projmatname,
                   option = c('linear', 'interaction', 'nonlinear'),
                   method = c(
                     'KennedyERC',
                     'spatialplus',
                     'bobb_exposurepenalized',
                     'spectral_discrete',
                     'keller_szpiro_selectingscale_outcome',
                     'keller_szpiro_selectingscale_preadjustment',
                     'unadjustedOLS',
                     'GPCERF',
                     'GPCERF_nn'#,
                     #'nonoracle'
                   ),
                   filename = NULL, 
                   adjmat = NULL) {
  # nsims is the number of simulations
  # a.vals is the vector of a exposure values to predict ERF for
  # projmat is the projection matrix of confounding
  # latnorm and longnorm are the normalized vectors of lat and long values
  # option is the form of the outcome model
  # method is the method used to predict ERF
  # filename is the name of the file to save the results of estimated ERF
  # adjmat is the adjacency matrix of the spatial locations
  
  # writes estimated ERFs to a csv file named filename
  
  option = match.arg(option)
  method = match.arg(method)

  # If filename is null, set filename to projmat+option+method.csv
  if (is.null(filename)) {
    filename = paste0('results/', projmatname, '_', option, '_', method, '.csv')
  }
  muests = matrix(NA, nrow = length(a.vals), ncol = nsims)

  # Estimate true ERF
  U = createU(latnorm, longnorm)
  mutrue = computemutrue(a.vals, U, latnorm, longnorm, 
                         option = option)
  for (sim in 1:nsims){
    # for every 200 simulations, print filename and simulation number as a vector
    if (sim %% 200 == 0) {
      print(c(filename, sim))
    }
    u = createU(latnorm, longnorm)
    x = createX(u, latnorm, longnorm)
    Adat = createA(u, x, projmat)
    a = Adat$A
    ac = Adat$Ac
    y = createY(u, x, a, option = 'linear')
    if (method == 'spatialplus'){
      erfest = spatialplus(
        y = y,
        a = a,
        x = x,
        latnorm = latnorm,
        longnorm = longnorm,
        a.vals = a.vals
      )
      muests[,sim] = erfest$erf
    }
    if (method == 'bobb_exposurepenalized'){
      erfest = bobb_exposurepenalized(
        y = y,
        a = a,
        x = x,
        latnorm = latnorm,
        longnorm = longnorm,
        a.vals = a.vals
      )
      muests[,sim] = erfest$erf
    }
    if (method == 'KennedyERC'){
      l = cbind(projmat %*% a, x)
      #colnames(l) = c('ac', 'x')
      erfest = npcausal::ctseff(y, 
                                a, 
                                x = l,
                                a.rng = c(min(a.vals), max(a.vals)),
                                n.pts = length(a.vals),
                                bw.seq = seq(0.2, 2, length.out = 50))
      muests[,sim] = erfest$res$est
    }
    if (method == 'spectral_discrete'){
      stopifnot(!is.null(adjmat))
      erfest = spectral_discrete(
        y = y,
        a = a,
        x = x,
        a.vals = a.vals,
        adj = adjmat
      )
      muests[,sim] = erfest$erf
    }
    if (method == 'keller_szpiro_selectingscale_outcome'){
      erfest = keller_szpiro_selectingscale(
        y = y,
        a = a,
        x = x,
        latnorm = latnorm,
        longnorm = longnorm,
        a.vals = a.vals,
        fitmethod = 'outcome',
        selection = 'AIC-NE'
      )
      muests[,sim] = erfest$erf
    }
    if (method == 'keller_szpiro_selectingscale_preadjustment'){
      erfest = keller_szpiro_selectingscale(
        y = y,
        a = a,
        x = x,
        latnorm = latnorm,
        longnorm = longnorm,
        a.vals = a.vals,
        fitmethod = 'preadjustment',
        selection = 'AIC-NE'
      )
      muests[,sim] = erfest$erf
    }
    if (method == 'GPCERF'){
      data = cbind.data.frame(Y = y, treat = a, X = x, Ac = projmat %*% a)
      gps_m = GPCERF::estimate_gps(cov_mt = data[,-(1:2)],
                           w_all = data$treat,
                           sl_lib = c("SL.xgboost"),
                           dnorm_log = FALSE)
      cerf_obj = GPCERF::estimate_cerf_gp(data,
                              w=a.vals,
                              gps_m,
                              params = list(alpha = c(0.1,1),
                                            beta = c(0.1,1),
                                            g_sigma = c(0.1,1),
                                            tune_app = "all"),
                              outcome_col = "Y",
                              treatment_col = "treat",
                              gps_col = "Ac",
                              gps_m = gps_m)
      muests[,sim] = cerf_obj$cerf
    }
    if (method == 'GPCERF_nn'){
      data = cbind.data.frame(Y = y, treat = a, X = x, Ac = projmat %*% a)
      gps_m = GPCERF::estimate_gps(cov_mt = data[,-(1:2)],
                           w_all = data$treat,
                           sl_lib = c("SL.xgboost"),
                           dnorm_log = FALSE)
      cerf_nngp_obj = GPCERF::estimate_cerf_nngp(data,
                                         w=a.vals,
                                         gps_m,
                                         params = list(alpha = c(0.1,1,5),
                                                       beta = c(0.1,1,5),
                                                       g_sigma = c(0.1,1,5),
                                                       tune_app = "all",
                                                       n_neighbor = 100,
                                                       block_size = 1e4),
                                         outcome_col = "Y",
                                         treatment_col = "treat",
                                         covariates_col = c('X', 'Ac'),
                                         nthread = 1)
      muests[,sim] = cerf_nngp_obj$posterior$mean
    }
    if (method == 'unadjustedOLS'){
      mod = lm(y ~ a + x)
      
      erf = rep(NA, length(a.vals))
      for (i in 1:length(a.vals)){
        preds = predict(mod, newdata = data.frame(a = a.vals[i], x = x))
        erf[i] = mean(preds)
      }
      muests[,sim] = erf
    }
  }
  # if (method == 'nonoracle'){
  #   # TODO
  #   return(NULL)
  # }
  
  # Create dataframe whose first column is a.vals, second is mutrue, and the rest are muests
  df = cbind(a.vals, mutrue, muests)

  # write results to file
  # Check if file exists
  if (file.exists(filename)){
    # write new sims to file as new columns
    olddf = read.csv(filename)
    newdf = cbind(olddf, muests)
    write.csv(newdf, filename, row.names = FALSE)
  }
  # if file for ERF ests does not exist create it and write results
  else{
    write.csv(df, filename, row.names = FALSE)
  }
  invisible(filename)
}
