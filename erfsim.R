library(ggplot2)
library(gridExtra)
library(mvtnorm)
erfsim = function(vvt, # projection matrix
                  beta1, # coeff of x
                  beta2, # coeff of x^2
                  betaxc, # coeff of x*cov
                  betac, # coeff of cov
                  sig, # gaussian error sd
                  x, # treatment
                  cov # covariate 
                  ) {
  vvtx = as.numeric(vvt %*% x)
  rx = x - vvtx
  n = nrow(vvt)
  err = rnorm(n, mean = 0, sd = sig)
  # Generate outcome
  y = beta1*x + beta2*x^2 + betaxc*x*cov + betac*cov + err
  # Estimate ERF with original data
  lm1 = lm(y ~ x + I(x^2) + x:cov + cov)
  lm3 = lm(y ~ x + I(x^2)) # excludes covariate
  # Predict for each value of x and average
  newx = x + 1
  erf1 = rep(NA, length(newx))
  erf3 = rep(NA, length(newx))
  for (i in 1:length(newx)){
    # ERF is average of predictions
    newdata1 = data.frame('x' = x, 'cov' = cov)
    newdata3 = data.frame('x' = x)
    newdata1$x = newx[i]
    newdata3$x = newx[i]
    res1 = predict(lm1, newdata =  newdata1)
    erf1[i] = mean(res1)
    res3 = predict(lm3, newdata =  newdata3)
    erf3[i] = mean(res3) # constant anyway
  }
  
  # Estimate ERF with projected exposure
  lm2 = lm(y ~ vvtx + I(vvtx^2) + cov + cov:vvtx + rx + I(rx^2) + cov:rx + rx:vvtx)
  lm4 = lm(y ~ vvtx + I(vvtx^2)+ rx + I(rx^2) + rx:vvtx) # excludes covariate

  # Predict for each value of x and average
  erf2 = rep(NA, length(newx))
  erf4 = rep(NA, length(newx))
  
  for (i in 1:length(newx)){
    newdata2 = data.frame('vvtx' = vvtx, 'cov' = cov, 'rx' = rx)
    newdata4 = data.frame('vvtx' = vvtx, 'rx' = rx)
    
    newdata2$vvtx = newx[i]-rx
    newdata4$vvtx = newx[i]-rx
    res2 = predict(lm2, newdata = newdata2)
    erf2[i] = mean(res2)
    res4 = predict(lm4, newdata = newdata4)
    erf4[i] = mean(res4)
  }
  df = cbind.data.frame(x,
                        erf1,
                        erf2,
                        erf3,
                        erf4)
  colnames(df) = c('predx', 'lm1_pred', 'lm2_pred', 'lm3_pred', 'lm4_pred')
  return(list('df' = df, 
              'vvtx' = data.frame('obs_vvtx' = vvtx))) # also save observed projected x
}

# fxn to create plots and run simulation
erfplot = function(n, # sample size
                   vvt, # projection matrix
                   nreps, # number of times that new y vector created
                   beta1 = 2, # coeff x
                   beta2 = 0, # coeff x^2
                   betaxc = 0, # coeff x*cov
                   betac = 0, # coeff cov
                   sig = 1, # gaussian error variance
                   realdat = F,
                   confounding = F) {
  gs = list()
  # stores predictions across simulations.
  allpred = data.frame(matrix(nrow = 4*n*nreps, ncol = 4))
  names(allpred) = c('predx', 'pred', 'model', 'sim')
  for (r in 1:nreps){
    if (!realdat){
      if (!confounding){
        x = rexp(n) -1
        cov = rgamma(n, shape = 2, rate = 1)
      }
      if (confounding){
        x = rt(n, df = 2)
        x_noconf = vvt %*% x
        x_conf = x - x_noconf
        cov_noconf = (diag(1,n)-vvt) %*% rnorm(n, mean = 2, sd = 1)
        cov_conf = x_conf
        cov = as.numeric(cov_conf + cov_noconf)
      }
    }
    if (realdat){
      return(NULL) # to do
    }
    # Calculate ERFs
    dfrep = erfsim(vvt=vvt, beta1=beta1, beta2=beta2, betaxc=betaxc, betac=betac, sig=sig, x=x, 
                   cov=cov)
    
    # There must be a better way to code this. This is terrible. 
    data1 = cbind.data.frame(dfrep$df$predx, dfrep$df$lm1_pred, rep('x', n), rep(r, n))
    data2 = cbind.data.frame(dfrep$df$predx, dfrep$df$lm2_pred, rep('proj_x', n), rep(r, n))
    data3 = cbind.data.frame(dfrep$df$predx, dfrep$df$lm3_pred, rep('x_nocov', n), rep(r, n))
    data4 = cbind.data.frame(dfrep$df$predx, dfrep$df$lm4_pred, rep('proj_x_nocov', n), rep(r, n))
    colnames(data1) = colnames(allpred)
    colnames(data2) = colnames(allpred)
    colnames(data3) = colnames(allpred)
    colnames(data4) = colnames(allpred)

    # Calculate the row indices
    start_idx = (4 * n * (r - 1) + 1)
    end_idx = (4 * r * n)

    # Update the rows in allpred
    allpred[start_idx:end_idx, ] = rbind(data1, data2, data3, data4)
    
    # Plot ERFs from a single run
    if (r <= 9){
      gs[[r]] = ggplot(dfrep$df, aes(x = predx)) + 
        geom_line(aes(y = lm1_pred, color = 'x ERF')) + 
        geom_line(aes(y = lm2_pred, color = 'proj x ERF')) + 
        geom_line(aes(y = lm3_pred, color = 'x ERF no cov')) + 
        geom_line(aes(y = lm4_pred, color = 'proj x ERF no cov')) + 
        geom_rug(data=dfrep$vvtx, aes(x = obs_vvtx), sides = 'b', inherit.aes = F) + 
        labs(y = 'prediction', x = 'x or proj x') + 
        scale_colour_manual("", 
                            breaks = c("x ERF", "proj x ERF", "x ERF no cov", "proj x ERF no cov"),
                            values = c("blue", "red", "lightblue", "pink")) +
        ylim(-10,10) +
        theme_minimal()
    }
  }
  # Plot ERFs across runs
  gs[[10]] = ggplot(allpred, aes(predx, pred, color = model)) +
    stat_smooth(aes(group = interaction(sim, model)), method = 'loess', se = F, lty = 3, size = 0.5) +
    stat_smooth(method = 'loess', se = F, lty = 1, size = 2) + 
    labs(title = 'All ERFs, with Loess') + 
    scale_colour_manual("", 
                        values = c(x = "blue", 
                                   proj_x = "red", 
                                   x_nocov = "lightblue", 
                                   proj_x_nocov = "pink")) +
    ylim(-10,10)
  return(gs)
}

# Synthetic data
n = 20
A = matrix(rexp(n^2), nrow = n, ncol = n)
v = svd(A)$u[,5:10] #291:300 # 1:10
vvt = v %*% t(v)
nreps = 100

filename = 'linearquadraticinteraction_confounding_byscale'
gs = erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps, filename=filename, 
             confounding = T)
png(paste('images/erfplots/', filename, '.jpeg', sep = ''), height = 1000, width = 1000, res = 100)
gs[[10]]
dev.off()

#filename = 'linearquadraticinteraction'
#erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps, filename=filename)

# grid.arrange(grobs = gs[1:9])
# 
# filename = 'linear'
# erfplot(n=n, vvt=vvt, nreps=nreps, filename=filename)
# 
# filename = 'quadratic'
# erfplot(n=n, vvt=vvt, beta2=-1, nreps=nreps, filename=filename)
# 
# filename = 'linearinteraction'
# erfplot(n=n, vvt=vvt, betaxc=1, betac = 1, nreps=nreps, filename=filename)
# 
# filename = 'linearquadraticinteraction'
# erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps, filename=filename)

# Real data
# set.seed(20)
# study = read.csv('/Users/sophie/Documents/SpatialConf/archived/Study_dataset_2010.csv')
# adj = read.csv("/Users/sophie/Documents/SpatialConf/archived/adjacency_matrix.csv",
#                header = F) # created from spacebench script
# x = study$qd_mean_pm25 - mean(study$qd_mean_pm25)
# cov = scale(study$gmet_mean_summer_rmn)
# n = 3109
# nreps = 20
# R = diag(rowSums(adj)) - adj # graph laplacian 
# E = eigen(R) # eigen component
# G = E$vectors

# v = G[,1:1000]
# vvt = v %*% t(v)
# filename = 'real_linearquadraticinteraction_spectral1'
# erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps,  
#         x = as.numeric(x), cov = as.numeric(cov), filename=filename)
# v = G[,1001:2000] 
# vvt = v %*% t(v)
# filename = 'real_linearquadraticinteraction_spectral2'
# erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps,  
#         x = as.numeric(x), cov = as.numeric(cov), filename=filename)
# v = G[,2001:3000]
# vvt = v %*% t(v)
# filename = 'real_linearquadraticinteraction_spectral3'
# erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps,  
#         x = as.numeric(x), cov = as.numeric(cov), filename=filename)
# v = G[,3099:3108] # 3099:3108 #1:10 
# vvt = v %*% t(v)
# filename = 'real_linearquadraticinteraction_spectral4'
# erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps,  
#         x = as.numeric(x), cov = as.numeric(cov), filename=filename)
# 
# # Nested
# groups = cbind(study$region, study$STATE)
# nest = nested_decomp_mats(groups)
# vvt = nest$decomp_mats[[1]]
# filename = 'real_linearquadraticinteraction_nested1'
# erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps,  
#         x = as.numeric(x), cov = as.numeric(cov), filename=filename)
# vvt = nest$decomp_mats[[2]]
# filename = 'real_linearquadraticinteraction_nested2'
# erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps,  
#         x = as.numeric(x), cov = as.numeric(cov), filename=filename)
# vvt = nest$decomp_mats[[3]]
# filename = 'real_linearquadraticinteraction_nested3'
# erfplot(n=n, vvt=vvt, betaxc=1, beta2=-2, betac=1, nreps=nreps,  
#         x = as.numeric(x), cov = as.numeric(cov), filename=filename)

# Simulation would be to generate x and cov independently from Gaussian Random Fields
# Replace the confounded part (v2v2t) cov with a highly correlated part of v2v2tx

