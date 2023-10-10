library(ggplot2)
library(gridExtra)

erfsim = function(v, # matrix used for projection vvt
                  beta1, # coeff of x
                  beta2, # coeff of x^2
                  betaxc, # coeff of x*cov
                  betac, # coeff of cov
                  sig, # gaussian error sd
                  x, # treatment
                  cov # covariate 
                  ) {
  vvt = v %*% t(v)
  vvtx = as.numeric(vvt %*% x)
  rx = x - vvtx
  n = nrow(v)
  err = rnorm(n, mean = 0, sd = sig)
  # Generate outcome
  y = beta1*x + beta2*x^2 + betaxc*x*cov + betac*cov + err
  # Estimate ERF with original data
  lm1 = lm(y ~ x + I(x^2) + x:cov + cov)
  # Predict for each value of x and average
  newx = seq(-1, 1 - 2 / n, by = 2 / n)
  erf1 = rep(NA, length(newx))
  for (i in 1:length(newx)){
    newdata = data.frame('x' = x, 'cov' = cov)
    newdata$x = newx[i]
    res = predict(lm1, newdata =  newdata)
    erf1[i] = mean(res)
  }
  
  # Estimate ERF with projected exposure
  lm2 = lm(y ~ vvtx + I(vvtx^2) + cov + cov:vvtx + rx + I(rx^2) + cov:rx + rx:vvtx)
  # Predict for each value of x and average
  newx = seq(-1, 1 - 2 / n, by = 2 / n)
  erf2 = rep(NA, length(newx))
  for (i in 1:length(newx)){
    newdata = data.frame('vvtx' = vvtx, 'cov' = cov, 'rx' = rx)
    newdata$vvtx = newx[i]
    #newdata$rx = rep(0, n)
    res = predict(lm2, newdata =  newdata)
    erf2[i] = mean(res)
  }
  df = cbind.data.frame(seq(-1, 1 - 2 / n, by = 2 / n),
                        erf1,
                        erf2)
  colnames(df) = c('predx', 'lm1_pred', 'lm2_pred')
  return(list('df' = df, 'vvtx' = data.frame('obs_vvtx' = vvtx)))
}

# fxn to create plots and run simulation
erfplot = function(n, # sample size
                   v, # matrix used for projection vvt
                   nreps, # number of times that new y vector created
                   beta1 = 2, # coeff x
                   beta2 = 0, # coeff x^2
                   betaxc = 0, # coeff x*cov
                   betac = 0, # coeff cov
                   sig = 1, # gaussian error variance
                   x = rexp(n) - 1, # exposure
                   cov = rgamma(n, shape = 2, rate = 1), # can I do this
                   filename) { # filename for jpeg output
  gs = list()
  # stores predictions across simulations.
  allpred = data.frame(matrix(nrow = 2*n*nreps, ncol = 4))
  names(allpred) = c('predx', 'pred', 'model', 'sim')
  for (r in 1:nreps){
    dfrep = erfsim(v=v, beta1=beta1, beta2=beta2, betaxc=betaxc, betac=betac, sig=sig, x=x, cov=cov)
    allpred[(2*n*(r-1)+1):((2*r-1)*n),] = cbind.data.frame(dfrep$df$predx, 
                                                           dfrep$df$lm1_pred, 
                                                           rep('x', n), 
                                                           rep(r, n))
    allpred[(((2*r-1)*n)+1):(2*r*n),] = cbind.data.frame(dfrep$df$predx, 
                                                         dfrep$df$lm2_pred, 
                                                         rep('proj_x', n), 
                                                         rep(r, n))
    # Plot data from a single run
    gs[[r]] = ggplot(dfrep$df, aes(x = predx)) + 
      geom_line(aes(y = lm1_pred, color = 'ERF from x')) + 
      geom_line(aes(y = lm2_pred, color = 'ERF from proj x')) + 
      geom_rug(data=dfrep$vvtx, aes(x = obs_vvtx), sides = 'b', inherit.aes = F) + 
      labs(y = 'prediction', x = 'x or proj x') + 
      scale_colour_manual("", 
                          breaks = c("ERF from x", "ERF from proj x"),
                          values = c("blue", "red")) +
      ylim(-4,4) +
      xlim(-1,1) +
      theme_minimal()
  }
  gs[[nreps + 1]] = ggplot(allpred, aes(predx, pred, color = model)) +
    stat_smooth(aes(group = interaction(sim, model)), method = 'loess', se = FALSE, lty = 3) +
    stat_smooth(method = 'loess', se = FALSE, lty = 1, size = 2) + 
    labs(title = 'All ERFs, with Loess')
  png(paste('images/erfplots/', filename, '.jpeg', sep = ''), height = 1024, width = 2000)
  # Just plot 9 of them
  lay <- rbind(c(1,2,3,10, 10, 10),
               c(4,5,6,10,10,10),
               c(7,8,9,10,10,10))
  # plot the 9 in a grid next to the plot of loess curves
  grid.arrange(grobs = gs[c(1:9, nreps+1)], layout_matrix = lay)
  dev.off()
}

# Synthetic data
n = 300
A = matrix(rexp(n^2), nrow = n, ncol = n)
v = svd(A)$u[,291:300]
nreps = 20
filename = 'linear'
erfplot(n=n, v=v, nreps=nreps, filename=filename)

filename = 'quadratic'
erfplot(n=n, v=v, beta2=-1, nreps=nreps, filename=filename)

filename = 'linearinteraction'
erfplot(n=n, v=v, betaxc=1, betac = 1, nreps=nreps, filename=filename)

filename = 'linearquadraticinteraction'
erfplot(n=n, v=v, betaxc=1, beta2=-2, betac=1, nreps=nreps, filename=filename)

# Real data
set.seed(20)
study = read.csv('/Users/sophie/Documents/SpatialConf/archived/Study_dataset_2010.csv')
adj = read.csv("/Users/sophie/Documents/SpatialConf/archived/adjacency_matrix.csv",
               header = F) # created from spacebench script
R = diag(rowSums(adj)) - adj # graph laplacian 
E = eigen(R) # eigen component
G = E$vectors
v = G[,1:10] # 3099:3108 #1:10 
x = study$qd_mean_pm25 - mean(study$qd_mean_pm25)
cov = scale(study$gmet_mean_summer_rmn)
n = 3109
nreps = 20
filename = 'real_linear'
erfplot(n=n, v=as.matrix(v), nreps=nreps, x = as.numeric(x), cov = as.numeric(cov), 
        filename=filename, sig = 1)
filename = 'real_linearquadraticinteraction'
erfplot(n=n, v=as.matrix(v), betaxc=1, beta2=-2, betac=1, nreps=nreps,  x = as.numeric(x), cov = as.numeric(cov), filename=filename)

