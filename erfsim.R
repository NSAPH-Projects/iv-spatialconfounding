library(ggplot2)
library(gridExtra)
library(mvtnorm)
erfsim = function(vvt, # projection matrix
                  beta1, # coeff of x
                  beta2, # coeff of x^2
                  beta3, # coeff of x^3
                  betaxc, # coeff of x*cov
                  betaxc2, # coeff of x*cov^2
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
  y = beta1*x + beta2*x^2 + beta3*x^3 + betaxc*x*cov + betaxc2*x^2*cov + betac*cov + err
  # Estimate ERF with original data
  lm1 = lm(y ~ I(x^2) + I(x^3) + x*cov + I(x^2):cov)
  lm3 = lm(y ~ x + I(x^2) + I(x^3)) # excludes covariates
  
  # Predict for each value of x and average
  # Choose newx range that is either shifted original or that vvtx covers well
  #newx = x + 1
  #newx = seq(-2.5, 2.5, 5/(n-1))
  newx = seq(-4,0.5, 4.5/(n-1))
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
  lm2 = lm(y ~ vvtx + I(vvtx^2) + I(vvtx^3) + 
             rx + I(rx^2) + I(rx^3) + 
             cov + cov:vvtx + cov:rx + rx:vvtx + 
             I(vvtx^2):rx + I(rx^2):vvtx + 
             I(vvtx^2):cov + I(rx^2):cov + vvtx:rx:cov)
  lm4 = lm(y ~ vvtx + I(vvtx^2) + I(vvtx^3) + 
             rx + I(rx^2) + I(rx^3) + 
             rx:vvtx + 
             I(vvtx^2):rx + I(rx^2):vvtx) # excludes confounder
  # Predict for each value of x and average
  erf2 = rep(NA, length(newx))
  erf4 = rep(NA, length(newx))
  
  for (i in 1:length(newx)){
    newdata2 = data.frame('vvtx' = vvtx, 'cov' = cov, 'rx' = rx)
    newdata4 = data.frame('vvtx' = vvtx, 'rx' = rx)
    
    newdata2$vvtx = newx[i]-rx
    newdata4$vvtx = newx[i]-rx
    res2 = predict(lm2, newdata = newdata2)
    erf2[i] = mean(res2) # technically should be doing differently (two means)
    # but ok for now because red curves look great
    res4 = predict(lm4, newdata = newdata4)
    erf4[i] = mean(res4)
  }
  df = cbind.data.frame(newx, #x
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
                   beta3 = 0, # coeff x^3
                   betaxc = 0, # coeff x*cov
                   betaxc2 = 0, # coeff x*cov^2
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
        x = rt(n, ncp = -1, df = 2)
        x_noconf = vvt %*% x
        x_conf = x - x_noconf
        cov_noconf = vvt %*% rnorm(n, mean = 2, sd = 1)
        cov_conf = 0.8*scale(x_conf) + sqrt(1-0.8^2)*(scale((diag(1,n)-vvt) %*% rnorm(n, mean = 0.5)))
        cov = as.numeric(cov_conf + cov_noconf)
      }
    }
    if (realdat){
      return(NULL) # to do
    }
    # Calculate ERFs
    dfrep = erfsim(vvt=vvt, beta1=beta1, beta2=beta2, beta3=beta3, betaxc=betaxc, 
                   betaxc2=betaxc2, betac=betac, sig=sig, x=x, cov=cov)
    
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
                            values = c("blue", "red", "darkturquoise", "coral")) +
        ylim(-10,10) +
        theme_minimal()
    }
  }
  # Plot ERFs across runs
  gs[[10]] = ggplot(allpred, aes(predx, pred, color = model)) +
    geom_line(aes(group = interaction(sim, model)), lty = 3, size = 0.5) +
    stat_smooth(se = F, lty = 1, size = 2) + 
    labs(title = 'All ERFs') + 
    scale_colour_manual("", 
                        values = c(x = "blue", 
                                   proj_x = "red", 
                                   x_nocov = "darkturquoise", 
                                   proj_x_nocov = "coral")) +
    ylim(-10,10)
  return(gs)
}

# Synthetic data
set.seed(15)
n = 50
A = matrix(rexp(n^2), nrow = n, ncol = n)
A[lower.tri(A)] = t(A)[lower.tri(A)]
E = eigen(A)
v = E$vectors[,1:30] 
#v = E$vectors[,20:49]
vvt = v %*% t(v)
nreps = 100

# Create plots
gs = erfplot(n=n, vvt=vvt,nreps=nreps, confounding =T)
g1 = gs[[10]] + labs(title = 'simple linear') 
gs = erfplot(n=n, vvt=vvt,nreps=nreps, confounding =T,
             beta1 = 4, beta2 = -2, beta3 = -1) 
g2 = gs[[10]] + labs(title = 'cubic') 
gs = erfplot(n=n, vvt=vvt,nreps=nreps, confounding =T,
             beta1 = 1, beta2 = -0.5, beta3 = -0.25, betac = -1)
g3 = gs[[10]] + labs(title = 'cubic x + additive confounder effect')  
gs = erfplot(n=n, vvt=vvt,nreps=nreps, confounding =T,
             beta1 = 1, beta2 = -0.5, beta3 = -0.25, betac = -1, betaxc = -1) 
g4 = gs[[10]] + labs(title = 'cubix x + interaction') 
gs = erfplot(n=n, vvt=vvt,nreps=nreps, confounding =T,
             betac = -5, betaxc = -1, betaxc2 = 2)
g5 = gs[[10]] + labs(title = 'quadratic confounder interaction')
gs = erfplot(n=n, vvt=vvt,nreps=nreps, confounding =T,
             beta2 = -2, beta3 = 2, betac = -1, betaxc = -1, betaxc2 = 2)
g6 = gs[[10]] + labs(title = 'cubic x + quadratic confounder interaction') + ylim(-15,15) 

filename = 'plots_combined_130_100_imperfect'
#filename = 'plots_combined_2049_100_imperfect'
png(paste('images/erfplots/', filename, '.jpeg', sep = ''), height = 2000, width = 2000, res = 100)
grid.arrange(grobs = list(g1, g2, g3, g4, g5, g6), 
             layout = matrix(1:6, nrow = 3, ncol = 2))
dev.off()

# Exploratory data analysis
n = 500
A = matrix(rexp(n^2), nrow = n, ncol = n)
A[lower.tri(A)] = t(A)[lower.tri(A)]
E = eigen(A)
#v = E$vectors[,1:30] 
v = E$vectors[,20:49]
vvt = v %*% t(v)
cors = matrix(NA, nrow = nreps, ncol = 6) 
# x cov, x1 cov1, x2 cov1, x2 co1, x2 cov2
for (i in 1:nreps){
  x = rt(n, ncp = -1, df = 2)
  x_noconf = vvt %*% x
  x_conf = x - x_noconf
  cov_noconf = vvt %*% rnorm(n, mean = 2, sd = 1)
  cov_conf = 0.8*scale(x_conf) + sqrt(1-0.8^2)*(scale((diag(1,n)-vvt) %*% rnorm(n, mean = 0.5)))
  cov = as.numeric(cov_conf + cov_noconf)
  cors[i,] = c(mean(x*cov), # 1-8
               mean(x_noconf*cov_noconf), # -.16-.16
               mean(x_conf*cov_noconf), # 0
               mean(x_noconf*cov_conf), # 0
               mean(x_conf*cov_conf), # 1-8
               mean(x_noconf*cov)) # -.16-.16
}
summary(cors)

