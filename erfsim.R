library(ggplot2)
library(gridExtra)

erfsim = function(v, beta1, beta2, betaxc, betac, sig){
  n = nrow(v)
  x = rexp(n)-1 # mean zero
  cov = rgamma(n, shape = 2, rate = 1)
  y = beta1*x + beta2*x^2+ betaxc*x*cov + betac*cov + rnorm(n, mean = 0, sd = sig)
  lm1 = lm(y~ x + I(x^2) + x:cov + cov)
  pred1 = predict(lm1, newdata = data.frame(x = seq(-1, 1-2/n, by = 2/n), 
                                            cov = rgamma(n, shape = 2, rate = 1)))

  # by construction orthonormal
  vvt = v%*%t(v)
  vvtx = as.numeric(vvt%*%x)
  print(mean(vvtx))
  rx = x - vvtx
  lm2 = lm(y ~ vvtx + I(vvtx^2) + cov + cov:vvtx + rx + I(rx^2) + cov:rx + rx:vvtx)
  pred2 = predict(lm2, newdata = data.frame(vvtx = seq(-1, 1-2/n, by = 2/n), rx = rep(0, n), 
                                            cov = rnorm(n, mean = 2, sd = 0.5)))
  pred3 = predict(lm1)
  pred4 = predict(lm2)
  df = cbind.data.frame(seq(-1, 1-2/n, by = 2/n),
                         pred1,
                         pred2, 
                        pred3, 
                        pred4)
  colnames(df) = c('predx', 'lm1_pred', 'lm2_pred', 'lm3_pred', 'lm4_pred')
  # average across values of x (binning) to create ERF: ok since cov indep of x?
  dffinal = aggregate(df, 
                      by=list(cut(df$predx,seq(-1,1,0.1))), 
                      mean) 
  return(list('df' = dffinal, 'vvtx' = data.frame('obs_vvtx' = vvtx)))
}

erfplot = function(n, A, v, nreps, beta1=2, beta2=0, betaxc=0, betac=0, sig=1, filename){
  gs = list()
  allpred = data.frame(matrix(nrow = 2*n*nreps, ncol = 4))
  names(allpred) = c('predx', 'pred', 'model', 'sim')
  for (r in 1:nreps){
    dfrep = erfsim(v=v, beta1, beta2, betaxc, betac, sig)
    allpred[(2*n*(r-1)+1):((2*r-1)*n),] = cbind.data.frame(dfrep$df$predx, dfrep$df$lm1_pred, 
                                                           rep('x', n), rep(r, n))
    allpred[(((2*r-1)*n)+1):(2*r*n),] = cbind.data.frame(dfrep$df$predx, dfrep$df$lm2_pred, 
                                                         rep('proj_x', n), rep(r, n))
    gs[[r]] = ggplot(dfrep$df, aes(x = predx)) + 
      geom_line(aes(y = lm1_pred, color = 'ERF from x')) + 
      geom_line(aes(y = lm2_pred, color = 'ERF from proj x')) + 
      geom_rug(data=dfrep$vvtx, aes(x = obs_vvtx), sides = 'b', inherit.aes = F) + 
      labs(y = 'prediction', x = 'x or proj x') + 
      scale_colour_manual("", 
                          breaks = c("ERF from x", "ERF from proj x"),
                          values = c("blue", "red")) +
      ylim(-4,4) +
      theme_minimal()
  }
  gs[[nreps + 1]] = ggplot(allpred, aes(predx, pred, color = model)) +
    stat_smooth(aes(group = interaction(sim, model)), method = 'loess', se = FALSE, lty = 3) +
    stat_smooth(method = 'loess', se = FALSE, lty = 1) + 
    labs(title = 'ERFs on a single plot')
  png(paste('images/erfplots/', filename, '.jpeg', sep = ''), height = 1024, width = 1200)
  # Just plot 8 of them
  do.call(grid.arrange,gs[c(1:8, nreps+1)])
  dev.off()
}

n = 300
A = matrix(rexp(n^2), nrow = n, ncol = n)
v = svd(A)$u[,1:10]
nreps = 20
filename = 'linear'
erfplot(n=n, A=A, v=v, nreps=nreps, filename=filename)

filename = 'quadratic'
erfplot(n=n, A=A, v=v, beta2=-1, nreps=nreps, filename=filename)

filename = 'linearinteraction'
erfplot(n=n, A=A, v=v, betaxc=1, betac = 1, nreps=nreps, filename=filename)

filename = 'linearquadraticinteraction'
erfplot(n=n, A=A, v=v, betaxc=1, beta2=-2, betac=1, nreps=nreps, filename=filename)

