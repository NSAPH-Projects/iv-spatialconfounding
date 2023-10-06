library(ggplot2)
library(gridExtra)

erfsim = function(n, beta1, beta2, betaxc, betac, sig){
  x = rexp(n)-1 # mean zero
  cov = rgamma(n, shape = 2, rate = 1)
  y = beta1*x + beta2*x^2+ betaxc*x*cov + betac*cov + rnorm(n, mean = 0, sd = sig)
  lm1 = lm(y~ x + I(x^2) + x:cov + cov)
  pred1 = predict(lm1, newdata = data.frame(x = seq(-1, 1-2/n, by = 2/n), 
                                            cov = rgamma(n, shape = 2, rate = 1)))

  A = matrix(rexp(n^2), nrow = n, ncol = n)
  v = svd(A)$u[,1:10] # by construction orthonormal
  vvt = v%*%t(v)
  vvtx = as.numeric(vvt%*%x)
  print(mean(vvtx))
  rx = x - vvtx
  lm2 = lm(y ~ vvtx + I(vvtx^2) + cov + cov:vvtx + rx + I(rx^2) + cov:rx + rx:I(vvtx^2))
  pred2 = predict(lm2, newdata = data.frame(vvtx = seq(-1, 1-2/n, by = 2/n), rx = rep(0, n), 
                                            cov = rnorm(n, mean = 2, sd = 0.5)))
  pred3 = predict(lm1)
  pred4 = predict(lm2)
  df = cbind.data.frame(x,
                         vvtx, 
                         seq(-1, 1-2/n, by = 2/n),
                         pred1,
                         pred2, 
                        pred3, 
                        pred4)
  colnames(df) = c('obs_x', 'obs_vvtx', 'predx', 'lm1_pred', 'lm2_pred', 'lm3_pred', 'lm4_pred')
  # average across values of x to create ERF
  #dffinal = aggregate(df, 
  #                    by=list(cut(df1$predx,seq(-1,1,0.1))), 
  #                    mean) 
  return(df)
}
# ERF: average over COVARIATE DISTRIBUTION
nreps = 9
gs = list()
for (rep in 1:nreps){
  dfrep = erfsim(n=500, beta1=2, beta2=0, betaxc=0, betac=0, sig=1)
  gs[[rep]] = ggplot(dfrep, aes(x = predx)) + 
    geom_line(aes(y = lm1_pred, color = 'ERF from x')) + 
    geom_line(aes(y = lm2_pred, color = 'ERF from proj x')) + 
    geom_rug(aes(x = obs_vvtx), sides = 'b') + 
    labs(title = 'ERF', y = 'prediction', x = 'x or proj x') + 
    scale_colour_manual("", 
                        breaks = c("ERF from x", "ERF from proj x"),
                        values = c("blue", "red")) + 
    theme_minimal()
}
png('images/erfs.jpeg', height = 1024 * 2, width = 1024 * 2)
do.call(grid.arrange,gs)
dev.off()