simfunc = function(nreps, dgm = c('linear', 'quadratic', 'interaction'),
                   interactionsd = 3){
  outcomemod = match.arg(dgm)
  xcoeffs = rep(NA, nreps)
  vvtxcoeffs = rep(NA, nreps)
  vvtyxcoeffs = rep(NA, nreps)
  vvtyvvtxcoeffs = rep(NA, nreps)
  for (i in 1:nreps){
    n = 10
    A = matrix(rexp(n^2), nrow = n, ncol = n)
    v = svd(A)$u[,c(1,2)] # by construction orthonormal
    vvt = v%*%t(v)
    x = rexp(n)
    if (outcomemod == 'linear'){
      y = 2*x
    }
    if (outcomemod == 'quadratic'){
      y = 2*x^2
    }
    if (outcomemod == 'interaction'){
      c = rnorm(n, mean = 0, sd = interactionsd) # mean 1
      y = 2*x + x*c
    }
    vvtx = vvt%*%x
    rx = x - vvtx
    vvty = vvt%*%y
    xcoeffs[i] = summary(lm(y~x))$coefficients[2,1]
    vvtxcoeffs[i] = summary(lm(y~vvtx + rx))$coefficients[2,1] 
    vvtyxcoeffs[i] = summary(lm(vvty ~ x))$coefficients[2,1]
    vvtyvvtxcoeffs[i] = summary(lm(vvty ~ vvtx + rx))$coefficients[2,1] 
  }
  df = data.frame(matrix(c(mean(xcoeffs),
                           mean(vvtxcoeffs),
                           mean(vvtyxcoeffs),
                           mean(vvtyvvtxcoeffs),
                           sd(xcoeffs),
                           sd(vvtxcoeffs),
                           sd(vvtyxcoeffs),
                           sd(vvtyvvtxcoeffs)), nrow = 4, ncol = 2, byrow = F), 
                  row.names = c('y~x', 'y~vvtx+rx', 'vvty~x', 'vvty ~ vvtx + rx'))
  colnames(df) = c('mean', 'sd')
  return(df)
}

set.seed(111)
simfunc(1000, dgm = 'linear')
simfunc(1000, dgm = 'interaction', interactionsd = 1)
simfunc(1000, dgm = 'interaction', interactionsd = 5)
simfunc(1000, dgm = 'interaction', interactionsd = 10)


