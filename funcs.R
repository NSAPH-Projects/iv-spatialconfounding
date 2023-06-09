library(MASS)
nested_decomp = function(vecs,group){
  # Create the decomposition matrices
  n = length(group)
  nst = as.data.frame(table(group))
  colnames(nst) = c('group', 'frequency')
  
  # create projection matrices
  G1 = matrix(NA, nrow = n, ncol = n)
  for (j in 1:nrow(G1)){
    st = group[j]
    stfac  = 1/nst$frequency[nst$group == st]
    for (k in 1:ncol(G1)){
      G1[j,k] = stfac*(st == group[k])
    }
  }
  
  A1 = G1
  A2 = diag(n)-G1 
  return(list('V1' = A1%*%vecs, 
              'V2' = A2%*%vecs))
}

spectral_decomp = function(vecs,A,inv = F){
  R = diag(rowSums(A))-A # precision of ICAR
  E = eigen(R) # eigen component
  D = E$val
  G = E$vec
  rm(E,R)
  if (inv){
    invG = solve(t(G)) 
    return(invG%*%vecs)
  }
  return(t(G)%*%vecs)
}

sim = function(n,# n^2 x n^2 is dim of adjacency matrix
               betax = 2,
               betaz = -1,
               betaxz = 0,
               sig = 1,
               rhox = c(0.7,0.1),
               outcome=c('linear', 'quadratic', 'interaction'), # outcome model
               decomposition = c('spectral', 'nested'), # data generating
               distribution = 'gaussian' # TO DO
               ){
  
  # Create coordinates
  df = as.data.frame(expand.grid(xcoord = 1:(n^2), ycoord = 1:(n^2)))
  df$state = as.numeric(paste((df$xcoord-1)%/%n,(df$ycoord-1)%/%n, sep = ''))
  df = df[order(df$state),]
  
  # Create adjacency matrix
  A = matrix(NA, nrow = nrow(df), ncol = nrow(df))
  for (i in 1:nrow(df)){
    coord = c(df$xcoord[i], df$ycoord[i])
    for (j in 1:i){
      A[i,j] = 1*(abs(df$xcoord[j]-coord[1]) + abs(df$ycoord[j]-coord[2]) == 1)
    }
  }
  for (i in 1:(nrow(df)-1)){
    for (j in (i+1):nrow(df)){
      A[i,j] = A[j,i]
    }
  }
  
  decomp = match.arg(decomposition)
  outcomemodel = match.arg(outcome)
  # Simulate X and Z
  if (decomp == 'nested'){
    stopifnot(length(rhox)==2)
    xz = mvrnorm(n^2, mu = c(0,0), 
                  Sigma = matrix(c(1,rhox[1],rhox[1],1), 
                                 nrow = 2, ncol = 2))
    X1 = xz[,1]
    Z1 = xz[,2]
    df$X1 = rep(X1,each = n^2)
    df$Z1 = rep(Z1, each = n^2)
    
    # add county-level variation
    xz = mvrnorm(n^4, mu = c(0,0), 
                  Sigma = matrix(c(1,rhox[2],rhox[2],1), 
                                 nrow = 2, ncol = 2))
    df$X2 = xz[,1]
    df$Z2 = xz[,2]
    df$X = df$X1 + df$X2
    df$Z = df$Z1 + df$Z2
  }
  
  if (decomp == 'spectral'){
    stopifnot(length(rhox)==n^4)
    Xstar = rep(NA, n^4)
    Zstar = rep(NA, n^4)
    for (i in 1:(n^4)){
      xz = mvrnorm(1, mu = c(0,0), 
                   Sigma = matrix(c(1,rhox[i],rhox[i],1), 
                                  nrow = 2, ncol = 2))
      Xstar[i] = xz[1]
      Zstar[i] = xz[2]
    }
    # Project into spatial domain
    spec = spectral_decomp(cbind(Xstar,Zstar),A,inv=T)
    df$X = spec[,1] # n^4 length vector
    df$Z = spec[,2]
  }
  # Simulate the outcome
  if (outcomemodel == 'linear'){
    df$Y = betax*df$X + betaz*df$Z + rnorm(length(df$X), mean = 0, 
                                  sd = sig)
  }
  if (outcomemodel == 'quadratic'){
    stopifnot(length(betax)>=2)
    df$Y = betax[1]*df$X + betax[2]*df$X^2 + betaz*df$Z + rnorm(length(df$X), mean = 0, 
                                                    sd = sig)
  }
  if (outcomemodel == 'interaction'){
    Y = betax*df$X + betaz*df$Z + betaxz*df$X*df$Z + rnorm(length(df$X), 
                                                           mean = 0, sd = sig)
  }
  return(list(
    'coord' = cbind(df$xcoord, df$ycoord),
    'A'=A,#adjacency mat
    'X'=df$X,#exposure
    'Y'=df$Y,#outcome
    'Z'=df$Z, #confounder
    'state' = df$state #nested group
    ))
}

analysis = function(A, # adjacency
                    X, # exposure
                    Y, # outcome
                    Z, # confounder
                    state, # nested group
                    decomposition = c('spectral', 'nested')
                    ){
  decomposition = match.arg(decomposition)
  if (decomposition == 'nested'){
    # Regress Y on X at each level
    nest = nested_decomp(cbind(X,Y), state)
    # CHECK because many repeated obs to a state
    # but maybe this makes sense since there are more 'counties' so adds weight
    X1 = nest$V1[,1]
    Y1 = nest$V1[,2]
    model1 = lm(Y1~X1)
    X2 = nest$V2[,1]
    Y2 = nest$V2[,2]
    model2 = lm(Y2~X2)
    betas = c(model2$coefficients[2], model1$coefficients[2])
  }
  if (decomposition == 'spectral'){
    spec = spectral_decomp(cbind(X,Y),A)
    Xstar = spec[,1]
    Ystar = spec[,2]
    # TO DO: WHAT IS THE BEST WAY TO Do THIS??
    # For now just split up the spectral r.v.s into 5 discrete scales
    # so that I have multiple observations per scale...
    betas = c()
    for (i in 1:5){ # 125
      model = lm(Ystar[(125*(i-1)):(125*i)] ~ Xstar[(125*(i-1)):(125*i)])
      betas = c(betas, model$coefficients[2])
    }
  }
  return(betas)
}


