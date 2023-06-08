
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
    R = diag(rowSums(A))-A # precision of ICAR
    E = eigen(R) # eigen component
    D = E$val
    G = E$vec
    invG = solve(t(G)) 
    rm(E,R)
    # Project into spatial domain
    df$X = invG%*%Xstar # n^4 length vector
    df$Z = invG%*%Zstar
  }
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
  return(list('A'=A,#adjacency mat
         'X'=df$X,#exposure
         'Y'=df$Y,#outcome
         'Z'=df$Z #confounder
         ))
}

# TO DO
analysis = function(A,
                    X,
                    Y,
                    Z,
                    decomposition = c('spectral', 'nested')
                    ){
  return
}