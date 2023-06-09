library(MASS)
library(stringr)
library(igraph)

nested_decomp = function(groups){
  # Create the decomposition matrices
  n = nrow(groups)
  l = ncol(groups) + 1

  # create projection matrices
  projs = list()
  for (i in 1:(l-1)){
    nst = as.data.frame(table(groups[,i]))
    colnames(nst) = c('group', 'frequency')
    G = matrix(NA, nrow = n, ncol = n)
    for (j in 1:nrow(G)){
      st = groups[j,i]
      stfac  = 1/nst$frequency[nst$group == st]
      for (k in 1:ncol(G)){
        G[j,k] = stfac*(st == groups[k,i])
      }
    }
    
    if (i == 1){
      projs[[i]] = G
    }
    else{
      projs[[i]] = G-projs[[i-1]]
    }
  }
  projs[[l]] = diag(n)-projs[[l-1]]
  
  return(projs)
}

spectral_decomp = function(A,inv = F){
  R = diag(rowSums(A))-A # precision of ICAR
  E = eigen(R) # eigen component
  D = E$val
  G = E$vec
  rm(E,R)
  if (inv){
    invG = solve(t(G)) 
    return(invG)
  }
  return(t(G))
}

sim = function(n,
               l=2, # levels of nested decomp
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
  print('Creating coordinates and groups')
  g = make_lattice(c(n^l,n^l))
  coords = layout_on_grid(g)
=  df = cbind.data.frame(xcoord = coords[,1], ycoord = coords[,2])
  for (i in 1:(l-1)){
    df[,(i+2)] = as.numeric(paste(str_pad((df$xcoord)%/%(n^(l-i)),2,pad = '0'),
                                  str_pad((df$ycoord)%/%(n^(l-i)),2,pad = '0'),
                                  sep = ''))
  }
  df = df[order(df[,(l+1)]),] # ordering by finest grid level should order by the coarser grids too

  # Create adjacency matrix
  print('Creating adjacency')
  A = as.matrix(as_adjacency_matrix(g))
  
  decomp = match.arg(decomposition)
  outcomemodel = match.arg(outcome)
  # Simulate X and Z
  print('Creating X and Z')
  if (decomp == 'nested'){ 
    stopifnot(length(rhox)==l)
    # there must be a faster way to do this
    Xmat = matrix(NA, ncol = l, nrow = n^(2*l))
    Zmat = matrix(NA, ncol = l, nrow = n^(2*l))
    for (i in 1:l){
      xz = mvrnorm(n^(2*i), mu = c(0,0), 
                   Sigma = matrix(c(1,rhox[i],rhox[i],1), 
                                  nrow = 2, ncol = 2))
      Xi = xz[,1]
      Zi = xz[,2]
      Xmat[,i] = rep(Xi,each = n^(2*(l-i)))
      Zmat[,i] = rep(Zi, each = n^(2*(l-i)))
    }
    # important here to check the order of the df
    df$X = rowSums(Xmat)
    df$Z = rowSums(Zmat)
  }
  
  if (decomp == 'spectral'){
    stopifnot(length(rhox)==n^(2*l))
    Xstar = rep(NA, n^(2*l))
    Zstar = rep(NA, n^(2*l))
    for (i in 1:(n^(2*l))){
      xz = mvrnorm(1, mu = c(0,0), 
                   Sigma = matrix(c(1,rhox[i],rhox[i],1), 
                                  nrow = 2, ncol = 2))
      Xstar[i] = xz[1]
      Zstar[i] = xz[2]
    }
    # Project into spatial domain
    spec = spectral_decomp(A,inv=T)
    df$X = spec %*% Xstar # n^4 length vector
    df$Z = spec %*% Zstar
  }
  # Simulate the outcome
  print('Simulating outcome')
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
    df$Y = betax*df$X + betaz*df$Z + betaxz*df$X*df$Z + rnorm(length(df$X), 
                                                           mean = 0, sd = sig)
  }
  return(list(
    'coord' = cbind(df$xcoord, df$ycoord),
    'A'=A,#adjacency mat
    'X'=df$X,#exposure
    'Y'=df$Y,#outcome
    'Z'=df$Z, #confounder
    'groups' = as.matrix(df[,3:(l+1)], nrow = nrow(df), ncol = ncol(df[,3:(l+1)])) #nested group
    ))
}

analysis = function(n, # subgroups in a group
                    A, # adjacency
                    X, # exposure
                    Y, # outcome
                    Z, # confounder
                    groups, # nested group
                    decomposition = c('spectral', 'nested')
                    ){
  decomposition = match.arg(decomposition)
  l = ncol(groups) + 1
  if (decomposition == 'nested'){
    # Regress Y on X at each level
    nest = nested_decomp(groups)
    # CHECK because many repeated obs to a state
    # but maybe this makes sense since there are more 'counties' so adds weight
    betas = rep(NA, l)
    for (i in 1:l){ 
      Xi = nest[[i]] %*% X
      Yi = nest[[i]] %*% Y
      modeli = lm(Yi~Xi)
      betas[i] = modeli$coefficients[2]
    }
    betas = rev(betas) # flip order since want small scale to big
  }
  if (decomposition == 'spectral'){
    spec = spectral_decomp(A)
    Xstar = spec %*% X
    Ystar = spec %*% Y
    # TO DO: WHAT IS THE BEST WAY TO Do THIS??
    # For now just split up the spectral r.v.s into discrete scales
    # so that I have multiple observations per scale...
    betas = c()
    num = n^(2*l-1)
    for (i in 1:n){ # 125
      model = lm(Ystar[(num*(i-1)):(num*i)] ~ Xstar[(num*(i-1)):(num*i)])
      betas = c(betas, model$coefficients[2])
    }
  }
  return(betas)
}


