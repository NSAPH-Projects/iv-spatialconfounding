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
      projs[[i]] = G-Reduce("+", projs[1:(i-1)]) # sum up the prev As
    }
  }
  projs[[l]] = diag(n)-Reduce("+", projs[1:(l-1)]) 
  
  return(projs)
}

spectral_decomp = function(A){
  R = diag(rowSums(A))-A # precision of ICAR
  E = eigen(R) # eigen component
  D = E$val
  G = E$vec
  rm(E,R)
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
               distribution = 'exponential', # distribution of X,
               truncate = NULL, # after this spatial level Z will not vary
               quiet = F,
               spec = NULL
               ){
  decomp = match.arg(decomposition)
  outcome = match.arg(outcome)

  # Create coordinates
  if (!quiet){
    print('Creating coordinates and groups')
  }
  g = make_lattice(c(n^l,n^l))
  coords = layout_on_grid(g)
  df = cbind.data.frame(xcoord = coords[,1], ycoord = coords[,2])
  for (i in 1:(l-1)){
    df[,(i+2)] = as.numeric(paste(str_pad((df$xcoord)%/%(n^(l-i)),2,pad = '0'),
                                  str_pad((df$ycoord)%/%(n^(l-i)),2,pad = '0'),
                                  sep = ''))
  }
  df = df[order(df[,(l+1)]),] # ordering by finest grid level should order by the coarser grids too

  # Create adjacency matrix
  if (!quiet){
    print('Creating adjacency')
  }
  A = as.matrix(as_adjacency_matrix(g))
  
  # Simulate X and Z
  if (!quiet){
    print('Creating X and Z')
  }
  if (decomp == 'nested'){ 
    if (is.null(truncate)){ # arg for stopping variation in Z after certain spatial scale
      truncate = l
    }
    stopifnot(length(rhox)>=min(l,truncate))
    # there must be a faster way to do this
    X = rep(0, n^(2*l))#matrix(NA, ncol = l, nrow = n^(2*l))
    Z = rep(0, n^(2*l))#matrix(NA, ncol = l, nrow = n^(2*l))
    
    for (i in 1:l){
      if (i > truncate){
        if (distribution == 'gaussian'){
          Xi = rnorm(n^(2*i), mean = 0, sd = 1)
        }
        if (distribution == 'exponential'){
          Xi = rexp(n^(2*i))-1
        }
        Zi = rep(0, n^(2*i)) 
      }
      else{
        if (distribution == 'gaussian'){
          xz = mvrnorm(n^(2*i), mu = c(0,0), 
                       Sigma = matrix(c(1,rhox[i],rhox[i],1), 
                                      nrow = 2, ncol = 2))
          Xi = xz[,1]
          Zi = xz[,2]
        }
        if (distribution == 'exponential'){
          Xi = rexp(n^(2*i))-1
          Zi = rhox[i]*Xi + sqrt(1-rhox[i]^2)*(rexp(n^(2*i))-1)
        }
      }
      X = X + rep(Xi,each = n^(2*(l-i)))
      Z = Z + rep(Zi, each = n^(2*(l-i)))
    }
    df$X = X
    df$Z = Z
  }
  
  if (decomp == 'spectral'){
    if (is.null(truncate)){ # arg for stopping variation in Z after certain spatial scale
      truncate = n^(2*l)
    }
    stopifnot(length(rhox)>=min(n^(2*l), truncate))
    Xstar = rep(NA, n^(2*l))
    Zstar = rep(NA, n^(2*l))
    if (distribution == 'gaussian'){
      for (i in 1:(n^(2*l))){
        if (i > truncate){
          if (distribution == 'gaussian'){
            Xstar[i] = rnorm(1, mean = 0, sd = 1)
          }
          if (distribution == 'exponential'){
            Xstar[i] = rexp(1)-1
          }
          Zstar[i] = 0
        }
        else{
          xz = mvrnorm(1, mu = c(0,0), 
                       Sigma = matrix(c(1,rhox[i],rhox[i],1), 
                                      nrow = 2, ncol = 2))
          Xstar[i] = xz[1]
          Zstar[i] = xz[2]
        }
      }
    }
    if (distribution == 'exponential'){
      Xstar = rexp(n^(2*l))-1
      Zstar = rhox*Xstar + sqrt(1-rhox^2)*(rexp(n^(2*l))-1)
    }
    # Project into spatial domain
    if (is.null(spec)){
      specinv = t(spectral_decomp(A))
    }
    else{
      specinv = t(spec)
    }
    df$X = specinv %*% Xstar # n^4 length vector
    df$Z = specinv %*% Zstar
  }
  # Simulate the outcome
  if (!quiet){
    print('Simulating outcome')
  }
  if (outcome == 'linear'){
    df$Y = betax*df$X + betaz*df$Z + rnorm(length(df$X), mean = 0, 
                                  sd = sig)
  }
  if (outcome == 'quadratic'){
    stopifnot(length(betax)>=2)
    df$Y = betax[1]*df$X + betax[2]*df$X^2 + betaz*df$Z + rnorm(length(df$X), mean = 0, 
                                                    sd = sig)
  }
  if (outcome == 'interaction'){
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
                    groups, # nested group
                    decomposition = c('spectral', 'nested'),
                    outcome = 'linear',
                    quiet = F,
                    nest = NULL,
                    spec = NULL,
                    return_decomps = F,
                    spectralmethod = 'bin'
                    ){
  decomposition = match.arg(decomposition)
  l = ncol(groups) + 1
  if (decomposition == 'nested'){
    # Regress Y on X at each level
    
    if (is.null(nest)){
      if (!quiet){
        print('Perform decomposition')
      }
      nest = nested_decomp(groups)
    }
    # CHECK because many repeated obs to a state
    # but maybe this makes sense since there are more 'counties' so adds weight
    if (!quiet){
      print('Calculate betahats')
    }
    betas = c() # CHANGED from rep(NA,)
    for (i in 1:l){ 
      Xi = nest[[i]] %*% X
      Yi = nest[[i]] %*% Y
      if (outcome == 'linear'){
        modeli = lm(Yi~Xi)
        betas = c(betas,modeli$coefficients[2])
      }
      if (outcome == 'quadratic'){
        modeli = lm(Yi~ poly(Xi,degree = 2))
        betas = cbind(betas, modeli$coefficients[2:3])
      }
    }
    #betas = rev(betas) # flip order since want small scale to big
    if (outcome == 'linear'){
      betas = rev(betas)
    }
    if (outcome == 'quadratic'){
      betas = betas[,ncol(betas):1]
    }
    if (return_decomps){
      return(list('betas' = betas, 'nest' = nest))
    }
  }
  if (decomposition == 'spectral'){
    if (is.null(spec)){
      if (!quiet){
        print('Perform decomposition')
      }
      spec = spectral_decomp(A)
    }
    Xstar = spec %*% X
    Ystar = spec %*% Y
    # TO DO: WHAT IS THE BEST WAY TO Do THIS??
    # For now just split up the spectral r.v.s into discrete scales
    # so that I have multiple observations per scale...
    if (!quiet){
      print('Calculate betahats')
    }
    betas = c()
    if (spectralmethod == 'bin'){
      num = n^(2*l-2) # -1
      for (i in 1:(n^2)){ # n
        if (outcome == 'linear'){
          model = lm(Ystar[(num*(i-1)):(num*i)] ~ Xstar[(num*(i-1)):(num*i)])
          betas = c(betas, model$coefficients[2])
        }
        if (outcome == 'quadratic'){
          model = lm(Ystar[(num*(i-1)):(num*i)] ~ poly(Xstar[(num*(i-1)):(num*i)], 2))
          betas = cbind(betas, model$coefficients[2:3])
        }
      }
    }
    if (spectralmethod == 'wls'){
      betas = c()
      for (i in 1:length(Xstar)){
        k = 10
        window = seq(max(1, i-k), min(i+k, length(Xstar)), by = 1)
        wt = exp(-0.1*abs(window-i))
        wt = wt/sum(wt)
        if (outcome == 'linear'){
          model = lm(Ystar[window]~Xstar[window], weights = wt)
          betas = c(betas, model$coefficients[2])
        }
        if (outcome == 'quadratic'){
          model = lm(Ystar[window]~poly(Xstar[window],2),weights = wt)
          betas = cbind(betas, model$coefficients[2:3])
        }
      }
    }
    
    if (return_decomps){
      return(list('betas' = betas, 'spec' = spec))
    }
  }
  return(betas)
}

coherence = function(n, # subgroups in a group
                    A, # adjacency
                    X, # exposure
                    Z, # confounder
                    groups, # nested group
                    decomposition = c('spectral', 'nested'),
                    quiet = F,
                    nest = NULL,
                    spec = NULL,
                    spectralmethod = 'bin'
  ){
  decomposition = match.arg(decomposition)
  l = ncol(groups) + 1
  if (decomposition == 'nested'){
    # Regress Y on X at each level
    if (is.null(nest)){
      if (!quiet){
        print('Perform decomposition')
      }
      nest = nested_decomp(groups)
    }
    
    # CHECK because many repeated obs to a state
    # but maybe this makes sense since there are more 'counties' so adds weight
    if (!quiet){
      print('Calculate correlations')
    }
    cors = c() # CHANGED from rep(NA,)
    for (i in 1:l){ 
      Xi = nest[[i]] %*% X
      Zi = nest[[i]] %*% Z
      cors = c(cors,cor(Xi,Zi))
    }
    cors = rev(cors)
  }
  if (decomposition == 'spectral'){
    if (is.null(spec)){
      if (!quiet){
        print('Perform decomposition')
      }
      spec = spectral_decomp(A)
    }
    Xstar = spec %*% X
    Zstar = spec %*% Z
    if (!quiet){
      print('Calculate correlations')
    }
    if (spectralmethod == 'bin'){
      cors = c()
      num = n^(2*l-2) # -1
      for (i in 1:(n^2)){ # n
        cors = c(cors, cor(Zstar[(num*(i-1)):(num*i)],Xstar[(num*(i-1)):(num*i)]))
      }
    }
    if (spectralmethod == 'wls'){
      # TO DO
      cors = c()
      k = 10
      for (i in 1:length(Xstar)){
        window = seq(max(1, i-k), min(i+k, length(Xstar)), by = 1)
        wt = exp(-0.1*abs(window-i))
        wt = wt/sum(wt)
        cors = c(cors, cov.wt(cbind(Xstar[window],Zstar[window]),wt=wt,cor = T)$cov[1,2])
      }
    }
    
  }
  return(cors)
}

simfunc = function(nsims=100, 
                   n=5, 
                   l=2,
                   outcome='linear', 
                   rhox, 
                   dgm, 
                   quiet=T,
                   objective='analysis',
                   nest=NULL, # assume given
                   spec=NULL, # assume given
                   truncate=NULL,
                   distribution='exponential',
                   betaxz=0,
                   spectralmethod = 'bin'
                   ){
  nestedmat = matrix(NA, nrow = nsims, ncol = l)
  if (spectralmethod == 'bin'){
    spectralmat = matrix(NA, nrow = nsims, ncol = n^2)
  }
  if (spectralmethod == 'wls'){
    spectralmat = matrix(NA, nrow = nsims, ncol = n^(2*l))
  }
  for (num in 1:nsims){
    outsim = sim(n = n, 
                 outcome = outcome, 
                 rhox = rhox,
                 decomposition = dgm, 
                 truncate = truncate,
                 quiet = quiet,
                 distribution = distribution,
                 betaxz=betaxz,
                 spec = spec)
    if (objective == 'analysis'){
      nestedmat[num,] = analysis(n = n, 
                                 A = outsim$A, 
                                 X = outsim$X, 
                                 Y = outsim$Y, 
                                 groups = outsim$groups,
                                 decomposition = 'nested', 
                                 quiet = quiet,
                                 nest = nest)
      spectralmat[num,] = analysis(n = n, 
                                   A = outsim$A, 
                                   X = outsim$X, 
                                   Y = outsim$Y, 
                                   groups = outsim$groups,
                                   decomposition = 'spectral', 
                                   quiet = quiet, 
                                   spec = spec,
                                   spectralmethod = spectralmethod)
    }
    if (objective == 'coherence'){
      nestedmat[num,] = coherence(n = n, 
                                  A = outsim$A, 
                                  X = outsim$X, 
                                  Z = outsim$Z, 
                                  groups = outsim$groups,
                                  decomposition = 'nested', 
                                  quiet = quiet, nest = nest)
      spectralmat[num,] = coherence(n = n, 
                                    A = outsim$A, 
                                    X = outsim$X, 
                                    Z = outsim$Z, 
                                    groups = outsim$groups,
                                    decomposition = 'spectral', 
                                    quiet = quiet, spec = spec,
                                    spectralmethod = spectralmethod)
    }
  }
  return(list('nestedmat' = nestedmat, 'spectralmat' = spectralmat))
}


plotfunc = function(n=5,
                    l=2,
                    nestedmat,
                    spectralmat,
                    ylab = 'betahat',
                    hline = 2,
                    ylim = c(0.5,3),
                    col='blue'
                    ){
  nested_df <- data.frame(Spatial_Scale = 1:ncol(nestedmat), cors = colMeans(nestedmat))
  spectral_df <- data.frame(Spatial_Scale = 1:ncol(spectralmat), cors = colMeans(spectralmat))
  
  # Plot for nestedmat
  plot_nested <- ggplot(nested_df, aes(x = Spatial_Scale, y = cors)) +
    geom_line(color = col) +
    geom_point(color = col) +
    xlab('Spatial Scale') +
    ylab(ylab) +
    ggtitle('Nested') +
    ylim(ylim[1],ylim[2]) +
    theme_minimal() +
    theme(legend.position = 'topright') +
    geom_hline(yintercept = hline, color = 'red', linetype = "dashed")
  
  # Plot for spectralmat
  plot_spectral <- ggplot(spectral_df, aes(x = Spatial_Scale, y = cors)) +
    geom_line(color = col) +
    xlab('Spatial Scale') +
    ylab(ylab) +
    ggtitle('Spectral') +
    ylim(ylim[1],ylim[2]) +
    theme_minimal() +
    theme(legend.position = 'topright') +
    geom_hline(yintercept = hline, color = 'red', linetype = "dashed")
  
  # Arrange and center the plots
  combined_plots = grid.arrange(plot_nested, plot_spectral, nrow = 1)
  
  return(combined_plots)
}

