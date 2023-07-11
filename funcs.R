library(MASS)
library(stringr)
library(igraph)
library(ggplot2)
library(gridExtra)
library(abind)

nested_decomp_mats <- function(groups, append_identity = TRUE) {
  if (!is.matrix(groups)) {
    groups <- as.matrix(groups, ncol = 1)
  }
  if (append_identity) {
    groups <- cbind(groups, 1:nrow(groups))
  }
  n <- nrow(groups)
  L <- ncol(groups)

  # create averaging matrices
  # TODO: currently using dense format
  #   but this is inefficient, more efficient to perform
  #   the decompositions by recursive groupby/tapply ops
  proj_mats <- list()
  decomp_mats <- list()
  for (l in 1:L) {
    P_l <- matrix(0, n, n)
    for (group in unique(groups[, l])) {
      ix <- which(groups[, l] == group)
      P_l[ix, ix] <- 1 / length(ix)
    }
    proj_mats[[l]] <- P_l
    if (l > 1) {
      decomp_mats[[l]] <- P_l - proj_mats[[l - 1]]
    } else {
      decomp_mats[[l]] <- P_l
    }
  }
  return(list(
    proj_mats = proj_mats,
    decomp_mats = decomp_mats
  ))
}

nested_decomp_fast <- function(x, groups) {
  if (!is.matrix(groups)) {
    groups <- as.matrix(groups, ncol = 1)
  }
  n <- nrow(groups)
  L <- ncol(groups)
  # first compute group averages using t apply
  avs <- matrix(0, n, L)
  decomp <- matrix(0, n, L)
  for (l in 1:L) {
    avs[, l] <- tapply(x, groups[, l], mean)[groups[, l]]
    if (l == 1) {
      decomp[, l] <- avs[, l]
    } else {
      decomp[, l] <- avs[, l] - avs[, l - 1]
    }
  }
  return(list(
    avs = avs,
    decomp = decomp
  ))
}


nested_decomp <- function(groups, append_identity = TRUE) {
  if (!is.matrix(groups)) {
    groups <- as.matrix(groups, ncol = 1)
  }
  if (append_identity) {
    groups <- cbind(groups, 1:nrow(groups))
  }
  n <- nrow(groups)
  L <- ncol(groups)

  # decomposition matrices
  decomp <- nested_decomp_mats(groups, append_identity = FALSE)

  # sequentially collect the orthogonal components
  # from the orthogonal decomposition matrices

  curr_ix <- 0
  prev_unique_values <- 0
  basis <- matrix(0, n, n)

  # list of unique values at each level
  num_unique_values <- numeric(L)
  for (l in seq_len(L - 1)) {
    num_unique_values[l] <- length(unique(groups[, l]))
  }
  num_unique_values[L] <- n

  for (l in seq_len(L)) {
    dmat_l <- decomp$decomp_mats[[l]]

    # get rank of orthogonal projection matrix
    rank <- num_unique_values[l] - prev_unique_values

    # obtain orthonormal basis of image
    eig <- eigen(dmat_l, symmetric = TRUE)
    basis_idx <- (curr_ix + 1):(curr_ix + rank)
    basis[, basis_idx] <- eig$vectors[, seq_len(rank)]

    # update recursion
    prev_unique_values <- num_unique_values[l]
    curr_ix <- curr_ix + rank
  }

  return(basis)
}

get_projection <- function(x, basis, k) {
  xstar <- t(basis) %*% x
  n <- length(xstar)
  if (k < n) {
    xstar[(k + 1):n] <- 0
  }
  xproj <- basis %*% xstar
  return(xproj)
}


spectral_decomp <- function(A, inv = FALSE) {
  R <- diag(rowSums(A)) - A # graph laplacian (ICAR precision)
  E <- eigen(R) # eigen component
  D <- E$val
  G <- E$vec

  if (inv) {
    return(G)
  }
  return(t(G))
}

make_coords_df = function(n, l, quiet=FALSE) {
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

  # Create adjacency matrix
  if (!quiet){
    print('Creating adjacency')
  }
  A = as.matrix(as_adjacency_matrix(g))
  return (list(
    df=df,
    adjacency_mat=A
  ))
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
               nest = NULL,
               spec = NULL,
               Z = NULL,
               G = NULL
){
  decomp = match.arg(decomposition)
  outcome = match.arg(outcome)
  
  if (is.null(G)){
    lattice = make_coords_df(n, l, quiet = quiet)
    df = lattice$df
    groups = as.matrix(df[,3:(l+1)], nrow = nrow(df), ncol = ncol(df[,3:(l+1)]))
    A = lattice$adjacency_mat
    N = n^(2*l)
  }
  
  else{
    # TO DO
    return(NULL)
  }
  
  # Simulate X and Z
  if (!quiet){
    print('Creating X and Z')
  }
  # Create Z in the spatial domain
  if (distribution == 'exponential'){
    if (is.null(Z)){
      Z = rexp(N)-1
    }
    else{
      stopifnot(length(Z) == N)
    }
    noise = rexp(N)-1
  }
  if (distribution == 'gaussian'){
    if (is.null(Z)){
      Z = rnorm(N)
    }
    else{
      stopifnot(length(Z) == N)
    }
    noise = rnorm(N)
  }
  # Project to nested domain and get X
  if (decomp == 'nested'){ 
    if (is.null(truncate)){
      truncate = l
    }
    if (is.null(nest)){
      nest = nested_decomp_mats(groups)
    }
    Zmat = matrix(NA, nrow = N, ncol = l)
    Xmat = matrix(NA, nrow = N, ncol = l)
    
    for (i in 1:l){ 
      if (i > truncate){
        Zmat[,i] = rep(0, N)
        Xmat[,i] = nest$decomp_mats[[i]]%*%noise
      }
      else{
        Zi = nest$decomp_mats[[i]] %*% Z
        Zmat[,i] = Zi
        Xmat[,i] = rhox[i]*Zi + sqrt(1-rhox[i]^2)*(nest$decomp_mats[[i]] %*% noise)
      }
    }
    df$Z = rowSums(Zmat) 
    df$X = rowSums(Xmat)
  }
  # Project to spectral domain and get X
  if (decomp == 'spectral'){
    if (is.null(spec)){
      spec = spectral_decomp(A)
    }
    Zstar = spec %*% Z
    Xstar = rhox*Zstar + sqrt(1-rhox^2)*(spec %*% noise)
    if (!is.null(truncate)){
      Zstar = c(Zstar[1:truncate], rep(0,N-truncate))
      Xstar = c(Xstar[1:truncate], (spec %*% noise)[(truncate+1):N])
    }
    df$Z = t(spec) %*% Zstar
    df$X = t(spec) %*% Xstar
  }

  # Simulate the outcome
  if (!quiet){
    print('Simulating outcome')
  }
  if (outcome == 'linear'){
    df$Y = betax*df$X + betaz*df$Z + rnorm(N, mean = 0, 
                                           sd = sig)
  }
  if (outcome == 'quadratic'){
    stopifnot(length(betax)>=2)
    df$Y = betax[1]*df$X + betax[2]*df$X^2 + betaz*df$Z + rnorm(N, mean = 0, 
                                                                sd = sig)
  }
  if (outcome == 'interaction'){
    df$Y = betax*df$X + betaz*df$Z + betaxz*df$X*df$Z + rnorm(N, 
                                                              mean = 0, sd = sig)
  }
  return(list(
    'coord' = cbind(df$xcoord, df$ycoord),
    'A'=A,#adjacency mat
    'X'=df$X,#exposure
    'Y'=df$Y,#outcome
    'Z'=df$Z, #confounder
    'groups' = groups #nested group
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
      nest = nested_decomp_mats(groups)
    }
    # CHECK because many repeated obs to a state
    # but maybe this makes sense since there are more 'counties' so adds weight
    if (!quiet){
      print('Calculate betahats')
    }
    betas = c() 
    for (i in 1:l){ 
      Xi = nest$decomp_mats[[i]] %*% X
      Yi = nest$decomp_mats[[i]] %*% Y
      if (outcome == 'linear'){
        modeli = lm(Yi~Xi)
        betas = c(betas,modeli$coefficients[2])
      }
      if (outcome == 'quadratic'){
        Xi2 = nest$decomp_mats[[i]] %*% (X^2)
        modeli = lm(Yi~ Xi + Xi2)
        betas = cbind(betas, modeli$coefficients[2:3])
      }
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
    Xstar2 = spec %*% (X^2) # only used if model is quadratic
    Ystar = spec %*% Y
    if (!quiet){
      print('Calculate betahats')
    }
    betas = c()
    if (spectralmethod == 'bin'){
      num = n^(2*l-2) # -1
      zeroeig = which(apply(spec, 1, function(x) length(unique(round(x,10)))) == 1)
      for (i in 1:(n^2)){ # n
        window = (num*(i-1) + 1):(num*i)
        # Take out observations corresponding to 0 eigenvalue
        window = setdiff(window, zeroeig)
        
        if (outcome == 'linear'){
          model = lm(Ystar[window] ~ Xstar[window])
          betas = c(betas, model$coefficients[2])
        }
        if (outcome == 'quadratic'){
          model = lm(Ystar[window] ~ Xstar[window] + Xstar2[window])
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
          model = lm(Ystar[window]~Xstar[window] + Xstar2[window],weights = wt)
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
      nest = nested_decomp_mats(groups)
    }
    
    # CHECK because many repeated obs to a state
    # but maybe this makes sense since there are more 'counties' so adds weight
    if (!quiet){
      print('Calculate correlations')
    }
    cors = c() # CHANGED from rep(NA,)
    for (i in 1:l){ 
      Xi = nest$decomp_mats[[i]] %*% X
      Zi = nest$decomp_mats[[i]] %*% Z
      cors = c(cors,cor(Xi,Zi))
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
    Zstar = spec %*% Z
    if (!quiet){
      print('Calculate correlations')
    }
    if (spectralmethod == 'bin'){
      cors = c()
      num = n^(2*l-2) # -1
      zeroeig = which(apply(spec, 1, function(x) length(unique(round(x,10)))) == 1)
      for (i in 1:(n^2)){ # n
        window = (num*(i-1)+1):(num*i)
        # Take out observations corresponding to eigenval 0
        window = setdiff(window, zeroeig)
        cors = c(cors, cor(Zstar[window],Xstar[window]))
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
        cors = c(cors, cov.wt(cbind(Xstar[window],Zstar[window]),wt=wt,cor = T)$cor[1,2])
      }
    }
    
  }
  return(cors)
}

simfunc = function(nsims=100, 
                   n=5, 
                   l=2,
                   betax = 2,
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
                   spectralmethod = 'bin',
                   Z = NULL
){
  if (outcome == 'quadratic'){ # can generalize to other models thru arg coeffs
    deg = 2
    nestedmats = array(NA, c(nsims, l, deg))
    if (spectralmethod == 'bin'){
      spectralmats = array(NA, c(nsims, n^2,deg))
    }
    if (spectralmethod == 'wls'){
      spectralmats = array(NA, c(nsims, n^(2*l),deg))
    }
    for (num in 1:nsims){
      outsim = sim(n = n, 
                   l=l,
                   betax = betax,
                   outcome = outcome, 
                   rhox = rhox,
                   decomposition = dgm, 
                   truncate = truncate,
                   quiet = quiet,
                   distribution = distribution,
                   betaxz=betaxz,
                   spec = spec,
                   Z = Z)
      nestedout = analysis(n = n, 
                           A = outsim$A, 
                           X = outsim$X, 
                           Y = outsim$Y, 
                           groups = outsim$groups,
                           decomposition = 'nested', 
                           outcome = outcome,
                           quiet = quiet,
                           nest = nest)
      spectralout = analysis(n = n, 
                             A = outsim$A, 
                             X = outsim$X, 
                             Y = outsim$Y, 
                             groups = outsim$groups,
                             decomposition = 'spectral', 
                             outcome = outcome,
                             quiet = quiet, 
                             spec = spec,
                             spectralmethod = spectralmethod)
      nestedmats[num,,] = t(nestedout)
      spectralmats[num,,] = t(spectralout)
    }
    return(list('nestedmats' = nestedmats, 'spectralmats' = spectralmats))
  }
  
  else{
    nestedmat = matrix(NA, nrow = nsims, ncol = l)
    if (spectralmethod == 'bin'){
      spectralmat = matrix(NA, nrow = nsims, ncol = n^2)
    }
    if (spectralmethod == 'wls'){
      spectralmat = matrix(NA, nrow = nsims, ncol = n^(2*l))
    }
    for (num in 1:nsims){
      outsim = sim(n = n, 
                   l=l,
                   betax = betax,
                   outcome = outcome, 
                   rhox = rhox,
                   decomposition = dgm, 
                   truncate = truncate,
                   quiet = quiet,
                   distribution = distribution,
                   betaxz=betaxz,
                   nest = nest,
                   spec = spec,
                   Z = Z)
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
}

bineigen = function(D, bins){
  speceigen = rep(NA, bins)
  num = floor(length(D)/bins)
  for (i in 1:bins){
    if (i == bins){
      window = (num*(i-1) + 1):length(D)
    }
    else{
      window = (num*(i-1) + 1):(num*i)
    }
    speceigen[i] = mean(D[window])
  }
  return(speceigen)
}

plotfunc = function(n=5,
                    l=2,
                    nestedmat,
                    spectralmat,
                    ylab = 'bias',
                    hline = 2,
                    ylim = c(-2,2),
                    col='blue',
                    mains = c('Nested', 'Spectral'),
                    D = NULL
                    ){
  nested_df <- data.frame(Spatial_Scale = 1:ncol(nestedmat), betas = colMeans(nestedmat))
  
  if (is.null(D)){
    speceigen = 1:ncol(spectralmat)
  }
  else{
    speceigen = bineigen(D, bins = ncol(spectralmat))
  }
  spectral_df <- data.frame(Spatial_Scale = speceigen, betas = colMeans(spectralmat))
  # Plot for nestedmat
  nested_df$Custom_Labels <- factor(nested_df$Spatial_Scale, levels = c(1, 2, 3),
                                    labels = c('9x9', '3x3', '1x1'))
  plot_nested <- ggplot(nested_df, aes(x = Custom_Labels, y = betas-hline)) + 
    geom_ribbon(aes(x = 1:3, ymin = apply(nestedmat, 2, quantile, probs = 0.25)-hline, 
                    ymax = apply(nestedmat, 2, quantile, probs = 0.75)-hline), 
                fill = "lightblue") +
    geom_line(aes(y = betas-hline, x = Spatial_Scale),color = col) +
    geom_point(color = col) +
    xlab('Level') +
    ylab(ylab) +
    ggtitle(mains[1]) +
    ylim(ylim[1],ylim[2]) +
    theme_minimal() +
    theme(legend.position = 'topright') +
    geom_hline(yintercept = 0, color = 'red', linetype = "dashed") + 
    scale_x_discrete(labels = c('9x9', '3x3', '1x1'))
  
  # Plot for spectralmat
  plot_spectral <- ggplot(spectral_df, aes(x = Spatial_Scale, y = betas-hline)) + 
    geom_ribbon(aes(ymin = apply(spectralmat, 2, quantile, probs = 0.25)-hline, 
                    ymax = apply(spectralmat, 2, quantile, probs = 0.75)-hline), 
                fill = "lightblue") +
    geom_line(color = col) +
    xlab('Eigenvalue of Laplacian') +
    ylab(ylab) +
    ggtitle(mains[2]) +
    ylim(ylim[1],ylim[2]) +
    theme_minimal() +
    theme(legend.position = 'topright') +
    geom_hline(yintercept = 0, color = 'red', linetype = "dashed") 
  
  # Arrange and center the plots
  combined_plots = grid.arrange(plot_nested, plot_spectral, nrow = 1)
  
  return(combined_plots)
}
