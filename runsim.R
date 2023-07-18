library(MASS)
library(stringr)
library(igraph)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(sp)
library(raster)
library(sf)
library(fields)

set.seed(2)

source('funcs.R')

# I first compute nest and spec in this chunk so they don't have to be recomputed
outsim = sim(
  n = 3,
  l = 3,
  outcome = 'linear',
  rhox = c(0.9, 0.5, 0.001),
  decomposition = 'nested',
  quiet = T
)
# Save one nested DGM sim for exploration
save(outsim, file = 'outsim.Rdata') 

analnested = analysis(
  n = 3,
  A = outsim$A,
  X = outsim$X,
  Y = outsim$Y,
  groups = outsim$groups,
  decomposition = 'nested',
  quiet = T,
  return_decomps = T
)
analspectral = analysis(
  n = 3,
  A = outsim$A,
  X = outsim$X,
  Y = outsim$Y,
  groups = outsim$groups,
  decomposition = 'spectral',
  quiet = T,
  return_decomps = T
)
nest = analnested$nest
spec = analspectral$spec

out = list()
# Sim 1
out[[1]] = simfunc(
  n = 3,
  l = 3,
  rhox = c(0.9, 0.5, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  spectralmethod = 'bin'
)

# Sim 2
out[[2]] = simfunc(
  n = 3,
  l = 3,
  rhox = c(0.9, 0.5, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  truncate = 1,
  spectralmethod = 'bin'
)

# Sim 3
out[[3]] = simfunc(
  n = 3,
  l = 3,
  rhox = seq(0, 1, by = 1 / (3 ^ 6 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  spectralmethod = 'bin'
)

# Sim 4
out[[4]] = simfunc(
  n = 3,
  l = 3,
  rhox = c(rep(0, 3^6 - 100), seq(1 / 100, 1, by = 1 / 100)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  spectralmethod = 'bin'
)

# Sim 5
out[[5]] = simfunc(
  n = 3,
  l = 3,
  rhox = c(0.9, 0.5, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  spectralmethod = 'bin'
)

# Sim 6
out[[6]] = simfunc(
  n = 3,
  l = 3,
  rhox = c(0.9, 0.5, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  truncate = 1,
  spectralmethod = 'bin'
)

# Sim 7
out[[7]] = simfunc(
  n = 3,
  l = 3,
  rhox = seq(0, 1, by = 1 / (3 ^ 6 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  spectralmethod = 'bin'
)

# Sim 8
out[[8]] = simfunc(
  n = 3,
  l = 3,
  rhox = seq(0, 1, by = 1 / (3 ^ 6 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  truncate = 312,
  spectralmethod = 'bin'
)

# Sim 9
out[[9]] = simfunc(
  n = 3,
  l = 3,
  rhox = c(0.001, 0.5, 0.9),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  spectralmethod = 'bin'
)

# Sim 10
out[[10]] = simfunc(
  n = 3,
  l = 3,
  rhox = seq(1, 0, by = -1 / (3 ^ 6 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  spectralmethod = 'bin'
)

# Sim 11
out[[11]] = simfunc(
  n = 3,
  l = 3,
  rhox = c(0.9, 0.5, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  outcome = 'interaction',
  betaxz = -1,
  spectralmethod = 'bin'
)

# Sim 12
out[[12]] = simfunc(
  n = 3,
  l = 3,
  rhox = seq(0, 1, by = 1 / (3 ^ 6 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  outcome = 'interaction',
  betaxz = -1,
  spectralmethod = 'bin'
)

# Sim 13
out[[13]] = simfunc(
  n = 3,
  l = 3,
  outcome = 'quadratic',
  betax = c(2, 1),
  rhox = c(0.9, 0.5, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  quiet = T,
  spectralmethod = 'bin'
)

# Sim 14
out[[14]] = simfunc(
  n = 3,
  l = 3,  
  outcome = 'quadratic', 
  betax = c(2,1),
  rhox = seq(0,1, by = 1/(3^6-1)),
  dgm = 'spectral', 
  nest = nest,
  spec = spec,
  quiet = T,
  spectralmethod = 'bin'
)

# Sim 8
out[[15]] = simfunc(
  n = 3,
  l = 3,
  rhox = c(rep(0, 3^6 - 100), seq(1 / 100, 1, by = 1 / 100)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  spectralmethod = 'bin'
)

# Real data simulations
study = read.csv('/Users/sophie/Documents/SpatialConf/archived/Study_dataset_2010.csv')
A = read.csv("/Users/sophie/Documents/SpatialConf/archived/adjacency_matrix.csv",
               header = F) # created from spacebench script
study$FIPS = str_pad(study$FIPS, 5, pad = '0')

# Create Graph Laplacian
R = diag(rowSums(A)) - A # graph laplacian 
E = eigen(R) # eigen component
D = E$val
G = E$vec

Z = study$gmet_mean_tmmn
groups = cbind(study$region, study$STATE)
nest = nested_decomp_mats(groups)
spec = t(G)
out[[16]] = simfunc(
  A = A,
  Z = Z,
  groups = groups,
  rhox = c(0.9, 0.5, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  spectralmethod = 'bin'
)

out[[17]] = simfunc(
  A = A,
  Z = Z,
  groups = groups,
  rhox = seq(0, 1, by = 1 / (length(Z) - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  spectralmethod = 'bin'
)

out[[18]] = simfunc(
  A = A,
  Z = Z,
  groups = groups,
  rhox = c(0.9, 0.5, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  spectralmethod = 'bin'
)

out[[19]] = simfunc(
  A = A,
  Z = Z,
  groups = groups,
  rhox = seq(0, 1, by = 1 / (length(Z) - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  spectralmethod = 'bin'
)

save(out, file = 'out.Rdata')
