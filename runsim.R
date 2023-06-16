library(MASS)
library(stringr)
library(igraph)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(dplyr)
library(scatterplot3d) # update later with scatter3d (interactive)
set.seed(2)

source('funcs.R')

# I first compute nest and spec in this chunk so they don't have to be recomputed
outsim = sim(
  n = 5,
  outcome = 'linear',
  rhox = c(0.9, 0.001),
  decomposition = 'nested',
  quiet = T
)
analnested = analysis(
  n = 5,
  A = outsim$A,
  X = outsim$X,
  Y = outsim$Y,
  groups = outsim$groups,
  decomposition = 'nested',
  quiet = T,
  return_decomps = T
)
analspectral = analysis(
  n = 5,
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
  rhox = c(0.9, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  spectralmethod = 'wls'
)

# Sim 2
out[[2]] = simfunc(
  rhox = c(0.9, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  truncate = 1,
  spectralmethod = 'wls'
)

# Sim 3
out[[3]] = simfunc(
  rhox = seq(0, 1, by = 1 / (5 ^ 4 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  spectralmethod = 'wls'
)

# Sim 4
out[[4]] = simfunc(
  rhox = c(rep(0, 525), seq(1 / 100, 1, by = 1 / 100)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  spectralmethod = 'wls'
)

# Sim 5
out[[5]] = simfunc(
  rhox = c(0.9, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  spectralmethod = 'wls'
)

# Sim 6
out[[6]] = simfunc(
  rhox = c(0.9, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  truncate = 1,
  spectralmethod = 'wls'
)

# Sim 7
out[[7]] = simfunc(
  rhox = seq(0, 1, by = 1 / (5 ^ 4 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  spectralmethod = 'wls'
)

# Sim 8
out[[8]] = simfunc(
  rhox = seq(0, 1, by = 1 / (5 ^ 4 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  objective = 'coherence',
  truncate = 312,
  spectralmethod = 'wls'
)

# Sim 9
out[[9]] = simfunc(
  rhox = c(0.001, 0.9),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  spectralmethod = 'wls'
)

# Sim 10
out[[10]] = simfunc(
  rhox = seq(1, 0, by = -1 / (5 ^ 4 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  spectralmethod = 'wls'
)

# Sim 11
out[[11]] = simfunc(
  rhox = c(0.9, 0.001),
  dgm = 'nested',
  nest = nest,
  spec = spec,
  outcome = 'interaction',
  betaxz = -1,
  spectralmethod = 'wls'
)

# Sim 12
out[[12]] = simfunc(
  rhox = seq(0, 1, by = 1 / (5 ^ 4 - 1)),
  dgm = 'spectral',
  nest = nest,
  spec = spec,
  outcome = 'interaction',
  betaxz = -1,
  spectralmethod = 'wls'
)

save(out, file = 'out.Rdata')