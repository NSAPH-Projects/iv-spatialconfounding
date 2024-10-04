# Load in all utility functions
source('simfuncs.R')

# Load in simulation data,
# a list of a.vals (list of vectors), latnorm (vector), longnorm (vector), adjmat (matrix)
load('sim.RData')

# Load subspaces of confounding
load('Vcs.RData')
load('Vucs.RData')

library(spconf)
library(npcausal)
library(mgcv)
library(splines)
library(Matrix)

set.seed(123)
args = commandArgs(trailingOnly = TRUE)
nsims = as.integer(args[1])
projmatname = args[2]
option = args[3]
method = args[4]
adjustmentmatname = args[5]

print(c(nsims, projmatname, option, method, adjustmentmatname))

Vc = Vcs[[projmatname]]
Vuc = Vucs[[projmatname]]
a.vals = simlist$a.vals[[projmatname]]
 # setwd('/n/home07/swoodward/simulation')
if (adjustmentmatname == '') {
  simfunc(
    nsims = nsims,
    a.vals = a.vals,
    latnorm = simlist$latnorm,
    longnorm = simlist$longnorm,
    projmatname = projmatname,
    option = option,
    method = method,
    adjmat = simlist$adjmat,
    Vc = Vc,
    Vuc = Vuc
  )
} else{
  # Load basis matrix for conf adjustment for sensitivity/nonoracle
  load('adjustmentmats.RData')
  adjustmentmat = adjustmentmats[[adjustmentmatname]]
  simfunc(
    nsims,
    a.vals = a.vals,
    latnorm = simlist$latnorm,
    longnorm = simlist$longnorm,
    option = option,
    projmatname = projmatname,
    method = method,
    adjustmentmatname = adjustmentmatname,
    adjustmentmat = adjustmentmat,
    Vc = Vc,
    Vuc = Vuc
  )
}

  


