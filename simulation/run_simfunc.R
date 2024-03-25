library(devtools) 
library(spconf) 
library(npcausal) 
library(mgcv) 
library(eCAR) 
library(INLA) 
library(GpGp) 
library(GPCERF) 

# Load in all utility functions
source('simfuncs.R')

# Load in simulation data, 
# a list of a.vals (vector), latnorm (vector), longnorm (vector), adjmat (matrix)
load('sim.RData')

# Load in list of projection matrices of confounding
load('projmats.RData')

set.seed(123)
args = commandArgs(trailingOnly = TRUE)
nsims = as.integer(args[1])
projmatname = args[2]
option = args[3]
method = args[4]

projmat = projmats[[projmatname]]

simfunc(nsims,
        a.vals = simlist$a.vals,
        latnorm = simlist$latnorm,
        longnorm = simlist$longnorm,
        projmat = projmat,
        projmatname = projmatname,
        option = option,
        method = method,
        adjmat = simlist$adjmat)

