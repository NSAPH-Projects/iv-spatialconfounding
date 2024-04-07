library(spconf) 
library(npcausal) 
library(mgcv) 
library(eCAR) 
#library(INLA) # ISSUE HERE
library(GpGp) 
#library(GPCERF) # ISSUE HERE

# Load in all utility functions
source('simfuncs.R')

# Load in simulation data, 
# a list of a.vals (list of vectors), latnorm (vector), longnorm (vector), adjmat (matrix)
load('sim.RData')

# Load in list of projection matrices of confounding
load('projmats.RData')

# Load basis matrix for conf adjustment for sensitivity/nonoracle
load('Esub.RData')

set.seed(123)
args = commandArgs(trailingOnly = TRUE)
nsims = as.integer(args[1])
projmatname = args[2]
option = args[3]
method = args[4]
sensitivity = args[5]

projmat = projmats[[projmatname]]
a.vals = simlist$a.vals[[projmatname]]

if (sensitivity == 'no'){
  simfunc(nsims,
          a.vals = a.vals,
          latnorm = simlist$latnorm,
          longnorm = simlist$longnorm,
          projmat = projmat,
          projmatname = projmatname,
          option = option,
          method = method,
          adjmat = simlist$adjmat)
}

if (sensitivity == 'yes'){
  simfunc_nonoracle(nsims,
                    a.vals = a.vals,
                    latnorm = simlist$latnorm,
                    longnorm = simlist$longnorm,
                    option = option,
                    projmat = projmat,
                    projmatname = projmatname,
                    Esub,
                    method = method)
}


