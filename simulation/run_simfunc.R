# Load in all utility functions
source('../funcs.R')
library(dplyr)
library(geosphere)
library(foreach)
library(doParallel)
library(Matrix)

# Load in simulation data: lat, lon, B_tps_full, B_gl_full, statemat, W
load('sim.RData')

set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
nsims                 <- as.integer(args[1])
confounding_mechanism <- args[2]
option                <- args[3]
select_basis          <- if (length(args) >= 4) as.logical(args[4]) else TRUE
n_cores               <- if (length(args) >= 5) as.integer(args[5]) else 1L

print(c(nsims, confounding_mechanism, option, select_basis, n_cores))

simfunc(nsims,
        simlist$lat,
        simlist$lon,
        confounding_mechanism = confounding_mechanism,
        option                = option,
        B_tps_full            = simlist$B_tps_full,
        B_gl_full             = simlist$B_gl_full,
        statemat              = simlist$statemat,
        W                     = simlist$W,
        select_basis          = select_basis,
        n_cores               = n_cores
        )


  


