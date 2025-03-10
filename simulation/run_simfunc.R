# Load in all utility functions
source('../funcs.R')
library(dplyr)
library(geosphere)

# Load in simulation data,
# a list of lat (vector), lon (vector), GFT (matrix), statemat (matrix)
load('sim.RData')

set.seed(123)
args <- commandArgs(trailingOnly = TRUE)
nsims <- as.integer(args[1])
rangeu <- args[2]
option <- args[3]

print(c(nsims, rangeu, option))

simfunc(nsims,
        simlist$lat,
        simlist$lon,
        rangeu = rangeu,
        option = option,
        GFT_conf = simlist$GFT_conf,
        statemat = simlist$statemat,
        within_state_GP = F, # change to T for confounding scenario 3
        bootstrap_conf = T
        )


  


