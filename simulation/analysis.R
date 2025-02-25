library(ggplot2)
library(sf)
library(tidyr)
library(mgcv)
library(dplyr)
library(utils)
library(xtable)
library(Matrix)
source('../funcs.R')
load('sim.RData')

# Calculate distance matrix
distmat <- geosphere::distm(cbind(simlist$lon, simlist$lat), 
                           fun = distHaversine)
# Standardize so approximately range (0,2)
distmat <- distmat/1000000

# read in results files

# Confounding scenarios 1-2
csvs_notwithinstate <- list.files('results_Mar1/', pattern = '.csv')
# Confounding scenario 3 (GP within state)
csvs_withinstate <- list.files('results_Mar1/within_state/', pattern = '.csv') 
csvs <- c(csvs_notwithinstate, csvs_withinstate)

# Create storage for metrics 
analysisdf <- data.frame(
  confounding_scenario = character(length(csvs)),
  rangeu = character(length(csvs)),
  option = character(length(csvs)),
  method = character(length(csvs)),
  bias = numeric(length(csvs)),
  RMSE = numeric(length(csvs)),
  se = numeric(length(csvs))
)

# Extract components from filenames
rangeu <- gsub("(_.*$)", "", csvs)
confounding_scenario <- c(ifelse(rangeu == '0.05', 1, 2)[1:length(csvs_notwithinstate)], 
                         rep(3, length(csvs_withinstate)))
option <- gsub("^(.*?_)(.*?)(_.*$)", "\\2", csvs)
method <- gsub(".*_(IV-[A-Za-z]+|[A-Za-z]+)\\.csv", "\\1", csvs)

method[method == 'spatialcoord'] <- 'Spatial coordinates'
method[method == 'baseline'] <- 'No confounding adjustment'

# Precompute true estimand for each outcome model and confounding mechanism
mutrues <- data.frame(expand.grid(rangeu = c(0.05, 0.1), 
                                 option = c('linear', 'nonlinear')))
mutrues$withinstate <- F
mutrues <- rbind(mutrues, data.frame(rangeu = 0.05, option = 'linear', withinstate = T))
mutrues <- rbind(mutrues, data.frame(rangeu = 0.05, option = 'nonlinear', withinstate = T))
mutrues$theta <- NA
mutrues$option <- as.character(mutrues$option)

for (i in 1:nrow(mutrues)){
  mutrues$theta[i] <- computemutrue(option = mutrues$option[i], 
                                     rangeu = mutrues$rangeu[i], 
                                    within_state_GP = mutrues$withinstate[i],
                                    distmat = distmat,
                                    statemat = simlist$statemat)
}
                
gs <- list()
# Loop through results to calculate ERF metrics and create plots.
for (i in 1:length(csvs)){
  filename <- csvs[i]
  print(filename)
  
  analysisdf$confounding_scenario[i] <- confounding_scenario[i]
  analysisdf$rangeu[i] <- rangeu[i]
  analysisdf$option[i] <- option[i]
  analysisdf$method[i] <- method[i]
  if (confounding_scenario[i] !=3){
    df_temp <- read.csv(file.path('results_Mar1/', filename))
  }
  else{
    df_temp <- read.csv(file.path('results_Mar1/within_state/', filename))
  }
  muests <- df_temp 
  # Convert muests to a vector, it's just a single column
  muests <- as.vector(as.matrix(muests))
  # Compute true ERF
  mutrue <- mutrues[mutrues$rangeu == rangeu[i] & 
                      mutrues$option == option[i] & 
                      mutrues$withinstate == ifelse(confounding_scenario[i] == 3, T, F),]$theta
  df_temp$mutrue <- mutrue
  
  # Save metrics in analysisdf
  analysisdf$bias[i] <- mean(muests, na.rm = T) - mutrue
  analysisdf$RMSE[i] <- mean((muests - mutrue)^2, na.rm = T)
  analysisdf$se[i] <- mean((muests - mutrue)^2, na.rm = T) / length(muests)
}

# Print bias table
analysisdf_bias <- analysisdf[,c(1,3:5)] %>% 
  pivot_wider(names_from = method, values_from = c(avgabsbias))
print(xtable(analysisdf_bias, digits = 3), include.rownames = F)

# Print RMSE table
analysisdf_RMSE <- analysisdf[,c(1,3:4, 6)] %>% 
  pivot_wider(names_from = method, values_from = c(avgRMSE))
print(xtable(analysisdf_RMSE, digits = 3), include.rownames = F)

# Plot the ERFs
png('images/ERFplots_tinyscalelinear.png', 
    width = 4500, height = 1200, res = 400)
gridExtra::grid.arrange(grobs = gs[1:4], ncol = 4)
dev.off()

png('images/ERFplots_tinyscalenonlinear.png', 
    width = 4500, height = 1200, res = 400)
gridExtra::grid.arrange(grobs = gs[5:8], ncol = 4)
dev.off()

png('images/ERFplots_smallscalelinear.png', 
    width = 4500, height = 1200, res = 400)
gridExtra::grid.arrange(grobs = gs[9:12], ncol = 4)
dev.off()

png('images/ERFplots_smallscalenonlinear.png',
    width = 4500, height = 1200, res = 400)
gridExtra::grid.arrange(grobs = gs[13:16], ncol = 4)
dev.off()

png('images/ERFplots_tinyscalelinear_withinstate.png', 
    width = 4500, height = 1200, res = 400)
gridExtra::grid.arrange(grobs = gs[17:20], ncol = 4)
dev.off()

png('images/ERFplots_tinyscalenonlinear_withinstate.png', 
    width = 4500, height = 1200, res = 400)
gridExtra::grid.arrange(grobs = gs[21:24], ncol = 4)
dev.off()