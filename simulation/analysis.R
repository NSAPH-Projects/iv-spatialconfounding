library(ggplot2)
library(sf)
library(tidyr)
library(mgcv)
library(dplyr)
library(utils)
library(xtable)
library(Matrix)
library(tidyverse)
library(parallel)
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
csvs_withinstate <- list.files('results_Mar1/within_state/', 
                               pattern = '^0\\.01.*\\.csv$')
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
confounding_scenario <- c(ifelse(rangeu == '0.01', 1, 2)[1:length(csvs_notwithinstate)], 
                         rep(3, length(csvs_withinstate)))
option <- gsub("^(.*?_)(.*?)(_.*$)", "\\2", csvs)
method <- gsub(".*_([^\\.]+)\\.csv", "\\1", csvs)

# Precompute true estimand for each outcome model and confounding mechanism
# mutrues <- data.frame(expand.grid(rangeu = c(0.01, 0.05), 
#                                  option = c('linear', 'nonlinear')))
# mutrues$withinstate <- F
# mutrues <- rbind(mutrues, data.frame(rangeu = 0.01, option = 'linear', withinstate = T))
# mutrues <- rbind(mutrues, data.frame(rangeu = 0.01, option = 'nonlinear', withinstate = T))
# mutrues$theta <- NA
# mutrues$option <- as.character(mutrues$option)
# 
# mutrues$theta <- unlist(mclapply(1:nrow(mutrues), function(i) {
#   computemutrue(option = mutrues$option[i], 
#                 rangeu = mutrues$rangeu[i], 
#                 within_state_GP = mutrues$withinstate[i],
#                 distmat = distmat,
#                 statemat = simlist$statemat)
# }, mc.cores = 2))  # Adjust the number of cores
# 
# mutrues$confounding_scenario <- c(1, 2, 1, 2, 3, 3)
# save(mutrues, file = 'results_Mar1/mutrues.RData')
load('results_Mar1/mutrues.RData')
                
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


# Suppose after pivoting, your table has these columns:
# "confounding_scenario", "rangeu", "option", "method", "bias", "RMSE", "se"
# Format the numeric columns (bias, RMSE, se) in scientific notation with 3 significant figures
analysisdf_bias <- analysisdf %>%
  mutate(
    bias = round(bias*1000, 3),
    RMSE = round(RMSE*1000, 3),
    se = round(se*1000, 3)
  )

# Pivot wider from the original analysisdf
analysisdf_bias <- analysisdf_bias[, c(1, 3:5)] %>% 
  pivot_wider(names_from = method, values_from = bias)
analysisdf_bias <- analysisdf_bias[,c("confounding_scenario", "option", "baseline", "spatialcoord", "IV-TPS", "IV-GraphLaplacian", "IV-TPS-spatialcoord", "IV-GraphLaplacian-spatialcoord")]

# Print using xtable and prevent xtable from reformatting the already-formatted text
print(xtable(analysisdf_bias), include.rownames = FALSE, sanitize.text.function = identity)

# Do the same with RMSE
analysisdf_RMSE <- analysisdf %>%
  mutate(
    bias = round(bias*1000, 3),
    RMSE = round(RMSE*1000, 3),
    se = round(se*1000, 3)
  )

analysisdf_RMSE <- analysisdf_RMSE[, c(1, 3:4, 6)] %>% 
  pivot_wider(names_from = method, values_from = RMSE)
analysisdf_RMSE <- analysisdf_RMSE[,c("confounding_scenario", "option", "baseline", "spatialcoord", "IV-TPS", "IV-GraphLaplacian", "IV-TPS-spatialcoord", "IV-GraphLaplacian-spatialcoord")]

print(xtable(analysisdf_RMSE), include.rownames = FALSE, sanitize.text.function = identity)

# Now create facet_wrap boxplots with ggplot2

folder <- "results_Mar1"

# List all CSV files in that folder (with full paths)
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

read_estimates <- function(i) {
  
  filename <- csvs[i]
  print(filename)
  
  confounding_scenario <- confounding_scenario[i]
  rangeu <- rangeu[i]
  option <- option[i]
  method <- method[i]
  if (confounding_scenario !=3){
    dat <- read.csv(file.path('results_Mar1/', filename))
  }
  else{
    dat <- read.csv(file.path('results_Mar1/within_state/', filename))
  }

  # If the CSV doesn't have a header and just one column, name it "estimate"
  if (!"estimate" %in% colnames(dat)) {
    names(dat)[1] <- "estimate"
  }
  # Add the new columns
  dat <- dat %>%
    mutate(confounding_scenario = confounding_scenario,
           option = option,
           method = method)
  return(dat)
}

# Read all files and combine into one data frame
df <- map_dfr(1:length(csvs), read_estimates)

# Ensure rangeu and option are factors in both data frames with the same levels:
df <- df %>% 
  mutate(confounding_mechanism = factor(confounding_scenario),
         option = factor(option, levels = c("linear", "nonlinear")))
mutrues <- mutrues %>% 
  mutate(confounding_mechanism = factor(confounding_scenario),
         option = factor(option, levels = c("linear", "nonlinear")))

# Create the boxplot with horizontal lines for theta
ggplot(df, aes(x = method, y = estimate))+#, fill = method)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(option ~ confounding_mechanism) +
  geom_hline(data = mutrues, aes(yintercept = theta), 
             color = "red", linetype = "dashed", size = 1) +
  labs(x = NULL, y = "Estimate") +                    # Remove x-axis title
  scale_fill_discrete(name = "Method") +              # Change legend title
  theme_bw() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(-0.5, 0.5)
