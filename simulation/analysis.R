library(ggplot2)
library(sf)
library(tidyr)
library(mgcv)
library(dplyr)
library(utils)
library(xtable)
library(Matrix)
library(tidyverse)
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
confounding_scenario <- c(ifelse(rangeu == '0.01', 1, 2)[1:length(csvs_notwithinstate)], 
                         rep(3, length(csvs_withinstate)))
option <- gsub("^(.*?_)(.*?)(_.*$)", "\\2", csvs)
method <- gsub(".*_(IV-[A-Za-z]+|[A-Za-z]+)\\.csv", "\\1", csvs)

method[method == 'spatialcoord'] <- 'Spatial coordinates'
method[method == 'baseline'] <- 'No confounding adjustment'

# Precompute true estimand for each outcome model and confounding mechanism
mutrues <- data.frame(expand.grid(rangeu = c(0.01, 0.05), 
                                 option = c('linear', 'nonlinear')))
mutrues$withinstate <- F
mutrues <- rbind(mutrues, data.frame(rangeu = 0.01, option = 'linear', withinstate = T))
mutrues <- rbind(mutrues, data.frame(rangeu = 0.01, option = 'nonlinear', withinstate = T))
mutrues$theta <- NA
mutrues$option <- as.character(mutrues$option)

for (i in 1:nrow(mutrues)){
  mutrues$theta[i] <- computemutrue(option = mutrues$option[i], 
                                     rangeu = mutrues$rangeu[i], 
                                    within_state_GP = mutrues$withinstate[i],
                                    distmat = distmat,
                                    statemat = simlist$statemat)
}
save(mutrues, file = 'results_Mar1/mutrues.RData')
                
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

# Now create facet_wrap boxplots with ggplot2

folder <- "results_Mar1"

# List all CSV files in that folder (with full paths)
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

read_estimates <- function(file) {
  # Extract parts from the filename, assuming they are separated by underscores.
  file_base <- basename(file)
  parts <- strsplit(file_base, "_")[[1]]
  rangeu <- parts[1]                     # e.g., "0.01" or "0.05"
  option <- parts[2]                   # "linear" or "nonlinear"
  # Remove ".csv" from the method part
  method <- sub(".csv", "", parts[3])
  
  # Read the CSV. Adjust header = TRUE/FALSE depending on your file.
  dat <- read.csv(file, header = TRUE)
  # If the CSV doesn't have a header and just one column, name it "estimate"
  if (!"estimate" %in% colnames(dat)) {
    names(dat)[1] <- "estimate"
  }
  # Add the new columns
  dat <- dat %>%
    mutate(rangeu = rangeu,
           option = option,
           method = method)
  return(dat)
}

# Read all files and combine into one data frame
df <- map_dfr(files, read_estimates)

# Remove last two rows of mutrue
mutrues <- mutrues[-c(5, 6),]

# Ensure rangeu and option are factors in both data frames with the same levels:
df <- df %>% 
  mutate(rangeu = factor(rangeu, levels = c("0.01", "0.05")),
         option = factor(option, levels = c("linear", "nonlinear")))
mutrues <- mutrues %>% 
  mutate(rangeu = factor(rangeu, levels = c("0.01", "0.05")),
         option = factor(option, levels = c("linear", "nonlinear")))

# Create the boxplot with horizontal lines for theta
ggplot(df, aes(x = method, y = estimate)) +
  geom_boxplot() +
  facet_grid(option ~ rangeu) +
  # Add horizontal lines: one line per facet matching on rangeu and option
  geom_hline(data = mutrues, aes(yintercept = theta), 
             color = "red", linetype = "dashed", size = 1) +
  labs(x = "Method", y = "Estimate") +
  theme_bw() +
  ylim(0,2.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
