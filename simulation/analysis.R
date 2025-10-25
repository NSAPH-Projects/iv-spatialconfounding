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
library(geosphere)
library(patchwork)

source('../funcs.R')
load('sim.RData')

# Calculate distance matrix
distmat <- distm(cbind(simlist$lon, simlist$lat), 
                           fun = distHaversine)
# Standardize so approximately range (0,2)
distmat <- distmat/1000000

csvs <- list.files('results_Sep6/', pattern = '.csv')


# Create storage for metrics 
analysisdf <- data.frame(
  confounding_mechanism = character(length(csvs)),
  option = character(length(csvs)),
  method = character(length(csvs)),
  bias = numeric(length(csvs)),
  RMSE = numeric(length(csvs)),
  se = numeric(length(csvs))
)

# Extract components from filenames
confounding_mechanism <- as.integer(sub("^conf(\\d+)_.*$", "\\1", csvs))
option <- sub("^conf\\d+_([^_]+)_.*$", "\\1", csvs)
method <- sub("^conf\\d+_[^_]+_([^_]+)\\.csv$", "\\1", csvs)

# Precompute true estimand for each outcome model and confounding mechanism
mutrues <- data.frame(expand.grid(
  confounding_mechanism = 1:7,
  option = c('linear', 'nonlinear')))
for (i in 1:nrow(mutrues)){
  oracle_results <- read.csv('results_Sep6/conf' %>% 
                               paste0(mutrues$confounding_mechanism[i], '_', 
                                      mutrues$option[i], '_oracle.csv'))
  mutrues$theta[i] = mean(as.vector(as.matrix(oracle_results)), na.rm = T)
}

# Loop through results to calculate metrics and create plots.
for (i in 1:length(csvs)){
  filename <- csvs[i]
  print(filename)
  
  analysisdf$confounding_mechanism[i] <- confounding_mechanism[i]
  analysisdf$option[i] <- option[i]
  analysisdf$method[i] <- method[i]
  df_temp <- read.csv(file.path('results_Sep6/', filename))
  
  muests <- df_temp 
  # Convert muests to a vector, it's just a single column
  muests <- as.vector(as.matrix(muests))
  
  # Compute true truncated exposure estimate
  mutrue <- mutrues[mutrues$confounding_mechanism == confounding_mechanism[i] & 
                      mutrues$option == option[i],]$theta 
  df_temp$mutrue <- mutrue
  
  # Save metrics in analysisdf
  analysisdf$bias[i] <- mean(muests, na.rm = T) - mutrue
  analysisdf$RMSE[i] <- sqrt(mean((muests - mutrue)^2, na.rm = T))
  analysisdf$se[i] <- sd(muests, na.rm = T)
}

# Format the numeric columns (bias, RMSE, se) in scientific notation with 3 decimals
analysisdf_bias <- analysisdf %>%
  mutate(
    bias = round(bias*100, 3),
    RMSE = round(RMSE*100, 3),
    se = round(se*100, 3)
  )

# Pivot wider from the original analysisdf
analysisdf_bias <- analysisdf_bias[, 1:4] %>% 
  pivot_wider(names_from = method, values_from = bias)
analysisdf_bias <- analysisdf_bias[,c("confounding_mechanism", 
                                      "option", 
                                      "oracle",
                                      "baseline", 
                                      "spatialcoord", 
                                      "trueIV",
                                      "IV-TPS", 
                                      "IV-GraphLaplacian",
                                      "trueIV-spatialcoord",
                                      "IV-TPS-spatialcoord",
                                      "IV-GraphLaplacian-spatialcoord"
                                      )]

# Print using xtable and prevent xtable from reformatting the already-formatted text
print(xtable(analysisdf_bias), 
      include.rownames = FALSE, sanitize.text.function = identity)
# Print absolute bias using xtable
print(xtable(analysisdf_bias %>%
               mutate(across(where(is.numeric), ~ abs(.)))), 
      include.rownames = FALSE, sanitize.text.function = identity)
# Print absolute bias with the reordered confounding scenarios
analysisdf_bias_reordered <- analysisdf_bias %>%
  mutate(confounding_mechanism = factor(confounding_mechanism, 
                                       levels = c("1", "7", "3", "5", "6", "2", "4"),
                                       labels = c("1", "2", "3", "4", "5", "6", "7"))) %>%
  arrange(confounding_mechanism, option)
print(xtable(analysisdf_bias_reordered %>%
               mutate(across(where(is.numeric), ~ abs(.)))), 
      include.rownames = FALSE, sanitize.text.function = identity)
# Do the same with RMSE
analysisdf_RMSE <- analysisdf %>%
  mutate(
    bias = round(bias*100, 3),
    RMSE = round(RMSE*100, 3),
    se = round(se*100, 3)
  )

analysisdf_RMSE <- analysisdf_RMSE[, c(1:3,5)] %>% 
  pivot_wider(names_from = method, values_from = RMSE)
analysisdf_RMSE <- analysisdf_RMSE[,c("confounding_mechanism", 
                                      "option", 
                                      "oracle",
                                      "baseline", 
                                      "spatialcoord",
                                      "trueIV",
                                      "IV-TPS", 
                                      "IV-GraphLaplacian",
                                      "trueIV-spatialcoord",
                                      "IV-TPS-spatialcoord",
                                      "IV-GraphLaplacian-spatialcoord"
                                      )]

print(xtable(analysisdf_RMSE), include.rownames = FALSE, sanitize.text.function = identity)
# Print the reordered RMSE
analysisdf_RMSE_reordered <- analysisdf_RMSE %>%
  mutate(confounding_mechanism = factor(confounding_mechanism, 
                                       levels = c("1", "7", "3", "5", "6", "2", "4"),
                                       labels = c("1", "2", "3", "4", "5", "6", "7"))) %>%
  arrange(confounding_mechanism, option)
print(xtable(analysisdf_RMSE_reordered), include.rownames = FALSE, sanitize.text.function = identity)

# Now create facet_wrap boxplots with ggplot2

folder <- "results_Sep6"

read_estimates <- function(i) {
  
  filename <- csvs[i]
  print(filename)
  
  confounding_mechanism <- confounding_mechanism[i]
  option <- option[i]
  method <- method[i]
  dat <- read.csv(file.path('results_Sep6/', filename))

  # If the CSV doesn't have a header and just one column, name it "estimate"
  if (!"estimate" %in% colnames(dat)) {
    names(dat)[1] <- "estimate"
  }
  # Add the new columns
  dat <- dat %>%
    mutate(confounding_mechanism = confounding_mechanism,
           option = option,
           method = method)
  return(dat)
}

# Read all files and combine into one data frame
df <- map_dfr(1:length(csvs), read_estimates)
# Rename the method "IV-GraphLaplacian" to "IV-GL" in df
df$method[df$method == 'IV-GraphLaplacian-spatialcoord'] = "IV-GL+spatialcoord"
df$method[df$method == 'IV-GraphLaplacian'] = "IV-GL"
df$method[df$method == 'IV-TPS-spatialcoord'] = "IV-TPS+spatialcoord"
df$method[df$method == 'trueIV-spatialcoord'] = "trueIV+spatialcoord"


desired_order <- c("oracle", "baseline", "spatialcoord", 
                   "trueIV", 
                   "IV-TPS", "IV-GL", 
                   "trueIV+spatialcoord", "IV-TPS+spatialcoord", "IV-GL+spatialcoord"
                   )
df$method <- factor(df$method, levels = desired_order)

df <- df %>% 
  mutate(confounding_mechanism = factor(confounding_mechanism),
         option = factor(option, levels = c("linear", "nonlinear")))
mutrues <- mutrues %>% 
  mutate(confounding_mechanism = factor(confounding_mechanism),
         option = factor(option, levels = c("linear", "nonlinear")))
# Reorder confounding mechanism: 1 (GP), 7 (Leroux), 3 (discrete), 5 (2conf) ,6 (flipped)
df$confounding_mechanism_reordered <- factor(df$confounding_mechanism, 
                                              levels = c("1", "7", "3", "5", "6", "2", "4"),
                                              labels = c("1", "2", "3", "4", "5", "6", "7"))
mutrues$confounding_mechanism_reordered <- factor(mutrues$confounding_mechanism, 
                                                 levels = c("1", "7", "3", "5", "6", "2", "4"),
                                                 labels = c("1", "2", "3", "4", "5", "6", "7"))
print(xtable(select(mutrues, confounding_mechanism_reordered, option, theta) %>% 
               arrange(confounding_mechanism_reordered), digits = 4), 
      include.rownames = FALSE)

# Create the boxplot with horizontal lines for theta
method_cols <- c(
  "oracle"                 = "gray",
  "baseline"               = "#D62728",
  "spatialcoord"           = "purple", 
  "trueIV"                 = "lightgreen", 
  "IV-TPS"                 = "#9ECAE1",
  "IV-GL"                  = "#FDAE6B", 
  "trueIV+spatialcoord"    = "#2E8B57", 
  "IV-TPS+spatialcoord"    = "#1F77B4", 
  "IV-GL+spatialcoord"     = "#E6550D" 
)
png("images/boxplot_Sep6.png", width = 2500, height = 1250, res = 200)
ggplot(df, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.5, outliers = F, staplewidth = 1) + #, draw_quantiles = c(0.5)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "blue") +
  ggh4x::facet_grid2(
    option ~ confounding_mechanism_reordered,
    scales = "free",           # allows different scales per row/col
    independent = "all"        # allows different scales **per panel**
  ) +
  #facet_grid(option ~ confounding_mechanism, scales = "free_y") +
  geom_hline(data = mutrues, aes(yintercept = theta), 
             color = "red", linetype = "twodash", size = 1) +
  #labs(x = NULL, y = "Truncated Exposure Effect Estimate") +   
  labs(
    x = "Confounding mechanism (1–7)",
    y = "Truncated Exposure Effect Estimate"
  ) +
  # Remove x-axis title
  scale_fill_manual(
    name = "Method",
    values = method_cols,
    breaks = names(method_cols) # matches your desired order
  ) +
  theme_bw() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") 
dev.off()

# Plot for the main text: Confounding scenarios 1,3,5,6
df_sub <- subset(df, confounding_mechanism_reordered %in% 1:5 & 
               method != 'trueIV' & method != 'trueIV+spatialcoord')
df_sub <- droplevels(df_sub)
mutrues_sub <- subset(mutrues, confounding_mechanism_reordered %in% 1:5)
mutrues_sub <- droplevels(mutrues_sub)

png("images/boxplot_Sep6_maintext.png", width = 2000, height = 1300, res = 200)
ggplot(df_sub, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.5, outliers = F, staplewidth = 1) + #, draw_quantiles = c(0.5)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "blue") +
  ggh4x::facet_grid2(
    option ~ confounding_mechanism_reordered,
    scales = "free",           # allows different scales per row/col
    independent = "all"        # allows different scales **per panel**
  ) +
  geom_hline(data = mutrues_sub, aes(yintercept = theta), 
             color = "red", linetype = "twodash", size = 1) +
  labs(
    x = "Confounding mechanism (1–5)",
    y = "Truncated Exposure Effect Estimate"
  ) +
  # Remove x-axis title
  scale_fill_manual(
    name = "Method",
    values = method_cols,
    breaks = names(method_cols) # matches your desired order
  ) +
  theme_bw() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") 
dev.off()
