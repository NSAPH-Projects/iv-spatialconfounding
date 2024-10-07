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
csvs_notwithinstate <- list.files('results_Oct1/', pattern = '.csv')
# Confounding scenario 3 (GP within state)
csvs_withinstate <- list.files('results_Oct1/within_state/', pattern = '.csv') 
csvs <- c(csvs_notwithinstate, csvs_withinstate)

# Create storage for metrics 
analysisdf <- data.frame(
  confounding_scenario = character(length(csvs)),
  rangeu = character(length(csvs)),
  option = character(length(csvs)),
  method = character(length(csvs)),
  avgabsbias = numeric(length(csvs)),
  avgRMSE = numeric(length(csvs)),
  avgse = numeric(length(csvs))
)

# Extract components from filenames
rangeu <- gsub("(_.*$)", "", csvs)
confounding_scenario <- c(ifelse(rangeu == '0.05', 1, 2)[1:length(csvs_notwithinstate)], 
                         rep(3, length(csvs_withinstate)))
option <- gsub("^(.*?_)(.*?)(_.*$)", "\\2", csvs)
method <- gsub(".*_(IV-[A-Za-z]+|[A-Za-z]+)\\.csv", "\\1", csvs)

method[method == 'spatialcoord'] <- 'Spatial coordinates'
method[method == 'baseline'] <- 'No confounding adjustment'

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
    df_temp <- read.csv(file.path('results_Oct1/', filename))
  }
  else{
    df_temp <- read.csv(file.path('results_Oct1/within_state/', filename))
  }
  muests <- df_temp[,-1] # remove a.vals
  avals <- df_temp$avals
  # Compute true ERF
  mutrue <- computemutrue(a.vals = avals,
                         option = option[i])
  df_temp$mutrue <- mutrue
  # Compute metrics
  met <- metrics(a.vals = avals, 
                muests,
                mutrue)
  # Save metrics in analysisdf
  analysisdf$avgabsbias[i] <- met$avgabsbias
  analysisdf$avgRMSE[i] <- met$avgRMSE
  analysisdf$avgse[i] <- met$avgse
  
  # Create plot of estimated ERFs and mutrue
  df_temp <- df_temp %>%
    pivot_longer(cols = -c(avals, mutrue), names_to = "sim", values_to = "value") %>%
    mutate(sim = factor(sim))

  gs[[i]] <- ggplot(df_temp, aes(x = avals, y = value)) +
    # plot the mean line for the simulations
    stat_summary(fun = mean, geom = "line", aes(color = "Mean"), size = 1.2) +

    # plot the 95% confidence interval
    stat_summary(fun.min = function(x) quantile(x, 0.025),
                 fun.max = function(x) quantile(x, 0.975),
                 geom = "ribbon", alpha = 0.2, fill = "blue") +

    # add mutrue as a dashed line
    geom_line(aes(x = avals, y = mutrue), color = "black", linetype = "dashed") +

    labs(x = "a", y = expression("E( Y"^"a" ~ ")")) +

    # scale colors, keeping the mean as a distinct color
    scale_color_manual(values = c("Mean" = "blue", "black")) +

    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    xlim(-2,2) +
    ylim(-7, 7) +
    ggtitle(method[i])
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