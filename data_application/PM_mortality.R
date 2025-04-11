source('../funcs.R')
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(igraph)
library(geosphere)
library(Matrix) 
library(stringr)
library(readr)
library(gridExtra)
library(xtable)
library(lmtest)
library(tidyfst)
library(fst)
library(spdep)
library(mgcv)
library(RSpectra)

set.seed(123)

# Outcome: zipcode-level all cause mortality 2011-2016
# Exposure: zipcode-level all source PM2.5 2001-2010
# Confounders: zipcode-level Census + gridMET + BRFSS 2000
# geography (shapefile): 2010

# Read in National Causal Dataset
aggregate_data <- read.csv('/n/dominici_nsaph_l3/Lab/projects/analytic/aggregated_2000-2016_medicare_mortality_pm25_zip/aggregate_data.csv')

# convert zip to string 
aggregate_data$zip <- str_pad(aggregate_data$zip, width = 5, 
                              side = 'left', pad = '0')

# Extract covariate data from 2000
covariate_data_sub <- aggregate_data[aggregate_data$year %in% 2000, ]

# Extract outcome data from 2014-2016
outcome_data_sub <- aggregate_data[aggregate_data$year %in% 2011:2016, ]

print(length(unique(covariate_data_sub$zip))) # 33833 zip codes
print(length(unique(outcome_data_sub$zip))) # 34464 zip codes

# Select columns of interest: for now omit individual-level stratification!
covariate_data_sub <- select(
  covariate_data_sub,
  zip,
  mean_bmi,
  smoke_rate,
  hispanic,
  pct_blk,
  medhouseholdincome,
  medianhousevalue,
  poverty,
  education,
  popdensity,
  pct_owner_occ,
  summer_tmmx,
  winter_tmmx,
  summer_rmax,
  winter_rmax
)

outcome_data_sub <- dplyr::select(
  outcome_data_sub,
  zip,
  dead,
  time_count
)

# Aggregate covariate and outcome data to the zip.
covariate_data_sub <- aggregate(.~zip, 
                                data = covariate_data_sub,
                                FUN = mean,
                                na.rm = T)
outcome_data_sub <- aggregate(.~zip,
                              data = outcome_data_sub,
                              FUN = sum,
                              na.rm = T)
# Subset outcome_data_sub to zip codes that had person years > 10
outcome_data_sub <- subset(outcome_data_sub, time_count > 10) # 33795/34464

# Compute deathrate from total who died in zip in 2011-2016 divided by total personyears of zip
outcome_data_sub$deathrate <- outcome_data_sub$dead / outcome_data_sub$time_count
gc()

# Initialize an empty list to store dataframes for each year
pm_data_list <- list()

# Read in Exposure data
# Loop through the years 2001 to 2010
for (year in 2001:2010) { 
  # Create the file path for each year
  file_path <- paste0("/n/dominici_nsaph_l3/Lab/exposure/pm25/PM25_v2/annual/", year, ".rds")
  
  # Read in the data and add a 'year' column
  pm_data <- read_rds(file_path) %>%
    as.data.frame() %>%
    mutate(year = year) # Add the year column
  
  # Append the dataframe to the list
  pm_data_list[[as.character(year)]] <- pm_data
}

# Combine all dataframes in the list into a single dataframe
combined_pm_data <- bind_rows(pm_data_list)

# Read in polygon geometry
zipcode_sf_polygon <- st_read('/n/dominici_nsaph_l3/Lab/data/shapefiles/zip_shape_files/2010/zip/polygon/ESRI10USZIP5_POLY_WGS84.shp')
length(unique(zipcode_sf_polygon$ZIP)) # 30425
nrow(zipcode_sf_polygon) # 30430
zipcode_sf_polygon = st_make_valid(zipcode_sf_polygon)

# Some zipcodes are repeated (cross state lines and recorded as two polygons)
# Merge those polygons
duplicate_zipcodes <- zipcode_sf_polygon %>%
  group_by(ZIP) %>%
  filter(n() > 1) %>%
  pull(ZIP) %>%
  unique()
duplicated_sf <- zipcode_sf_polygon %>%
  filter(ZIP %in% duplicate_zipcodes)
merged_duplicated_sf <- duplicated_sf %>%
  group_by(ZIP) %>%
  summarise(geometry = st_union(geometry))
non_duplicated_sf <- zipcode_sf_polygon %>%
  filter(!ZIP %in% duplicate_zipcodes)
zipcode_sf_polygon <- bind_rows(non_duplicated_sf, merged_duplicated_sf)

# Read in point geometry
zipcode_sf_point <- st_read('/n/dominici_nsaph_l3/Lab/data/shapefiles/zip_shape_files/2010/zip/point/ESRI10USZIP5_POINT_WGS84_POBOX.shp')
length(unique(zipcode_sf_point$ZIP)) # 10945
zipcode_sf_point <- st_as_sf(zipcode_sf_point)

# average over years before 2014 to get PM exposure. 
total_pm <- aggregate(pm25~ZIP, data = combined_pm_data, 
                      FUN = mean,
                      na.rm = T)
mean(total_pm$ZIP %in% c(zipcode_sf_polygon$ZIP, zipcode_sf_point$ZIP)) # 98 percent.

polygon_data <- merge(zipcode_sf_polygon, total_pm, by = 'ZIP')
point_data <- merge(zipcode_sf_point, total_pm, by = 'ZIP')

pm25_limits <- range(c(polygon_data$pm25, point_data$pm25), na.rm = TRUE)

# Plot exposure (average PM 2001-2010) on the ZCTAs only
png(
  filename = paste0('images/PMaverage_2001_2010.png'), 
  width = 2000,
  height = 1000,
  res = 200
)
ggplot() +
  xlim(-125,-65) +
  ylim(25, 50) +
  # Plot polygons, color representing PM2.5 exposure
  geom_sf(data = polygon_data, aes(fill = pm25), color = NA) +
  scale_fill_viridis_c(name = "ug/m^3", limits = pm25_limits) +
  labs(title = "PM2.5 Exposure averaged over 2001-2010 across US zip codes") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),     
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))
dev.off()


# Now calculate centroids
polygon_data$centroid <- st_centroid(polygon_data)
polygon_data$centroid_coords <- st_coordinates(polygon_data$centroid)
polygon_data$x <- polygon_data$centroid_coords[, "X"]
polygon_data$y <-polygon_data$centroid_coords[, "Y"]

# Extract coordinates for point_data
point_coords <- st_coordinates(point_data)
point_data$x <- point_coords[, "X"]
point_data$y <- point_coords[, "Y"]

# Select necessary columns and bind the two datasets together
polygon_data_clean <- polygon_data %>% 
  dplyr::select(ZIP, pm25, x, y)

point_data_clean <- point_data %>% 
  dplyr::select(ZIP, pm25, x, y)

# Bind the two datasets
combined_data <- rbind(polygon_data_clean, point_data_clean)

mod = mgcv::gam(combined_data$pm25 ~ s(combined_data$x,combined_data$y,k=4,fx=T)) # unpenalized

combined_data$Ac_TPS = predict(mod)
combined_data$Auc_TPS = combined_data$pm25-combined_data$Ac_TPS
print(var(combined_data$Ac_TPS)/var(combined_data$pm25)) # 0.22

adj <- st_intersects(combined_data, sparse = T) # captures both polygon-polygon adj and point-in-polygon
adjacency_matrix <- sparseMatrix(
  i = rep(seq_along(adj), lengths(adj)),  # Row indices based on the list
  j = unlist(adj),                        # Flatten the list into column indices
  x = 1,                                  # Value for adjacency (1)
  dims = c(length(adj), length(adj))      # Dimensions of the sparse matrix
)
diag(adjacency_matrix) <- 0

# Graph Laplacian
D <- Diagonal(x = rowSums(adjacency_matrix))

# Subtract the adjacency matrix from the sparse diagonal matrix
L <- D - adjacency_matrix

start_time = Sys.time()
eig <- eigs_sym(L, k = 4, which = "SM", opts = list(tol = 1e-3, maxitr = 10000))
print(Sys.time() - start_time) # 33 sec
gc()

mod = lm(combined_data$pm25 ~ eig$vectors[,2:4]) # %*% t(eig$vectors) %*% combined_data$pm25
combined_data$Ac_GraphLaplacian = predict(mod)
combined_data$Auc_GraphLaplacian = combined_data$pm25-combined_data$Ac_GraphLaplacian
print(var(combined_data$Ac_GraphLaplacian)/var(combined_data$pm25)) # 0.20

g1 <- ggplot() +
  xlim(-125, -65) +
  ylim(25, 50) +
  
  # Plot polygons and points, both using the same color scale for PM2.5 exposure
  geom_sf(data = combined_data, aes(fill = pm25, color = pm25), size = 0.5) +
  
  scale_fill_viridis_c(name = "ug/m^3") +  # For both polygons and points
  scale_color_viridis_c(guide = "none") +  # Hides the redundant color legend
  
  labs(title = "PM2.5 Exposure (A) averaged over 2001-2010 across US zip codes") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))

g2 <- ggplot() +
  xlim(-125, -65) +
  ylim(25, 50) +
  
  # Plot polygons and points, both using the same color scale for Auc_TPS
  geom_sf(data = combined_data, aes(fill = Auc_TPS, color = Auc_TPS), size = 0.5) +
  
  scale_fill_viridis_c(name = "ug/m^3") +  # For both polygons and points
  scale_color_viridis_c(guide = "none") +  # Hides the redundant color legend
  
  labs(title = "Estimate of Auc using a thin plate spline") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))

g3 <- ggplot() +
  xlim(-125, -65) +
  ylim(25, 50) +
  
  # Plot polygons and points, both using the same color scale for Auc_GraphLaplacian
  geom_sf(data = combined_data, aes(fill = Auc_GraphLaplacian, color = Auc_GraphLaplacian), size = 0.5) +
  
  scale_fill_viridis_c(name = "ug/m^3") +  # For both polygons and points
  scale_color_viridis_c(guide = "none") +  # Hides the redundant color legend
  
  labs(title = "Estimate of Auc using the Graph Laplacian") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))

g4 <- ggplot() +
  xlim(-125, -65) +
  ylim(25, 50) +
  
  # Plot polygons and points, both using the same color scale for Auc_TPS
  geom_sf(data = combined_data, aes(fill = Ac_TPS, color = Ac_TPS), size = 0.5) +
  
  scale_fill_viridis_c(name = "ug/m^3") +  # For both polygons and points
  scale_color_viridis_c(guide = "none") +  # Hides the redundant color legend
  
  labs(title = "Estimate of Ac using a thin plate spline") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))

g5 <- ggplot() +
  xlim(-125, -65) +
  ylim(25, 50) +
  
  # Plot polygons and points, both using the same color scale for Auc_TPS
  geom_sf(data = combined_data, aes(fill = Ac_GraphLaplacian, color = Ac_GraphLaplacian), size = 0.5) +
  
  scale_fill_viridis_c(name = "ug/m^3") +  # For both polygons and points
  scale_color_viridis_c(guide = "none") +  # Hides the redundant color legend
  
  labs(title = "Estimate of Ac using the Graph Laplacian") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))

layout_matrix <- rbind(c(NA,1,1,NA),
                       c(2,2,4,4),
                       c(3,3,5,5))  

png(
  filename = paste0('images/PM_ivs_20pct_2010.png'),
  width = 2500,
  height = 3000,
  res = 300
)
grid.arrange(g1, g2, g3, g4, g5, layout_matrix = layout_matrix)
dev.off()
gc()

# Merge covariate_data_sub with combined_data
combined_data_covariates <- merge(combined_data, covariate_data_sub, 
                                  by.x = 'ZIP',
                                  by.y = 'zip') #  33829
combined_data_covariates_outcome <- merge(combined_data_covariates,
                                          outcome_data_sub,
                                          by.x = 'ZIP',
                                          by.y = 'zip')
nrow(combined_data_covariates_outcome) # 33255

covs <- st_drop_geometry(
  dplyr::select(
    combined_data_covariates_outcome,
    mean_bmi,
    smoke_rate,
    hispanic,
    pct_blk,
    medhouseholdincome,
    medianhousevalue,
    poverty,
    education,
    popdensity,
    pct_owner_occ,
    summer_tmmx,
    winter_tmmx,
    summer_rmax,
    winter_rmax
  )
)

# Data characteristics table
round(cbind(apply(covs,2,mean),apply(covs,2,sd)),3)
round(c(mean(combined_data_covariates_outcome$pm25), 
        sd(combined_data_covariates_outcome$pm25)),3)
round(c(mean(combined_data_covariates_outcome$deathrate), 
        sd(combined_data_covariates_outcome$deathrate)),3)

# Confounder matrix
covs <- as.data.frame(scale(covs))

# Confounder matrix excluding spatial confounders
xmat <- covs %>%
  dplyr::select(
    -summer_tmmx,
    -winter_tmmx,
    -summer_rmax,
    -winter_rmax
  )

############################################### ESTIMATE TRUNCATED EXPOSURE EFFECT FOR EACH CUTOFF WITH EACH METHOD ####################################

cutoffs <- seq(6,12,1) 
delta <- 1
y <- combined_data_covariates_outcome$deathrate 
a <- combined_data_covariates_outcome$pm25 

for (cutoff in cutoffs){
  # Baseline
  xsub <- xmat[a > cutoff - delta,]
  colnames(xsub) <- colnames(xmat)
  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction"),
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  baseline_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
                      mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  baseline_ci_norm <- c(baseline_est - 1.96*sqrt(vhat/length(y)),
                        baseline_est + 1.96*sqrt(vhat/length(y)))
  print(paste0("Baseline: ", c(baseline_est, baseline_ci_norm)))
  save(baseline_ci_norm, file = paste0('results_Mar24/baseline_ci_', cutoff, '.RData'))
  save(baseline_est, file = paste0('results_Mar24/baseline_est_', cutoff, '.RData'))
  
  # Spatialcoord
  xsub <- cbind(xmat,
                combined_data_covariates_outcome$x,
                combined_data_covariates_outcome$y)[which(a>cutoff-delta),]
  colnames(xsub) <- c(colnames(xmat),'lon','lat')
  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.glm", "SL.glm.interaction", "SL.mean", "SL.gam"),
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  spatialcoord_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
                          mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  spatialcoord_ci_norm <- c(spatialcoord_est - 1.96*sqrt(vhat/length(y)),
                            spatialcoord_est + 1.96*sqrt(vhat/length(y)))
  print(c(spatialcoord_est, spatialcoord_ci_norm))
  save(spatialcoord_ci_norm, file = paste0('results_Mar24/spatialcoord_ci_', cutoff, '.RData'))
  save(spatialcoord_est, file = paste0('results_Mar24/spatialcoord_est_', cutoff, '.RData'))
  
  # IV-TPS
  xsub <- cbind(xmat,
                combined_data_covariates_outcome$Ac_TPS)[which(a>cutoff-delta),]
  colnames(xsub) <- c(colnames(xmat),'Ac_TPS')

  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  IV_TPS_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
                    mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  IV_TPS_ci_norm <- c(IV_TPS_est - 1.96*sqrt(vhat/length(y)),
                      IV_TPS_est + 1.96*sqrt(vhat/length(y)))
  print(c(IV_TPS_est, IV_TPS_ci_norm))
  save(IV_TPS_ci_norm, file = paste0('results_Mar24/IV_TPS_ci_', cutoff, '.RData'))
  save(IV_TPS_est, file = paste0('results_Mar24/IV_TPS_est_', cutoff, '.RData'))
  
  # IV-GraphLaplacian
  xsub <- cbind(xmat,
                combined_data_covariates_outcome$Ac_GraphLaplacian)[which(a>cutoff-delta),]
  colnames(xsub) <- c(colnames(xmat),'Ac_GraphLaplacian')
 
  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  IV_GraphLaplacian_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
                               mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  IV_GraphLaplacian_ci_norm <- c(IV_GraphLaplacian_est - 1.96*sqrt(vhat/length(y)),
                                 IV_GraphLaplacian_est + 1.96*sqrt(vhat/length(y)))
  print(c(IV_GraphLaplacian_est, IV_GraphLaplacian_ci_norm))
  save(IV_GraphLaplacian_ci_norm,
       file = paste0('results_Mar24/IV_GraphLaplacian_ci_', cutoff, '.RData'))
  save(IV_GraphLaplacian_est,
       file = paste0('results_Mar24/IV_GraphLaplacian_est_', cutoff, '.RData'))
  
  # IV_TPS-spatialcoord
  xsub <- cbind(xmat,
                combined_data_covariates_outcome$Ac_TPS,
                combined_data_covariates_outcome$x,
                combined_data_covariates_outcome$y)[which(a>cutoff-delta),]
  colnames(xsub) <- c(colnames(xmat),'Ac_TPS', 'lon', 'lat')
  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  IV_TPS_spatialcoord_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
                                 mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  IV_TPS_spatialcoord_ci_norm <- c(IV_TPS_spatialcoord_est - 1.96*sqrt(vhat/length(y)),
                                   IV_TPS_spatialcoord_est + 1.96*sqrt(vhat/length(y)))
  print(c(IV_TPS_spatialcoord_est, IV_TPS_spatialcoord_ci_norm))
  save(IV_TPS_spatialcoord_ci_norm,
       file = paste0('results_Mar24/IV_TPS_spatialcoord_ci_', cutoff, '.RData'))
  save(IV_TPS_spatialcoord_est,
       file = paste0('results_Mar24/IV_TPS_spatialcoord_est_', cutoff, '.RData'))
  
  # IV_GraphLaplacian-spatialcoord
  xsub <- cbind(xmat,
                combined_data_covariates_outcome$Ac_GraphLaplacian,
                combined_data_covariates_outcome$x,
                combined_data_covariates_outcome$y)[which(a>cutoff-delta),]
  colnames(xsub) <- c(colnames(xmat),'Ac_GraphLaplacian', 'lon', 'lat')
  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  IV_GraphLaplacian_spatialcoord_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
                                            mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  IV_GraphLaplacian_spatialcoord_ci_norm <- c(IV_GraphLaplacian_spatialcoord_est - 1.96*sqrt(vhat/length(y)),
                                              IV_GraphLaplacian_spatialcoord_est + 1.96*sqrt(vhat/length(y)))
  print(c(IV_GraphLaplacian_spatialcoord_est, IV_GraphLaplacian_spatialcoord_ci_norm))
  save(IV_GraphLaplacian_spatialcoord_ci_norm,
       file = paste0('results_Mar24/IV_GraphLaplacian_spatialcoord_ci_', cutoff, '.RData'))
  save(IV_GraphLaplacian_spatialcoord_est,
       file = paste0('results_Mar24/IV_GraphLaplacian_spatialcoord_est_', cutoff, '.RData'))
  
  # Oracle
  xsub <- covs[a > cutoff - delta,]
  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction"),
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  oracle_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
                    mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  print((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
           mean(y[a<=cutoff])*mean(a<=cutoff)))
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  oracle_ci_norm <- c(oracle_est - 1.96*sqrt(vhat/length(y)),
                      oracle_est + 1.96*sqrt(vhat/length(y)))
  print(cutoff)
  print(paste0("Oracle: ", c(oracle_est, oracle_ci_norm)))
  save(oracle_ci_norm,
       file = paste0('results_Mar24/oracle_ci_', cutoff, '.RData'))
  save(oracle_est,
       file = paste0('results_Mar24/oracle_est_', cutoff, '.RData'))
  gc()
}

################################################# Exposure Response Curves ########################################################

qs <- quantile(a, probs = c(0.025, 0.975))

# BASELINE
start_time <- Sys.time()
erfest <- ctseff(
  y = y,
  a = a,
  x = xmat,
  n.pts = 100,
  a.rng = qs,
  sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
  bw.seq = seq(sd(a)/10, sd(a), length.out = 100),
  savephi = F
)
print(Sys.time() - start_time) 
save(erfest, file = paste0('results_Mar24/baseline_erf_.RData'))

# ORACLE
start_time <- Sys.time()
erfest <- ctseff(
  y = y,
  a = a,
  x = covs,
  n.pts = 100,
  a.rng = qs,
  sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
  bw.seq = seq(sd(a)/10, sd(a), length.out = 100),
  savephi = F
)
print(Sys.time() - start_time) 
save(erfest, file = paste0('results_Mar24/oracle_erf_.RData'))

# SPATIALCOORD
xsub <- cbind(xmat,
              combined_data_covariates_outcome$x,
              combined_data_covariates_outcome$y)
colnames(xsub) <- c(colnames(xmat),'lon', 'lat')
start_time <- Sys.time()
erfest <- ctseff(
  y = y,
  a = a,
  x = xsub,
  n.pts = 100,
  a.rng = qs,
  sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
  bw.seq = seq(sd(a)/10, sd(a), length.out = 100),
  savephi = F
)
print(Sys.time() - start_time) 
save(erfest, file = paste0('results_Mar24/spatialcoord_erf_.RData'))

# IV_TPS
xsub <- cbind(xmat,
              combined_data_covariates_outcome$Ac_TPS)
colnames(xsub) <- c(colnames(xmat),'Ac_TPS')
start_time <- Sys.time()
erfest <- ctseff(
  y = y,
  a = a,
  x = xsub,
  n.pts = 100,
  a.rng = qs,
  sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
  bw.seq = seq(sd(a)/10, sd(a), length.out = 100),
  savephi = F
)
print(Sys.time() - start_time) 
save(erfest, file = paste0('results_Mar24/IV_TPS_erf_.RData'))

# IV_GraphLaplacian
xsub <- cbind(xmat,
              combined_data_covariates_outcome$Ac_GraphLaplacian)
colnames(xsub) <- c(colnames(xmat),'Ac_GraphLaplacian')
start_time <- Sys.time()
erfest <- ctseff(
  y = y,
  a = a,
  x = xsub,
  n.pts = 100,
  a.rng = qs,
  sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
  bw.seq = seq(sd(a)/10, sd(a), length.out = 100),
  savephi = F
)
print(Sys.time() - start_time) 
save(erfest, file = paste0('results_Mar24/IV_GraphLaplacian_erf_.RData'))

# IV_TPS_spatialcoord
xsub <- cbind(xmat,
              combined_data_covariates_outcome$Ac_TPS,
              combined_data_covariates_outcome$x,
              combined_data_covariates_outcome$y)
colnames(xsub) <- c(colnames(xmat),'Ac_TPS', 'lon', 'lat')
start_time <- Sys.time()
erfest <- ctseff(
  y = y,
  a = a,
  x = xsub,
  n.pts = 100,
  a.rng = qs,
  sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
  bw.seq = seq(sd(a)/10, sd(a), length.out = 100),
  savephi = F
)
print(Sys.time() - start_time) 
save(erfest, file = paste0('results_Mar24/IV_TPS_spatialcoord_erf_.RData'))

# IV_GraphLaplacian_spatialcoord
xsub <- cbind(xmat,
              combined_data_covariates_outcome$Ac_GraphLaplacian,
              combined_data_covariates_outcome$x,
              combined_data_covariates_outcome$y)
colnames(xsub) <- c(colnames(xmat),'Ac_GraphLaplacian', 'lon', 'lat')
start_time <- Sys.time()
erfest <- ctseff(
  y = y,
  a = a,
  x = xsub,
  n.pts = 100,
  a.rng = qs,
  sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), 
  bw.seq = seq(sd(a)/10, sd(a), length.out = 100),
  savephi = F
)
print(Sys.time() - start_time) 
save(erfest, 
     file = paste0('results_Mar24/IV_GraphLaplacian_spatialcoord_erf_.RData'))


########################################## PLOTS OF TRUNCATED EXPOSURE EFFECT ESTIMATES ##############################################

# Create dataframe to store estimates
cutoffs <- 6:12
df <- data.frame(
  cutoff = numeric(0),
  method = character(0),
  point_est = numeric(0),
  lower = numeric(0),
  upper = numeric(0)
)
for (j in 1:length(cutoffs)){
  cutoff <- cutoffs[j]
  load(paste0('results_Mar24/baseline_ci_', cutoff, '.RData'))
  load(paste0('results_Mar24/baseline_est_', cutoff, '.RData'))
  load(paste0('results_Mar24/oracle_ci_', cutoff, '.RData'))
  load(paste0('results_Mar24/oracle_est_', cutoff, '.RData'))
  load(paste0('results_Mar24/spatialcoord_ci_', cutoff, '.RData'))
  load(paste0('results_Mar24/spatialcoord_est_', cutoff, '.RData'))
  load(paste0('results_Mar24/IV_TPS_ci_', cutoff, '.RData'))
  load(paste0('results_Mar24/IV_TPS_est_', cutoff, '.RData'))
  load(paste0('results_Mar24/IV_GraphLaplacian_ci_', cutoff, '.RData'))
  load(paste0('results_Mar24/IV_GraphLaplacian_est_', cutoff, '.RData'))
  load(paste0('results_Mar24/IV_TPS_spatialcoord_ci_', cutoff, '.RData'))
  load(paste0('results_Mar24/IV_TPS_spatialcoord_est_', cutoff, '.RData'))
  load(paste0('results_Mar24/IV_GraphLaplacian_spatialcoord_ci_', cutoff, '.RData'))
  load(paste0('results_Mar24/IV_GraphLaplacian_spatialcoord_est_', cutoff, '.RData'))
  # Combine the estimates and intervals into a single data frame
  newdf <- data.frame(
    cutoff = rep(cutoff, 7),
    method = factor(c("oracle", "baseline", "spatialcoord", "IV-TPS", "IV-GraphLaplacian",
                      "IV-TPS+spatialcoord", "IV-GraphLaplacian+spatialcoord"),
                    levels = c("oracle", "baseline", "spatialcoord", "IV-TPS", "IV-GraphLaplacian",
                               "IV-TPS+spatialcoord", "IV-GraphLaplacian+spatialcoord")),
    point_est = c(oracle_est, baseline_est, spatialcoord_est, IV_TPS_est,
                  IV_GraphLaplacian_est, IV_TPS_spatialcoord_est, IV_GraphLaplacian_spatialcoord_est),
    lower = c(oracle_ci_norm[1], baseline_ci_norm[1], spatialcoord_ci_norm[1], IV_TPS_ci_norm[1],
              IV_GraphLaplacian_ci_norm[1], IV_TPS_spatialcoord_ci_norm[1], IV_GraphLaplacian_spatialcoord_ci_norm[1]),
    upper = c(oracle_ci_norm[2], baseline_ci_norm[2], spatialcoord_ci_norm[2], IV_TPS_ci_norm[2],
              IV_GraphLaplacian_ci_norm[2], IV_TPS_spatialcoord_ci_norm[2], IV_GraphLaplacian_spatialcoord_ci_norm[2])
  )
  df <- rbind.data.frame(df, newdf)
}

df$cutoff_factor <- as.factor(df$cutoff)

# PLOT RESULTS

# Define custom facet labels
cutoff_labels <- c("6" = "cutoff (µg/m³): 6",
                   "7" = "cutoff (µg/m³): 7",
                   "8" = "cutoff (µg/m³): 8",
                   "9" = "cutoff (µg/m³): 9", 
                   "10" = "cutoff (µg/m³): 10", 
                   "11" = "cutoff (µg/m³): 11",
                   "12" = "cutoff (µg/m³): 12",
                   "13" = "cutoff (µg/m³): 13",
                   "14" = "cutoff (µg/m³): 14",
                   "15" = "cutoff (µg/m³): 15")
df$point_est <- 100*(df$point_est-1)
df$lower <- 100*(df$lower-1)
df$upper <- 100*(df$upper-1)

png(
  filename = paste0('images/estimates_cis_2010_apr3_full.png'),
  width = 4100,
  height = 2000,
  res = 407
)
ggplot(df, aes(x = method, y = point_est, color = method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = method), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ cutoff_factor, nrow = 1, labeller = as_labeller(cutoff_labels)) +
  scale_color_manual(values = c("oracle" = "forestgreen", 
                                "baseline" = "red", 
                                "spatialcoord" = "black", 
                                "IV-TPS" = "black", 
                                "IV-GraphLaplacian" = "black", 
                                "IV-TPS+spatialcoord" = "black", 
                                "IV-GraphLaplacian+spatialcoord" = "black")) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Method", 
       y = "Reduction in Mortality Rate (%)",
       title = expression("Impact of enforcing average " * PM[2.5] * " (2001-2010) below a cutoff on mortality rate in the Medicare population (2011–2016)")) +
  ylim(-8.25, 0.5)
dev.off()

# One number summary of distance to oracle in terms of Hausdorff distance
cutoffs <- 6:12
methods <- c('oracle',
             'baseline', 
             'spatialcoord', 
             'IV-TPS', 
             'IV-GraphLaplacian', 
             'IV-TPS+spatialcoord',
             'IV-GraphLaplacian+spatialcoord')
df$hd <- NA
for (cutoff in cutoffs){
  interval1 <- c(df$lower[df$method == 'oracle' & df$cutoff == cutoff], 
                 df$upper[df$method == 'oracle' & df$cutoff == cutoff])
  for (method in methods){
    interval2 <- c(df$lower[df$method == method & df$cutoff == cutoff], 
                   df$upper[df$method == method & df$cutoff == cutoff])
    df$hd[df$method == method & df$cutoff == cutoff] <- hausdorff_distance(interval1, interval2)
  }
}
aggregate(.~method, 
          data = df,
          FUN = mean,
          na.rm = T)
df_hd <- data.frame(
  hd6 = df$hd[df$cutoff == 6],
  hd7 = df$hd[df$cutoff == 7],
  hd8 = df$hd[df$cutoff == 8],
  hd9 = df$hd[df$cutoff == 9],
  hd10 = df$hd[df$cutoff == 10],
  hd11 = df$hd[df$cutoff == 11],
  hd12 = df$hd[df$cutoff == 12],
  average = aggregate(.~method, 
                      data = df,
                      FUN = mean,
                      na.rm = T)$hd)
df_hd_t <- t(df_hd)
df_hd_t <- round(df_hd_t,3)
df_hd_t <- as.data.frame(df_hd_t)
colnames(df_hd_t) <- df$method[df$cutoff == 6]
print(xtable(df_hd_t))

########################################## PLOTS FOR THE EXPOSURE RESPONSE CURVE #########################################

methods <- c('baseline', 
             'oracle', 
             'spatialcoord', 
             'IV_TPS', 
             'IV_GraphLaplacian',
             'IV_TPS_spatialcoord',
             'IV_GraphLaplacian_spatialcoord'
)
df <- data.frame(
  a.vals = seq(qs[1], qs[2], length.out = 100)
)
cils <- data.frame(
  a.vals = seq(qs[1], qs[2], length.out = 100)
)
cius <- data.frame(
  a.vals = seq(qs[1], qs[2], length.out = 100)
)
for (j in 1:length(methods)){
  method <- methods[j]
  load(paste0('results_Mar24/', method, '_erf_.RData'))
  df[j+1] = 100*erfest$res$est
  cils[j+1] = 100*erfest$res$ci.ll
  cius[j+1] = 100*erfest$res$ci.ul
}
colnames(df)[2:ncol(df)] <- methods
colnames(cils)[2:ncol(df)] <- methods
colnames(cius)[2:ncol(df)] <- methods

# Convert data to long format for ggplot
data_long <- df %>%
  pivot_longer(cols = -a.vals, names_to = "Method", values_to = "Value")

gmain <- ggplot(data_long, aes(x = a.vals, y = Value, color = Method)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("baseline" = "red", 
                                "oracle" = "black", 
                                "spatialcoord" = "purple", 
                                "IV_TPS" = "darkgreen", 
                                "IV_TPS_spatialcoord" = "seagreen2",
                                "IV_GraphLaplacian" = "blue",
                                "IV_GraphLaplacian_spatialcoord" = "lightskyblue")) +
  labs(x = "PM Exposure Averaged over 2001-2010, µg/m³", 
       y = "Mortality Rate 1/(100 person-years)", 
       title = "ERC between air pollution exposure \n and mortality rate in US zipcodes",
       color = "Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(3.9, 5.1)

data <- data.frame(
  a.vals = df$a.vals,
  Baseline_est = df$baseline,
  Baseline_ci_ll = cils$baseline,
  Baseline_ci_ul = cius$baseline,
  Oracle_est = df$oracle,
  Oracle_ci_ll = cils$oracle,
  Oracle_ci_ul = cius$oracle,
  SpatialCoordinates_est = df$spatialcoord,
  SpatialCoordinates_ci_ll = cils$spatialcoord,
  SpatialCoordinates_ci_ul = cius$spatialcoord,
  IV_TPS_est = df$IV_TPS,
  IV_TPS_ci_ll = cils$IV_TPS,
  IV_TPS_ci_ul = cius$IV_TPS,
  IV_GraphLaplacian_est = df$IV_GraphLaplacian,
  IV_GraphLaplacian_ci_ll = cils$IV_GraphLaplacian,
  IV_GraphLaplacian_ci_ul = cius$IV_GraphLaplacian,
  IV_TPS_spatialcoord_est = df$IV_TPS_spatialcoord,
  IV_TPS_spatialcoord_ci_ll = cils$IV_TPS_spatialcoord,
  IV_TPS_spatialcoord_ci_ul = cius$IV_TPS_spatialcoord,
  IV_GraphLaplacian_spatialcoord_est = df$IV_GraphLaplacian_spatialcoord,
  IV_GraphLaplacian_spatialcoord_ci_ll = cils$IV_GraphLaplacian_spatialcoord,
  IV_GraphLaplacian_spatialcoord_ci_ul = cius$IV_GraphLaplacian_spatialcoord
)

data_filtered <-  data
# Plot00: Oracle only
plot00 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # Baseline with confidence intervals
  geom_ribbon(aes(ymin = Oracle_ci_ll, ymax = Oracle_ci_ul), fill = "black", alpha = 0.2) +
  labs(x = "PM Exposure Averaged over 2001-2010, µg/m³", 
       y = "Mortality Rate 1/(100 person-years)", 
       title = "Oracle",
       color = "Method") +
  scale_color_manual(values = c("Oracle" = "black")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(3.9, 5.1)

# Plot 0: Baseline and Oracle
plot0 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # Baseline with confidence intervals
  geom_ribbon(aes(ymin = Baseline_ci_ll, ymax = Baseline_ci_ul), fill = "red", alpha = 0.2) +
  labs(x = "PM Exposure Averaged over 2001-2010, µg/m³", 
       y = "Mortality Rate 1/(100 person-years)", 
       title = "Baseline",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(3.9, 5.1)

# Plot 1: Baseline, Oracle, and Spatial Coordinates (with Confidence Bands)
plot1 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # Spatial Coordinates with confidence intervals
  geom_ribbon(aes(ymin = SpatialCoordinates_ci_ll, ymax = SpatialCoordinates_ci_ul), fill = "purple", alpha = 0.2) +
  geom_line(aes(y = SpatialCoordinates_est, color = "SpatialCoordinates"), size = 1) +
  labs(x = "PM Exposure Averaged over 2001-2010, µg/m³", 
       y = "Mortality Rate 1/(100 person-years)", 
       title = "Spatial Coordinates",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black", 
                                "SpatialCoordinates" = "purple")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(3.9, 5.1)

# Plot 2: Baseline, Oracle, and IV-TPS (with Confidence Bands)
plot2 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # IV-TPS with confidence intervals
  geom_ribbon(aes(ymin = IV_TPS_ci_ll, ymax = IV_TPS_ci_ul), fill = "darkgreen", alpha = 0.1) +
  geom_line(aes(y = IV_TPS_est, color = "IV_TPS"), size = 1) +
  labs(x = "PM Exposure Averaged over 2001-2010, µg/m³", 
       y = "Mortality Rate 1/(100 person-years)", 
       title = "IV-TPS",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black", 
                                "IV_TPS" = "darkgreen")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+ 
  ylim(3.9, 5.1)

# Plot 3: Baseline, Oracle, and IV-GraphLaplacian (with Confidence Bands)
plot3 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # IV-GraphLaplacian with confidence intervals
  geom_ribbon(aes(ymin = IV_GraphLaplacian_ci_ll, ymax = IV_GraphLaplacian_ci_ul), fill = "blue", alpha = 0.1) +
  geom_line(aes(y = IV_GraphLaplacian_est, color = "IV_GraphLaplacian"), size = 1) +
  labs(x = "PM Exposure Averaged over 2001-2010, µg/m³", 
       y = "Mortality Rate 1/(100 person-years)", 
       title = "IV-GraphLaplacian",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black", 
                                "IV_GraphLaplacian" = "blue")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+ 
  ylim(3.9, 5.1)

# Plot 4: Baseline, Oracle, and IV-TPS_spatialcoord (with Confidence Bands)
plot4 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # IV-TPS with confidence intervals
  geom_ribbon(aes(ymin = IV_TPS_spatialcoord_ci_ll, ymax = IV_TPS_spatialcoord_ci_ul), fill = "seagreen2", alpha = 0.3) +
  geom_line(aes(y = IV_TPS_spatialcoord_est, color = "IV_TPS_spatialcoord"), size = 1) +
  labs(x = "PM Exposure Averaged over 2001-2010, µg/m³",
       y = "Mortality Rate 1/(100 person-years)",
       title = "IV-TPS+spatialcoord",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red",
                                "Oracle" = "black",
                                "IV_TPS_spatialcoord" = "seagreen2")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  ylim(3.9, 5.1)

# Plot 5: Baseline, Oracle, and IV-GraphLaplacian_spatialcoord (with Confidence Bands)
plot5 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # IV-GraphLaplacian with confidence intervals
  geom_ribbon(aes(ymin = IV_GraphLaplacian_spatialcoord_ci_ll, ymax = IV_GraphLaplacian_spatialcoord_ci_ul), fill = "lightskyblue", alpha = 0.3) +
  geom_line(aes(y = IV_GraphLaplacian_spatialcoord_est, color = "IV_GraphLaplacian_spatialcoord"), size = 1) +
  labs(x = "PM Exposure Averaged over 2001-2010, µg/m³",
       y = "Mortality Rate 1/(100 person-years)",
       title = "IV-GraphLaplacian+spatialcoord",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red",
                                "Oracle" = "black",
                                "IV_GraphLaplacian_spatialcoord" = "lightskyblue")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  ylim(3.9, 5.1)
# Create the plot
png(
  filename = paste0('images/erfs_Apr2.png'),
  width = 10000,
  height = 10000,
  res = 800
)
# Define a layout matrix
layout_matrix <- rbind(c(1,2),
                       c(3,4),
                       c(5,6),
                       c(7,8)
)  

# Arrange the plots according to the layout
grid.arrange(gmain, plot00, plot0, plot1, plot2, plot3, plot4, plot5,
             layout_matrix = layout_matrix)
dev.off()


############################################# SENSITIVITY ANALYSIS FOR DF in TWO DECOMPOSITIONS ######################################

cutoff <- 9
ks <- seq(4,8,by = 1) # Vary df in spline for TPS. 4 is minimum
for (k in ks){
  print(k)
  mod = mgcv::gam(combined_data$pm25 ~ s(combined_data$x,combined_data$y,k=k,fx=T)) # unpenalized
  combined_data[[paste0('Ac_TPS_', k)]] = predict(mod)
  combined_data[[paste0('Auc_TPS_', k)]] = combined_data$pm25-combined_data$Ac_TPS
}
gc()

adj <- st_intersects(combined_data, sparse = T) # captures both polygon-polygon adj and point-in-polygon
adjacency_matrix <- sparseMatrix(
  i = rep(seq_along(adj), lengths(adj)),  # Row indices based on the list
  j = unlist(adj),                        # Flatten the list into column indices
  x = 1,                                  # Value for adjacency (1)
  dims = c(length(adj), length(adj))      # Dimensions of the sparse matrix
)
diag(adjacency_matrix) <- 0

# Graph Laplacian
D <- Diagonal(x = rowSums(adjacency_matrix))

# Subtract the adjacency matrix from the sparse diagonal matrix
L <- D - adjacency_matrix

start_time = Sys.time()
eig <- eigs_sym(L, k = 7, which = "SM", opts = list(tol = 1e-3, maxitr = 10000))
print(Sys.time() - start_time) # 5 min
gc()

ks <- seq(1,7,by = 1) # Vary number of eigenvectors for IV-GraphLaplacian.
for (k in ks){
  print(k)
  mod = lm(combined_data$pm25 ~ eig$vectors[,(ncol(eig$vectors)-k+1):ncol(eig$vectors)]) 
  combined_data[[paste0('Ac_GraphLaplacian_', k)]] = predict(mod)
  combined_data[[paste0('Auc_GraphLaplacian_', k)]] = combined_data$pm25-combined_data$Ac_GraphLaplacian
}
combined_data_covariates <- merge(combined_data, covariate_data_sub, 
                                  by.x = 'ZIP',
                                  by.y = 'zip') #  33829
combined_data_covariates_outcome <- merge(combined_data_covariates,
                                          outcome_data_sub,
                                          by.x = 'ZIP',
                                          by.y = 'zip')
nrow(combined_data_covariates_outcome) # 33255
gc()

xmat <- data.frame('Intercept' = rep(1, nrow(combined_data_covariates_outcome)))

# IV_TPS
delta <- 1
y <- combined_data_covariates_outcome$deathrate
a <- combined_data_covariates_outcome$pm25
for (k in seq(4,8,by = 1)){
  print(k)
  xsub <- cbind(xmat,
                combined_data_covariates_outcome[[paste0('Ac_TPS_', k)]])[which(a>cutoff-delta),]
  colnames(xsub) <- c(colnames(xmat),'Ac_TPS')
  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), # exclude RFs # "SL.gam",
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  IV_TPS_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) +
                    mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  IV_TPS_ci_norm <- c(IV_TPS_est - 1.96*sqrt(vhat/length(y)),
                      IV_TPS_est + 1.96*sqrt(vhat/length(y)))
  print(c(IV_TPS_est, IV_TPS_ci_norm))
  save(IV_TPS_ci_norm, file = paste0('results_Mar24/sensitivity/IV_TPS_ci_', k, '.RData'))
  save(IV_TPS_est, file = paste0('results_Mar24/sensitivity/IV_TPS_est_', k, '.RData'))
}

gc()

# IV-GraphLaplacian
for (k in seq(1,7,by = 1)){
  print(k)
  xsub <- cbind(xmat,
                combined_data_covariates_outcome[[paste0('Ac_GraphLaplacian_', k)]])[which(a>cutoff-delta),]
  colnames(xsub) <- c(colnames(xmat),'Ac_GraphLaplacian')
  start_time <- Sys.time()
  erfest <- ctseff(
    y = y[a > cutoff - delta],
    a = a[a > cutoff - delta],
    x = xsub,
    n.pts = 5,
    a.rng = c(cutoff - delta, cutoff + delta),
    sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), # exclude RFs
    bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
  )
  print(Sys.time() - start_time) 
  IV_GraphLaplacian_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) + 
                               mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
  vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
  IV_GraphLaplacian_ci_norm <- c(IV_GraphLaplacian_est - 1.96*sqrt(vhat/length(y)), 
                                 IV_GraphLaplacian_est + 1.96*sqrt(vhat/length(y))) 
  print(c(IV_GraphLaplacian_est, IV_GraphLaplacian_ci_norm))
  save(IV_GraphLaplacian_ci_norm, file = paste0('results_Mar24/sensitivity/IV_GraphLaplacian_ci_', k, '.RData'))
  save(IV_GraphLaplacian_est, file = paste0('results_Mar24/sensitivity/IV_GraphLaplacian_est_', k, '.RData'))
}

load(paste0('results_Mar24/oracle_ci_', cutoff, '.RData'))
load(paste0('results_Mar24/oracle_est_', cutoff, '.RData'))

# Baseline
xsub <- xmat[which(a>cutoff-delta),]
xsub <- as.matrix(xsub, ncol = 1)
colnames(xsub) <- colnames(xmat)
start_time <- Sys.time()
erfest <- ctseff(
  y = y[a > cutoff - delta],
  a = a[a > cutoff - delta],
  x = xsub,
  n.pts = 5,
  a.rng = c(cutoff - delta, cutoff + delta),
  sl.lib = c("SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"), # exclude RFs
  bw.seq = seq(sd(a)/10, sd(a), length.out = 100)
)
print(Sys.time() - start_time) 
baseline_est <- ((erfest$res$est[erfest$res$a.vals == cutoff]*mean(a>cutoff) + 
                    mean(y[a<=cutoff])*mean(a<=cutoff)))/mean(y)
vhat <- asymptotic_variance_delta(y, a, erfest, cutoff, delta)
baseline_ci_norm <- c(baseline_est - 1.96*sqrt(vhat/length(y)), 
                      baseline_est + 1.96*sqrt(vhat/length(y))) 

# PLOT IV TPS for each k by k, as well as oracle and baseline.
ks <- 4:8 
ks <- as.factor(ks)
df <- data.frame(
  k = numeric(0),
  point_est = numeric(0),
  lower = numeric(0),
  upper = numeric(0)
)

for (j in 1:length(ks)){
  k <- ks[j]
  load(paste0('results_Mar24/sensitivity/IV_TPS_ci_', k, '.RData'))
  load(paste0('results_Mar24/sensitivity/IV_TPS_est_', k, '.RData'))
  # Combine the estimates and intervals into a single data frame
  newdf <- data.frame(
    k = k,
    point_est = IV_TPS_est,
    lower = IV_TPS_ci_norm[1],
    upper = IV_TPS_ci_norm[2]
  )
  df <- rbind.data.frame(df, newdf)
}
df <- rbind.data.frame(data.frame(k = 'baseline', 
                                  point_est = baseline_est,
                                  lower = baseline_ci_norm[1],
                                  upper = baseline_ci_norm[2]),
                       df)
df <- rbind.data.frame(data.frame(k = 'oracle', 
                                  point_est = oracle_est,
                                  lower = oracle_ci_norm[1],
                                  upper = oracle_ci_norm[2]),
                       df)
df$k <- factor(df$k, levels = c("oracle", "baseline", as.character(4:8)))

df$point_est <- 100*(df$point_est-1)
df$lower <- 100*(df$lower-1)
df$upper <- 100*(df$upper-1)

png(
  filename = paste0('images/estimates_cis_sensitivity_TPS_nocov.png'),
  width = 2000,
  height = 1500,
  res = 300
)
ggplot(df, aes(x = k, y = point_est, color = k)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  # Add a line only through the 4-8 points; group=1 ensures they are connected in order
  geom_line(data = subset(df, !(k %in% c("oracle", "baseline"))), aes(group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  scale_color_manual(values = c("oracle" = "forestgreen", 
                                "baseline" = "red", 
                                "4" = "black", 
                                "5" = "black", 
                                "6" = "black", 
                                "7" = "black", 
                                "8" = "black")) +
  labs(x = "Dimension of Thin Plate Spline basis", 
       y = "Reduction in Mortality Rate (%)",
       title = expression(atop("Impact of enforcing average " * PM[2.5] * " (2001-2010) below 9 µg/m³ on",
                               "mortality rate in the Medicare population (2011-2016)")))

dev.off()

# PLOT IV Graph Laplacian for each k by k, as well as oracle and baseline.
ks <- 3:7
ks <- as.factor(ks)
df <- data.frame(
  k = numeric(0),
  point_est = numeric(0),
  lower = numeric(0),
  upper = numeric(0)
)

for (j in 1:length(ks)){
  k <- ks[j]
  load(paste0('results_Mar24/sensitivity/IV_GraphLaplacian_ci_', k, '.RData'))
  load(paste0('results_Mar24/sensitivity/IV_GraphLaplacian_est_', k, '.RData'))
  # Combine the estimates and intervals into a single data frame
  newdf <- data.frame(
    k = k,
    point_est = IV_GraphLaplacian_est,
    lower = IV_GraphLaplacian_ci_norm[1],
    upper = IV_GraphLaplacian_ci_norm[2]
  )
  df <- rbind.data.frame(df, newdf)
}
df <- rbind.data.frame(data.frame(k = 'baseline', 
                                  point_est = baseline_est,
                                  lower = baseline_ci_norm[1],
                                  upper = baseline_ci_norm[2]),
                       df)
df <- rbind.data.frame(data.frame(k = 'oracle', 
                                  point_est = oracle_est,
                                  lower = oracle_ci_norm[1],
                                  upper = oracle_ci_norm[2]),
                       df)
df$k <- factor(df$k, levels = c("oracle", "baseline", as.character(ks)))
df$point_est <- 100*(df$point_est-1)
df$lower <- 100*(df$lower-1)
df$upper <- 100*(df$upper-1)

png(
  filename = paste0('images/estimates_cis_sensitivity_GraphLaplacian_nocov.png'),
  width = 2000,
  height = 1500,
  res = 300
)
ggplot(df, aes(x = k, y = point_est, color = k)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  # Add a line only through the 4-8 points; group=1 ensures they are connected in order
  geom_line(data = subset(df, !(k %in% c("oracle", "baseline"))), aes(group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("oracle" = "forestgreen", 
                                "baseline" = "red", 
                                "3" = "black", 
                                "4" = "black", 
                                "5" = "black",
                                "6" = "black",
                                "7" = "black")) +
  labs(x = "Dimension of Eigenbasis of Graph Laplacian", 
       y = "Reduction in Mortality Rate (%)",
       title = expression(atop("Impact of enforcing average " * PM[2.5] * " (2001-2010) below 9 µg/m³ on",
                               "mortality rate in the Medicare population (2011-2016)")))
dev.off()