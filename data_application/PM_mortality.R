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

# Outcome: zipcode-level all cause mortality 2014-2016
# Exposure: zipcode-level all source PM2.5 2001-2013
# Confounders: zipcode-level Census + gridMET + BRFSS 2000
# geography (shapefile): 2013

# Read in National Causal Dataset
aggregate_data <- read.csv('/n/dominici_nsaph_l3/Lab/projects/analytic/aggregated_2000-2016_medicare_mortality_pm25_zip/aggregate_data.csv')

# convert zip to string 
aggregate_data$zip <- str_pad(aggregate_data$zip, width = 5, 
                              side = 'left', pad = '0')

# Extract covariate data from 2000
covariate_data_sub <- aggregate_data[aggregate_data$year %in% 2000, ]

# Extract outcome data from 2014-2016
outcome_data_sub <- aggregate_data[aggregate_data$year %in% 2014:2016, ]

print(length(unique(covariate_data_sub$zip))) # 33833 zip codes
print(length(unique(outcome_data_sub$zip))) # 34238 zip codes

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
  winter_rmax,
  region
)

# remove region string and replace it with indicators
covariate_data_sub$region <- as.factor(covariate_data_sub$region)

# Create design matrix (indicator variables) for 'region'
region_design_matrix <- model.matrix(~ region - 1, data=covariate_data_sub)
covariate_data_sub <- cbind.data.frame(covariate_data_sub %>% select(-region), 
                                      region_design_matrix)

outcome_data_sub <- select(
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

# Compute deathrate from total who died in zip in 2014-2016 divided by total personyears of zip
outcome_data_sub$deathrate <- outcome_data_sub$dead / outcome_data_sub$time_count

# Initialize an empty list to store dataframes for each year
pm_data_list <- list()

# Read in Exposure data
# Loop through the years 2001 to 2013
for (year in 2001:2013) {
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
zipcode_sf_polygon <- st_read('/n/dominici_nsaph_l3/Lab/data/shapefiles/zip_shape_files/2013/zip/polygon/ESRI13USZIP5_POLY_WGS84.shp')
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
zipcode_sf_point <- st_read('/n/dominici_nsaph_l3/Lab/data/shapefiles/zip_shape_files/2013/zip/point/ESRI13USZIP5_POINT_WGS84_POBOX.shp')
length(unique(zipcode_sf_point$ZIP)) # 10945
zipcode_sf_point <- st_as_sf(zipcode_sf_point)

# average over years before 2014 to get PM exposure. 
total_pm <- aggregate(pm25~ZIP, data = combined_pm_data, 
                     FUN = mean,
                     na.rm = T)
mean(total_pm$ZIP %in% c(zipcode_sf_polygon$ZIP, zipcode_sf_point$ZIP)) # 96.8 percent.

polygon_data <- merge(zipcode_sf_polygon, total_pm, by = 'ZIP')
point_data <- merge(zipcode_sf_point, total_pm, by = 'ZIP')

pm25_limits <- range(c(polygon_data$pm25, point_data$pm25), na.rm = TRUE)

# Plot exposure (average PM 2001-2013) on the ZCTAs only
png(
  filename = paste0('images/PMaverage_2001_2013.png'),
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
  labs(title = "PM2.5 Exposure averaged over 2001-2013 across US zip codes") +
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
  select(ZIP, pm25, x, y)

point_data_clean <- point_data %>% 
  select(ZIP, pm25, x, y)

# Bind the two datasets
combined_data <- rbind(polygon_data_clean, point_data_clean)

# Check the combined data
head(combined_data)
nrow(combined_data) # 40657

maxdf <- 10
for (k in 2:(maxdf + 1)){
  mod <- mgcv::gam(combined_data$pm25 ~ s(combined_data$x,combined_data$y,k=k,fx=T)) # unpenalized
  combined_data[[paste0("Ac_TPS_", k)]] <- predict(mod)
  combined_data[[paste0("Auc_TPS_", k)]] <- combined_data$pm25-combined_data[[paste0("Ac_TPS_", k)]]
  print(var(predict(mod), na.rm = T)/var(combined_data$pm25, na.rm = T))
}
# mod <- mgcv::gam(combined_data$pm25 ~ s(combined_data$x,combined_data$y,k=6,fx=T)) # unpenalized
# combined_data$Ac_TPS <- predict(mod)
# combined_data$Auc_TPS <- combined_data$pm25-combined_data$Ac_TPS

adj <- st_intersects(combined_data, sparse = T) # captures both polygon-polygon adj and point-in-polygon
adjacency_matrix <- sparseMatrix(
  i = rep(seq_along(adj), lengths(adj)),  # Row indices based on the list
  j = unlist(adj),                        # Flatten the list into column indices
  x = 1,                                  # Value for adjacency (1)
  dims = c(length(adj), length(adj))      # Dimensions of the sparse matrix
)

# Graph Laplacian
D <- Diagonal(x = rowSums(adjacency_matrix))

# Subtract the adjacency matrix from the sparse diagonal matrix
L <- D - adjacency_matrix

start_time <- Sys.time()
maxdf <- 8
eig <- eigs_sym(L, k = 1+maxdf, which = "SM", opts = list(tol = 1e-3, maxitr = 3000)) 
print(eig$values)
print(Sys.time() - start_time) # 33 sec
for (k in 2:(1+maxdf)){
  mod <- lm(combined_data$pm25 ~ eig$vectors[,1:k])
  combined_data[[paste0("Ac_GraphLaplacian_",k)]] <- predict(mod)
  combined_data[[paste0("Auc_GraphLaplacian_",k)]] <- combined_data$pm25-combined_data[[paste0("Ac_GraphLaplacian_",k)]]
  print(var(predict(mod), na.rm = T)/var(combined_data$pm25, na.rm = T))
}
colnames(combined_data)[which(colnames(combined_data)== "Ac_GraphLaplacian_6")] <- 'Ac_GraphLaplacian'
colnames(combined_data)[which(colnames(combined_data)== "Auc_GraphLaplacian_6")] <- 'Auc_GraphLaplacian'

# mod <- lm(combined_data$pm25 ~ eig$vectors) 
# combined_data$Ac_GraphLaplacian <- predict(mod)
# combined_data$Auc_GraphLaplacian <- combined_data$pm25-combined_data$Ac_GraphLaplacian

g1 <- ggplot() +
  xlim(-125, -65) +
  ylim(25, 50) +
  
  # Plot polygons and points, both using the same color scale for PM2.5 exposure
  geom_sf(data = combined_data, aes(fill = pm25, color = pm25), size = 0.5) +
  
  scale_fill_viridis_c(name = "ug/m^3") +  # For both polygons and points
  scale_color_viridis_c(guide = "none") +  # Hides the redundant color legend
  
  labs(title = "PM2.5 Exposure (A) averaged over 2001-2013 across US zip codes") +
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
  
  labs(title = "Estimate of the instrumental variable (Auc) using a thin plate spline") +
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
  
  labs(title = "Estimate of the instrumental variable (Auc) using the Graph Laplacian") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(50, "points"))

png(
  filename = paste0('images/PM_ivs_5.png'),
  width = 2000,
  height = 3000,
  res = 200
)
grid.arrange(g1, g2, g3, nrow = 3)
dev.off()

# Merge covariate_data_sub with combined_data
combined_data_covariates <- merge(combined_data, covariate_data_sub, 
                                 by.x = 'ZIP',
                                 by.y = 'zip') #  33829
combined_data_covariates_outcome <- merge(combined_data_covariates,
                                         outcome_data_sub,
                                         by.x = 'ZIP',
                                         by.y = 'zip')
nrow(combined_data_covariates_outcome) # 33454

covs <- st_drop_geometry(
  select(
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
    winter_rmax,
    regionMIDWEST,
    regionNORTHEAST,
    regionSOUTH,
    regionWEST
  )
)

# Data characteristics Table
round(cbind(apply(covs, 2, mean), apply(covs, 2, sd)), 3)
round(c(
  mean(combined_data_covariates_outcome$pm25),
  sd(combined_data_covariates_outcome$pm25)
), 3)
round(c(
  mean(combined_data_covariates_outcome$deathrate),
  sd(combined_data_covariates_outcome$deathrate)
), 3)

# Transform outcome variable to make model fitting easier
combined_data_covariates_outcome$logdeathrate <- log(combined_data_covariates_outcome$deathrate + 0.01)

# Calculate ERF for 100 exposure values in the 10%-90% percentile range of exposure
qs <- quantile(combined_data_covariates_outcome$pm25, probs = c(0.1,0.9))
a.vals <- seq(qs[1], qs[2], length.out = 100)

# Used for bandwidth selection
sdpm <- sd(combined_data_covariates_outcome$pm25)

# 'Nonspatial' confounders
xmat <- covs %>%
  select(
    -regionMIDWEST,
    -regionNORTHEAST,
    -regionSOUTH,
    -regionWEST,
    -summer_tmmx,
    -winter_tmmx,
    -summer_rmax,
    -winter_rmax
  )

# Baseline: adjust for nonspatial confounders only
erfest_baseline <- ctseff(y = combined_data_covariates_outcome$logdeathrate,
                              a = combined_data_covariates_outcome$pm25,
                              x = xmat,
                          bw.seq = seq(sdpm/2,sdpm,length.out = 50),
                          a.rng = c(min(a.vals), max(a.vals)),
                          sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean")) # exclude RFs

# Spatial coordinates: adjust for nonspatial confounders and spatial coordinates
erfest_spatialcoord <-  ctseff(y = combined_data_covariates_outcome$logdeathrate,
                               a = combined_data_covariates_outcome$pm25,
                               x = cbind(xmat,
                                         combined_data_covariates_outcome$x,
                                         combined_data_covariates_outcome$y
                               ),
                               bw.seq = seq(sdpm/2,sdpm,length.out = 50),
                               a.rng = c(min(a.vals), max(a.vals)),
                               sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean")) # exclude RFs

# IV-TPS: adjust for nonspatial confounders and predicted values from 5df TPS
erfest_IV_TPS <- ctseff(y = combined_data_covariates_outcome$logdeathrate,
                            a = combined_data_covariates_outcome$pm25,
                            x = cbind(xmat,
                                      combined_data_covariates_outcome$Ac_TPS
                              ),
                        bw.seq = seq(sdpm/2,sdpm,length.out = 50),
                        a.rng = c(min(a.vals), max(a.vals)),
                        sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"))

# IV-GraphLaplacian: adjust for nonspatial confounders and projection of exposure onto 5 largest scale eigenvec of GL
erfest_IV_GraphLaplacian <- ctseff(y = combined_data_covariates_outcome$logdeathrate,
                             a = combined_data_covariates_outcome$pm25,
                             x = cbind(xmat,
                                       combined_data_covariates_outcome$Ac_GraphLaplacian
                                       ),
                             bw.seq = seq(sdpm/2,sdpm,length.out = 50),
                             a.rng = c(min(a.vals), max(a.vals)),
                             sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"))

# Oracle: adjust for all confounders
erfest_oracle <- ctseff(y = combined_data_covariates_outcome$logdeathrate,
                               a = combined_data_covariates_outcome$pm25,
                               x = covs,
                        bw.seq = seq(sdpm/2,sdpm,length.out = 50),
                        a.rng = c(min(a.vals), max(a.vals)),
                        sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"))

# Compute Causal Risk Ratio from the oracle
(exp(erfest_oracle$res$est[which.min(abs(a.vals - 12))])-0.01)/(exp(erfest_oracle$res$est[which.min(abs(a.vals - 9))])-0.01)

# Plots
data <- data.frame(
  a.vals = a.vals,
  Baseline = 100*(exp(erfest_baseline$res$est) - 0.01),
  Oracle = 100*(exp(erfest_oracle$res$est) - 0.01),
  SpatialCoordinates = 100*(exp(erfest_spatialcoord$res$est) - 0.01),
  IV_TPS = 100*(exp(erfest_IV_TPS$res$est) - 0.01),
  IV_GraphLaplacian = 100*(exp(erfest_IV_GraphLaplacian$res$est) - 0.01)
)

# Convert data to long format for ggplot
data_long <- data %>%
  pivot_longer(cols = -a.vals, names_to = "Method", values_to = "Value")

gmain <- ggplot(data_long, aes(x = a.vals, y = Value, color = Method)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black", 
                                "SpatialCoordinates" = "purple", 
                                "IV_TPS" = "darkgreen", 
                                "IV_GraphLaplacian" = "blue")) +
  labs(x = "PM Exposure Averaged over 2001-2013, µg/m³", 
       y = "Mortality Rate in 2014-2016, 1/(100 person-years)", 
       title = "ERC between air pollution exposure \n and mortality rate in US zipcodes",
       color = "Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(3,5.5)

data <- data.frame(
  a.vals = a.vals,
  Baseline_est = 100*(exp(erfest_baseline$res$est) - 0.01),
  Baseline_ci_ll = 100*(exp(erfest_baseline$res$ci.ll) - 0.01),
  Baseline_ci_ul = 100*(exp(erfest_baseline$res$ci.ul) - 0.01),
  Oracle_est = 100*(exp(erfest_oracle$res$est) - 0.01),
  Oracle_ci_ll = 100*(exp(erfest_oracle$res$ci.ll) - 0.01),
  Oracle_ci_ul = 100*(exp(erfest_oracle$res$ci.ul) - 0.01),
  SpatialCoordinates_est = 100*(exp(erfest_spatialcoord$res$est) - 0.01),
  SpatialCoordinates_ci_ll = 100*(exp(erfest_spatialcoord$res$ci.ll) - 0.01),
  SpatialCoordinates_ci_ul = 100*(exp(erfest_spatialcoord$res$ci.ul) - 0.01),
  IV_TPS_est = 100*(exp(erfest_IV_TPS$res$est) - 0.01),
  IV_TPS_ci_ll = 100*(exp(erfest_IV_TPS$res$ci.ll) - 0.01),
  IV_TPS_ci_ul = 100*(exp(erfest_IV_TPS$res$ci.ul) - 0.01),
  IV_GraphLaplacian_est = 100*(exp(erfest_IV_GraphLaplacian$res$est) - 0.01),
  IV_GraphLaplacian_ci_ll = 100*(exp(erfest_IV_GraphLaplacian$res$ci.ll) - 0.01),
  IV_GraphLaplacian_ci_ul = 100*(exp(erfest_IV_GraphLaplacian$res$ci.ul) - 0.01)
)

data_filtered <-  data
# Plot00: Oracle only
plot00 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # Baseline with confidence intervals
  geom_ribbon(aes(ymin = Oracle_ci_ll, ymax = Oracle_ci_ul), fill = "black", alpha = 0.2) +
  labs(x = "PM Exposure Averaged over 2001-2013, µg/m³", 
       y = "Mortality Rate, 1/(100 person-years)", 
       title = "Oracle",
       color = "Method") +
  scale_color_manual(values = c("Oracle" = "black")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(3,5.5)

# Plot 0: Baseline and Oracle
plot0 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # Baseline with confidence intervals
  geom_ribbon(aes(ymin = Baseline_ci_ll, ymax = Baseline_ci_ul), fill = "red", alpha = 0.2) +
  labs(x = "PM Exposure Averaged over 2001-2013, µg/m³", 
       y = "Mortality Rate, 1/(100 person-years)", 
       title = "Baseline",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(3,5.5)

# Plot 1: Baseline, Oracle, and Spatial Coordinates (with Confidence Bands)
plot1 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # Spatial Coordinates with confidence intervals
  geom_ribbon(aes(ymin = SpatialCoordinates_ci_ll, ymax = SpatialCoordinates_ci_ul), fill = "purple", alpha = 0.2) +
  geom_line(aes(y = SpatialCoordinates_est, color = "SpatialCoordinates"), size = 1) +
  labs(x = "PM Exposure Averaged over 2001-2013, µg/m³", 
       y = "Mortality Rate, 1/(100 person-years)", 
       title = "Spatial Coordinates",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black", 
                                "SpatialCoordinates" = "purple")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(3,5.5)

# Plot 2: Baseline, Oracle, and IV-TPS (with Confidence Bands)
plot2 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # IV-TPS with confidence intervals
  geom_ribbon(aes(ymin = IV_TPS_ci_ll, ymax = IV_TPS_ci_ul), fill = "green", alpha = 0.2) +
  geom_line(aes(y = IV_TPS_est, color = "IV_TPS"), size = 1) +
  labs(x = "PM Exposure Averaged over 2001-2013, µg/m³", 
       y = "Mortality Rate, 1/(100 person-years)", 
       title = "IV-TPS",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black", 
                                "IV_TPS" = "darkgreen")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+ 
  ylim(3,5.5)

# Plot 3: Baseline, Oracle, and IV-GraphLaplacian (with Confidence Bands)
plot3 <- ggplot(data_filtered, aes(x = a.vals)) +
  # Baseline
  geom_line(aes(y = Baseline_est, color = "Baseline"), size = 1) +
  # Oracle
  geom_line(aes(y = Oracle_est, color = "Oracle"), size = 1) +
  # IV-GraphLaplacian with confidence intervals
  geom_ribbon(aes(ymin = IV_GraphLaplacian_ci_ll, ymax = IV_GraphLaplacian_ci_ul), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = IV_GraphLaplacian_est, color = "IV_GraphLaplacian"), size = 1) +
  labs(x = "PM Exposure Averaged over 2001-2013, µg/m³", 
       y = "Mortality Rate, 1/(100 person-years)", 
       title = "IV-GraphLaplacian",
       color = "Method") +
  scale_color_manual(values = c("Baseline" = "red", 
                                "Oracle" = "black", 
                                "IV_GraphLaplacian" = "blue")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+ 
  ylim(3,5.5)

# Create the plot
png(
  filename = paste0('images/erfs.png'),
  width = 10000,
  height = 10000,
  res = 800
)
# Define a layout matrix
layout_matrix <- rbind(c(1,2),
                       c(3,4),
                       c(5,6))  

# Arrange the plots according to the layout
grid.arrange(gmain, plot00, plot0, plot1, plot2, plot3, 
             layout_matrix = layout_matrix)
dev.off()

###################################### SENSITIVITY ANALYSIS #################################
Ac_cols <- grep("^Ac", colnames(combined_data_covariates_outcome), value = TRUE)
# For each candidate, estimate the ERC
erfs_sens <- list()
for (Ac_col in Ac_cols){
  erfs_sens <- ctseff(y = combined_data_covariates_outcome$logdeathrate,
                      a = combined_data_covariates_outcome$pm25,
                      x = cbind(xmat, combined_data_covariates_outcome[[Ac_col]]),
                      bw.seq = seq(sdpm/2,sdpm,length.out = 50),
                      a.rng = c(min(a.vals), max(a.vals)),
                      sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.glm.interaction", "SL.mean"))
  erfs_sens[[Ac_col]] <- out
}
