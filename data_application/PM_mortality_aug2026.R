source('../funcs.R')
source('../R/basis_selection.R')
source('../R/truncated_effect.R')
source('../R/outer_crossfit.R')
source('../R/data_application_basis.R')
library(dplyr)
library(ggplot2)
library(sf)
library(Matrix)
library(stringr)
library(readr)
library(xtable)
library(mgcv)
library(RSpectra)

set.seed(123)

HIDDEN_PATH <- '' # Fill in lab path

# Outcome: zipcode-level all cause mortality 2011-2016
# Exposure: zipcode-level all source PM2.5 2001-2010
# Confounders: zipcode-level Census + gridMET + BRFSS 2000
# geography (shapefile): 2010

# Read in National Causal Dataset
aggregate_data <- read.csv(paste0(HIDDEN_PATH, 'projects/analytic/aggregated_2000-2016_medicare_mortality_pm25_zip/aggregate_data.csv'))

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
  file_path <- paste0(HIDDEN_PATH, "exposure/pm25/PM25_v2/annual/", year, ".rds")
  
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
zipcode_sf_polygon <- st_read(paste0(HIDDEN_PATH, 'data/shapefiles/zip_shape_files/2010/zip/polygon/ESRI10USZIP5_POLY_WGS84.shp'))
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
zipcode_sf_point <- st_read(paste0(HIDDEN_PATH, 'data/shapefiles/zip_shape_files/2010/zip/point/ESRI10USZIP5_POINT_WGS84.shp'))
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
#
# Estimation of truncated exposure effect, selecting scale of conf: truncated
# spatial basis (basis_k_max columns) + sequential stability criterion for
# basis selection, done separately per cutoff and per outer fold; nuisance
# estimation on the full A >= cutoff population;
# 5-fold cross-fitting; local-quadratic boundary regression (local_degree =
# 2). 
#
# basis_k_max = 50 and the n_uc candidate range mirror what's already written
# in the main text. NEITHER IS TIMING-VALIDATED AT THIS SCALE (n ~ 33,255).
# Run data_application/timing_test.R first and adjust basis_k_max / step /
# alpha below if a full run is not tractable: cost scales as candidates x
# 5 outer folds x 7 cutoffs x 2 IV methods, and each candidate itself needs an
# inner 2-fold SuperLearner fit. Increasing `step` (e.g. to 5, giving 4
# candidates instead of 20) is the first knob to pull if timing is too slow.

cutoffs <- seq(6, 12, 1)
y <- combined_data_covariates_outcome$deathrate
a <- combined_data_covariates_outcome$pm25

sl_library <- c("SL.gam", "SL.glm", "SL.mean", "SL.glm.interaction")
n_grid <- 25L
alpha <- 1
basis_k_max <- 50L

# Build the truncated bases once, on the full pre-merge combined_data extent
# (same extent L/adjacency_matrix above were built on), then row-align to the
# final n = 33,255 analysis population by ZIP -- the same alignment the
# fixed-k Ac_TPS/Ac_GraphLaplacian columns above already rely on via merge().
B_tps_full <- build_tps_basis_truncated(
  lat = combined_data$y, lon = combined_data$x, k_max = basis_k_max
)
B_gl_full <- build_gl_basis_truncated(L, k_max = basis_k_max)
zip_idx <- match(combined_data_covariates_outcome$ZIP, combined_data$ZIP)
stopifnot(!anyNA(zip_idx))
B_tps <- B_tps_full[zip_idx, , drop = FALSE]
B_gl  <- B_gl_full[zip_idx, , drop = FALSE]

# Full sweep n_uc in {1, ..., basis_k_max - 1}. select_contiguous_candidate()
# tests every candidate against n_uc_cands[1] and walks toward larger n_uc
#
# Cost: 49 candidates x 5 outer folds x 7 cutoffs x 2 IV methods (each
# candidate needs an inner 2-fold SuperLearner fit). If timing_test.R shows
# this isn't tractable, raise `step` first (e.g. step = 2L halves the
# candidate count) before shrinking basis_k_max itself.
n_uc_cands <- make_candidate_grid(
  m = basis_k_max, step = 1L,
  n_uc_star = 1L, min_confounded = 1L
)

lon_n <- (combined_data_covariates_outcome$x - min(combined_data_covariates_outcome$x)) /
  (max(combined_data_covariates_outcome$x) - min(combined_data_covariates_outcome$x))
lat_n <- (combined_data_covariates_outcome$y - min(combined_data_covariates_outcome$y)) /
  (max(combined_data_covariates_outcome$y) - min(combined_data_covariates_outcome$y))

set.seed(20260827)
folds <- make_random_folds(length(y))

fit_method <- function(method, cutoff) {
  t0 <- Sys.time()
  fit <- tryCatch({
    if (method == "oracle") {
      estimate_crossfit_truncated_effect(
        y = y, a = a, w = covs, cutoff = cutoff, outer_folds = folds,
        final_nuisance_args = list(sl_library = sl_library),
        n_grid = n_grid, constrain = TRUE, local_degree = 2L
      )
    } else if (method == "baseline") {
      estimate_crossfit_truncated_effect(
        y = y, a = a, w = xmat, cutoff = cutoff, outer_folds = folds,
        final_nuisance_args = list(sl_library = sl_library),
        n_grid = n_grid, constrain = TRUE, local_degree = 2L
      )
    } else if (method == "spatialcoord") {
      estimate_crossfit_truncated_effect(
        y = y, a = a, w = cbind(xmat, Longitude = lon_n, Latitude = lat_n),
        cutoff = cutoff, outer_folds = folds,
        final_nuisance_args = list(sl_library = sl_library),
        n_grid = n_grid, constrain = TRUE, local_degree = 2L
      )
    } else {
      B_method <- if (method == "IV-TPS") B_tps else B_gl
      estimate_selected_truncated_effect(
        y = y, a = a, x = xmat, B = B_method, cutoff = cutoff,
        n_uc_cands = n_uc_cands, outer_folds = folds, alpha = alpha,
        instrument_scale = "small",
        candidate_args = list(
          nuisance_args = list(sl_library = sl_library),
          n_grid = n_grid, constrain = TRUE, local_degree = 2L
        ),
        final_nuisance_args = list(sl_library = sl_library),
        n_grid = n_grid, constrain = TRUE, local_degree = 2L
      )
    }
  }, error = function(e) {
    message(sprintf("%s @ cutoff=%d failed: %s", method, cutoff, e$message))
    NULL
  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  print(sprintf("%s @ cutoff=%d: %.1f sec", method, cutoff, elapsed))
  if (is.null(fit)) {
    return(data.frame(
      method = method, cutoff = cutoff, psi = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      n_uc_mean = NA_real_, time_s = elapsed
    ))
  }
  n_uc_mean <- if (!is.null(fit$selected_n_uc)) mean(fit$selected_n_uc) else NA_real_
  data.frame(
    method = method, cutoff = cutoff, psi = fit$psi,
    ci_lower = fit$confidence_interval[1], ci_upper = fit$confidence_interval[2],
    n_uc_mean = n_uc_mean, time_s = elapsed
  )
}

methods_data_app <- c("oracle", "baseline", "spatialcoord", "IV-TPS", "IV-GraphLaplacian")
results_truncated_effect <- do.call(rbind, lapply(cutoffs, function(cutoff) {
  out <- do.call(rbind, lapply(methods_data_app, fit_method, cutoff = cutoff))
  gc()
  out
}))

dir.create("results_aug2026", showWarnings = FALSE)
write.csv(
  results_truncated_effect,
  file = "results_aug2026/truncated_effect_estimates.csv",
  row.names = FALSE
)
print(results_truncated_effect)


########################################## PLOTS OF TRUNCATED EXPOSURE EFFECT ESTIMATES ##############################################

# Build the plotting dataframe directly from the new pipeline's output
# (results_aug2026/truncated_effect_estimates.csv, written above), rather
# than the legacy per-method-per-cutoff .RData files under results_Mar24/
# that the old ctseff()-based pipeline produced. Reading from the CSV
# (rather than reusing results_truncated_effect in memory) keeps this block
# runnable on its own, and keeps every result under results_aug2026/.
df <- read.csv("results_aug2026/truncated_effect_estimates.csv") %>%
  rename(point_est = psi, lower = ci_lower, upper = ci_upper) %>%
  mutate(method = factor(method, levels = methods_data_app))

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
  filename = paste0('images/estimates_cis_2010_aug2026.png'),
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
                                "IV-GraphLaplacian" = "black")) +
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
methods <- methods_data_app
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


############################################# SENSITIVITY ANALYSIS: EFFECT ESTIMATE VS n_uc #############################################
#
# This holds n_uc FIXED at the same value across
# every outer fold (rather than the per-fold sequential-stability choice) and
# traces out how the truncated-effect estimate and its CI move as that fixed
# n_uc varies -- using the same estimate_crossfit_truncated_effect() call
# already used for oracle/baseline/spatialcoord above, with a fixed Ac
# appended as the extra adjustment column via project_Ac() (the same
# projection the selection procedure uses internally, here applied on the
# full sample since n_uc is fixed rather than data-adaptively chosen per
# fold). Swept across cutoffs 6-9 and with/without the 10 measured
# covariates. Oracle and baseline(fullcov) are reused from
# results_truncated_effect above rather than recomputed.

sens_cutoffs <- 6:9
n_uc_grid <- seq(5, 45, by = 5)
xmat_nocov <- data.frame(Intercept = rep(1, nrow(xmat)))
covset_list <- list(fullcov = xmat, nocov = xmat_nocov)

fit_fixed <- function(w, cutoff) {
  tryCatch(
    estimate_crossfit_truncated_effect(
      y = y, a = a, w = w, cutoff = cutoff, outer_folds = folds,
      final_nuisance_args = list(sl_library = sl_library),
      n_grid = n_grid, constrain = TRUE, local_degree = 2L
    ),
    error = function(e) {
      message(sprintf("cutoff=%d failed: %s", cutoff, e$message))
      NULL
    }
  )
}

record <- function(cutoff, covset, method, n_uc, psi, ci_lower, ci_upper) {
  data.frame(
    cutoff = cutoff, covset = covset, method = method, n_uc = n_uc,
    psi = psi, ci_lower = ci_lower, ci_upper = ci_upper
  )
}

record_fit <- function(cutoff, covset, method, n_uc, fit) {
  if (is.null(fit)) {
    record(cutoff, covset, method, n_uc, NA_real_, NA_real_, NA_real_)
  } else {
    record(cutoff, covset, method, n_uc, fit$psi,
           fit$confidence_interval[1], fit$confidence_interval[2])
  }
}

sens_rows <- list()
for (cutoff in sens_cutoffs) {
  oracle_row <- results_truncated_effect[
    results_truncated_effect$method == "oracle" &
      results_truncated_effect$cutoff == cutoff, ]
  baseline_fullcov_row <- results_truncated_effect[
    results_truncated_effect$method == "baseline" &
      results_truncated_effect$cutoff == cutoff, ]

  sens_rows[[length(sens_rows) + 1]] <- record(
    cutoff, "fullcov", "oracle", NA_integer_,
    oracle_row$psi, oracle_row$ci_lower, oracle_row$ci_upper
  )
  sens_rows[[length(sens_rows) + 1]] <- record(
    cutoff, "nocov", "oracle", NA_integer_,
    oracle_row$psi, oracle_row$ci_lower, oracle_row$ci_upper
  )
  sens_rows[[length(sens_rows) + 1]] <- record(
    cutoff, "fullcov", "baseline", NA_integer_,
    baseline_fullcov_row$psi, baseline_fullcov_row$ci_lower, baseline_fullcov_row$ci_upper
  )

  baseline_nocov_fit <- fit_fixed(xmat_nocov, cutoff)
  sens_rows[[length(sens_rows) + 1]] <- record_fit(
    cutoff, "nocov", "baseline", NA_integer_, baseline_nocov_fit
  )

  for (covset_name in names(covset_list)) {
    w_base <- covset_list[[covset_name]]

    for (n_uc_val in n_uc_grid) {
      for (method in c("IV-TPS", "IV-GraphLaplacian")) {
        B_method <- if (method == "IV-TPS") B_tps else B_gl
        Ac <- project_Ac(
          B = B_method, A = a, n_uc = n_uc_val,
          train_idx = seq_along(a), predict_idx = seq_along(a),
          instrument_scale = "small"
        )$Ac
        fit <- fit_fixed(cbind(w_base, Ac = Ac), cutoff)
        sens_rows[[length(sens_rows) + 1]] <- record_fit(
          cutoff, covset_name, method, n_uc_val, fit
        )
      }
    }
  }
  gc()
}

sens_df <- do.call(rbind, sens_rows)
write.csv(sens_df, "results_aug2026/sensitivity_nuc_estimates.csv", row.names = FALSE)

# ---- Plots: estimate + CI vs n_uc, one panel per cutoff x covariate set ----
sens_df_pct <- sens_df %>%
  mutate(point_est = 100 * (psi - 1), lower = 100 * (ci_lower - 1), upper = 100 * (ci_upper - 1))

sens_lines <- sens_df_pct %>% filter(method %in% c("IV-TPS", "IV-GraphLaplacian"))
sens_refs  <- sens_df_pct %>% filter(method %in% c("oracle", "baseline"))

dir.create("images/sensitivity_nuc", showWarnings = FALSE, recursive = TRUE)

for (cutoff in sens_cutoffs) {
  for (covset_name in names(covset_list)) {
    df_lines <- sens_lines %>% filter(cutoff == !!cutoff, covset == covset_name)
    df_refs  <- sens_refs  %>% filter(cutoff == !!cutoff, covset == covset_name)

    covset_title <- ifelse(
      covset_name == "fullcov",
      "adjusting for measured covariates",
      "without measured covariates"
    )

    g <- ggplot() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_line(
        data = df_lines,
        aes(x = n_uc, y = point_est, color = method), linewidth = 0.8
      ) +
      geom_point(
        data = df_lines,
        aes(x = n_uc, y = point_est, color = method), size = 2
      ) +
      geom_errorbar(
        data = df_lines,
        aes(x = n_uc, ymin = lower, ymax = upper, color = method), width = 1
      ) +
      geom_hline(
        data = df_refs,
        aes(yintercept = point_est, linetype = method), color = "gray35", linewidth = 0.7
      ) +
      scale_x_continuous(breaks = n_uc_grid) +
      scale_color_manual(
        name = NULL,
        values = c("IV-TPS" = "#1b9e77", "IV-GraphLaplacian" = "#d95f02")
      ) +
      scale_linetype_manual(
        name = NULL,
        values = c("baseline" = "solid", "oracle" = "dotdash")
      ) +
      labs(
        x = expression(n[uc]),
        y = "Reduction in Mortality Rate (%)",
        title = bquote(atop(
          "Impact of enforcing average " * PM[2.5] * " (2001-2010) below " *
            .(cutoff) ~ mu * g / m^3 ~ " on",
          "mortality rate in the Medicare population (2011-2016), " * .(covset_title)
        ))
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "top"
      )

    ggsave(
      g,
      file = sprintf("images/sensitivity_nuc/estimates_cis_nuc_cutoff_%d_%s.png", cutoff, covset_name),
      width = 1500, height = 1000, units = "px", dpi = 240
    )
  }
}
