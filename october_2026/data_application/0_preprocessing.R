library(dplyr)
library(ggplot2)
library(sf)
library(Matrix)
library(stringr)
library(readr)
library(xtable)
library(mgcv)
library(RSpectra)
library(here)
library(gridExtra)

set.seed(123)

HIDDEN_PATH <- '' # Fill in lab path
base_path <- paste0(HIDDEN_PATH, 'projects/wildfirePM-mortality-spatialconf_IV/spatial-scale-conf/iv-spatialconfounding-sophi890-issue54')
# Run from october_2026/data_application/; shared modules live in ../R/.
# base_path is only used for lab-server data and large outputs.
source('../R/basis_selection.R')
source('../R/truncated_effect.R')
source('../R/basis_construction.R')
results_dir <- file.path(base_path, "results_oct2026")

basis_k_max <- 50L
basis_run_tag <- paste0("k", basis_k_max, "_allzip_knn")

# Use the smallest value in this sequence that produces one connected graph for GL
knn_k_candidates <- seq_len(30L)

# for example plots later
n_c_candidate_grid <- c(
  49L, 45L, 40L, 35L, 30L, 25L,
  20L, 15L, 10L, 6L, 2L
)

n_uc_cands <- sort(as.integer(basis_k_max - n_c_candidate_grid))

# These diagnostics are inexpensive projections; evaluate more dimensions than
# will necessarily be used in the outcome analysis.
n_c_diagnostic_grid <- sort(unique(c(
  6L, 10L, 15L, 20L, 25L, 30L, 40L, 49L
)))

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

# Extract outcome data from 2011-2016
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
rm(aggregate_data)
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
length(unique(zipcode_sf_polygon$ZIP)) # 30377
nrow(zipcode_sf_polygon) # 30377
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
length(unique(zipcode_sf_point$ZIP)) # 41033
zipcode_sf_point <- st_as_sf(zipcode_sf_point)

# average over years before 2010 to get PM exposure. 
total_pm <- aggregate(pm25~ZIP, data = combined_pm_data, 
                      FUN = mean,
                      na.rm = T)
mean(total_pm$ZIP %in% c(zipcode_sf_polygon$ZIP, zipcode_sf_point$ZIP)) # 97 percent.
rm(combined_pm_data, pm_data_list, pm_data)

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

polygon_data_clean <- polygon_data %>%
  dplyr::select(ZIP, pm25)

point_data_clean <- point_data %>%
  dplyr::select(ZIP, pm25)

combined_data <- rbind(
  polygon_data_clean,
  point_data_clean
)

combined_data <- combined_data %>%
  group_by(ZIP) %>%
  filter(
    !(
      st_geometry_type(geometry) %in% c("POINT", "MULTIPOINT") &
        any(
          st_geometry_type(geometry) %in%
            c("POLYGON", "MULTIPOLYGON")
        )
    )
  ) %>%
  ungroup()

stopifnot(!anyDuplicated(combined_data$ZIP))

nrow(combined_data) # 40471

# -----------------------------------------------------------------------------
# Merge exposure, covariates, and outcome
# -----------------------------------------------------------------------------

analysis_sf <- merge(
  combined_data,
  covariate_data_sub,
  by.x = "ZIP",
  by.y = "zip"
)

analysis_sf <- merge(
  analysis_sf,
  outcome_data_sub,
  by.x = "ZIP",
  by.y = "zip"
)

message("Rows after the data merge: ", nrow(analysis_sf))
# 33254

covariate_names <- c(
  "mean_bmi",
  "smoke_rate",
  "hispanic",
  "pct_blk",
  "medhouseholdincome",
  "medianhousevalue",
  "poverty",
  "education",
  "popdensity",
  "pct_owner_occ",
  "summer_tmmx",
  "winter_tmmx",
  "summer_rmax",
  "winter_rmax"
)
# Define the complete analysis population before constructing either basis.
required_numeric <- c(
  "pm25",
  "deathrate",
  covariate_names
)

analysis_values <- sf::st_drop_geometry(analysis_sf)[
  , required_numeric, drop = FALSE
]


rm(
  analysis_values,
  covariate_data_sub,
  outcome_data_sub,
  zipcode_sf_polygon,
  zipcode_sf_point,
  duplicate_zipcodes,
  duplicated_sf,
  merged_duplicated_sf,
  non_duplicated_sf,
  polygon_data,
  point_data,
  polygon_data_clean,
  point_data_clean,
  combined_data,
  total_pm,
  pm25_limits
)

gc()

# Use a continental-US projected coordinate system for spatial operations and
# the TPS coordinates. Coordinates extracted below are converted to kms
analysis_sf <- sf::st_transform(analysis_sf, crs = 5070)

# -----------------------------------------------------------------------------
# Sparse connected components
# -----------------------------------------------------------------------------

find_sparse_components <- function(adj) {
  n <- length(adj)
  membership <- integer(n)
  component_sizes <- integer(0)
  queue <- integer(n)
  component_id <- 0L
  
  for (start in seq_len(n)) {
    if (membership[start] != 0L) {
      next
    }
    
    component_id <- component_id + 1L
    queue_head <- 1L
    queue_tail <- 1L
    queue[1L] <- start
    membership[start] <- component_id
    component_size <- 0L
    
    while (queue_head <= queue_tail) {
      vertex <- queue[queue_head]
      queue_head <- queue_head + 1L
      component_size <- component_size + 1L
      
      neighbors <- adj[[vertex]]
      new_neighbors <- neighbors[membership[neighbors] == 0L]
      
      if (length(new_neighbors) > 0L) {
        membership[new_neighbors] <- component_id
        
        new_positions <- seq.int(
          queue_tail + 1L,
          queue_tail + length(new_neighbors)
        )
        
        queue[new_positions] <- new_neighbors
        queue_tail <- queue_tail + length(new_neighbors)
      }
    }
    
    component_sizes[component_id] <- component_size
  }
  
  list(
    membership = membership,
    csize = component_sizes,
    no = component_id
  )
}

# Convert a symmetric sparse matrix into the adjacency-list representation used
# by find_sparse_components().
sparse_matrix_to_adjacency <- function(W) {
  # Expand a symmetric sparse representation before reading its compressed
  # columns; otherwise a dsCMatrix would expose only its stored triangle.
  W_c <- methods::as(W, "generalMatrix")
  W_c <- methods::as(W_c, "CsparseMatrix")
  n <- nrow(W_c)
  
  lapply(seq_len(n), function(vertex) {
    first <- W_c@p[vertex] + 1L
    last <- W_c@p[vertex + 1L]
    
    if (first > last) {
      integer(0)
    } else {
      W_c@i[seq.int(first, last)] + 1L
    }
  })
}

# -----------------------------------------------------------------------------
# One proximity graph
# -----------------------------------------------------------------------------

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop(
    "The RANN package is required to construct the all-ZIP kNN graph.",
    call. = FALSE
  )
}

graph_points <- sf::st_point_on_surface(analysis_sf)
graph_coordinates <- sf::st_coordinates(graph_points)
n_graph <- nrow(graph_coordinates)

stopifnot(
  n_graph == nrow(analysis_sf),
  ncol(graph_coordinates) >= 2L,
  all(is.finite(graph_coordinates[, 1:2, drop = FALSE]))
)

if (n_graph <= max(knn_k_candidates)) {
  stop(
    "The analysis sample is too small for knn_k_candidates.",
    call. = FALSE
  )
}

duplicate_coordinate_count <- sum(
  duplicated(graph_coordinates[, 1:2, drop = FALSE])
)

if (duplicate_coordinate_count > 0L) {
  warning(
    paste0(
      duplicate_coordinate_count,
      " ZIPs share a representative-point coordinate with an earlier ZIP."
    ),
    call. = FALSE
  )
} # 374 ZIPs share a representative-point coordinate with an earlier ZIP. 

nearest <- RANN::nn2(
  data = graph_coordinates,
  query = graph_coordinates,
  k = max(knn_k_candidates) + 1L
)

# Remove each observation from its own neighbor order. This is robust to
# duplicate coordinates, for which the self-match need not be column 1.
neighbor_matrix <- t(vapply(seq_len(n_graph), function(i) {
  candidates <- nearest$nn.idx[i, ]
  candidates <- unique(candidates[candidates != i])
  if (length(candidates) < max(knn_k_candidates)) {
    stop(
      "A ZIP does not have enough distinct neighbors for knn_k_candidates.",
      call. = FALSE
    )
  }
  
  as.integer(head(candidates, max(knn_k_candidates)))
}, integer(max(knn_k_candidates))))

connectivity_records <- vector("list", length(knn_k_candidates))
W <- NULL
selected_k <- NA_integer_
selected_edge_from <- selected_edge_to <- NULL
component_info <- NULL

# Here we increase knn until GL graph is connected
for (candidate_index in seq_along(knn_k_candidates)) {
  candidate_k <- knn_k_candidates[candidate_index]
  edge_from_candidate <- rep(seq_len(n_graph), each = candidate_k)
  edge_to_candidate <- as.vector(
    t(neighbor_matrix[, seq_len(candidate_k), drop = FALSE])
  )
  
  W_directed <- Matrix::sparseMatrix(
    i = edge_from_candidate,
    j = edge_to_candidate,
    x = rep.int(1, length(edge_from_candidate)),
    dims = c(n_graph, n_graph)
  )
  
  W_candidate <- 1 * ((W_directed + Matrix::t(W_directed)) > 0)
  diag(W_candidate) <- 0
  W_candidate <- Matrix::drop0(W_candidate)
  
  candidate_components <- find_sparse_components(
    sparse_matrix_to_adjacency(W_candidate)
  )
  
  connectivity_records[[candidate_index]] <- data.frame(
    k = candidate_k,
    n_components = candidate_components$no,
    largest_component = max(candidate_components$csize)
  )
  
  message(
    "Symmetric kNN graph with k = ",
    candidate_k,
    ": ",
    candidate_components$no,
    " connected component(s)."
  )
  
  if (candidate_components$no == 1L) {
    W <- W_candidate
    selected_k <- candidate_k
    selected_edge_from <- edge_from_candidate
    selected_edge_to <- edge_to_candidate
    component_info <- candidate_components
    break
  }
}

connectivity_diagnostics <- do.call(
  rbind,
  connectivity_records[!vapply(connectivity_records, is.null, logical(1))]
)

if (is.null(W) || !is.finite(selected_k)) {
  stop(
    paste0(
      "No connected all-ZIP graph was obtained. Increase ",
      "knn_k_candidates only after inspecting the geography."
    ),
    call. = FALSE
  )
}

degree <- Matrix::rowSums(W)
edge_distance_km <- sqrt(
  rowSums(
    (
      graph_coordinates[selected_edge_from, 1:2, drop = FALSE] -
        graph_coordinates[selected_edge_to, 1:2, drop = FALSE]
    )^2
  )
) / 1000

component_summary <- list(
  n_before = n_graph,
  n_after = n_graph,
  n_removed = 0L,
  n_components = component_info$no,
  largest_component_sizes = component_info$csize
)

graph_summary <- list(
  graph_type = "symmetric_binary_knn",
  coordinate_crs = 5070L,
  n_vertices = n_graph,
  n_undirected_edges = as.numeric(sum(W) / 2),
  k_candidates = knn_k_candidates,
  selected_k = selected_k,
  duplicate_coordinate_count = duplicate_coordinate_count,
  connectivity_diagnostics = connectivity_diagnostics,
  degree_quantiles = stats::quantile(
    degree,
    c(0, 0.25, 0.5, 0.75, 1),
    names = TRUE
  ),
  directed_neighbor_distance_km_quantiles = stats::quantile(
    edge_distance_km,
    c(0, 0.5, 0.9, 0.99, 1),
    names = TRUE
  )
)

message(
  "Selected connected all-ZIP graph with k = ",
  selected_k,
  "; degree quantiles: ",
  paste(signif(graph_summary$degree_quantiles, 3L), collapse = ", ")
)

message(
  "Directed-neighbor distance quantiles in km: ",
  paste(
    signif(
      graph_summary$directed_neighbor_distance_km_quantiles,
      3L
    ),
    collapse = ", "
  )
)

if (max(edge_distance_km) > 500) {
  warning(
    paste0(
      "The selected graph contains a directed kNN edge longer than 500 km. ",
      "Confirm that the analysis data contain only the intended contiguous-US ",
      "ZIP codes before fitting outcome models."
    ),
    call. = FALSE
  )
}

stopifnot(
  nrow(W) == nrow(analysis_sf),
  isTRUE(isSymmetric(W)),
  component_info$no == 1L,
  all(degree >= selected_k),
  component_summary$n_removed == 0L
)

L <- Matrix::Diagonal(
  x = degree
) - W

rm(
  W_directed,
  W_candidate,
  graph_points,
  graph_coordinates,
  nearest,
  neighbor_matrix,
  edge_from_candidate,
  edge_to_candidate,
  selected_edge_from,
  selected_edge_to,
  edge_distance_km,
  degree,
  candidate_components,
  connectivity_records
)

gc()

# -----------------------------------------------------------------------------
# Projected coordinates and analysis matrices
# -----------------------------------------------------------------------------

analysis_points <- sf::st_point_on_surface(analysis_sf)
basis_coordinates <- sf::st_coordinates(analysis_points)

analysis_sf$basis_x <- basis_coordinates[, "X"] / 1000
analysis_sf$basis_y <- basis_coordinates[, "Y"] / 1000

rm(analysis_points, basis_coordinates)
gc()

covs_raw <- sf::st_drop_geometry(
  dplyr::select(analysis_sf, all_of(covariate_names))
)

covs_scaled <- scale(covs_raw)
covariate_center <- attr(covs_scaled, "scaled:center")
covariate_scale <- attr(covs_scaled, "scaled:scale")
covs <- as.data.frame(covs_scaled)

if (!all(is.finite(as.matrix(covs)))) {
  stop(
    "The scaled covariate matrix contains non-finite values. Check for a zero-variance covariate.",
    call. = FALSE
  )
}

# Exclude the four measured spatial meteorological confounders from xmat.
xmat <- covs %>%
  select(
    -summer_tmmx,
    -winter_tmmx,
    -summer_rmax,
    -winter_rmax
  )

y <- analysis_sf$deathrate
a <- analysis_sf$pm25
zip <- as.character(analysis_sf$ZIP)

scale_01 <- function(z) {
  z_range <- range(z)
  
  if (!is.finite(diff(z_range)) || diff(z_range) <= 0) {
    stop("Cannot rescale a constant or non-finite coordinate.", call. = FALSE)
  }
  
  (z - z_range[1L]) / diff(z_range)
}

lon_n <- scale_01(analysis_sf$basis_x)
lat_n <- scale_01(analysis_sf$basis_y)

# -----------------------------------------------------------------------------
# Construct the two truncated bases
# -----------------------------------------------------------------------------

message(
  "Constructing the ",
  basis_k_max,
  "-column truncated TPS basis."
)

B_tps <- build_tps_basis_truncated(
  lat = analysis_sf$basis_y,
  lon = analysis_sf$basis_x,
  k_max = basis_k_max
)

gc()

message(
  "Constructing the ",
  basis_k_max,
  "-column truncated graph-Laplacian basis."
)

B_gl <- build_gl_basis_truncated(
  L = L,
  k_max = basis_k_max
)

gc()

# -----------------------------------------------------------------------------
# Validate dimensions, orthogonality, and smooth-to-rough ordering
# -----------------------------------------------------------------------------

stopifnot(
  identical(dim(B_tps), c(length(y), basis_k_max)),
  identical(dim(B_gl), c(length(y), basis_k_max))
)

roughness_tps <- attr(B_tps, "roughness")

if (is.null(roughness_tps)) {
  stop(
    "B_tps has no roughness attribute. Check data_application_basis.R.",
    call. = FALSE
  )
}

roughness_gl <- colSums(B_gl * (L %*% B_gl)) /
  colSums(B_gl^2)

ordering_tolerance <- 1e-8

stopifnot(
  length(roughness_tps) == basis_k_max,
  length(roughness_gl) == basis_k_max,
  all(diff(roughness_tps) >= -ordering_tolerance),
  all(diff(roughness_gl) >= -ordering_tolerance)
)

tps_orthogonality_error <- max(
  abs(crossprod(B_tps) - diag(basis_k_max))
)

gl_orthogonality_error <- max(
  abs(crossprod(B_gl) - diag(basis_k_max))
)

message(
  "Maximum TPS orthogonality error: ",
  signif(tps_orthogonality_error, 4L)
)

message(
  "Maximum graph-Laplacian orthogonality error: ",
  signif(gl_orthogonality_error, 4L)
)

if (tps_orthogonality_error > 1e-6) {
  warning("The TPS basis is not numerically orthonormal.", call. = FALSE)
}

if (gl_orthogonality_error > 1e-6) {
  warning(
    "The graph-Laplacian basis is not numerically orthonormal.",
    call. = FALSE
  )
}

# The graph constant is column 1. Because n_uc_cands is
# increasing, the selection procedure starts with 49 confounded directions
# (n_uc = 1) and sequentially permits more retained directions to enter B^uc
stopifnot(
  all(n_uc_cands >= 1L),
  all(n_uc_cands < basis_k_max),
  !is.unsorted(n_uc_cands, strictly = TRUE),
  identical(
    sort(as.integer(basis_k_max - n_uc_cands), decreasing = TRUE),
    as.integer(n_c_candidate_grid)
  )
)

# -----------------------------------------------------------------------------
# Diagnose residual association with the omitted confounders
# -----------------------------------------------------------------------------

omitted_weather_names <- c(
  "summer_tmmx",
  "winter_tmmx",
  "summer_rmax",
  "winter_rmax"
)

X_diagnostic <- cbind(
  intercept = 1,
  as.matrix(xmat)
)
X_diagnostic_qr <- qr(X_diagnostic)

weather_residuals <- qr.resid(
  X_diagnostic_qr,
  as.matrix(covs[, omitted_weather_names, drop = FALSE])
)

diagnose_basis_leakage <- function(B, basis_name, n_c_values) {
  do.call(rbind, lapply(n_c_values, function(n_c) {
    if (n_c < 1L || n_c >= ncol(B)) {
      stop("Every diagnostic n_c must be between 1 and ncol(B) - 1.")
    }
    
    Bc <- B[, seq_len(n_c), drop = FALSE]
    Ac <- qr.fitted(qr(Bc), a)
    Auc <- a - Ac
    Auc_residual <- qr.resid(X_diagnostic_qr, Auc)
    
    partial_correlations <- drop(
      stats::cor(Auc_residual, weather_residuals)
    )
    
    names(partial_correlations) <- paste0(
      "partial_cor_",
      omitted_weather_names
    )
    
    weather_design <- cbind(
      intercept = 1,
      weather_residuals
    )
    Auc_weather_residual <- qr.resid(qr(weather_design), Auc_residual)
    Auc_residual_ss <- sum(Auc_residual^2)
    
    joint_weather_R2 <- if (Auc_residual_ss > 0) {
      1 - sum(Auc_weather_residual^2) / Auc_residual_ss
    } else {
      NA_real_
    }
    
    data.frame(
      basis = basis_name,
      n_c = as.integer(n_c),
      n_uc_internal = as.integer(ncol(B) - n_c),
      exposure_R2 = 1 -
        sum(Auc^2) / sum((a - mean(a))^2),
      Auc_sd = stats::sd(Auc),
      Auc_residual_sd = stats::sd(Auc_residual),
      joint_weather_partial_R2 = joint_weather_R2,
      max_abs_weather_partial_cor = max(abs(partial_correlations)),
      as.list(partial_correlations),
      check.names = FALSE
    )
  }))
}

basis_leakage_diagnostics <- dplyr::bind_rows(
  diagnose_basis_leakage(
    B = B_tps,
    basis_name = "TPS",
    n_c_values = n_c_diagnostic_grid
  ),
  diagnose_basis_leakage(
    B = B_gl,
    basis_name = "GraphLaplacian",
    n_c_values = n_c_diagnostic_grid
  )
) %>%
  dplyr::arrange(n_c, basis)

print(basis_leakage_diagnostics)

leakage_diagnostics_path <- file.path(
  results_dir,
  paste0(
    "basis_confounder_leakage_diagnostics_",
    basis_run_tag,
    ".csv"
  )
)

utils::write.csv(
  basis_leakage_diagnostics,
  leakage_diagnostics_path,
  row.names = FALSE
)

rm(
  X_diagnostic,
  X_diagnostic_qr,
  weather_residuals
)
gc()

# The graph matrices are unnecessary downstream once B_gl constructed
rm(W, L)
gc()

# -----------------------------------------------------------------------------
# Save preprocessing outputs
# -----------------------------------------------------------------------------

analysis_data <- sf::st_drop_geometry(analysis_sf)
analysis_data$row_id <- seq_len(nrow(analysis_data))

basis_map_sf <- dplyr::select(analysis_sf, ZIP)
basis_map_sf$row_id <- seq_len(nrow(basis_map_sf))
basis_map_sf <- dplyr::select(basis_map_sf, row_id, ZIP)

stopifnot(
  identical(zip, as.character(basis_map_sf$ZIP)),
  nrow(analysis_data) == nrow(B_tps),
  nrow(analysis_data) == nrow(B_gl)
)

preprocessed <- list(
  analysis_data = analysis_data,
  zip = zip,
  y = y,
  a = a,
  covs = covs,
  xmat = xmat,
  covariate_names = covariate_names,
  covariate_center = covariate_center,
  covariate_scale = covariate_scale,
  lon_n = lon_n,
  lat_n = lat_n,
  B_tps = B_tps,
  B_gl = B_gl,
  roughness_tps = roughness_tps,
  roughness_gl = roughness_gl,
  basis_k_max = basis_k_max,
  n_uc_cands = n_uc_cands,
  n_c_candidate_grid = n_c_candidate_grid,
  n_c_diagnostic_grid = n_c_diagnostic_grid,
  basis_leakage_diagnostics = basis_leakage_diagnostics,
  component_summary = component_summary,
  graph_summary = graph_summary,
  basis_run_tag = basis_run_tag,
  coordinate_crs = 5070L,
  tps_orthogonality_error = tps_orthogonality_error,
  gl_orthogonality_error = gl_orthogonality_error
)

preprocessed_path <- file.path(
  results_dir,
  paste0(
    "preprocessed_analysis_and_bases_",
    basis_run_tag,
    ".rds"
  )
)

map_path <- file.path(
  results_dir,
  paste0(
    "basis_map_geometry_",
    basis_run_tag,
    ".rds"
  )
)

saveRDS(
  preprocessed,
  file = preprocessed_path,
  compress = "gzip"
)

saveRDS(
  basis_map_sf,
  file = map_path,
  compress = "gzip"
)

message(
  "Saved ",
  length(y),
  " observations and two ",
  basis_k_max,
  "-column spatial bases."
)

message("Modeling bundle: ", preprocessed_path)
message("Map geometry: ", map_path)
message("Leakage diagnostics: ", leakage_diagnostics_path)

# -----------------------------------------------------------------------------
# Low-memory plots of large- and small-scale retained basis functions
# -----------------------------------------------------------------------------

basis_plot_dir <- file.path(
  results_dir,
  paste0("basis_checks_", basis_run_tag)
)

dir.create(
  basis_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

plot_basis_points <- function(
    values,
    basis_x,
    basis_y,
    title,
    subtitle,
    filename
) {
  stopifnot(
    length(values) == length(basis_x),
    length(values) == length(basis_y),
    all(is.finite(values))
  )
  
  # Use the same magnitude for the positive and negative color limits.
  z_limit <- max(abs(values))
  
  if (!is.finite(z_limit) || z_limit == 0) {
    z_limit <- 1
  }
  
  palette <- grDevices::hcl.colors(
    n = 101L,
    palette = "Blue-Red 3",
    rev = TRUE
  )
  
  color_breaks <- seq(
    -z_limit,
    z_limit,
    length.out = length(palette) + 1L
  )
  
  color_index <- findInterval(
    values,
    color_breaks,
    all.inside = TRUE
  )
  
  grDevices::png(
    filename = filename,
    width = 1800,
    height = 1100,
    res = 200,
    bg = "white"
  )
  
  # Ensure the graphics device closes even if plotting fails.
  on.exit({
    grDevices::dev.off()
    invisible(gc())
  }, add = TRUE)
  
  graphics::par(
    mar = c(0.5, 0.5, 3.5, 0.5),
    xaxs = "i",
    yaxs = "i"
  )
  
  graphics::plot(
    basis_x,
    basis_y,
    col = palette[color_index],
    pch = 15,
    cex = 0.35,
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = title,
    sub = subtitle
  )
  
  legend_positions <- c(1L, 51L, 101L)
  legend_values <- signif(c(-z_limit, 0, z_limit), 3L)
  
  graphics::legend(
    "bottomleft",
    legend = legend_values,
    fill = palette[legend_positions],
    title = "Basis value",
    horiz = TRUE,
    bty = "n",
    cex = 0.8
  )
  
  invisible(NULL)
}

# TPS columns 1--3 are the affine null space. Column 4 is therefore the
# smoothest non-affine TPS direction.
plot_basis_points(
  values = B_tps[, 4L],
  basis_x = analysis_data$basis_x,
  basis_y = analysis_data$basis_y,
  title = "TPS: large-scale basis function",
  subtitle = paste(
    "Column 4; roughness =",
    signif(roughness_tps[4L], 4L)
  ),
  filename = file.path(
    basis_plot_dir,
    "TPS_large_scale_column_04.png"
  )
)

# The final column is the smallest-scale direction among the retained TPS
# columns.
plot_basis_points(
  values = B_tps[, basis_k_max],
  basis_x = analysis_data$basis_x,
  basis_y = analysis_data$basis_y,
  title = "TPS: smallest-scale retained basis function",
  subtitle = paste(
    "Column",
    basis_k_max,
    "; roughness =",
    signif(roughness_tps[basis_k_max], 4L)
  ),
  filename = file.path(
    basis_plot_dir,
    sprintf(
      "TPS_small_scale_column_%02d.png",
      basis_k_max
    )
  )
)

# For the connected graph, column 1 is the constant zero-eigenvalue direction.
# Column 2 is the smoothest nonconstant graph-Laplacian direction.
plot_basis_points(
  values = B_gl[, 2L],
  basis_x = analysis_data$basis_x,
  basis_y = analysis_data$basis_y,
  title = "Graph Laplacian: large-scale basis function",
  subtitle = paste(
    "Column 2; roughness =",
    signif(roughness_gl[2L], 4L)
  ),
  filename = file.path(
    basis_plot_dir,
    "GL_large_scale_column_02.png"
  )
)

# The final column is the smallest-scale direction among the retained GL
# columns.
plot_basis_points(
  values = B_gl[, basis_k_max],
  basis_x = analysis_data$basis_x,
  basis_y = analysis_data$basis_y,
  title = "Graph Laplacian: smallest-scale retained basis function",
  subtitle = paste(
    "Column",
    basis_k_max,
    "; roughness =",
    signif(roughness_gl[basis_k_max], 4L)
  ),
  filename = file.path(
    basis_plot_dir,
    sprintf(
      "GL_small_scale_column_%02d.png",
      basis_k_max
    )
  )
)

message("Basis plots saved in: ", basis_plot_dir)

# -----------------------------------------------------------------------------
# Recreate the exposure / IV decomposition maps (g1-g5)
# -----------------------------------------------------------------------------
project_onto_basis <- function(B, n_cols, a) {
  Bc <- B[, seq_len(n_cols), drop = FALSE]
  Ac <- as.numeric(Bc %*% crossprod(Bc, a))
  list(Ac = Ac, Auc = a - Ac)
}
tps_decomp <- project_onto_basis(B_tps, 2L, a)
gl_decomp  <- project_onto_basis(B_gl, 6L, a)

analysis_sf$Ac_TPS             <- tps_decomp$Ac
analysis_sf$Auc_TPS            <- tps_decomp$Auc
analysis_sf$Ac_GraphLaplacian  <- gl_decomp$Ac
analysis_sf$Auc_GraphLaplacian <- gl_decomp$Auc

message(
  "TPS confounded share of exposure variance: ",
  signif(var(analysis_sf$Ac_TPS) / var(a), 3L)
)
message(
  "Graph-Laplacian confounded share of exposure variance: ",
  signif(var(analysis_sf$Ac_GraphLaplacian) / var(a), 3L)
)

# analysis_sf is in crs 5070 (meters) by this point in the script, but the
# xlim/ylim below are degree bounds -- reproject a plotting-only copy back to
# lon/lat so the maps aren't clipped to nothing.
plot_sf <- sf::st_transform(analysis_sf, crs = 4326)

ug_m3_label <- expression(mu * g / m^3)

map_theme <- theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(30, "points"),
        plot.title = element_text(hjust = 0.5))

g1 <- ggplot() +
  xlim(-125, -65) + ylim(25, 50) +
  geom_sf(data = plot_sf, aes(fill = pm25, color = pm25), size = 0.5) +
  scale_fill_viridis_c(name = ug_m3_label) +
  scale_color_viridis_c(guide = "none") +
  labs(title = expression(
    paste(PM[2.5], " exposure (A) averaged over 2001-2010 across US zip codes"
  ))) +
  map_theme

g2 <- ggplot() +
  xlim(-125, -65) + ylim(25, 50) +
  geom_sf(data = plot_sf, aes(fill = Auc_TPS, color = Auc_TPS), size = 0.5) +
  scale_fill_viridis_c(name = ug_m3_label) +
  scale_color_viridis_c(guide = "none") +
  labs(title = expression(paste("Estimate of ", A^{"uc"}, " using a thin plate spline basis"))) +
  map_theme

g3 <- ggplot() +
  xlim(-125, -65) + ylim(25, 50) +
  geom_sf(data = plot_sf, aes(fill = Auc_GraphLaplacian, color = Auc_GraphLaplacian), size = 0.5) +
  scale_fill_viridis_c(name = ug_m3_label) +
  scale_color_viridis_c(guide = "none") +
  labs(title = expression(paste("Estimate of ", A^{"uc"}, " using a graph Laplacian basis"))) +
  map_theme

g4 <- ggplot() +
  xlim(-125, -65) + ylim(25, 50) +
  geom_sf(data = plot_sf, aes(fill = Ac_TPS, color = Ac_TPS), size = 0.5) +
  scale_fill_viridis_c(name = ug_m3_label) +
  scale_color_viridis_c(guide = "none") +
  labs(title = expression(paste("Estimate of ", A^{"c"}, " using a thin plate spline basis"))) +
  map_theme

g5 <- ggplot() +
  xlim(-125, -65) + ylim(25, 50) +
  geom_sf(data = plot_sf, aes(fill = Ac_GraphLaplacian, color = Ac_GraphLaplacian), size = 0.5) +
  scale_fill_viridis_c(name = ug_m3_label) +
  scale_color_viridis_c(guide = "none") +
  labs(title = expression(paste("Estimate of ", A^{"c"}, " using a graph Laplacian basis"))) +
  map_theme

layout_matrix <- rbind(
  c(NA, 1, 1, NA),
  c(2, 2, 4, 4),
  c(3, 3, 5, 5)
)

png(
  filename = "images/PM_ivs_20pct_2010.png",
  width = 2500,
  height = 3000,
  res = 300
)
grid.arrange(g1, g2, g3, g4, g5, layout_matrix = layout_matrix)
dev.off()
gc()