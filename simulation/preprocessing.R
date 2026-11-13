# Fixed inputs for the simulation, saved to sim.rds.
#
# Builds, once:
#   - the contiguous-US outline and the 48 contiguous states plus DC (Census
#     2018 cartographic boundaries), projected to CONUS Albers equal-area (km);
#   - 2e6 uniform reference points over the US, used to standardize every
#     surface to mean 0 and variance 1 and to compute the true psi by Monte
#     Carlo;
#   - the fixed spatial surfaces, each a sum of 30 cosine waves with
#     wavelengths in one band (funcs.R, wavelength_bands):
#       large    1200-2500 km  G in LS and stateJump, U1 in twoU   (seed 101)
#       large2   1200-2500 km  U2 in twoU
#       medium   800-1500 km   G in MS                             (seed 202)
#       x        500-1500 km   F, the spatial part of X
#   - the coordinate function of coordnoise (center and scale);
#   - for stateJump, one fixed offset b_r ~ N(0, 1) per state (seed 601), the
#     state of every reference point, and the center and scale of
#     G_LS + b_state over the reference points.
#
# Large-scale band: 1200-2500 km rather than longer wavelengths. A 2500-4500
# km surface is almost a low-order smooth function of the coordinates (the
# spatial coordinates method's Super Learner fits ~98% of its variance from
# (x, y) at n = 1000, vs ~70% at 1200-2500 km), so LS would test only the
# Gilbert et al. setting, which coordnoise already covers.
#
# Seeds: over the US, independent large-scale surfaces can correlate strongly
# by chance, so x and large2 each take the first seed, from a fixed start,
# whose |r| with every surface it appears with is below 0.1:
#   x: from 203, against large, medium, coord and a small-scale surface
#      (200-500 km, seed 301; kept only for this check so the seed choice is
#      unchanged from the run reported in the paper);
#   large2: from 102, against large and x.
#
# Usage (from simulation/): Rscript preprocessing.R
# Downloads the state shapefile; peak memory ~1.5 GB; a few minutes (state
# lookup for the 2e6 reference points).

suppressMessages({
  library(sf)
})
source("funcs.R")

set.seed(20261013)

# ---- Outline and states ----
temp_zip <- tempfile(fileext = ".zip")
download.file("https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_20m.zip",
              temp_zip, mode = "wb")
unzip(temp_zip, exdir = tempdir())
states_raw <- st_read(file.path(tempdir(), "cb_2018_us_state_20m.shp"), quiet = TRUE)
unlink(temp_zip)

# Drop Alaska, Hawaii, and Puerto Rico; project to km.
states_raw <- states_raw[!states_raw$STATEFP %in% c("02", "15", "72"), ]
outline <- st_make_valid(st_union(st_transform(states_raw, us_crs_km)))
outline <- st_sf(id = 1, geometry = outline)

# ---- Reference points ----
reference_points <- sample_us_points(outline, n_truth_points)
sim <- list(outline = outline, reference_points = reference_points)

# Coordinate function of coordnoise, standardized over the US.
coord_bbox <- unclass(st_bbox(outline))[c("xmin", "ymin", "xmax", "ymax")]
coord_raw <- coord_surface_raw(reference_points, coord_bbox)
sim$coord_u <- list(bbox = coord_bbox, center = mean(coord_raw), scale = sd(coord_raw))

# ---- Surfaces ----
check_points <- reference_points[1:2e5, ]
values_at <- function(surface) eval_surface(surface, check_points)

# First seed >= start whose surface has |r| < tol with every vector in others.
first_uncorrelated <- function(band, start, others, tol = 0.1) {
  for (seed in start + 0:99) {
    surface <- make_band_surface(band, reference_points, seed = seed)
    v <- values_at(surface)
    if (all(abs(vapply(others, function(o) cor(v, o), numeric(1))) < tol)) {
      message(sprintf("band %g-%g km: seed %d", band[1], band[2], seed))
      return(surface)
    }
  }
  stop("No seed found")
}

surfaces <- list(
  medium_conf = make_band_surface(wavelength_bands$medium, reference_points, seed = 202),
  large = make_band_surface(wavelength_bands$large, reference_points, seed = 101)
)
small <- make_band_surface(wavelength_bands$small, reference_points, seed = 301)
sim$surfaces <- surfaces
g_values <- list(large = values_at(surfaces$large), medium_conf = values_at(surfaces$medium_conf),
                 small = values_at(small), coord = coord_surface(sim, check_points))
surfaces$x <- first_uncorrelated(wavelength_bands$x, 203, g_values)
surfaces$large2 <- first_uncorrelated(wavelength_bands$large, 102,
                                      list(g_values$large, values_at(surfaces$x)))
sim$surfaces <- surfaces

# ---- States (mechanism stateJump) ----
states <- states_raw[, c("STATEFP", "NAME")]
states <- st_make_valid(st_transform(states, us_crs_km))
states <- states[order(states$STATEFP), ]
sim$states <- states
set.seed(601)
offsets <- setNames(rnorm(nrow(states)), states$NAME)
sim$jump <- list(offsets = unname(offsets), ref_state = NULL)
chunks <- split(seq_len(nrow(reference_points)), ceiling(seq_len(nrow(reference_points)) / 2e5))
ref_state <- unlist(lapply(chunks, function(i) state_of(sim, reference_points[i, , drop = FALSE])),
                    use.names = FALSE)
raw <- eval_surface(surfaces$large, reference_points) + sim$jump$offsets[ref_state]
sim$jump$ref_state <- ref_state
sim$jump$center <- mean(raw)
sim$jump$scale <- sd(raw)
message(sprintf("stateJump: %d states; offsets account for %.0f%% of Var(G_LS + b)",
                nrow(states), 100 * var(sim$jump$offsets[ref_state]) / var(raw)))

vals <- sapply(c(names(surfaces), "coord", "jump"), function(k) eval_named_surface(sim, k, check_points))
print(round(cor(vals), 2))

saveRDS(sim, "sim.rds")
