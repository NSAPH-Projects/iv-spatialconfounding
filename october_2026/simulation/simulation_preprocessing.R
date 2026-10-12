library(sf)
library(mgcv)
library(dplyr)
library(utils)
library(gridExtra)
library(ggplot2)
library(geoR)
library(mvtnorm)
library(geosphere)
library(Matrix)
source('funcs.R')

# Import geographies
temp_dir <- tempdir()
temp_zip <- tempfile(tmpdir = temp_dir, fileext = ".zip")
download.file(
  'https://www2.census.gov/geo/tiger/TIGER2010/COUNTY/2010/tl_2010_us_county10.zip',
  temp_zip,
  mode = "wb"
)
unzip(temp_zip, exdir = temp_dir)
uscounties <- st_read(file.path(temp_dir, 'tl_2010_us_county10.shp'))
unlink(temp_zip)

# Exclude Hawaii, Alaska, PR (disconnected from graph)
uscounties <- uscounties %>% filter(!STATEFP10 %in% c('02', '15', '72'))

# Add EPA regions
# source https://www.epa.gov/aboutepa/regional-and-geographic-offices
states <- c('CT', 'ME', 'MA', 'NH', 'RI', 'VT', 
            'NY', 'NJ', 'PR', 'VI', 
            'DE', 'DC', 'MD', 'PA', 'VA', 'WV', 
            'AL', 'FL', 'GA', 'KY', 'MS', 'NC', 'SC', 'TN', 
            'IL', 'IN', 'MI', 'MN', 'OH', 'WI', 
            'AR', 'LA', 'NM', 'OK', 'TX', 
            'IA', 'KS', 'MO', 'NE', 
            'CO', 'MT', 'ND', 'SD', 'UT', 'WY', 
            'AZ', 'CA', 'HI', 'NV', 
            'AK', 'ID', 'OR', 'WA')
statefps  <- c('09', '23', '25', '33', '44', '50', 
               '36', '34', '72', '78', 
               '10', '11', '24', '42', '51', '54', 
               '01', '12', '13', '21', '28', '37', '45', '47', 
               '17', '18', '26', '27', '39', '55', 
               '05', '22', '35', '40', '48', 
               '19', '20', '29', '31', 
               '08', '30', '38', '46', '56', '49', 
               '04', '06', '15', '32', 
               '02', '16', '41', '53')
regions <- c(rep(1,6),
             rep(2,4),
             rep(3,6),
             rep(4,8),
             rep(5,6),
             rep(6,5),
             rep(7,4),
             rep(8,6),
             rep(9,4),
             rep(10,4))
region_data <- data.frame(State = states, STATEFP = statefps, Region = regions)
uscounties <- left_join(uscounties, region_data, by = c('STATEFP10' = 'STATEFP'))

# Subset to just EPA region 6 for computation time
# this is nice because ~500 counties
uscounties <- subset(uscounties, Region %in% 6:7)

# yields 915 counties
n <- nrow(uscounties)

lat <- as.numeric(uscounties$INTPTLAT10)
lon <- as.numeric(uscounties$INTPTLON10)
distmat <- geosphere::distm(cbind(lon, lat), 
                            fun = distHaversine)
# scale distance matrix so range is between (0,2)
distmat <- distmat/1000000

# Extract adjacency matrix
adjacency_list <- st_touches(uscounties)
adjmat <- sapply(adjacency_list, function(row) {
  binary_row = rep(0, nrow(uscounties))
  binary_row[row] = 1
  binary_row
})
adjmat <- (adjmat + t(adjmat)) > 0

# Graph Laplacian
L <- diag(rowSums(adjmat)) - adjmat

# sim.RData holds the bases used by every simulation run. It is rebuilt only
# when rebuild_sim_data is TRUE (or the file is missing); otherwise the existing
# file is loaded, so this script can be rerun just to redraw the figure below
# without changing the inputs of runs already completed.
rebuild_sim_data <- FALSE

if (rebuild_sim_data || !file.exists('sim.RData')) {
  # Build ordered GL and TPS bases (largest scale in col 1, smallest in col n).
  # These are passed to simfunc() and used for basis projection.
  source('../R/basis_construction.R')
  B_gl_full  <- build_gl_basis_truncated(L, k_max = nrow(adjmat)-1)
  B_tps_full <- build_tps_basis_truncated(lat, lon, k_max = nrow(adjmat)-1)

  # Indicator matrix for states
  statemat <- model.matrix(~-1 + State, data = uscounties)

  # Save simulation data
  simlist <- list(
    'lat'        = lat,
    'lon'        = lon,
    'B_tps_full' = B_tps_full,
    'B_gl_full'  = B_gl_full,
    'statemat'   = statemat,
    'W'          = adjmat
  )
  save(simlist,
       file = 'sim.RData')
} else {
  load('sim.RData')
  # The shapefile must reproduce the counties (and their order) in sim.RData.
  stopifnot(
    length(simlist$lat) == n,
    isTRUE(all.equal(simlist$lat, lat)),
    isTRUE(all.equal(simlist$lon, lon)),
    all(simlist$W == adjmat)
  )
  B_gl_full  <- simlist$B_gl_full
  B_tps_full <- simlist$B_tps_full
}

####################### Plot data on maps for each confounding mechanism. #########################
# Mechanisms are plotted in the same order as the confounding_mechanism switch in funcs.R (simfunc()).

set.seed(33)

# Helper: extract vector column from nsims=1 output
col1 <- function(x) x[, 1]

# Helper: drop the per-panel color-scale legend (32 panels with individual
# legends is too cluttered; the diverging scale already conveys sign/shape).
strip_legend <- function(gs) lapply(gs, function(g) g + theme(legend.position = "none"))

# MECHANISM 1: TPS-basis confounding (IV-TPS recovers Ac exactly)
dat <- compute_data_TPS_basis(B_tps_full = B_tps_full, nsims = 1)
uscounties$Auc <- col1(dat$Auc)
uscounties$Ac  <- col1(dat$Ac)
uscounties$U   <- col1(dat$U)
uscounties$A   <- uscounties$Ac + uscounties$Auc
gs1 <- strip_legend(plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c(
    expression(paste("Exposure ", A)),
    expression(paste("Instrumental variation ", A^{uc})),
    expression(paste("Non-instrumental variation ", A^c)),
    expression(paste("Confounder ", U))
  ),
  xlimits = c(-115,-85), ylimits = c(25,40)
))

# MECHANISM 2: GL-basis confounding (IV-GraphLaplacian recovers Ac exactly)
dat <- compute_data_GL_basis(B_gl_full = B_gl_full, nsims = 1)
uscounties$Auc <- col1(dat$Auc)
uscounties$Ac  <- col1(dat$Ac)
uscounties$U   <- col1(dat$U)
uscounties$A   <- uscounties$Ac + uscounties$Auc
gs2 <- strip_legend(plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c('', '', '', ''),
  xlimits = c(-115,-85), ylimits = c(25,40)
))

# MECHANISM 3: Two-confounder TPS-basis
dat <- compute_data_TPS_2U(B_tps_full = B_tps_full, nsims = 1)
uscounties$Auc <- col1(dat$Auc)
uscounties$Ac  <- col1(dat$Ac)
uscounties$U1  <- col1(dat$U1)
uscounties$U2  <- col1(dat$U2)
uscounties$A   <- uscounties$Ac + uscounties$Auc
gs3 <- strip_legend(plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U1', 'U2'),
  c('', '', '', '', ''),
  xlimits = c(-115,-85), ylimits = c(25,40)
))

# MECHANISM 4: Reversed-scale TPS (confounding in high-frequency TPS components)
dat <- compute_data_TPS_reversed(B_tps_full = B_tps_full, nsims = 1)
uscounties$Auc <- col1(dat$Auc)
uscounties$Ac  <- col1(dat$Ac)
uscounties$U   <- col1(dat$U)
uscounties$A   <- uscounties$Ac + uscounties$Auc
gs4 <- strip_legend(plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c('', '', '', ''),
  xlimits = c(-115,-85), ylimits = c(25,40)
))

# MECHANISM 5: Bivariate Leroux CAR
dat <- compute_data_leroux(W = adjmat, nsims = 1)
uscounties$Auc <- col1(dat$Auc)
uscounties$Ac  <- col1(dat$Ac)
uscounties$U   <- col1(dat$U)
uscounties$A   <- uscounties$Ac + uscounties$Auc
gs5 <- strip_legend(plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c('', '', '', ''),
  xlimits = c(-115,-85), ylimits = c(25,40)
))

# MECHANISM 6: Coordinate-based nonlinear
dat <- compute_data_spatialcoord(lat = simlist$lat, long = simlist$lon,
                                 nsims = 1, distmat = distmat)
uscounties$Auc <- col1(dat$Auc)
uscounties$Ac  <- col1(dat$Ac)
uscounties$U   <- col1(dat$U)
uscounties$A   <- uscounties$Ac + uscounties$Auc
gs6 <- strip_legend(plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c('', '', '', ''),
  xlimits = c(-115,-85), ylimits = c(25,40)
))

# MECHANISM 7: Two-confounder GL-basis
dat <- compute_data_GL_2U(B_gl_full = B_gl_full, nsims = 1)
uscounties$Auc <- col1(dat$Auc)
uscounties$Ac  <- col1(dat$Ac)
uscounties$U1  <- col1(dat$U1)
uscounties$U2  <- col1(dat$U2)
uscounties$A   <- uscounties$Ac + uscounties$Auc
gs7 <- strip_legend(plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U1', 'U2'),
  c('', '', '', '', ''),
  xlimits = c(-115,-85), ylimits = c(25,40)
))

# MECHANISM 8: Reversed-scale GL (confounding in high-frequency GL components)
dat <- compute_data_GL_reversed(B_gl_full = B_gl_full, nsims = 1)
uscounties$Auc <- col1(dat$Auc)
uscounties$Ac  <- col1(dat$Ac)
uscounties$U   <- col1(dat$U)
uscounties$A   <- uscounties$Ac + uscounties$Auc
gs8 <- strip_legend(plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c('', '', '', ''),
  xlimits = c(-115,-85), ylimits = c(25,40)
))

row_label <- function(label) {
  grid::textGrob(
    label = label,
    x = grid::unit(1, "npc"),
    y = grid::unit(0.5, "npc"),
    just = c("right", "center"),
    gp = grid::gpar(
      fontsize = 40, 
      fontface = "bold"
    )
  )
}

dir.create("images", showWarnings = FALSE)
png(
  "images/all-decomp.jpeg",
  height = 6000,
  width = 5500,
  res = 140
)

ggpubr::ggarrange(
  row_label("Confounding\n mechanism 1"), gs1[[1]], gs1[[2]], gs1[[3]],
  ggpubr::ggarrange(gs1[[4]], ncol = 1),
  
  row_label("Confounding\n mechanism 2"), gs2[[1]], gs2[[2]], gs2[[3]],
  ggpubr::ggarrange(gs2[[4]], ncol = 1),
  
  row_label("Confounding\n mechanism 3"), gs3[[1]], gs3[[2]], gs3[[3]],
  ggpubr::ggarrange(gs3[[4]], gs3[[5]], ncol = 2),
  
  row_label("Confounding\n mechanism 4"), gs4[[1]], gs4[[2]], gs4[[3]],
  ggpubr::ggarrange(gs4[[4]], ncol = 1),
  
  row_label("Confounding\n mechanism 5"), gs5[[1]], gs5[[2]], gs5[[3]],
  ggpubr::ggarrange(gs5[[4]], ncol = 1),
  
  row_label("Confounding\n mechanism 6"), gs6[[1]], gs6[[2]], gs6[[3]],
  ggpubr::ggarrange(gs6[[4]], ncol = 1),
  
  row_label("Confounding\n mechanism 7"), gs7[[1]], gs7[[2]], gs7[[3]],
  ggpubr::ggarrange(gs7[[4]], gs7[[5]], ncol = 2),
  
  row_label("Confounding\n mechanism 8"), gs8[[1]], gs8[[2]], gs8[[3]],
  ggpubr::ggarrange(gs8[[4]], ncol = 1),
  
  nrow = 8,
  ncol = 5,
  widths = c(0.40, 1, 1, 1, 1)
)

dev.off()
