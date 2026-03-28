library(sf)
library(mgcv)
library(dplyr)
library(utils)
library(gridExtra)
library(ggplot2)
library(geoR)
library(mvtnorm)
library(geosphere)
source('../funcs.R')

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
uscounties <- subset(uscounties, Region %in% 6)

# yields 503 counties
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

# Build ordered GL and TPS bases (largest scale in col 1, smallest in col n).
# These are passed to simfunc() and used for basis projection.
B_gl_full  <- build_gl_basis(adjmat)
B_tps_full <- build_tps_basis(lat, lon)

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

####################### Plot data on maps from 3 confounding mechanisms. #########################

set.seed(33)

# FIRST Confounding mechanism 
Sigma_GP <- compute_Sigma_GP(distmat = distmat,
                    rangeu = 0.01, 
                    rangec = 0.5)
dat <- compute_data_GP(n = 1, Sigma_GP = Sigma_GP)
Ac <- dat$Ac 
Auc <- dat$Auc
U <- dat$U
A <- Ac + Auc
uscounties$A <- A
uscounties$Ac <- Ac
uscounties$Auc <- Auc
uscounties$U <- U
gs1 <- plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c(
    'Exposure A',
    'Unconfounded Exposure Auc',
    'Confounded Exposure Ac',
    'Unmeasured Confounder U'
  ),
  xlimits = c(-115,-85),
  ylimits = c(25,40)
)

# SECOND confounding mechanism
Sigma_GP <- compute_Sigma_GP(distmat = distmat,
                            rangeu = 0.05, 
                            rangec = 0.5)
dat <- compute_data_GP(n = 1, Sigma_GP = Sigma_GP)
Ac <- dat$Ac 
Auc <- dat$Auc
U <- dat$U
A <- Ac + Auc
uscounties$A <- A
uscounties$Ac <- Ac
uscounties$Auc <- Auc
uscounties$U <- U
gs2 <- plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c(
    '',
    '',
    '',
    ''
  ),
  xlimits = c(-115,-85),
  ylimits = c(25,40)
)

# THIRD confounding mechanism
dat <- compute_data_GP_state(distmat = distmat,
                            rangeu = 0.01, 
                            rangec = 0.5,
                            n = 1,
                            statemat = statemat)
Ac <- dat$Ac 
Auc <- dat$Auc
U <- dat$U
A <- Ac + Auc
uscounties$A <- A
uscounties$Ac <- Ac
uscounties$Auc <- Auc
uscounties$U <- U
gs3 <- plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c(
    '',
    '',
    '',
    ''
  ),
  xlimits = c(-115,-85),
  ylimits = c(25,40)
)

# FOURTH confounding mechanism 
dat <- compute_data_spatialcoord(lat = simlist$lat, long = simlist$lon, nsims = 1)
Ac <- dat$Ac 
Auc <- dat$Auc
U <- dat$U
A <- Ac + Auc
uscounties$A <- A
uscounties$Ac <- Ac
uscounties$Auc <- Auc
uscounties$U <- U
gs4 <- plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c(
    '',
    '',
    '',
    ''
  ),
  xlimits = c(-115,-85),
  ylimits = c(25,40)
)

# FIFTH confounding mechanism
rangeu <- 0.01
Sigma_GP <- compute_Sigma_GP_2U(distmat = distmat,
                                kappa = 2,
                                rangeu = rangeu,
                                rangec = 0.5,
                                rangez1 = 0.5,
                                rangez2 = 0.3)
dat <- compute_data_GP_2U(n = 1, Sigma_GP = Sigma_GP)
Ac <- dat$Ac 
Auc <- dat$Auc
U1 <- dat$U1
U2 <- dat$U2
A <- Ac + Auc
uscounties$A <- A
uscounties$Ac <- Ac
uscounties$Auc <- Auc
uscounties$U1 <- U1
uscounties$U2 <- U2
gs5 <- plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U1', 'U2'),
  c(
    '',
    '',
    '',
    '',
    ''
  ),
  xlimits = c(-115,-85),
  ylimits = c(25,40)
)

# SIXTH confounding mechanism
rangeu <- 0.1
Sigma_GP <- compute_Sigma_GP(distmat = distmat,
                             rangeu = rangeu, 
                             rangec = 0.01)
dat <- compute_data_GP(n = 1, Sigma_GP = Sigma_GP)
Ac <- dat$Ac 
Auc <- dat$Auc
U <- dat$U
A <- Ac + Auc
uscounties$A <- A
uscounties$Ac <- Ac
uscounties$Auc <- Auc
uscounties$U <- U
gs6 <- plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c(
    '',
    '',
    '',
    ''
  ),
  xlimits = c(-115,-85),
  ylimits = c(25,40)
)

# SEVENTH confounding mechanism
W <- adjmat
dat <- compute_data_leroux(W = W, n = 1)
Ac <- dat$Ac 
Auc <- dat$Auc
U <- dat$U
A <- Ac + Auc
uscounties$A <- A
uscounties$Ac <- Ac
uscounties$Auc <- Auc
uscounties$U <- U
gs7 <- plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac', 'U'),
  c(
    '',
    '',
    '',
    ''
  ),
  xlimits = c(-115,-85),
  ylimits = c(25,40)
)

png(
  'images/all-three-decomp.jpeg',
  height = 5400,
  width = 5600,
  res = 140
)
ggpubr::ggarrange(gs1[[1]], gs1[[2]], gs1[[3]], gs1[[4]],  # First row (gs1)
                  gs7[[1]], gs7[[2]], gs7[[3]], gs7[[4]], # Seventh row (gs7)
                  gs3[[1]], gs3[[2]], gs3[[3]], gs3[[4]], # Third row (gs3)
                  gs5[[1]], gs5[[2]], gs5[[3]], ggpubr::ggarrange(gs5[[4]], gs5[[5]], ncol = 2), # Fifth row (gs5)
                  gs6[[1]], gs6[[2]], gs6[[3]], gs6[[4]], # Sixth row (gs6)
                  gs2[[1]], gs2[[2]], gs2[[3]], gs2[[4]], # Second row (gs2)
                  gs4[[1]], gs4[[2]], gs4[[3]], gs4[[4]], # Fourth row (gs4)
                  nrow = 7, ncol = 4)
dev.off()

