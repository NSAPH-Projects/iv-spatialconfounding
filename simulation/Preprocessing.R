library(sf)
library(mgcv)
library(dplyr)
library(utils)
source('simfuncs.R')

# Import geographies
temp_dir = tempdir()
temp_zip = tempfile(tmpdir = temp_dir, fileext = ".zip")
download.file(
  'https://www2.census.gov/geo/tiger/TIGER2010/COUNTY/2010/tl_2010_us_county10.zip',
  temp_zip,
  mode = "wb"
)
unzip(temp_zip, exdir = temp_dir)
uscounties = st_read(file.path(temp_dir, 'tl_2010_us_county10.shp'))
unlink(temp_zip)

# Exclude Hawaii, Alaska, PR (disconnected from graph)
uscounties = uscounties %>% filter(!STATEFP10 %in% c('02', '15', '72'))
# yields 3109 counties
n = nrow(uscounties)

# Extract normalized lat + long
uscounties$longnorm = scale(as.numeric(uscounties$INTPTLON10))
uscounties$latnorm = scale(as.numeric(uscounties$INTPTLAT10))

# Extract adjacency matrix
adjacency_list = st_touches(uscounties)
adjmat = sapply(adjacency_list, function(row) {
  binary_row = rep(0, nrow(uscounties))
  binary_row[row] = 1
  binary_row
})
adjmat = (adjmat + t(adjmat)) > 0

# Create sim.RData
# a list of a.vals (vector), latnorm (vector), longnorm (vector), adjmat (matrix)
a.vals = seq(-3.5,4.5,length.out = 100)
latnorm = uscounties$latnorm
longnorm = uscounties$longnorm

simlist = list(
  'a.vals' = a.vals,
  'latnorm' = latnorm,
  'longnorm' = longnorm,
  'adjmat' = adjmat
)
save(simlist,
     file = 'simulation/sim.RData')

# Create projmats.RData
# a list of projection matrices of confounding

# Thin plate spline
yfake = rnorm(3109)
gam_model = mgcv::gam(yfake ~ s(latnorm, longnorm, bs="tp", k=11), method="REML")
b = predict(gam_model, type="lpmatrix")
b = b[,-1] # remove intercept
b15 = b[,1:5]
b610 = b[,6:10]
Vtps15 = b15 %*% solve(t(b15) %*% b15) %*% t(b15)
Vtps610 = b610 %*% solve(t(b610) %*% b610) %*% t(b610)

# Graph Fourier
L = diag(rowSums(adjmat)) - adjmat
E = eigen(L)
Vgft15 = E$vectors[,(n-5):(n-1)]
Vgft610 = E$vectors[,(n-10):(n-6)]

# Nested 
groups = uscounties$STATEFP10
nest = nested_decomp_mats(groups)
Vstate = nest$decomp_mats[[1]]

projmats = list(
  'TPS15' = Vtps15,
  'TPS610' = Vtps610,
  'GFT15' = Vgft15,
  'GFT610' = Vgft610,
  'nestedstate' = Vstate
)
save(projmats,
     file = 'simulation/projmats.RData')