library(sf)
library(mgcv)
library(dplyr)
library(utils)
library(gridExtra)
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

latnorm = uscounties$latnorm
longnorm = uscounties$longnorm

# Extract adjacency matrix
adjacency_list = st_touches(uscounties)
adjmat = sapply(adjacency_list, function(row) {
  binary_row = rep(0, nrow(uscounties))
  binary_row[row] = 1
  binary_row
})
adjmat = (adjmat + t(adjmat)) > 0


# Create projmats.RData
# a list of projection matrices of confounding

# # Thin plate spline
# yfake = rnorm(3109)
# gam_model = mgcv::gam(yfake ~ s(latnorm, longnorm, bs="tp", k=11), method="REML")
# b = predict(gam_model, type="lpmatrix")
# b = b[,-1] # remove intercept
# b15 = b[,1:5]
# b610 = b[,6:10]
# Vtps15 = b15 %*% solve(t(b15) %*% b15) %*% t(b15)
# Vtps610 = b610 %*% solve(t(b610) %*% b610) %*% t(b610)

# Graph Fourier
L = diag(rowSums(adjmat)) - adjmat
E = eigen(L)
Vgft150 = E$vectors[,(n-50):(n-1)]
Vgft51100 = E$vectors[,(n-100):(n-51)]
GFT150 = Vgft150 %*% t(Vgft150)
GFT51100 = Vgft51100 %*% t(Vgft51100)

# Nested s
groups = uscounties$STATEFP10
nest = nested_decomp_mats(groups)
Vstate = nest$decomp_mats[[1]]

projmats = list(
  # 'TPS15' = Vtps15,
  # 'TPS610' = Vtps610,
  'GFT150' = GFT150,
  'GFT51100' = GFT51100,
  'nestedstate' = Vstate
)
save(projmats,
     file = 'projmats.RData')

# Create sim.RData
# a list latnorm (vector), longnorm (vector), adjmat (matrix), and a.vals (list of vectors)
U = createU(latnorm, longnorm)
X = createX(U, latnorm, longnorm)
a.vals = list()
a.vals[['GFT150']] = seq(-10,15,length.out = 100)
a.vals[['GFT51100']] = seq(-10,15,length.out = 100)
# a.vals[['TPS15']] = seq(-10,15,length.out = 100)
# a.vals[['TPS610']] = seq(-10,15,length.out = 100)
a.vals[['nestedstate']] = seq(-10,25,length.out = 100)

simlist = list(
  'a.vals' = a.vals,
  'latnorm' = latnorm,
  'longnorm' = longnorm,
  'adjmat' = adjmat
)
save(simlist,
     file = 'sim.RData')

## Create example plot of data
uscounties$U = createU(latnorm, longnorm)
uscounties$X = createX(uscounties$U, latnorm, longnorm)
Adat = createA(uscounties$U, uscounties$X, projmat = projmats$GFT150)
uscounties$A = Adat$A
uscounties$Auc = Adat$Auc
uscounties$Ac = Adat$Ac
uscounties$Y = createY(uscounties$U, uscounties$X, uscounties$A, 'nonlinear')

gs = plotfunc(
  uscounties,
  c('U', 'X', 'A', 'Y', 'Auc', 'Ac'),
  c(
    'Unmeasured Confounder U',
    'Measured Confounder X',
    'Exposure A',
    'Outcome Y',
    'Unconfounded Exposure Auc',
    'Confounded Exposure Ac'
  )
)
png('images/exampledata_gft150.jpeg', height = 5000, width = 4000, res = 200)
grid.arrange(grobs = gs, 
             layout = matrix(1:6, nrow = 3, ncol = 2))
dev.off()

# Plot the first 20 eigenvectors on the usmap
n = ncol(E$vectors)
for (i in 101:130) {
  uscounties[[paste0('eig', i)]] = E$vectors[,n-i+1]
}
gs = plotfunc(
  uscounties,
  paste0('eig', 101:130),
  paste0('Eigenvector ', 101:130)
)
png('images/eigenvectors100130.jpeg', height = 5000, width = 10000, res = 200)
grid.arrange(grobs = gs, 
             layout = matrix(1:30, nrow = 5, ncol = 6))
dev.off()