library(sf)
library(mgcv)
library(dplyr)
library(utils)
library(gridExtra)
library(ggplot2)
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

# Add EPA regions
# source https://www.epa.gov/aboutepa/regional-and-geographic-offices
states = c('CT', 'ME', 'MA', 'NH', 'RI', 'VT', 
           'NY', 'NJ', 'PR', 'VI', 
           'DE', 'DC', 'MD', 'PA', 'VA', 'WV', 
           'AL', 'FL', 'GA', 'KY', 'MS', 'NC', 'SC', 'TN', 
           'IL', 'IN', 'MI', 'MN', 'OH', 'WI', 
           'AR', 'LA', 'NM', 'OK', 'TX', 
           'IA', 'KS', 'MO', 'NE', 
           'CO', 'MT', 'ND', 'SD', 'UT', 'WY', 
           'AZ', 'CA', 'HI', 'NV', 
           'AK', 'ID', 'OR', 'WA')
statefps  = c('09', '23', '25', '33', '44', '50', 
              '36', '34', '72', '78', 
              '10', '11', '24', '42', '51', '54', 
              '01', '12', '13', '21', '28', '37', '45', '47', 
              '17', '18', '26', '27', '39', '55', 
              '05', '22', '35', '40', '48', 
              '19', '20', '29', '31', 
              '08', '30', '38', '46', '56', '49', 
              '04', '06', '15', '32', 
              '02', '16', '41', '53')
regions = c(rep(1,6),
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
uscounties = left_join(uscounties, region_data, by = c('STATEFP10' = 'STATEFP'))

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

# Graph Fourier
L = diag(rowSums(adjmat)) - adjmat
E = eigen(L)
Vcgft150 = E$vectors[,(n-50):(n-1)]
Vucgft150 = Vuc = E$vectors[,1:(n-51)]
Vcgft51100 = E$vectors[,(n-100):(n-51)]
Vucgft51100 = cbind(E$vectors[,1:(n-101)], E$vectors[, (n-50):(n-1)])

# Nested decomposition
groups = cbind(uscounties$Region, uscounties$STATEFP10)
#nest = nested_decomp_mats(groups)
nest = nested_decomp(groups)
nestedbasis = nest$basis
Vregion = nest$basis[,1:nest$ranks[1]]
Vstate = nest$basis[,(nest$ranks[1]+1):(nest$ranks[1] + nest$ranks[2])]
Vcounty = nest$basis[,(nest$ranks[1] + nest$ranks[2] + 1):ncol(nest$basis)]
#Vregion = nest$decomp_mats[[1]]
#Vstate = nest$decomp_mats[[2]]

Vcs = list(
  'GFT150' = Vcgft150,
  'GFT51100' = Vcgft51100,
  'nestedstate' = Vstate,
  'nestedregion' = Vregion
)
save(Vcs, file = 'Vcs.RData')

Vucs = list(
  'GFT150' = Vucgft150,
  'GFT51100' = Vucgft51100,
  'nestedstate' = cbind(Vcounty, Vregion),
  'nestedregion' = cbind(Vcounty, Vstate)
)
save(Vucs, file = 'Vucs.RData')

# Create sim.RData
# a list latnorm (vector), longnorm (vector), adjmat (matrix), and a.vals (list of vectors)
# Ranges from a single creation of A, 5 and 95 percent quantiles
set.seed(33)
U = createU(latnorm, longnorm)
X = createX(U, latnorm, longnorm)
a.vals = list()
for (projmatname in names(Vcs)){
  A = createA_fromV(U, X, Vc = Vcs[[projmatname]], Vuc = Vucs[[projmatname]])$A
  a.vals[[projmatname]] = seq(quantile(A, probs = 0.05), 
                              quantile(A, probs = 0.95),
                              length.out = 100)
}

simlist = list(
  'a.vals' = a.vals,
  'latnorm' = latnorm,
  'longnorm' = longnorm,
  'adjmat' = adjma
)
save(simlist,
     file = 'sim.RData')

### USED FOR SENSITIVITY ####
# Create adjustmentmats.RData
adjustmentmats = list()
Esub = E$vectors[,(n-200):(n-1)]
# reverse column order of Esub so it's in increasing order of eigenvalues
Esub = Esub[,ncol(Esub):1]
# bin the columns of Esub into 10 bins
binsizes = c(20,50,100)
for (binsize in binsizes){
  for (j in 1:(200/binsize)){
    print(c(((j-1)*binsize+1),(binsize*j)))
    binmat = Esub[,((j-1)*binsize+1):(binsize*j)]
    adjustmentmats[[paste0('GFT',((j-1)*binsize+1),(binsize*j))]] = binmat %*% t(binmat)
  }
}
adjustmentmats[['nestedstate']] = Vstate
adjustmentmats[['nestedregion']] = Vregion
save(adjustmentmats, file = 'adjustmentmats.RData')

# Create adjustmentmats2.RData 
adjustmentmats = list()
Esub = E$vectors[,(n-1):(n-200)] 
# reverse column order of Esub so it's in increasing order of eigenvalues
#Esub = Esub[,ncol(Esub):1]
# bin the columns of Esub into 10 bins
for (j in 1:10){
  print(c(1,(10*j)))
  binmat = Esub[,1:(10*j)]
  adjustmentmats[[paste0('GFT',1,(10*j))]] = binmat %*% t(binmat)
}
save(adjustmentmats, file = 'adjustmentmats2.RData')

################# Create example plots of data ####################

projmatname = 'nestedregion'
uscounties$U = createU(latnorm, longnorm)
uscounties$X = createX(uscounties$U, latnorm, longnorm)
Adat = createA_fromV(uscounties$U, uscounties$X, Vc = Vcs[[projmatname]], Vuc = Vucs[[projmatname]])
uscounties$A = Adat$A
uscounties$Auc = Adat$Auc
uscounties$Ac = Adat$Ac
uscounties$Y = createY(uscounties$U, uscounties$X, uscounties$A, 'nonlinear')

gs = plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac'),
  c(
    'Exposure A',
    'Unconfounded Exposure Auc',
    'Confounded Exposure Ac'
  )
)
png(
  paste0('images/exampledata_', projmatname, '_smoothed.jpeg'),
  height = 1000,
  width = 4000,
  res = 100
)
# plot using grid.arrange in a single row of three columns
grid.arrange(grobs = gs, nrow = 1)
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


##### May 5: smoother exposure ######
#Vc = E$vectors[,(n-100):(n-51)]
#Vuc = cbind(E$vectors[,1:(n-101)], E$vectors[,(n-50):(n-1)])
Vc = E$vectors[,(n-50):(n-1)]
Vuc = E$vectors[,1:(n-51)]
save(Vc, file = 'Vc.RData')
save(Vuc, file = 'Vuc.RData')
uscounties$U = createU(latnorm, longnorm)
uscounties$X = createX(uscounties$U, latnorm, longnorm)
Adat = createA_fromV(
  uscounties$U,
  uscounties$X,
  Vc = Vc,
  Vuc = Vuc,
  trunc = 1500
)
# eliminate 1000 high freq eigenvec variation ^^
uscounties$A = Adat$A
uscounties$Auc = Adat$Auc
uscounties$Ac = Adat$Ac

gs = plotfunc(
  uscounties,
  c('A', 'Auc', 'Ac'),
  c(
    'Exposure A',
    'Unconfounded Exposure Auc',
    'Confounded Exposure Ac'
  )
)
png(
  'images/exampledata_GFT51100_smoothed.jpeg',
  height = 1000,
  width = 4000,
  res = 100
)
# plot using grid.arrange in a single row of three columns
grid.arrange(grobs = gs, nrow = 1)
dev.off()