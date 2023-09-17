library(spData)
library(sf)
library(spdep)
library(dplyr)
library(ggpubr)
library(gridExtra)

spscale = function(x, nbs, E){ # make sure order of $X$ matches A
  normx = x/sqrt(sum(x^2)) # normalize via Euclidean norm
  totaldist = 0
  for (node in 1:length(x)){
    totaldist = totaldist + sum((normx[node] - normx[nbs[[node]]])^2, na.rm = T)
  }
  return(totaldist/(2*E))
}

spscale2 = function(v, L, E){
  normv = v/sqrt(sum(v^2))
  return(sum(normv * (L %*% normv))/E)
}

boston.tr <- sf::st_read(system.file("shapes/boston_tracts.shp",
                                     package="spData")[1])

g1 = ggplot(boston.tr) +
  geom_sf(aes(fill = DIS), col = NA, size = 0.005) + 
  theme_minimal() + 
  scale_fill_viridis_c() + 
  theme(plot.title = element_text(size = 24 * 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))
g1

boston.tr$logCRIM = log(boston.tr$CRIM)
g2 = ggplot(boston.tr) +
  geom_sf(aes(fill = logCRIM), col = NA, size = 0.005) + 
  theme_minimal() + 
  scale_fill_viridis_c() + 
  theme(plot.title = element_text(size = 24 * 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))
g2

# Arrange g1 and g2 in a single plot and post to overleaf. 
png('images/casestudya.jpeg', height = 500, width = 1000)
grid.arrange(g1, g2, ncol = 2)
dev.off()

nri = read.csv('National_Risk_Index_Census_Tracts.csv')
#nrow(nri) # 85154 census tracts
cal = nri[nri$STATE == 'California',] # 9106
caltr = st_read('tl_2021_06_tract/tl_2021_06_tract.shp')
#nrow(caltr) # 9129
#mean(cal$TRACTFIPS %in% as.numeric(caltr$GEOID)) # all covered
caltr$TRACTFIPS = as.numeric(caltr$GEOID)

calmerged = merge(caltr, cal, by = 'TRACTFIPS')
#nrow(calmerged) # 9106

g1 = ggplot(calmerged) +
  geom_sf(aes(fill = SOVI_SCORE), col = NA, size = 0.005) + 
  theme_minimal() + 
  scale_fill_viridis_c() + 
  theme(plot.title = element_text(size = 24 * 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))

g2 = ggplot(calmerged) +
  geom_sf(aes(fill = HWAV_AFREQ), col = NA, size = 0.005) + 
  theme_minimal() + 
  scale_fill_viridis_c() + 
  theme(plot.title = element_text(size = 24 * 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))

png('images/casestudyb.jpeg', height = 500, width = 1000)
grid.arrange(g1, g2, ncol = 2)
dev.off()

#### Calculate spatial scale of all variables ####
nbsbos = poly2nb(boston.tr)
1/spscale(boston.tr$DIS, nbs = nbsbos)
1/spscale(boston.tr$logCRIM, nbs = nbsbos)

nbscal = poly2nb(calmerged)
1/spscale(calmerged$SOVI_SCORE, nbs = nbscal)
1/spscale(calmerged$HWAV_AFREQ, nbs = nbscal)

## Decompose using Graph Laplacian
source('funcs.R')
# Create adjacency matrices
bosadj = nb2mat(neighbours = nbsbos, style = 'B', zero.policy = T)
bosspec = spectral_decomp(bosadj)

caladj = nb2mat(neighbours = nbscal, style = 'B', zero.policy = T)
calspec = spectral_decomp(caladj)




# Boston map
bosmap = cbind(boston.tr, t(bosspec))
for (j in 38:543){
  names(bosmap)[j] = paste('eig', j-37, sep = '')
}
gs = list()
idxs = c(50, 100, 150, 200, 250, 300, 350, 400, 450, 475, 500, 505)
eigens = names(bosmap)[37 + idxs]

# compute spatial scale for eigenvector
D = rep(NA, length(idxs))
for (i in 1:length(idxs)){
  id = idxs[i]
  D[i] = spscale(bosspec[id,], nbs = nbsbos)
}

for (i in 1:length(eigens)){
  gs[[i]] = ggplot(bosmap) +
    geom_sf(aes_string(fill = eigens[i]), color = NA) +
    scale_fill_viridis_c() + 
    labs(title = paste('spscale =', round(D[i],3))) + 
    theme_minimal() +
    theme(plot.title = element_text(size = 24 * 2,hjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          line = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_line(colour = "transparent"))
}

png('images/bosmap.jpeg', height = 1024 * 0.6 * 2, width = 1024 * 2)
do.call(grid.arrange,gs)
dev.off()

# Plot of eigenvalue vs spatial scale for boston. 
L = diag(rowSums(bosadj)) - bosadj # graph laplacian 
E = eigen(L) # eigen component
D = E$val
G = E$vec
edges = sum(bosadj)
spscales = rep(NA, length(D))
spscales2 = rep(NA, length(D))
for (i in 1:length(spscales)){
  spscales[i] = spscale(G[,i], nbs = nbsbos, E = edges)
  spscales2[i] = spscale2(G[,i], L = L, E = edges)
}

plot(D, spscales, type = 'l', xlab = 'eigenvalue', ylab = 'spscale')
lines(D, spscales2, type = 'l', col = 'red')

## What if subspace is spanned by multiple vectors
sims = 30
spscales2 = rep(NA, length(D)/11)
eigenavg = rep(NA, length(D)/11)
for (i in 1:length(spscales)){
  idxs = (11*(i-1)+1):(11*i)
  eigenavg[i] = mean(D[idxs])
  spscales_sims = rep(NA, nsims)
  for (sim in 1:nsims){
    c1s = runif(10, min = 0, max = 0.4)
    c2 = sqrt(1-sum(c1s^2))
    coeffs = c(c1s, c2)
    vec = rep(0, nrow(G))
    for (j in 1:length(coeffs)){
      vec = vec + G[,idxs[j]]*coeffs[j]
    }
    spscales_sims[sim] = spscale2(vec, L = L, E = edges)
  }
  spscales2[i] = max(spscales_sims)
}

#plot(D, spscales, type = 'l', xlab = 'eigenvalue', ylab = 'spscale')
plot(eigenavg, spscales2*edges, col = 'red', pch = 19)
abline(a = 0, b = 1, col = 'black')


# California map
calmap = cbind(calmerged, t(calspec))
for (j in 483:9588){
  names(calmap)[j] = paste('eig', j-482, sep = '')
}
gs = list()
idxs = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 9050, 9100, 9104)
eigens = names(calmap)[482 + idxs]

# compute spatial scale for eigenvector
D = rep(NA, length(idxs))
for (i in 1:length(idxs)){
  id = idxs[i]
  D[i] = spscale(calspec[id,], nbs = nbscal)
}

for (i in 1:length(eigens)){
  gs[[i]] = ggplot(calmap) +
    geom_sf(aes_string(fill = eigens[i]), color = NA) +
    scale_fill_viridis_c() + 
    labs(title = paste('spscale =', round(D[i],4))) + 
    theme_minimal() +
    theme(plot.title = element_text(size = 24 * 2,hjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          line = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_line(colour = "transparent"))
}

png('images/calmap.jpeg', height = 1024 * 0.6 * 2, width = 1024 * 2)
do.call(grid.arrange,gs)
dev.off()

# Plot of eigenvalue vs spatial scale for california. 
L = diag(rowSums(caladj)) - caladj # graph laplacian 
E = eigen(L) # eigen component
D = E$val
G = E$vec
spscales = rep(NA, length(D))
for (i in 1:length(spscales)){
  spscales[i] = spscale(G[,i], nbs = nbscal)
}

plot(D, spscales, type = 'l', xlab = 'eigenvalue', ylab = 'spscale')





## California vulnerability using nested decomposition
# https://wifire-data.sdsc.edu/dataset/counties-in-california/resource/248bb029-e66b-4cc1-b551-ec4b3642ea3f?inner_span=True
regcal = read.csv('California_County_Boundaries.csv')
calmerged$County_FIPS_ID = as.numeric(calmerged$COUNTYFP)
regcal = regcal[,c(2,3,6,7)]
cal_reg = merge(x=calmerged,y=regcal, 
             by="County_FIPS_ID", all.x=TRUE)
groups = cbind(cal_reg$AdminRegion, cal_reg$CountyName)
nest = nested_decomp_mats(groups)
nestedvuln = matrix(data = NA, nrow = nrow(cal_reg), ncol = 3)
for (i in 1:3){
  nestedvuln[,i] = nest$decomp_mats[[i]] %*% cal_reg$SOVI_SCORE
}
labels = c('Region', 'County', 'Tract')
colnames(nestedvuln) = labels
cal_reg = cbind(cal_reg, nestedvuln)

gs = list()
for (i in 1:3){
  gs[[i]] = ggplot(cal_reg) +
    geom_sf(aes_string(fill = labels[i]), col = NA) + 
    theme_minimal() + 
    scale_fill_viridis_c() + 
    labs(title = labels[i]) + 
    theme(plot.title = element_text(size = 24 * 2, hjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          line = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_line(colour = "transparent"))
}
gs[[4]] = ggplot(cal_reg) +
  geom_sf(aes(fill = SOVI_SCORE), col = NA) + 
  theme_minimal() + 
  scale_fill_viridis_c() +   
  labs(title = 'Total') + 
  theme(plot.title = element_text(size = 24 * 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))

png('images/nested_vuln_cal.jpeg', height = 1024 * 0.6, width = 1024 * 2)
grid.arrange(grobs = gs, ncol = 4)
dev.off()




#### California vulnerability using spectral decomposition
Xstar = calspec %*% scale(calmerged$SOVI_SCORE)
idxs = c(9080, 9094, 9100)
eigenvecs = t(calspec)[,idxs]
for (col in 1:ncol(eigenvecs)){
  eigenvecs[,col] = Xstar[idxs[col]]*eigenvecs[,col]
}
calmap = cbind(calmerged, eigenvecs)
for (j in 483:485){
  names(calmap)[j] = paste('eig', idxs[j-482], sep = '')
}
gs = list()
eigens = names(calmap)[483:485]

for (i in 1:length(eigens)){
  gs[[i]] = ggplot(calmap) +
    geom_sf(aes_string(fill = eigens[i]), color = NA) +
    scale_fill_viridis_c() + 
    ggtitle(paste("\u03BB = ", round(D[idxs[i]],4))) + 
    theme_minimal() +
    theme(plot.title = element_text(size = 24 * 2, hjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          line = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_line(colour = "transparent"))
}
gs[[4]] = ggplot(calmerged) +
  geom_sf(aes(fill = SOVI_SCORE), col = NA) + 
  theme_minimal() + 
  scale_fill_viridis_c() +   
  labs(title = 'Total') + 
  theme(plot.title = element_text(size = 24 * 2, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))

png('images/spectral_vuln_cal.jpeg', height = 1024*0.6, width = 1024 * 2)
grid.arrange(grobs = gs, ncol = 4)
dev.off()

edges = sum(caladj)

set.seed(23)
nsims = 30
spscales_nested = rep(NA, 3)
bdedges = rep(NA, 3)
for (i in 1:3){
  # Generate coefficients
  spscales_sims = rep(NA, nsims)
  for (sim in 1:nsims){
    coeffs = runif(9106, min = 0, max = 0.1)
    vec1 = nest$decomp_mats[[i]] %*% coeffs
    spscales_sims[sim] = spscale(vec1, nbscal, edges)
  }
  spscales_nested[i] = max(spscales_sims)
  # calculate number of boundary edges
  pij = nest$proj_mats[[i]]
  bdes = 0
  for (k in 1:nrow(caladj)){
    idxs = 1*(pij[k,] == 0)
    bdes = bdes + sum(caladj[k,idxs], na.rm = T)
  }
  bdedges[i] = bdes
}

plot(bdedges/edges, spscales_nested, type = 'l')


