library(spData)
library(sf)
library(spdep)
library(dplyr)
library(ggpubr)
library(gridExtra)

spscale = function(x, nbs){ # make sure order of $X$ matches A
  avgdists = rep(NA, length(x))
  for (node in 1:length(x)){
    avgdists[node] = mean(abs(x[node] - x[nbs[[node]]]), na.rm = T)
  }
  return(avgdists)
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
mean(spscale(scale(boston.tr$DIS), nbs = nbsbos))
mean(spscale(scale(boston.tr$logCRIM), nbs = nbsbos))

nbscal = poly2nb(calmerged)
mean(spscale(scale(calmerged$SOVI_SCORE), nbs = nbscal), na.rm = T) # 2 NAs
mean(spscale(scale(calmerged$HWAV_AFREQ), nbs = nbscal), na.rm = T)
