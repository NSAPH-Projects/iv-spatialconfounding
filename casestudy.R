library(spData)
library(sf)
library(spdep)
library(dplyr)
library(ggpubr)
library(gridExtra)

scale = function(x, nbs){ # make sure order of $X$ matches A
  avgdists = rep(NA, length(x))
  for (node in 1:length(x)){
    avgdists[node] = mean(abs(x[node] - x[nbs[[node]]]), na.rm = T)
  }
  return(avgdists)
}

#data(boston)
#str(boston.c)
#boston_sf <- boston.c %>% st_as_sf(., coords = c("LAT","LON"))
#plot(boston_sf["CRIM"])
#plot(boston_sf["MEDV"], logz = T)

boston.tr <- sf::st_read(system.file("shapes/boston_tracts.shp",
                                     package="spData")[1])
#boston.tr = boston.tr %>% mutate(across(where(is.numeric), scale))
#bostonmerged <- merge(boston.tr, boston.c,by = "TOWN")
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
caltr = st_read('/Users/sophie/Downloads/tl_2021_06_tract/tl_2021_06_tract.shp')
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

#### Calculate scale of all variables ####
nbs = poly2nb(boston.tr)
mean(scale(boston.tr$DIS, nbs = nbs))
mean(scale(boston.tr$logCRIM, nbs = nbs))

nbs = poly2nb(cal)
mean(scale(calmerged$SOVI_SCORE, nbs = nbs))
mean(scale(calmerged$HWAV_AFREQ, nbs = nbs))
