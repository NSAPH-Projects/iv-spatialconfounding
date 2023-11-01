``` r
library(gasper)
library(spData)
```

    ## To access larger datasets in this package, install the spDataLarge
    ## package with: `install.packages('spDataLarge',
    ## repos='https://nowosad.github.io/drat/', type='source')`

``` r
library(sf)
```

    ## Linking to GEOS 3.11.0, GDAL 3.5.3, PROJ 9.1.0; sf_use_s2() is TRUE

``` r
library(spdep)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(ggplot2)
source('funcs.R')
```

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

## Import and Process Data

``` r
boston.tr <- sf::st_read(system.file("shapes/boston_tracts.shp",
                                     package="spData")[1])
```

    ## Reading layer `boston_tracts' from data source 
    ##   `/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/spData/shapes/boston_tracts.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 506 features and 36 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -71.52311 ymin: 42.00305 xmax: -70.63823 ymax: 42.67307
    ## Geodetic CRS:  NAD27

``` r
boston_counties <- st_read("townssurvey_shp/TOWNSSURVEY_POLY.shp")
```

    ## Reading layer `TOWNSSURVEY_POLY' from data source 
    ##   `/Users/sophie/Documents/SpatialConf/townssurvey_shp/TOWNSSURVEY_POLY.shp' 
    ##   using driver `ESRI Shapefile'
    ## Simple feature collection with 1239 features and 22 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 33863.73 ymin: 777606.4 xmax: 330837 ymax: 959743
    ## Projected CRS: NAD83 / Massachusetts Mainland

``` r
# sourced from https://www.mass.gov/info-details/massgis-data-municipalities#downloads-

boston.tr$TOWN[1:132] = 'Boston'
boston.tr$TOWN[186:189] = 'Saugus'
boston.tr$TOWN = toupper(boston.tr$TOWN)
#boston.tr$TOWN %in% boston_counties$TOWN 
bosmerged = merge(boston.tr, cbind.data.frame('TOWN' = boston_counties$TOWN, 
                                              'FIPS' = boston_counties$FIPS_COUNT),
                  by = 'TOWN')
bosmerged = distinct(bosmerged)
nrow(bosmerged)
```

    ## [1] 506

``` r
length(unique(bosmerged$TRACT))
```

    ## [1] 506

``` r
length(unique(bosmerged$TOWN))
```

    ## [1] 78

``` r
length(unique(bosmerged$FIPS))
```

    ## [1] 5

## Crime and Distance: Exploratory Analysis

``` r
bosmerged$logCRIM = log(bosmerged$CRIM)
bosmerged$logDIS = log(bosmerged$DIS)
par(mfrow = c(1,2))
hist(bosmerged$logCRIM)
hist(bosmerged$logDIS)
```

![](boston_spatial_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
g1 = ggplot(bosmerged) +
  geom_sf(aes(fill = logCRIM), color = NA) +
  scale_fill_viridis_c() + 
  labs(title = 'Log Crime') + 
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))

g2 = ggplot(bosmerged) +
  geom_sf(aes(fill = logDIS), color = NA) +
  scale_fill_viridis_c() + 
  labs(title = 'Log Distance to City') + 
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))

grid.arrange(grobs = list(g1, g2))
```

![](boston_spatial_files/figure-markdown_github/unnamed-chunk-3-2.png)

## Nested Filtering

``` r
groups = cbind(bosmerged$FIPS, bosmerged$TOWN)
nest = nested_decomp_mats(groups)
bosmerged$nested_county = nest$decomp_mats[[1]] %*% bosmerged$logCRIM
bosmerged$nested_town = nest$decomp_mats[[2]] %*% bosmerged$logCRIM
bosmerged$nested_tract = nest$decomp_mats[[3]] %*% bosmerged$logCRIM

g3 = ggplot(bosmerged) +
  geom_sf(aes(fill = nested_county), color = NA) +
  scale_fill_viridis_c() + 
  labs(title = 'Log Distance to City') + 
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))
g4 = ggplot(bosmerged) +
  geom_sf(aes(fill = nested_town), color = NA) +
  scale_fill_viridis_c() + 
  labs(title = 'Log Distance to City') + 
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))
g5 = ggplot(bosmerged) +
  geom_sf(aes(fill = nested_tract), color = NA) +
  scale_fill_viridis_c() + 
  labs(title = 'Log Distance to City') + 
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "transparent"))

grid.arrange(grobs = list(g3,g4,g5), ncol = 3)
```

![](boston_spatial_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
corcounty = cor(bosmerged$nested_county, bosmerged$logDIS)
cortown = cor(bosmerged$nested_town, bosmerged$logDIS)
cortract= cor(bosmerged$nested_tract, bosmerged$logDIS)

datnested = data.frame(Level = c('County', 'Town', 'Tract'),
                   Correlation = c(corcounty, cortown, cortract))

# Plot correlations
ggplot(datnested, aes(x = Level, y = Correlation)) +
  geom_point(shape = 19) +
  geom_hline(yintercept = cor(bosmerged$logCRIM, bosmerged$logDIS), 
             linetype = "dashed", color = "red") +
  labs(title = "Nested Filtering",
       y = "Cor(VVtX, U)") +
  theme_minimal() +
  scale_x_discrete(labels = c('County', 'Town', 'Tract'))
```

![](boston_spatial_files/figure-markdown_github/unnamed-chunk-4-2.png)

## Fourier Filtering

``` r
nbsbos = poly2nb(bosmerged) 
bosadj = nb2mat(neighbours = nbsbos, style = 'B', zero.policy = T)
L = diag(rowSums(bosadj)) - bosadj
E = eigen(L)
evalues = E$values
evectors = E$vectors

gs = list()
maxevals = rep(NA, 10)
cors = rep(NA, 10)
for (i in 1:10){
  ixs = ((i-1)*50+6):(50*i+5)
  V = evectors[,ixs]
  maxeval = evalues[(i-1)*50+6]
  maxevals[i] = maxeval
  name = paste('eigen', round(maxeval,2), sep = '')
  print(name)
  eigenpart = V%*%t(V) %*% bosmerged$logCRIM
  bosmerged[[name]] = eigenpart
  gs[[i]] = ggplot(bosmerged) +
    geom_sf(aes_string(fill = name), color = NA) +
    scale_fill_viridis_c() + 
    labs(title = paste('eigenval =', round(maxeval,2))) + 
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          line = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_line(colour = "transparent"))
  cors[i] = cor(bosmerged$logDIS, eigenpart)
}
```

    ## [1] "eigen11.8"

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## [1] "eigen9.35"
    ## [1] "eigen8.35"
    ## [1] "eigen7.47"
    ## [1] "eigen6.72"
    ## [1] "eigen5.95"
    ## [1] "eigen5.04"
    ## [1] "eigen3.9"
    ## [1] "eigen2.88"
    ## [1] "eigen1.47"

``` r
grid.arrange(grobs = gs)
```

![](boston_spatial_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# Plot correlations
ggplot(data.frame(maxevals = maxevals, cors = cors), aes(x = maxevals, y = cors)) +
  geom_point(shape = 19) +
  geom_hline(yintercept = cor(bosmerged$logCRIM, bosmerged$logDIS), 
             linetype = "dashed", color = "red") +
  labs(x = "lambda", y = "Corr(U, VVtX)", title = "Fourier Filtering") +
  theme_minimal()
```

![](boston_spatial_files/figure-markdown_github/unnamed-chunk-5-2.png)
