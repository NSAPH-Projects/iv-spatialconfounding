library(ggplot2)
library(sf)
library(mgcv)
library(dplyr)
library(utils)
source('simfuncs.R')

load('sim.RData')
load('projmats.RData')

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

projmat = projmats[['nestedstate']]
U = createU(uscounties$latnorm, uscounties$longnorm)
X = createX(U, uscounties$latnorm, uscounties$longnorm)
Adat = createA(U,X,projmat)
A = Adat$A

# extract all filenames of csvs in clusterresults
csvs = list.files('clusterresults', pattern = '.csv')

# Initialize dataframe that will store bias + MSEs. 
# this dataframe has columns outcomemod, projmat, method, bias, mse
# and will have length(csvs) rows
analysisdf = data.frame(
  projmat = character(length(csvs)),
  outcomemod = character(length(csvs)),
  method = character(length(csvs)),
  bias = numeric(length(csvs)),
  mse = numeric(length(csvs))
)

for (i in 1:length(csvs)){
  filename = csvs[i]
  # read in the csv
  df_temp = read.csv(file.path('clusterresults', filename))
  # remove .csv from filename
  filename = gsub('.csv', '', filename)
  # split filename by _
  filename = unlist(strsplit(filename, '_'))
  analysisdf$projmat[i] = filename[1]
  analysisdf$outcomemod[i] = filename[2]
  analysisdf$method[i] = filename[3]
  muests = df_temp[, -(1:2)]
  met = metrics(a.vals,
          muests,
          df_temp$mutrue,
          A)
  analysisdf$bias[i] = met$avgabsbias
  analysisdf$mse[i] = met$avgRMSE
  # Use matplot to plot
  matplot(
    x = a.vals,
    y = cbind(df_temp$mutrue, muests),
    type = "l",
    xlab = "Treatment level A=a",
    ylab = expression("E( Y"^"a" ~ ")"),
    col = c('black', rep('red', ncol(muests))),
    # Title is the filename
    main = paste(filename)
  )
}


