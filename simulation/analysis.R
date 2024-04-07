library(ggplot2)
library(sf)
library(mgcv)
library(dplyr)
library(utils)
library(xtable)
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

# Initialize dataframe that will store bias + ses. 
# this dataframe has columns outcomemod, projmat, method, bias, se
# and will have length(csvs) rows
analysisdf = data.frame(
  outcomemod = character(length(csvs)),
  projmat = character(length(csvs)),
  method = character(length(csvs)),
  absolutebias = numeric(length(csvs)),
  se = numeric(length(csvs))
)

for (i in 1:length(csvs)){
  if (i %% 15 == 1) { # If i is 1 or a multiple of 13 (e.g., 1, 13, 25,...), open a new png device
    png(
      filename = sprintf("images/plots_group_%d.png", ceiling(i / 15)),
      width = 1200,
      height = 800,
      res = 120
    )
    par(mfrow = c(3, 5)) # Set up the plot layout for 12 plots per page
  }
  filename = csvs[i]
  # read in the csv
  df_temp = read.csv(file.path('clusterresults', filename))
  # remove .csv from filename
  filename = gsub('.csv', '', filename)
  # split filename by _
  filename = unlist(strsplit(filename, '_'))
  analysisdf$projmat[i] = filename[1]
  analysisdf$outcomemod[i] = filename[2]
  if (filename[3] == 'keller'){
    analysisdf$method[i] = paste(filename[3], filename[6])
  }
  else{
    analysisdf$method[i] = filename[3]
  }
  muests = df_temp[, -(1:2)]
  met = metrics(simlist$a.vals,
          muests,
          df_temp$mutrue,
          A)
  analysisdf$absolutebias[i] = met$avgabsbias
  analysisdf$se[i] = met$avgse
  # Use matplot to plot
  matplot(
    x = simlist$a.vals,
    y = cbind(df_temp$mutrue, muests),
    type = "l",
    xlab = "Treatment level A=a",
    ylab = expression("E( Y"^"a" ~ ")"),
    col = c('black', rep('red', ncol(muests))),
    # Title is the filename
    main = paste(analysisdf$projmat[i], analysisdf$outcomemod[i], analysisdf$method[i])
  )
  if (i %% 15 == 0 || i == length(csvs)) { # If i is a multiple of 12 or the last plot, close the png device
    dev.off()
  }
}
# sort analysis df so that outcome mod is in order of linear, interaction, nonlinear
analysisdf$outcomemod = factor(analysisdf$outcomemod, levels = c('linear', 'interaction', 'nonlinear'))
analysisdf = analysisdf[order(analysisdf$outcomemod),]
analysisdflinear = analysisdf[analysisdf$outcomemod == 'linear',]
analysisdfinteraction = analysisdf[analysisdf$outcomemod == 'interaction',]
analysisdfnonlinear = analysisdf[analysisdf$outcomemod == 'nonlinear',]

# print latex table without rownames
print(xtable(analysisdflinear[,-1]), include.rownames=FALSE)
print(xtable(analysisdfinteraction[,-1]), include.rownames=FALSE)
print(xtable(analysisdfnonlinear[,-1]), include.rownames=FALSE)

