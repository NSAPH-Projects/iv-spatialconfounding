library(ggplot2)
library(sf)
library(tidyr)
library(mgcv)
library(dplyr)
library(utils)
library(xtable)
library(Matrix)
source('simfuncs.R')

load('sim.RData')
load('Vcs.RData')
load('Vucs.RData')
load('adjustmentmats.RData')

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

#U = createU(uscounties$latnorm, uscounties$longnorm)
#X = createX(U, uscounties$latnorm, uscounties$longnorm)

# extract all filenames of csvs in clusterresults
csvs = list.files('clusterresults_apr27', pattern = 'linear_KennedyERC_.csv') # _.csv
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

# Extract components from filenames
projmat = gsub("(_.*$)", "", csvs)
option = gsub("^(.*?_)(.*?)(_.*$)", "\\2", csvs)
method = gsub("^(.*?_.*?_)(.*)(\\_.csv$)", "\\2", csvs)

# Define custom orders
projmat_order = factor(projmat, levels = c("nestedregion", "nestedstate", "GFT150", "GFT51100"))
option_order = factor(option, levels = c("linear", "interaction", "nonlinear"))
method_order = factor(method, levels = c("unadjustedOLS", "bobb_exposurepenalized", "keller_szpiro_selectingscale_outcome","spatialplus", "KennedyERC", "spatialcoord", "oracleU"))

# Create an ordered dataframe
ordered_csvs = data.frame(csvs, projmat_order, option_order, method_order)

# Order dataframe by defined factors
ordered_csvs = ordered_csvs[order(ordered_csvs$projmat_order, ordered_csvs$option_order, ordered_csvs$method_order), "csvs"]

set.seed(123)
for (i in 1:length(ordered_csvs)){
  # if (i %% 12 == 1) { # If i is 1 plus a multiple of 6  open a new png device
  #   png(
  #     filename = sprintf("images/plots_group_%d_apr20.png", ceiling(i / 12)),
  #     width = 2000,
  #     height = 1200,
  #     res = 150
  #   )
  #   par(mfrow = c(3, 4)) # Set up the plot layout for 12 plots per page
  # }
  filename = ordered_csvs[i]
  # read in the csv
  df_temp = read.csv(file.path('clusterresults_apr27', filename))
  # remove .csv from filename
  filename = gsub('_.csv', '', filename)
  # split filename by _
  filename = unlist(strsplit(filename, '_'))
  analysisdf$projmat[i] = filename[1]
  analysisdf$outcomemod[i] = filename[2]
  # if (filename[3] == 'keller'){
  #   analysisdf$method[i] = paste(filename[3], filename[6])
  # }
  # else{
  analysisdf$method[i] = filename[3]
  #}
  #Vc = Vcs[[filename[1]]]
  #Vuc = Vucs[[filename[1]]]
  #Adat = createA_fromV(U,X,Vc = Vc, Vuc = Vuc)
  #A = Adat$A
  if (filename[2] == 'linear'){
    minmax = c(-5,15)
  }
  if (filename[2] == 'interaction'){
    minmax = c(-20,20)
  }
  if (filename[2] == 'nonlinear'){
    minmax = c(1,8)
  }
  # muests is last 100 columns of df_temp
  df_temp = df_temp[,c(1,(ncol(df_temp)-99):(ncol(df_temp)))]
  #muests = df_temp[, -1]
  # a.vals = simlist$a.vals[[filename[1]]]

  df_temp$mutrue = computemutrue(df_temp$a.vals, latnorm, longnorm, 
                         option = filename[2])
  #df_temp$a.vals = a.vals#simlist$a.vals[[filename[1]]]
  # met = metrics(a.vals = simlist$a.vals[[filename[1]]],
  #         muests = muests,
  #         mutrue = mutrue,
  #         A = A)
  #analysisdf$absolutebias[i] = met$avgabsbias
  #analysisdf$se[i] = met$avgse
  # Use matplot to plot
  df_temp = df_temp %>% 
    pivot_longer(cols = -c(a.vals, mutrue), names_to = "sim", values_to = "value") %>% 
    mutate(sim = factor(sim))
  g = ggplot(df_temp, aes(x = a.vals, y = value, color = sim)) + 
    geom_line() + 
    # add mutrue as a line
    geom_line(aes(x = a.vals, y = mutrue), color = "black", linetype = "dashed") +
    labs(x = "a", y = expression("E( Y"^"a" ~ ")")) + 
    scale_color_manual(values = c(rep(rgb(1, 0, 0, 0.2), 100), 'black')) + 
    theme_minimal() + 
    theme(
      legend.position = "none"
    ) +
    ylim(minmax)
  
  png(
    filename = sprintf("images/erfs/apr27/%s_%s_%s_2.png", analysisdf$projmat[i], 
                       analysisdf$outcomemod[i], analysisdf$method[i]),
    width = 1000,
    height = 1000,
    res = 250
  )
  # matplot(
  #   x = simlist$a.vals[[filename[1]]],
  #   y = cbind(muests, mutrue),
  #   type = "l",
  #   xlab = "Treatment level A=a",
  #   ylab = expression("E( Y"^"a" ~ ")"),
  #   col = c(rep(rgb(1, 0, 0, 0.2), ncol(muests)), 'black'),  # red with 5% transparency
  #   lty = c(rep(1, ncol(muests)), 1),
  #   lwd = c(rep(1, ncol(muests)), 2),
  #   ylim = minmax
  #   )
  print(g)
  dev.off()
  # if (i %% 12 == 0 || i == length(csvs)) { # If i is a multiple of 12 or the last plot, close the png device
  #   dev.off()
  # }
}
analysisdf[,4:5] = round(analysisdf[,4:5], 3)
analysisdf

# print latex table without rownames
# print(xtable(analysisdflinear[,-1]), include.rownames=FALSE)
# print(xtable(analysisdfinteraction[,-1]), include.rownames=FALSE)
# print(xtable(analysisdfnonlinear[,-1]), include.rownames=FALSE)

## Non oracle method!
csvs = list.files('clusterresults_may1', pattern = '.csv')

# Initialize dataframe that will store bias + ses. 
# this dataframe has columns outcomemod, projmat, method, bias, se
# and will have length(csvs) rows
analysisdf_sensitivity = data.frame(
  outcomemod = character(length(csvs)),
  projmat = character(length(csvs)),
  adjustmentbasis = character(length(csvs)),
  absolutebias = numeric(length(csvs)),
  se = numeric(length(csvs))
)

# Extract components from filenames
projmat = gsub("(_.*$)", "", csvs)
option = gsub("^(.*?_)(.*?)(_.*$)", "\\2", csvs)

# Define custom orders
projmat_order = factor(projmat, levels = c("nestedregion", "nestedstate", "GFT150", "GFT51100"))
option_order = factor(option, levels = c("linear", "interaction", "nonlinear"))

# Create an ordered dataframe
ordered_csvs = data.frame(csvs, projmat_order, option_order)

# Order dataframe by defined factors
ordered_csvs = ordered_csvs[order(ordered_csvs$projmat_order, ordered_csvs$option_order), "csvs"]

# for adjustmentmat in adjustmentmats calculate spscale
spscales = rep(NA,length(adjustmentmats))
for (i in 1:length(adjustmentmats)){
  print(i)
  adjustmentmat = adjustmentmats[[i]]
  spscales[i] = spscale(adjustmentmat, adjmat, nsims = 100)
}
names(spscales) = names(adjustmentmats)
# sort spscales
round(sort(spscales),6)
save(spscales, file = 'spscales.RData')

# Alternatively if just GFT
E = eigen(diag(rowSums(adjmat))-adjmat)
Esub = E$vectors[,(n-1):(n-200)] 
spscales = rep(NA,length(adjustmentmats))
for (j in 1:10){
  print(c(1,(10*j)))
  spscales[j] = mean(E$values[(n-10*j):(n-1)])
}
names(spscales) = names(adjustmentmats)
# sort spscales
round(sort(spscales),6)


bin_to_yval <- list(
  `120` = 20, `2140` = 20, `4160` = 20, `6180` = 20, `81100` = 20,
  `101120` = 20, `121140` = 20, `141160` = 20, `161180` = 20, `181200` = 20,
  `150` = 50, `51100` = 50, `101150` = 50, `151200` = 50,
  `1100` = 100, `101200` = 100,
  `110` = 10, `130` = 30, `140` = 40, `160` = 60,
  `170` = 70, `180` = 80, `190` = 90
)

set.seed(123)
gs = list()
xvals = rep(NA, length(ordered_csvs))
yvals = rep(NA, length(ordered_csvs))
for (i in 1:length(ordered_csvs)){
  
  filename = ordered_csvs[i]
  df_temp = read.csv(file.path('clusterresults_may1', filename))
  filename = gsub('.csv', '', filename)
  # split filename by _
  filename = unlist(strsplit(filename, '_'))
  analysisdf_sensitivity$projmat[i] = filename[1]
  analysisdf_sensitivity$outcomemod[i] = filename[2]
  analysisdf_sensitivity$adjustmentbasis[i] = filename[4]
  df_temp$mutrue = computemutrue(df_temp$a.vals, latnorm, longnorm, option = filename[2])
  xvals[i] = spscales[[filename[4]]]
  # if filename[4] starts with nested, yval is region or state
  if (grepl('nested', filename[4])){
    yvals[i] = rankMatrix(adjustmentmats[[filename[4]]]) # nested
  } else {
    # remove GFT from filename[4]
    bin = gsub('GFT', '', filename[4])
    yvals[i] = bin_to_yval[[bin]]
  }
  # if filename
  projmat = projmats[[filename[1]]]
  if (filename[2] == 'linear'){
    minmax = c(0,5)
  }
  if (filename[2] == 'interaction'){
    minmax = c(-10,10)
  }
  if (filename[2] == 'nonlinear'){
    minmax = c(1,5)
  }
  
  muests = df_temp[, -(1:2)]
  # Use matplot to plot
  df_temp = df_temp %>% 
    pivot_longer(cols = -c(a.vals, mutrue), names_to = "sim", values_to = "value") %>% 
    mutate(sim = factor(sim))
  
  gs[[i]] = ggplot(df_temp, aes(x = a.vals, y = value, color = sim)) + 
    geom_line() + 
    # add mutrue as a line
    geom_line(aes(x = a.vals, y = mutrue), color = "black", linetype = "dashed") +
    labs(x = "a", y = expression("E( Y"^"a" ~ ")")) + 
    scale_color_manual(values = c(rep(rgb(1, 0, 0, 0.2), ncol(muests)), 'black')) + 
    theme_minimal() + 
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ylim(minmax) + # add title
    ggtitle(analysisdf_sensitivity$adjustmentbasis[i])
}
analysisdf_sensitivity[,4:5] = round(analysisdf_sensitivity[,4:5], 3)
analysisdf_sensitivity

# print latex table without rownames
# print(xtable(analysisdflinear[,-1]), include.rownames=FALSE)
# print(xtable(analysisdfinteraction[,-1]), include.rownames=FALSE)
# print(xtable(analysisdfnonlinear[,-1]), include.rownames=FALSE)


# Base plot and positions
# nestedstate_idxs are the files that begin with nestedstate
# nestedregion_idxs are the files that begin with nestedregion
combs = expand.grid(unique(projmat_order), unique(option_order))
colnames(combs) = c('projmat', 'outcomemod')

for (c in 1:nrow(combs)){
  projmatname = combs$projmat[c]
  outcomemodname = combs$outcomemod[c]
  filename = paste0(projmatname, '_', outcomemodname)
  idxs = which(analysisdf_sensitivity$projmat == projmatname &
                 analysisdf_sensitivity$outcomemod == outcomemodname)
  base_plot = ggplot() + xlim(0, 0.00017) + ylim(0, 110)
  positions = data.frame(xmin = xvals[idxs] - 0.000005,
                          xmax = xvals[idxs] + 0.000005,
                          ymin = yvals[idxs] - 7,
                          ymax = yvals[idxs]  + 7)
  
  # Add plots to the base plot
  plots = gs[idxs]
  
  for (i in seq_along(plots)) {
    p = ggplotGrob(plots[[i]])
    base_plot = base_plot +
      annotation_custom(grob = p, xmin = positions$xmin[i], xmax = positions$xmax[i],
                        ymin = positions$ymin[i], ymax = positions$ymax[i])
  }
  
  # Highlight the second plot
  highlight_index = which(analysisdf_sensitivity[idxs,]$adjustmentbasis == projmatname)
  base_plot2 = base_plot +
    geom_rect(aes(xmin = positions$xmin[highlight_index], 
                  xmax = positions$xmax[highlight_index]+0.000004,
                  ymin = positions$ymin[highlight_index], 
                  ymax = positions$ymax[highlight_index]),
              color = "green", fill = NA, size = 2, linetype = "solid") + 
    #ggtitle(paste("Sensitivity Analysis: ", projmatname, outcomemodname)) + 
    xlab(expression("Spatial Frequency of Confounding Adjustment Subspace V, " * E[v %~% "Unif"(V) * ", " * "|" * v * "|" == 1] * 
                      "mean"[A[ij] == 1] * "|v"[i] - "v"[j] * "|"^{2})) +
    ylab('Dimension of Basis') + 
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5), 
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), 
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), 
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.line.x = element_line(color = 'blue',
                                 arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both"))
    ) + # relabel ticks on y axis
    # scale_y_continuous(breaks = 0:5,
    #                    labels=c("",
    #                             "nested",
    #                             "GFT binsize 100",
    #                             "GFT binsize 50",
    #                             "GFT binsize 20",
    #                             ""),
    #                    limits = c(0.5, 4.5)) + 
    annotate(
      "text",
      x = 0.000005,
      y = 0.55,
      label = "Large spatial scale",
      size = 5,
      color = 'blue',
      fontface = "italic"
    ) +
    annotate(
      "text",
      x = 0.00016,
      y = 0.55,
      label = "Small spatial scale",
      size = 5,
      color = 'blue',
      fontface = "italic"
    )
  
  png(
    filename = paste0("images/erfs/sensitivity/may1/", filename, "3.png"),
    width = 3000,
    height = 1800,
    res = 156
  )
  print(base_plot2)
  dev.off()
}


####################################################################
nestedregion_idxs = which(
  analysisdf_sensitivity$projmat == 'nestedregion' &
    analysisdf_sensitivity$outcomemod == 'nonlinear'
)
# convert yvals to numeric vector 1:4
base_plot = ggplot() + xlim(0, 0.00017) + ylim(0, 5)
positions = data.frame(xmin = xvals[nestedregion_idxs] - 0.0000065,
                        xmax = xvals[nestedregion_idxs] + 0.0000065,
                        ymin = yvals[nestedregion_idxs] - 0.45,
                        ymax = yvals[nestedregion_idxs]  + 0.45)

# Add plots to the base plot
plots = gs[nestedregion_idxs]

for (i in seq_along(plots)) {
  p = ggplotGrob(plots[[i]])
  base_plot = base_plot +
    annotation_custom(grob = p, xmin = positions$xmin[i], xmax = positions$xmax[i],
                      ymin = positions$ymin[i], ymax = positions$ymax[i])
}

# Highlight the second plot
highlight_index = 17
base_plot2 = base_plot +
  geom_rect(
    aes(
      xmin = positions$xmin[highlight_index],
      xmax = positions$xmax[highlight_index] + 0.000003,
      ymin = positions$ymin[highlight_index],
      ymax = positions$ymax[highlight_index]
    ),
    color = "green",
    fill = NA,
    size = 2,
    linetype = "solid"
  ) +
  #ggtitle("Sensitivity Analysis: Nested Region Nonlinear") + 
  xlab(expression("Spatial Frequency of Confounding Adjustment Subspace V, " * E[v %~% "Unif"(V) * ", " * "|" * v * "|" == 1] * 
                    "mean"[A[ij] == 1] * "|v"[i] - "v"[j] * "|^2")) +
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5), 
    axis.title.x = element_text(size = 20,
                                face = "bold",
                                margin = margin(
                                  t = 20,
                                  r = 0,
                                  b = 0,
                                  l = 0
                                )),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.line.x = element_line(color = 'blue',
                               arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both"))
  ) + # relabel ticks on y axis
  scale_y_continuous(breaks = 0:5,
                   labels=c("",
                            "nested",
                            "GFT binsize 100",
                            "GFT binsize 50",
                            "GFT binsize 20",
                            ""),
                   limits = c(0.5, 4.5)) + 
  annotate("text", x = 0, y = 0.55, label = "Large spatial scale", size = 5, color = 'blue', fontface = "italic") +
  annotate("text", x = 0.00017, y = 0.55, label = "Small spatial scale", size = 5, color = 'blue', fontface = "italic")

png(
  filename = "images/test.png",
  width = 2000,
  height = 1200,
  res = 100
)
base_plot2
dev.off()

############### MAY 7 ############

combs = expand.grid(unique(projmat_order), unique(option_order))
colnames(combs) = c('projmat', 'outcomemod')

for (c in 1:nrow(combs)){
  projmatname = combs$projmat[c]
  outcomemodname = combs$outcomemod[c]
  filename = paste0(projmatname, '_', outcomemodname)
  idxs = which(analysisdf_sensitivity$projmat == projmatname &
                 analysisdf_sensitivity$outcomemod == outcomemodname)
  base_plot = ggplot() + xlim(0, 0.33) + ylim(0, 110)
  positions = data.frame(xmin = xvals[idxs] - 0.02,
                         xmax = xvals[idxs] + 0.02,
                         ymin = yvals[idxs] - 10,
                         ymax = yvals[idxs]  + 10)
  
  # Add plots to the base plot
  plots = gs[idxs]
  
  for (i in seq_along(plots)) {
    p = ggplotGrob(plots[[i]])
    base_plot = base_plot +
      annotation_custom(grob = p, xmin = positions$xmin[i], xmax = positions$xmax[i],
                        ymin = positions$ymin[i], ymax = positions$ymax[i])
  }
  
  # Highlight the second plot
  highlight_index = which(analysisdf_sensitivity[idxs,]$adjustmentbasis == projmatname)
  base_plot2 = base_plot +
    geom_rect(aes(xmin = positions$xmin[highlight_index], 
                  xmax = positions$xmax[highlight_index]+0.000003,
                  ymin = positions$ymin[highlight_index], 
                  ymax = positions$ymax[highlight_index]),
              color = "green", fill = NA, size = 2, linetype = "solid") + 
    #ggtitle(paste("Sensitivity Analysis: ", projmatname, outcomemodname)) + 
    xlab(expression("Spatial Frequency of Confounding Adjustment Subspace V, " * E[v %~% "Unif"(V) * ", " * "|" * v * "|" == 1] * 
                      "mean"[A[ij] == 1] * "|v"[i] - "v"[j] * "|"^{2})) +
    ylab('Dimension of Basis') + 
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 35, face = "bold", hjust = 0.5), 
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), 
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(t = 20, r = 0, b = 0, l = 0)), 
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.line.x = element_line(color = 'blue',
                                 arrow = grid::arrow(length = unit(0.3, "cm"), ends = "both"))
    ) + # relabel ticks on y axis
    # scale_y_continuous(breaks = 0:5,
    #                    labels=c("",
    #                             "nested",
    #                             "GFT binsize 100",
    #                             "GFT binsize 50",
    #                             "GFT binsize 20",
    #                             ""),
    #                    limits = c(0.5, 4.5)) + 
    annotate(
      "text",
      x = 0,
      y = 0.55,
      label = "Large spatial scale",
      size = 5,
      color = 'blue',
      fontface = "italic"
    ) +
    annotate(
      "text",
      x = 0.35,
      y = 0.55,
      label = "Small spatial scale",
      size = 5,
      color = 'blue',
      fontface = "italic"
    )
  
  png(
    filename = paste0("images/erfs/sensitivity/may7/", filename, ".png"),
    width = 2000,
    height = 1200,
    res = 100
  )
  print(base_plot2)
  dev.off()
}


