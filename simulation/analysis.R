library(ggplot2)
library(sf)
library(tidyr)
library(mgcv)
library(dplyr)
library(utils)
library(xtable)
library(Matrix)
library(tidyverse)
library(parallel)
library(geosphere)
library(patchwork)

source('../funcs.R')
load('sim.RData')

RESULTS_DIR <- "results_Mar27/"

csvs <- list.files(RESULTS_DIR, pattern = '\\.csv$')
csvs <- csvs[!grepl('_(ci_lower|ci_upper|time|n_uc)\\.csv$', csvs)]


# Create storage for metrics
analysisdf <- data.frame(
  confounding_mechanism = character(length(csvs)),
  option = character(length(csvs)),
  method = character(length(csvs)),
  bias = numeric(length(csvs)),
  RMSE = numeric(length(csvs)),
  se = numeric(length(csvs)),
  coverage = numeric(length(csvs))
)

# Extract components from filenames
confounding_mechanism <- as.integer(sub("^conf(\\d+)_.*$", "\\1", csvs))
option <- sub("^conf\\d+_([^_]+)_.*$", "\\1", csvs)
method <- sub("^conf\\d+_[^_]+_([^_]+)\\.csv$", "\\1", csvs)

# Precompute true estimand for each outcome model and confounding mechanism
mutrues <- data.frame(expand.grid(
  confounding_mechanism = 1:7,
  option = c('linear', 'nonlinear')))
for (i in 1:nrow(mutrues)){
  fname <- paste0(RESULTS_DIR, 'conf', mutrues$confounding_mechanism[i],
                  '_', mutrues$option[i], '_oracle.csv')
  mutrues$theta[i] <- if (file.exists(fname))
    mean(as.vector(as.matrix(read.csv(fname))), na.rm = TRUE)
  else
    NA_real_
}

# Loop through results to calculate metrics and create plots.
for (i in 1:length(csvs)){
  filename <- csvs[i]
  print(filename)
  
  analysisdf$confounding_mechanism[i] <- confounding_mechanism[i]
  analysisdf$option[i] <- option[i]
  analysisdf$method[i] <- method[i]
  df_temp <- read.csv(file.path(RESULTS_DIR, filename))
  
  muests <- df_temp 
  # Convert muests to a vector, it's just a single column
  muests <- as.vector(as.matrix(muests))
  
  # Compute true truncated exposure estimate
  mutrue <- mutrues[mutrues$confounding_mechanism == confounding_mechanism[i] & 
                      mutrues$option == option[i],]$theta 
  df_temp$mutrue <- mutrue
  
  # Save metrics in analysisdf
  analysisdf$bias[i] <- mean(muests, na.rm = T) - mutrue
  analysisdf$RMSE[i] <- sqrt(mean((muests - mutrue)^2, na.rm = T))
  analysisdf$se[i] <- sd(muests, na.rm = T)

  # Compute coverage from CI files if they exist
  fname_lower <- file.path(RESULTS_DIR,
                           sub('\\.csv$', '_ci_lower.csv', filename))
  fname_upper <- file.path(RESULTS_DIR,
                           sub('\\.csv$', '_ci_upper.csv', filename))
  if (file.exists(fname_lower) && file.exists(fname_upper)) {
    lower <- as.vector(as.matrix(read.csv(fname_lower)))
    upper <- as.vector(as.matrix(read.csv(fname_upper)))
    analysisdf$coverage[i] <- mean(lower <= mutrue & upper >= mutrue, na.rm = TRUE)
  } else {
    analysisdf$coverage[i] <- NA_real_
  }
}

# Ensure all CM × option × method combinations appear even if CSVs are missing
all_methods <- c("oracle", "baseline", "spatialcoord", "trueIV",
                 "IV-TPS", "IV-GraphLaplacian",
                 "trueIV-spatialcoord", "IV-TPS-spatialcoord", "IV-GraphLaplacian-spatialcoord")
full_grid <- expand.grid(
  confounding_mechanism = as.character(1:7),
  option = c("linear", "nonlinear"),
  method = all_methods,
  stringsAsFactors = FALSE
)
analysisdf$confounding_mechanism <- as.character(analysisdf$confounding_mechanism)
analysisdf <- full_grid %>%
  left_join(analysisdf, by = c("confounding_mechanism", "option", "method"))

# Format the numeric columns (bias, RMSE, se) in scientific notation with 3 decimals
analysisdf_bias <- analysisdf %>%
  mutate(
    bias = round(bias*100, 3),
    RMSE = round(RMSE*100, 3),
    se = round(se*100, 3)
  )

# Pivot wider from the original analysisdf
desired_method_cols <- c("confounding_mechanism", "option",
                         "oracle", "baseline", "spatialcoord", "trueIV",
                         "IV-TPS", "IV-GraphLaplacian",
                         "trueIV-spatialcoord", "IV-TPS-spatialcoord",
                         "IV-GraphLaplacian-spatialcoord")
analysisdf_bias <- analysisdf_bias[, 1:4] %>%
  pivot_wider(names_from = method, values_from = bias) %>%
  select(any_of(desired_method_cols))

# Print using xtable and prevent xtable from reformatting the already-formatted text
print(xtable(analysisdf_bias), 
      include.rownames = FALSE, sanitize.text.function = identity)
# Print absolute bias using xtable
print(xtable(analysisdf_bias %>%
               mutate(across(where(is.numeric), ~ abs(.)))), 
      include.rownames = FALSE, sanitize.text.function = identity)
# Print absolute bias with the reordered confounding scenarios
analysisdf_bias_reordered <- analysisdf_bias %>%
  mutate(confounding_mechanism = factor(confounding_mechanism,
                                       levels = c("1", "2", "3", "4", "5", "6", "7"))) %>%
  arrange(confounding_mechanism, option)
print(xtable(analysisdf_bias_reordered %>%
               mutate(across(where(is.numeric), ~ abs(.)))), 
      include.rownames = FALSE, sanitize.text.function = identity)
# Do the same with RMSE
analysisdf_RMSE <- analysisdf %>%
  mutate(
    bias = round(bias*100, 3),
    RMSE = round(RMSE*100, 3),
    se = round(se*100, 3)
  )

analysisdf_RMSE <- analysisdf_RMSE[, c(1:3,5)] %>%
  pivot_wider(names_from = method, values_from = RMSE) %>%
  select(any_of(desired_method_cols))

print(xtable(analysisdf_RMSE), include.rownames = FALSE, sanitize.text.function = identity)
# Print the reordered RMSE
analysisdf_RMSE_reordered <- analysisdf_RMSE %>%
  mutate(confounding_mechanism = factor(confounding_mechanism,
                                       levels = c("1", "2", "3", "4", "5", "6", "7"))) %>%
  arrange(confounding_mechanism, option)
print(xtable(analysisdf_RMSE_reordered), include.rownames = FALSE, sanitize.text.function = identity)

# Now create facet_wrap boxplots with ggplot2

read_estimates <- function(i) {
  tryCatch({
    filename <- csvs[i]
    print(filename)

    cm_i     <- confounding_mechanism[i]
    option_i <- option[i]
    method_i <- method[i]
    dat <- read.csv(file.path(RESULTS_DIR, filename))

    # If the CSV doesn't have a header and just one column, name it "estimate"
    if (!"estimate" %in% colnames(dat)) {
      names(dat)[1] <- "estimate"
    }
    # Add the new columns
    dat <- dat %>%
      mutate(confounding_mechanism = cm_i,
             option = option_i,
             method = method_i)
    return(dat)
  }, error = function(e) {
    warning(paste("Skipping", csvs[i], ":", conditionMessage(e)))
    NULL
  })
}

# Read all files and combine into one data frame
df <- map_dfr(1:length(csvs), read_estimates)
# Rename the method "IV-GraphLaplacian" to "IV-GL" in df
df$method[df$method == 'IV-GraphLaplacian-spatialcoord'] = "IV-GL+spatialcoord"
df$method[df$method == 'IV-GraphLaplacian'] = "IV-GL"
df$method[df$method == 'IV-TPS-spatialcoord'] = "IV-TPS+spatialcoord"
df$method[df$method == 'trueIV-spatialcoord'] = "trueIV+spatialcoord"


desired_order <- c("oracle", "baseline", "spatialcoord", 
                   "trueIV", 
                   "IV-TPS", "IV-GL", 
                   "trueIV+spatialcoord", "IV-TPS+spatialcoord", "IV-GL+spatialcoord"
                   )
df$method <- factor(df$method, levels = desired_order)

df <- df %>%
  mutate(confounding_mechanism = factor(confounding_mechanism,
                                        levels = as.character(1:7)),
         option = factor(option, levels = c("linear", "nonlinear")))
mutrues <- mutrues %>% 
  mutate(confounding_mechanism = factor(confounding_mechanism),
         option = factor(option, levels = c("linear", "nonlinear")))
# Confounding mechanisms 1-7 are already in the correct display order in the CSVs.
df$confounding_mechanism_reordered <- factor(df$confounding_mechanism,
                                              levels = c("1", "2", "3", "4", "5", "6", "7"))
mutrues$confounding_mechanism_reordered <- factor(mutrues$confounding_mechanism,
                                                 levels = c("1", "2", "3", "4", "5", "6", "7"))
print(xtable(select(mutrues, confounding_mechanism_reordered, option, theta) %>% 
               arrange(confounding_mechanism_reordered), digits = 4), 
      include.rownames = FALSE)

# Create the boxplot with horizontal lines for theta
method_cols <- c(
  "oracle"                 = "gray",
  "baseline"               = "#D62728",
  "spatialcoord"           = "purple", 
  "trueIV"                 = "lightgreen", 
  "IV-TPS"                 = "#9ECAE1",
  "IV-GL"                  = "#FDAE6B", 
  "trueIV+spatialcoord"    = "#2E8B57", 
  "IV-TPS+spatialcoord"    = "#1F77B4", 
  "IV-GL+spatialcoord"     = "#E6550D" 
)
png("images/boxplot_Sep6.png", width = 2500, height = 1250, res = 200)
ggplot(df, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.5, outliers = F, staplewidth = 1) + #, draw_quantiles = c(0.5)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "blue") +
  ggh4x::facet_grid2(
    option ~ confounding_mechanism_reordered,
    scales = "free",           # allows different scales per row/col
    independent = "all"        # allows different scales **per panel**
  ) +
  #facet_grid(option ~ confounding_mechanism, scales = "free_y") +
  geom_hline(data = mutrues, aes(yintercept = theta), 
             color = "red", linetype = "twodash", size = 1) +
  #labs(x = NULL, y = "Truncated Exposure Effect Estimate") +   
  labs(
    x = "Confounding mechanism (1–7)",
    y = "Truncated Exposure Effect Estimate"
  ) +
  # Remove x-axis title
  scale_fill_manual(
    name = "Method",
    values = method_cols,
    breaks = names(method_cols) # matches your desired order
  ) +
  theme_bw() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") 
dev.off()

################# COVERAGE TABLES AND PLOTS #################

# Coverage table (pivot wider, same pattern as bias/RMSE)
analysisdf_cov <- analysisdf %>%
  mutate(coverage = round(coverage, 3))

analysisdf_cov <- analysisdf_cov[, c(1:3, 7)] %>%
  pivot_wider(names_from = method, values_from = coverage) %>%
  select(any_of(desired_method_cols))

analysisdf_cov_reordered <- analysisdf_cov %>%
  mutate(confounding_mechanism = factor(confounding_mechanism,
                                        levels = c("1", "2", "3", "4", "5", "6", "7"))) %>%
  arrange(confounding_mechanism, option)
print(xtable(analysisdf_cov_reordered, digits = 3),
      include.rownames = FALSE, sanitize.text.function = identity)

# Long-format coverage data for plotting
# Only methods that produce CIs (non-NA coverage)
cov_long <- analysisdf %>%
  filter(!is.na(coverage)) %>%
  mutate(
    confounding_mechanism = factor(confounding_mechanism,
                                   levels = as.character(1:7)),
    option = factor(option, levels = c("linear", "nonlinear")),
    method = recode(method,
                    "IV-GraphLaplacian-spatialcoord" = "IV-GL+spatialcoord",
                    "IV-GraphLaplacian"              = "IV-GL",
                    "IV-TPS-spatialcoord"            = "IV-TPS+spatialcoord",
                    "trueIV-spatialcoord"            = "trueIV+spatialcoord"),
    method = factor(method, levels = desired_order)
  )

# Coverage point plot: one point per method per panel, reference line at 0.95
png("images/coverage_Mar27.png", width = 2500, height = 1250, res = 200)
ggplot(cov_long, aes(x = method, y = coverage, color = method)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "twodash", linewidth = 0.8) +
  ggh4x::facet_grid2(
    option ~ confounding_mechanism,
    scales = "free_x",
    independent = "none"
  ) +
  scale_color_manual(name = "Method", values = method_cols,
                     breaks = names(method_cols)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95, 1)) +
  labs(x = "Confounding mechanism (1–7)", y = "Coverage (nominal 0.95)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
dev.off()

# Main-text version: CMs 1–5, drop trueIV variants
cov_long_sub <- subset(cov_long, confounding_mechanism %in% 1:5 &
                         !method %in% c("trueIV", "trueIV+spatialcoord"))
cov_long_sub <- droplevels(cov_long_sub)

png("images/coverage_Mar27_maintext.png", width = 2000, height = 1300, res = 200)
ggplot(cov_long_sub, aes(x = method, y = coverage, color = method)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "twodash", linewidth = 0.8) +
  ggh4x::facet_grid2(
    option ~ confounding_mechanism,
    scales = "free_x",
    independent = "none"
  ) +
  scale_color_manual(name = "Method", values = method_cols,
                     breaks = names(method_cols)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 0.95, 1)) +
  labs(x = "Confounding mechanism (1–5)", y = "Coverage (nominal 0.95)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
dev.off()

# Plot for the main text: Confounding scenarios 1,3,5,6
df_sub <- subset(df, confounding_mechanism_reordered %in% 1:5 & 
               method != 'trueIV' & method != 'trueIV+spatialcoord')
df_sub <- droplevels(df_sub)
mutrues_sub <- subset(mutrues, confounding_mechanism_reordered %in% 1:5)
mutrues_sub <- droplevels(mutrues_sub)

################# n_uc SELECTION HISTOGRAMS #################

n_uc_csvs <- list.files(RESULTS_DIR, pattern = '_n_uc\\.csv$')

read_n_uc <- function(fname) {
  tryCatch({
    m <- regmatches(fname,
                    regexec('^conf(\\d+)_([^_]+)_(.+)_n_uc\\.csv$', fname,
                            perl = TRUE))[[1]]
    dat <- read.csv(file.path(RESULTS_DIR, fname))
    data.frame(
      confounding_mechanism = m[2],
      option                = m[3],
      method                = m[4],
      n_uc                  = as.vector(as.matrix(dat))
    )
  }, error = function(e) NULL)
}

df_n_uc <- map_dfr(n_uc_csvs, read_n_uc)

if (nrow(df_n_uc) > 0) {
  df_n_uc <- df_n_uc %>%
    mutate(
      confounding_mechanism = factor(confounding_mechanism,
                                     levels = as.character(1:7)),
      option = factor(option, levels = c("linear", "nonlinear")),
      method = recode(method,
                      "IV-GraphLaplacian-spatialcoord" = "IV-GL+spatialcoord",
                      "IV-GraphLaplacian"              = "IV-GL",
                      "IV-TPS-spatialcoord"            = "IV-TPS+spatialcoord",
                      "trueIV-spatialcoord"            = "trueIV+spatialcoord"),
      method = factor(method, levels = desired_order)
    ) %>%
    filter(!is.na(n_uc))

  n_uc_method_cols <- method_cols[levels(droplevels(df_n_uc$method))]

  png("images/n_uc_hist_Mar27.png", width = 3000, height = 1800, res = 200)
  print(
    ggplot(df_n_uc, aes(x = n_uc, fill = option)) +
      geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
      facet_grid(method ~ confounding_mechanism,
                 scales = "free_y") +
      scale_fill_manual(values = c("linear" = "#1F77B4", "nonlinear" = "#E6550D")) +
      labs(x = expression(n[uc]~"(selected)"),
           y = "Count",
           fill = "Outcome model",
           title = "Selected number of unconfounded basis components") +
      theme_bw() +
      theme(legend.position = "top",
            strip.text = element_text(size = 8))
  )
  dev.off()
}

png("images/boxplot_Sep6_maintext.png", width = 2000, height = 1300, res = 200)
ggplot(df_sub, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.5, outliers = F, staplewidth = 1) + #, draw_quantiles = c(0.5)) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "blue") +
  ggh4x::facet_grid2(
    option ~ confounding_mechanism_reordered,
    scales = "free",           # allows different scales per row/col
    independent = "all"        # allows different scales **per panel**
  ) +
  geom_hline(data = mutrues_sub, aes(yintercept = theta), 
             color = "red", linetype = "twodash", size = 1) +
  labs(
    x = "Confounding mechanism (1–5)",
    y = "Truncated Exposure Effect Estimate"
  ) +
  # Remove x-axis title
  scale_fill_manual(
    name = "Method",
    values = method_cols,
    breaks = names(method_cols) # matches your desired order
  ) +
  theme_bw() +   
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") 
dev.off()
