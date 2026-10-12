# Full-sample spatial-basis selection for the PM2.5 analysis
#
# For each cutoff, adjustment set, and spatial basis, apply the sequential
# stability rule in ../R/basis_selection.R (shared with ../simulation):
#   1. reference: the smallest B^c explaining at least r2_reference of
#      Var(PM2.5);
#   2. move one basis function at a time from B^c to B^uc;
#   3. stop at the first Bonferroni-significant departure from the reference
#      and keep the previous dimension (n_ac = 1 if nothing departs).
# The unadjusted baseline and the oracle are fit separately for comparison.

library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

rm(list = ls())
gc()
HIDDEN_PATH <- ""
base_path <- file.path(
  HIDDEN_PATH,
  "projects",
  "wildfirePM-mortality-spatialconf_IV",
  "spatial-scale-conf",
  "iv-spatialconfounding-sophi890-issue54"
)

# Run from october_2026/data_application/; shared modules live in ../R/.
# base_path is only used for lab-server data and large outputs.
source(file.path("..", "R", "basis_construction.R"))
source(file.path("..", "R", "truncated_effect.R"))
source(file.path("..", "R", "basis_selection.R"))

results_dir <- file.path(base_path, "results_oct2026")
images_dir <- file.path(base_path, "images_oct2026")
output_dir <- file.path(results_dir, "full_sample_basis_selection")
# Same folder 0_preprocessing.R writes the preprocessed bundle to.
data_dir <- results_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(images_dir, recursive = TRUE, showWarnings = FALSE)

expected_basis_run_tag <- "k50_allzip_knn"
preprocessed_path <- file.path(
  data_dir,
  paste0(
    "preprocessed_analysis_and_bases_",
    expected_basis_run_tag,
    ".rds"
  )
)

# -----------------------------------------------------------------------------
# Analysis settings
# -----------------------------------------------------------------------------

cutoffs <- 6:12
n_grid <- 25L
density_trim <- 0.01
constrain <- TRUE

# The selection rule (same settings as ../simulation).
r2_reference <- 0.7
family_alpha <- 0.05
selection_seed <- 12001L

sl_library <- c(
  "SL.gam",
  "SL.glm",
  "SL.mean",
  "SL.glm.interaction"
)

# -----------------------------------------------------------------------------
# Load the analysis data and bases
# -----------------------------------------------------------------------------

prep <- readRDS(preprocessed_path)

y <- as.numeric(prep$y)
a <- as.numeric(prep$a)
xmat <- as.matrix(prep$xmat)
oracle_covariates <- as.matrix(prep$covs)
B_tps <- as.matrix(prep$B_tps)
B_gl <- as.matrix(prep$B_gl)

validate_ordered_basis(B_tps, prep$roughness_tps)
validate_ordered_basis(B_gl, prep$roughness_gl)

n <- length(y)
stopifnot(
  length(a) == n,
  nrow(xmat) == n,
  nrow(oracle_covariates) == n,
  nrow(B_tps) == n,
  nrow(B_gl) == n
)

# n_ac counts the constant basis column, so n_ac = 1 is no spatial adjustment.
basis_specifications <- list(
  "IV-TPS" = list(
    B = B_tps,
    label = "Thin plate spline"
  ),
  "IV-GraphLaplacian" = list(
    B = B_gl,
    label = "Graph Laplacian"
  )
)

adjustment_matrices <- list(
  # without_covariates = matrix(
  #   1,
  #   nrow = n,
  #   ncol = 1L,
  #   dimnames = list(NULL, "Intercept")
  # ),
  with_covariates = xmat
)

exposure_sd <- stats::sd(a)
bandwidth_candidates <- seq(
  exposure_sd / 10,
  exposure_sd,
  length.out = 100L # same grid as ../simulation (R/truncated_effect.R default)
)

effect_args <- list(
  nuisance_args = list(sl_library = sl_library),
  n_grid = n_grid,
  bw_seq = bandwidth_candidates,
  density_trim = density_trim,
  constrain = constrain
)

# -----------------------------------------------------------------------------
# Fit baselines, oracle comparisons, and spatial candidates
# -----------------------------------------------------------------------------

baseline_fits <- list()
oracle_fits <- list()
selection_objects <- list()
diagnostic_tables <- list()
selected_tables <- list()

for (cutoff in cutoffs) {
  message("Cutoff = ", cutoff)
  
  # Fit oracle
  set.seed(selection_seed)
  oracle_fit <- do.call(
    estimate_truncated_effect_full,
    c(
      list(
        y = y,
        a = a,
        w = oracle_covariates,
        cutoff = cutoff
      ),
      effect_args
    )
  )
  oracle_fits[[as.character(cutoff)]] <- oracle_fit
  
  for (adjustment_name in names(adjustment_matrices)) {
    x <- adjustment_matrices[[adjustment_name]]
    baseline_key <- paste(cutoff, adjustment_name, sep = "|")
    
    # Fit baseline
    set.seed(selection_seed)
    baseline_fit <- do.call(
      estimate_truncated_effect_full,
      c(
        list(y = y, a = a, w = x, cutoff = cutoff),
        effect_args
      )
    )
    baseline_fits[[baseline_key]] <- baseline_fit
    
    for (basis_name in names(basis_specifications)) {
      specification <- basis_specifications[[basis_name]]
      
      message("Selecting ", basis_name, " with ", adjustment_name)
      
      selection <- select_basis_dimension(
        y = y,
        a = a,
        x = x,
        B = specification$B,
        cutoff = cutoff,
        effect_estimator = estimate_truncated_effect_full,
        effect_args = effect_args,
        r2_reference = r2_reference,
        family_alpha = family_alpha,
        seed = selection_seed
      )
      
      selection_key <- paste(
        cutoff,
        adjustment_name,
        basis_name,
        sep = "|"
      )
      selection_objects[[selection_key]] <- selection
      
      diagnostics <- selection$diagnostics %>%
        mutate(
          cutoff = cutoff,
          adjustment = adjustment_name,
          basis = basis_name,
          .before = 1L
        )
      diagnostic_tables[[selection_key]] <- diagnostics
      
      selected_tables[[selection_key]] <- diagnostics %>%
        filter(selected) %>%
        transmute(
          cutoff,
          adjustment,
          basis,
          selected_n_ac = n_ac,
          selected_n_uc = n_uc,
          exposure_R2,
          psi,
          ci_lower,
          ci_upper,
          standard_error
        )
      
      message(
        "  selected q = ", selection$selected_n_ac,
        "; n_uc = ", selection$selected_n_uc
      )
    }
    gc()
  }
}

selection_diagnostics <- bind_rows(diagnostic_tables)
selected_results <- bind_rows(selected_tables)

oracle_results <- bind_rows(lapply(names(oracle_fits), function(cutoff) {
  fit <- oracle_fits[[cutoff]]
  data.frame(
    cutoff = as.numeric(cutoff),
    psi = fit$psi,
    ci_lower = fit$confidence_interval[1L],
    ci_upper = fit$confidence_interval[2L],
    standard_error = fit$standard_error
  )
}))

write.csv(
  selection_diagnostics,
  file.path(output_dir, "selection_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  selected_results,
  file.path(output_dir, "selected_results.csv"),
  row.names = FALSE
)
write.csv(
  oracle_results,
  file.path(output_dir, "oracle_results.csv"),
  row.names = FALSE
)

# Unadjusted baseline (measured covariates only, no spatial adjustment).
baseline_results <- bind_rows(lapply(names(baseline_fits), function(key) {
  fit <- baseline_fits[[key]]
  parts <- strsplit(key, "|", fixed = TRUE)[[1L]]
  data.frame(
    cutoff = as.numeric(parts[1L]),
    adjustment = parts[2L],
    psi = fit$psi,
    ci_lower = fit$confidence_interval[1L],
    ci_upper = fit$confidence_interval[2L],
    standard_error = fit$standard_error
  )
}))
write.csv(
  baseline_results,
  file.path(output_dir, "baseline_results.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    selections = selection_objects,
    baselines = baseline_fits,
    oracles = oracle_fits,
    settings = list(
      cutoffs = cutoffs,
      r2_reference = r2_reference,
      family_alpha = family_alpha,
      density_trim = density_trim,
      constrain = constrain,
      n_grid = n_grid
    )
  ),
  file.path(output_dir, "full_sample_selection_fits.rds")
)

selected_results <- read.csv(file.path(output_dir, 'selected_results.csv'))
print(selected_results)

# -----------------------------------------------------------------------------
# Plot estimates and selected dimensions
# -----------------------------------------------------------------------------
selection_diagnostics <- read.csv(file.path(output_dir, "selection_diagnostics.csv"))
oracle_results <- read.csv(file.path(output_dir, "oracle_results.csv"))
baseline_results <- read.csv(file.path(output_dir, "baseline_results.csv"))
basis_labels <- vapply(
  basis_specifications,
  function(specification) specification$label,
  character(1)
)

adjustment_labels <- c(
  #without_covariates = "Without measured covariates",
  with_covariates = "With measured covariates"
)

for (cutoff in cutoffs) {
  plot_data <- selection_diagnostics %>%
    filter(.data$cutoff == .env$cutoff) %>%
    mutate(
      basis_label = factor(
        basis_labels[basis],
        levels = unname(basis_labels)
      ),
      adjustment_label = factor(
        adjustment_labels[adjustment],
        levels = unname(adjustment_labels)
      ),
      percent_change = 100 * (psi - 1),
      percent_lower = 100 * (ci_lower - 1),
      percent_upper = 100 * (ci_upper - 1)
    )
  
  selected_lines <- plot_data %>%
    filter(selected) %>%
    select(basis_label, adjustment_label, selected_n_ac = n_ac)
  
  baseline_lines <- plot_data %>%
    distinct(basis_label, adjustment, adjustment_label) %>%
    inner_join(
      baseline_results %>% filter(.data$cutoff == .env$cutoff),
      by = "adjustment"
    ) %>%
    transmute(basis_label, adjustment_label, baseline = 100 * (psi - 1))
  
  oracle_row <- oracle_results %>%
    filter(.data$cutoff == .env$cutoff)
  
  oracle_reference <- plot_data %>%
    distinct(basis_label, adjustment_label) %>%
    mutate(
      oracle = 100 * (oracle_row$psi[1] - 1),
      oracle_lower = 100 * (oracle_row$ci_lower[1] - 1),
      oracle_upper = 100 * (oracle_row$ci_upper[1] - 1)
    )
  
  selection_plot <- ggplot(
    plot_data,
    aes(x = n_ac, y = percent_change)
  ) +
    geom_hline(yintercept = 0, color = "gray75") +
    geom_rect(
      data = oracle_reference,
      aes(ymin = oracle_lower, ymax = oracle_upper),
      xmin = -Inf,
      xmax = Inf,
      inherit.aes = FALSE,
      fill = "forestgreen",
      alpha = 0.10
    ) +
    geom_hline(
      data = oracle_reference,
      aes(yintercept = oracle),
      inherit.aes = FALSE,
      color = "forestgreen",
      linewidth = 0.8
    ) +
    geom_hline(
      data = baseline_lines,
      aes(yintercept = baseline),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "gray35"
    ) +
    geom_line(color = "#2C7FB8", linewidth = 0.8) +
    geom_errorbar(
      aes(ymin = percent_lower, ymax = percent_upper),
      width = 0.15,
      color = "gray40"
    ) +
    geom_point(
      aes(fill = selected),
      shape = 21,
      size = 3,
      color = "#2C7FB8"
    ) +
    geom_vline(
      data = selected_lines,
      aes(xintercept = selected_n_ac),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "#D7301F",
      linewidth = 0.9
    ) +
    scale_fill_manual(
      values = c(`FALSE` = "white", `TRUE` = "#D7301F"),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = sort(unique(plot_data$n_ac))
    ) +
    facet_grid(adjustment_label ~ basis_label) +
    labs(
      x = expression("Dimension of spatial basis of confounding adjustment, " * n[c]),
      y = "Reduction in mortality rate (%)",
      title = bquote(
        "Sensitivity to spatial basis dimension at cutoff " *
          .(cutoff) ~ mu * plain(g) / plain(m)^3
      )#,
      # caption = paste(
      #   "Error bars are pointwise 95% confidence intervals.",
      #   "The green line and band show the oracle estimate and its 95% confidence interval.",
      #   "The dotted line shows the baseline estimate.",
      #   "The dashed vertical line marks the selected dimension: the last one",
      #   "before the first Bonferroni-adjusted departure from the reference."
      # )
    ) +
    theme_classic(base_size = 15) +
    theme(
      strip.background = element_blank(),
      legend.position = "none",
      plot.caption = element_text(hjust = 0)
    ) 
  
  print(selection_plot)
  
  ggsave(
    file.path(
      images_dir,
      paste0("full_sample_basis_selection_cutoff", cutoff, ".png")
    ),
    plot = selection_plot,
    width = 11,
    height = 7,
    dpi = 300
  )
}

# -----------------------------------------------------------------------------
# Plot selected estimates across cutoffs
# -----------------------------------------------------------------------------

baseline_plot <- baseline_results %>%
  transmute(
    cutoff,
    adjustment,
    method = "Baseline",
    psi,
    ci_lower,
    ci_upper
  )

selected_plot <- selected_results %>%
  transmute(
    cutoff,
    adjustment,
    method = basis,
    psi,
    ci_lower,
    ci_upper
  )

# Repeat each oracle estimate for each adjustment specification.
oracle_plot <- bind_rows(lapply(
  unique(selected_results$adjustment),
  function(adjustment_value) {
    oracle_results %>%
      transmute(
        cutoff,
        adjustment = adjustment_value,
        method = "Oracle",
        psi,
        ci_lower,
        ci_upper
      )
  }
))

# SMW 9/11 adds spatial coordinates results (previously run)
previous_results <- read.csv(file.path(base_path, "results_sep2026/truncated_effect_estimates_full_nuisance_k50_allzip_knn_K5_median_nuc41to48_alpha1_trim001_grid25_constrained_LL.csv"))
spatialcoord <- previous_results[previous_results$method == 'spatialcoord',]
spatialcoord_plot <- spatialcoord %>%
  filter(
    cutoff %in% cutoffs,
    is.finite(psi),
    is.finite(ci_lower),
    is.finite(ci_upper)
  ) %>%
  transmute(
    cutoff,
    adjustment = "with_covariates",
    method = "Spatial coordinates",
    psi,
    ci_lower,
    ci_upper
  )

across_cutoff_data <- bind_rows(
  oracle_plot,
  baseline_plot,
  spatialcoord_plot,
  selected_plot
) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Oracle",
        "Baseline",
        "Spatial coordinates",
        "IV-TPS",
        "IV-GraphLaplacian"
      ),
      labels = c(
        "Oracle",
        "Baseline",
        "Spatial coordinates",
        "IV-TPS",
        "IV-GL"
      ),
    ),
    adjustment_label = factor(
      adjustment,
      levels = c(
        "without_covariates",
        "with_covariates"
      ),
      labels = c(
        "Without measured covariates",
        "With measured covariates"
      )
    ),
    percent_change = 100 * (psi - 1),
    percent_lower = 100 * (ci_lower - 1),
    percent_upper = 100 * (ci_upper - 1)
  )

selected_estimate_plot <- ggplot(
  across_cutoff_data,
  aes(
    x = method,
    y = percent_change,
    ymin = percent_lower,
    ymax = percent_upper,
    color = method
  )
) +
  geom_hline(
    yintercept = 0,
    color = "grey75",
    linewidth = 0.4
  ) +
  geom_errorbar(
    width = 0.15,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  facet_grid(
    adjustment_label ~ cutoff,
    labeller = label_bquote(
      cols = .(cutoff) ~ mu * plain(g) / plain(m)^3
    )
  ) +
  scale_color_manual(
    values = c(
      "Oracle" = "forestgreen",
      "Baseline" = "red3",
      "Spatial coordinates" = "black",
      "IV-TPS" = "black",
      "IV-GL" = "black"
    ),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = "Reduction in mortality rate (%)",
    title = expression(
      "Impact of enforcing average " * plain(PM)[2.5] *
        " (2001-2010) below a cutoff on mortality rate in the Medicare population"
    )
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

print(selected_estimate_plot)

ggsave(
  file.path(
    images_dir,
    "selected_estimates_across_cutoffs.png"
  ),
  plot = selected_estimate_plot,
  width = 15,
  height = 6,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Hausdorff distance between each method CI and the oracle CI
# -----------------------------------------------------------------------------

oracle_intervals <- across_cutoff_data %>%
  filter(method == "Oracle") %>%
  transmute(
    cutoff,
    adjustment,
    oracle_lower = ci_lower,
    oracle_upper = ci_upper
  )

hausdorff_results <- across_cutoff_data %>%
  filter(method != "Oracle") %>%
  left_join(
    oracle_intervals,
    by = c("cutoff", "adjustment")
  ) %>%
  mutate(
    hausdorff = pmax(
      abs(ci_lower - oracle_lower),
      abs(ci_upper - oracle_upper)
    )
  ) %>%
  arrange(adjustment, method, cutoff)

# Distances at each cutoff
print(
  hausdorff_results %>%
    select(
      adjustment,
      method,
      cutoff,
      hausdorff
    ) %>%
    mutate(hausdorff = round(hausdorff, 3))
)

# Average distance across cutoffs
hausdorff_summary <- hausdorff_results %>%
  group_by(adjustment, method) %>%
  summarise(
    mean_hausdorff = mean(hausdorff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(adjustment, mean_hausdorff)

print(
  hausdorff_summary %>%
    mutate(mean_hausdorff = round(mean_hausdorff, 3))
)

write.csv(
  hausdorff_results,
  file.path(output_dir, "hausdorff_distance_from_oracle.csv"),
  row.names = FALSE
)

write.csv(
  hausdorff_summary,
  file.path(output_dir, "mean_hausdorff_distance_from_oracle.csv"),
  row.names = FALSE
)
