library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

rm(list = ls())
gc()
# Run from data_application/, after 0_preprocessing.R.
base_path <- ".." # Path to this repository

source(file.path(base_path, "R", "basis_construction.R"))
source(file.path(base_path, "R", "truncated_effect.R"))
source(file.path(base_path, "R", "basis_selection.R"))

results_dir <- file.path(base_path, "data_application", "results")
images_dir <- file.path(base_path, "data_application", "images")
output_dir <- file.path(results_dir, "full_sample_basis_selection")
data_dir <- file.path(base_path, "data_application", "preprocessed")

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
cutoffs <- 6:12
n_grid <- 25L
density_trim <- 0.01
constrain <- TRUE

# The selection rule.
r2_reference <- 0.5 # 0.25 worked well, 0.7 IV-TPS too noisy at heavy spatial confounding adjustments
family_alpha <- 0.05
selection_seed <- 12001L

sl_library <- c(
  "SL.gam",
  "SL.glm",
  "SL.mean",
  "SL.glm.interaction",
  "SL.earth"
)

prep <- readRDS(preprocessed_path)

y <- as.numeric(prep$y)
a <- as.numeric(prep$a)
xmat <- as.matrix(prep$xmat)
oracle_covariates <- as.matrix(prep$covs)
B_tps <- as.matrix(prep$B_tps)
B_gl <- as.matrix(prep$B_gl)

# Coordinates rescaled to [0, 1] (0_preprocessing.R), for the spatial
# coordinates method.
spatial_coordinates <- cbind(lon = as.numeric(prep$lon_n), lat = as.numeric(prep$lat_n))

validate_ordered_basis(B_tps, prep$roughness_tps)
validate_ordered_basis(B_gl, prep$roughness_gl)

n <- length(y)

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
  without_covariates = matrix(
    1,
    nrow = n,
    ncol = 1L,
    dimnames = list(NULL, "Intercept")
  ),
  with_covariates = xmat
)

exposure_sd <- stats::sd(a)
bandwidth_candidates <- seq(
  exposure_sd / 10,
  exposure_sd,
  length.out = 100L
)

effect_args <- list(
  nuisance_args = list(sl_library = sl_library),
  n_grid = n_grid,
  bw_seq = bandwidth_candidates,
  density_trim = density_trim,
  constrain = constrain
)

baseline_fits <- list()
spatialcoord_fits <- list()
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
    
    # Fit spatial coordinates: the adjustment covariates plus the coordinates
    # rescaled to [0, 1], as in the simulation (method_covariates()). The
    # constant intercept column is dropped, so without covariates the
    # adjustment set is the coordinates alone.
    set.seed(selection_seed)
    x_nonconstant <- x[, apply(x, 2L, stats::var) > 0, drop = FALSE]
    spatialcoord_fits[[baseline_key]] <- do.call(
      estimate_truncated_effect_full,
      c(
        list(y = y, a = a, w = cbind(x_nonconstant, spatial_coordinates), cutoff = cutoff),
        effect_args
      )
    )
    
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
        seed = selection_seed,
        full_trajectory = T
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
baseline_results <- do.call(rbind, lapply(names(baseline_fits), function(k) {
  f <- baseline_fits[[k]]
  p <- strsplit(k, "|", fixed = TRUE)[[1]]
  data.frame(cutoff = as.numeric(p[1]), adjustment = p[2], psi = f$psi,
             ci_lower = f$confidence_interval[1], ci_upper = f$confidence_interval[2])
}))

write.csv(baseline_results, file.path(output_dir, "baseline_results.csv"), row.names = FALSE)
spatialcoord_results <- do.call(rbind, lapply(names(spatialcoord_fits), function(k) {
  f <- spatialcoord_fits[[k]]
  p <- strsplit(k, "|", fixed = TRUE)[[1]]
  data.frame(cutoff = as.numeric(p[1]), adjustment = p[2], psi = f$psi,
             ci_lower = f$confidence_interval[1], ci_upper = f$confidence_interval[2])
}))
write.csv(spatialcoord_results, file.path(output_dir, "spatialcoord_results.csv"), row.names = FALSE)

saveRDS(
  list(
    selections = selection_objects,
    baselines = baseline_fits,
    spatialcoords = spatialcoord_fits,
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

selection_diagnostics <- read.csv(file.path(output_dir, "selection_diagnostics.csv"))
oracle_results <- read.csv(file.path(output_dir, "oracle_results.csv"))
baseline_results <- read.csv(file.path(output_dir, "baseline_results.csv"))
basis_labels <- vapply(
  basis_specifications,
  function(specification) specification$label,
  character(1)
)

adjustment_labels <- c(
  without_covariates = "Without measured covariates",
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
    select(basis_label, adjustment_label, selected_n_uc = n_uc)
  
  baseline_lines <- plot_data %>%
    distinct(basis_label, adjustment, adjustment_label) %>%
    inner_join(filter(baseline_results, .data$cutoff == .env$cutoff), by = "adjustment") %>%
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
    aes(x = n_uc, y = percent_change)
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
      aes(xintercept = selected_n_uc),
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
      breaks = sort(unique(plot_data$n_uc))
    ) +
    facet_grid(adjustment_label ~ basis_label) +
    labs(
      x = expression("Number of spatial basis functions used to approximate the instrumental variation, " * n[uc]),
      y = "Reduction in mortality rate (%)",
      title = bquote(
        "Cutoff " *
          .(cutoff) ~ mu * plain(g) / plain(m)^3
      )#,
      # caption = paste(
      #   "Error bars are pointwise 95% confidence intervals.",
      #   "The green line and band show the oracle estimate and its 95% confidence interval.",
      #   "The dotted line shows the baseline estimate.",
      #   "The dashed vertical line marks the first Bonferroni-adjusted",
      #   "paired departure from baseline."
      # )
    ) +
    theme_classic(base_size = 15) +
    theme(
      strip.background = element_blank(),
      legend.position = "none",
      plot.caption = element_text(hjust = 0),
      plot.title = element_text(hjust = 0.5)
    ) +
    # coord_cartesian zooms without dropping points (ylim() would remove
    # out-of-range estimates and break the lines and error bars).
    coord_cartesian(ylim = c(-8.25, 0.5))
  
  print(selection_plot)
  
  ggsave(
    file.path(
      images_dir,
      paste0("full_sample_basis_selection_cutoff", cutoff, ".png")
    ),
    plot = selection_plot,
    width = 11,
    height = 11,
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

# Spatial coordinates, fit above with the same estimator settings.
spatialcoord_results <- read.csv(file.path(output_dir, "spatialcoord_results.csv"))
spatialcoord_plot <- spatialcoord_results %>%
  filter(
    is.finite(psi),
    is.finite(ci_lower),
    is.finite(ci_upper)
  ) %>%
  transmute(
    cutoff,
    adjustment,
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
    percent_change = 100 * (psi - 1),
    percent_lower = 100 * (ci_lower - 1),
    percent_upper = 100 * (ci_upper - 1)
  )

# One plot per adjustment specification.
for (adjustment_name in names(adjustment_labels)) {
  selected_estimate_plot <- ggplot(
    filter(across_cutoff_data, .data$adjustment == .env$adjustment_name),
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
      ~ cutoff,
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
      subtitle = adjustment_labels[[adjustment_name]],
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
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  print(selected_estimate_plot)
  ggsave(
    file.path(
      images_dir,
      paste0("selected_estimates_across_cutoffs_", adjustment_name, ".png")
    ),
    plot = selected_estimate_plot,
    width = 15,
    height = 6,
    dpi = 300
  )
}

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
