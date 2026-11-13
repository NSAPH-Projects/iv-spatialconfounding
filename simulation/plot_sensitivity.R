# Plot the sensitivity trajectories saved by sensitivity.R.
#
# Usage (from simulation/):
#   Rscript plot_sensitivity.R [rds[,rds2,...]] [png] [scenarios]
#   scenarios: comma-separated subset to plot (default: the supplement set
#              Mechanisms 1-5, LS,MS,twoU,coordnoise,stateJump)
#   defaults: images/sensitivity.rds -> .png
# Rows are mechanisms, columns are method x outcome model (IV-TPS linear,
# IV-GL linear, IV-TPS nonlinear, IV-GL nonlinear); each panel has its own x
# axis. Each trajectory shows n_uc^0, where the walk starts (open triangle),
# and the selected n_uc (filled diamond); the dashed red line is the dimension
# capturing 95% of the variance of A_c's spatial confounded part U, median
# over the replicates; the psi plot adds the true psi (dotted grey). The x
# axis is n_uc = m - n_ac, so the walk runs left to right.
# Also writes <png>_correlation: corr(estimated A_c(n_ac), U) along the same
# trajectories, with the same layout and legend.
suppressMessages({
  library(ggplot2)
})

# sensitivity.R sets rds_path / png_path before sourcing this file.
if (!exists("rds_path")) {
  args <- commandArgs(trailingOnly = TRUE)
  rds_path <- if (length(args) >= 1L) args[[1]] else "images/sensitivity.rds"
  png_path <- if (length(args) >= 2L) args[[2]] else sub("\\.rds$", ".png", strsplit(rds_path, ",")[[1]][1])
}

# rds_path may list several files (comma-separated). A scenario in a later
# file replaces the same scenario in earlier files, so a rerun of updated
# scenarios can be combined with an older full run.
saved <- NULL
for (path in strsplit(rds_path, ",")[[1]]) {
  new <- readRDS(path)
  if (is.null(saved)) { saved <- new; next }
  replaced <- unique(new$trajectories$scenario)
  saved$trajectories <- rbind(saved$trajectories[!saved$trajectories$scenario %in% replaced, ],
                              new$trajectories[, names(saved$trajectories)])
  saved$truth <- rbind(saved$truth[!saved$truth$scenario %in% replaced, ], new$truth)
}
if (!exists("plot_scenarios")) {
  args <- commandArgs(trailingOnly = TRUE)
  plot_scenarios <- if (length(args) >= 3L) strsplit(args[[3]], ",")[[1]] else
    c("LS", "MS", "twoU", "coordnoise", "stateJump")
}
res <- saved$trajectories[saved$trajectories$scenario %in% plot_scenarios, ]
truth <- saved$truth[saved$truth$scenario %in% plot_scenarios, ]

# Supplement numbering (same as analysis.R); rows follow this order.
scenario_labels <- c(LS = "Mechanism 1", MS = "Mechanism 2", twoU = "Mechanism 3",
                     coordnoise = "Mechanism 4", stateJump = "Mechanism 5")
option_labels <- c(linear = "linear", tanh = "nonlinear")
bases <- intersect(c("TPS", "GL"), unique(res$basis))
options_ <- intersect(names(option_labels), unique(res$option))
scenarios_ <- intersect(names(scenario_labels), unique(res$scenario))

label_rows <- function(df) {
  df$scenario <- factor(scenario_labels[df$scenario], levels = scenario_labels[scenarios_])
  df$column <- factor(paste0("IV-", df$basis, ", ", option_labels[df$option]),
                      levels = c(outer(paste0("IV-", bases), option_labels[options_], paste, sep = ", ")))
  df
}

res$rep <- factor(res$rep)
res$n_uc <- res$m - res$n_ac
# Median over replicates of the dimension capturing 95% of Var(A^c).
true_dim <- unique(res[, c("scenario", "option", "basis", "rep", "m", "n_ac_95")])
true_dim$n_uc <- true_dim$m - true_dim$n_ac_95
true_dim <- label_rows(aggregate(n_uc ~ scenario + option + basis, true_dim, stats::median))
truth_lines <- label_rows(merge(truth[truth$scenario %in% scenarios_ & truth$option %in% options_, ],
                                data.frame(basis = bases)))
cor_res <- res  # unlabelled copy for the correlation plot
res <- label_rows(res)
selected_pts <- subset(res, n_ac == selected)
reference_pts <- subset(res, n_ac == reference)

# ---- Shared legend (both plots) ----
# Markers: open triangle = n_uc^0 (reference, where the walk starts); filled
# diamond = selected n_uc. Lines: dashed red = dimension capturing 95% of
# Var(A^c); the psi plot adds the true psi (dotted grey). Each legend entry is
# drawn by one layer only (show.legend), so keys do not overlap.
marker_scale <- scale_shape_manual(
  name = NULL, values = c(reference = 2, selected = 18),
  labels = c(reference = expression(n[uc]^0), selected = expression(paste("Selected ", n[uc]))))
marker_guide <- guide_legend(order = 2, override.aes = list(colour = "black", size = c(4.5, 6.5)))
line_labels <- c(true_dim = expression(paste("Dimension capturing 95% of Var(", U, ") (median over replicates)")),
                 truth = expression(paste("True ", psi)))
line_scale <- function(keys) scale_linetype_manual(
  name = NULL, values = c(true_dim = "dashed", truth = "dotted")[keys], labels = line_labels[keys])
line_guide <- function(keys) guide_legend(order = 3, override.aes = list(colour = c(true_dim = "red", truth = "grey40")[keys]))
legend_theme <- theme(legend.position = "top", legend.key.width = unit(2, "lines"))

p <- ggplot(res, aes(n_uc, psi, colour = rep, group = rep)) +
  geom_hline(data = truth_lines, aes(yintercept = psi_true, linetype = "truth"), colour = "grey40",
             key_glyph = "path", inherit.aes = FALSE) +
  geom_vline(data = true_dim, aes(xintercept = n_uc, linetype = "true_dim"), colour = "red",
             key_glyph = "path", inherit.aes = FALSE) +
  geom_line(alpha = 0.8) +
  geom_point(size = 0.4) +
  geom_point(data = reference_pts, aes(shape = "reference"), size = 4, stroke = 1.2,
             show.legend = c(shape = TRUE, colour = FALSE)) +
  geom_point(data = selected_pts, aes(shape = "selected"), size = 6, show.legend = FALSE) +
  ggh4x::facet_grid2(scenario ~ column, scales = "free", independent = "x") +
  scale_colour_brewer(palette = "Set1", name = "Replicate") +
  marker_scale + line_scale(c("true_dim", "truth")) +
  guides(colour = guide_legend(order = 1), shape = marker_guide, linetype = line_guide(c("true_dim", "truth"))) +
  labs(x = expression(paste("Number of basis functions ", n[uc],
                            " approximating the instrumental variation")),
       y = "Truncated exposure effect estimate") +
  theme_minimal(base_size = 22) +
  theme(strip.text = element_text(face = "bold"), panel.spacing.x = unit(3, "lines")) +
  legend_theme

ggsave(png_path, p, width = 4 + 4 * length(levels(res$column)), height = max(10, 3.2 * length(scenarios_) + 3), dpi = 200,
       bg = "white", limitsize = FALSE)
message("Wrote ", png_path)

# ---- Correlation of the estimated A^c with the true A^c along the trajectories ----
# Same layout as the psi plot (rows: mechanisms; columns: method x outcome).
# The projection depends on A and the basis only, so the curve is the same for
# both outcome models; the selected dimension (diamond) can differ.
# The estimated A_c also contains projected X and noise, so the correlation
# stays below 1 even at the right dimension.
if ("cor_ac" %in% names(cor_res)) {
  cor_res$n_uc <- cor_res$m - cor_res$n_ac
  curves <- label_rows(cor_res)
  sel_cor <- subset(curves, n_ac == selected)
  ref_cor <- subset(curves, n_ac == reference)
  p_cor <- ggplot(curves, aes(n_uc, cor_ac, colour = rep, group = rep)) +
    geom_vline(data = true_dim, aes(xintercept = n_uc, linetype = "true_dim"), colour = "red",
               key_glyph = "path", inherit.aes = FALSE) +
    geom_line(alpha = 0.8) +
    geom_point(data = ref_cor, aes(shape = "reference"), size = 2.5, stroke = 0.9,
               show.legend = c(shape = TRUE, colour = FALSE)) +
    geom_point(data = sel_cor, aes(shape = "selected"), size = 4, show.legend = FALSE) +
    ggh4x::facet_grid2(scenario ~ column, scales = "free", independent = "x") +
    scale_colour_brewer(palette = "Set1", name = "Replicate") +
    marker_scale + line_scale("true_dim") +
    guides(colour = guide_legend(order = 1), shape = marker_guide, linetype = line_guide("true_dim")) +
    labs(x = expression(paste("Number of basis functions ", n[uc],
                              " approximating the instrumental variation")),
         y = expression(paste("Cor(", widehat(A[c]), ", ", U, ")"))) +
    theme_minimal(base_size = 22) +
    theme(strip.text = element_text(face = "bold"), panel.spacing.x = unit(3, "lines")) +
    legend_theme

  cor_path <- sub("\\.png$", "_correlation.png", png_path)
  ggsave(cor_path, p_cor, width = 4 + 4 * length(levels(curves$column)),
         height = max(10, 3.2 * length(scenarios_) + 3), dpi = 200, bg = "white", limitsize = FALSE)
  message("Wrote ", cor_path)
}
