# Combined sensitivity / selection figure.
#
# Layout as plot_sensitivity.R: rows are Mechanisms 1-5, columns are method x
# outcome model (IV-TPS linear, IV-GL linear, IV-TPS nonlinear, IV-GL
# nonlinear), each panel with its own x axis (n_uc = m - n_ac). Each panel
# shows
#   - the trajectories of the estimate in the five sensitivity replicates,
#     with the selected n_uc (filled diamond);
#   - the true psi (dotted grey) and the mean baseline estimate over all
#     simulated datasets (dashed red), as horizontal lines;
#   - along its bottom, a histogram of the selected n_uc over all simulated
#     datasets (from results_dir), scaled to a band below the trajectories.
#
# Usage (from simulation/):
#   Rscript plot_combined.R [rds[,rds2,...]] [results_dir] [png] [scenarios]
#   defaults: images/sensitivity.rds, results_final,
#             images/sensitivity_selection.png, LS,MS,twoU,coordnoise,stateJump
#   results_dir needs summary.csv (analysis.R) and the IV result files.

suppressMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
rds_path <- if (length(args) >= 1L) args[[1]] else "images/sensitivity.rds"
results_dir <- if (length(args) >= 2L) args[[2]] else "results_final"
png_path <- if (length(args) >= 3L) args[[3]] else "images/sensitivity_selection.png"
plot_scenarios <- if (length(args) >= 4L) strsplit(args[[4]], ",")[[1]] else
  c("LS", "MS", "twoU", "coordnoise", "stateJump")

scenario_labels <- c(LS = "Mechanism 1", MS = "Mechanism 2", twoU = "Mechanism 3",
                     coordnoise = "Mechanism 4", stateJump = "Mechanism 5")
option_labels <- c(linear = "linear", tanh = "nonlinear")
method_cols <- c("IV-TPS" = "#9ECAE1", "IV-GL" = "#FDAE6B")

# ---- Trajectories (later files replace earlier ones scenario by scenario) ----
traj <- NULL
truth <- NULL
for (path in strsplit(rds_path, ",")[[1]]) {
  new <- readRDS(path)
  if (!is.null(traj)) {
    replaced <- unique(new$trajectories$scenario)
    traj <- traj[!traj$scenario %in% replaced, ]
    truth <- truth[!truth$scenario %in% replaced, ]
  }
  traj <- rbind(traj, new$trajectories[, c("scenario", "option", "basis", "rep", "m", "n_ac", "psi", "selected")])
  truth <- rbind(truth, new$truth)
}
traj <- traj[traj$scenario %in% plot_scenarios, ]
traj$method <- paste0("IV-", traj$basis)
traj$n_uc <- traj$m - traj$n_ac

# ---- Selections over all datasets, and the mean baseline estimate ----
files <- list.files(results_dir, pattern = "^[A-Za-z0-9]+_[a-z]+_iv(tps|gl)(_r[0-9]+)?\\.csv$", full.names = TRUE)
sel <- do.call(rbind, lapply(files, function(f) {
  d <- read.csv(f)
  d[!is.na(d$n_ac), c("scenario", "option", "method", "n_ac")]
}))
sel <- sel[sel$scenario %in% plot_scenarios, ]
m_basis <- 999L  # n - 1 columns
sel$n_uc <- m_basis - sel$n_ac
summ <- read.csv(file.path(results_dir, "summary.csv"))
base <- summ[summ$method == "baseline" & summ$scenario %in% plot_scenarios, ]
base$baseline_mean <- base$psi_true + base$bias
hlines <- rbind(
  data.frame(scenario = truth$scenario, option = truth$option, y = truth$psi_true, line = "truth"),
  data.frame(scenario = base$scenario, option = base$option, y = base$baseline_mean, line = "baseline")
)
hlines <- hlines[hlines$scenario %in% plot_scenarios, ]
hlines <- merge(hlines, data.frame(method = c("IV-TPS", "IV-GL")))

# ---- Histogram band: below the trajectories, per mechanism row ----
# Rows share a y axis; the band spans the bottom `band` share of the row's
# range, separated from the trajectories by a small gap.
band <- 0.22
gap <- 0.03
row_range <- do.call(rbind, lapply(split(rbind(
  data.frame(scenario = traj$scenario, y = traj$psi),
  data.frame(scenario = hlines$scenario, y = hlines$y)), ~scenario), function(d) {
  r <- range(d$y)
  h <- diff(r) / (1 - band - gap)  # total height including the band
  data.frame(scenario = d$scenario[1], top = r[2], data_bottom = r[1],
             band_top = r[1] - gap * h, band_bottom = r[1] - (gap + band) * h)
}))

# Bins per panel over the panel's x range (trajectories and selections).
hist_rects <- do.call(rbind, lapply(split(sel, list(sel$scenario, sel$option, sel$method), drop = TRUE), function(d) {
  tr <- traj[traj$scenario == d$scenario[1] & traj$option == d$option[1] & traj$method == d$method[1], ]
  xr <- range(c(d$n_uc, tr$n_uc))
  width <- max(1, diff(xr) / 40)
  breaks <- seq(xr[1] - width / 2, xr[2] + width, by = width)
  counts <- hist(d$n_uc, breaks = breaks, plot = FALSE)$counts
  rr <- row_range[row_range$scenario == d$scenario[1], ]
  keep <- counts > 0
  data.frame(scenario = d$scenario[1], option = d$option[1], method = d$method[1],
             xmin = head(breaks, -1)[keep], xmax = tail(breaks, -1)[keep],
             ymin = rr$band_bottom,
             ymax = rr$band_bottom + (rr$band_top - rr$band_bottom) * counts[keep] / max(counts))
}))
baselines <- merge(unique(hist_rects[, c("scenario", "option", "method")]), row_range)

# ---- Labels ----
column_levels <- c("IV-TPS, linear", "IV-GL, linear", "IV-TPS, nonlinear", "IV-GL, nonlinear")
relabel <- function(df) {
  df$column <- factor(paste0(df$method, ", ", option_labels[df$option]), levels = column_levels)
  df$scenario <- factor(scenario_labels[df$scenario], levels = scenario_labels[plot_scenarios])
  df
}
traj <- relabel(traj)
traj$rep <- factor(traj$rep)
selected_pts <- subset(traj, n_ac == selected)
hlines <- relabel(hlines)
hist_rects <- relabel(hist_rects)
baselines <- relabel(baselines)

line_labels <- c(truth = expression(paste("True ", psi)), baseline = "Mean baseline estimate")

p <- ggplot() +
  geom_rect(data = hist_rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = method),
            alpha = 0.9, inherit.aes = FALSE) +
  geom_hline(data = baselines, aes(yintercept = band_bottom), colour = "grey60", linewidth = 0.3) +
  geom_hline(data = hlines, aes(yintercept = y, linetype = line, colour = line), linewidth = 0.8,
             key_glyph = "path") +
  geom_line(data = traj, aes(n_uc, psi, group = rep, colour = rep), alpha = 0.85) +
  geom_point(data = selected_pts, aes(n_uc, psi, colour = rep, shape = "selected"), size = 5) +
  ggh4x::facet_grid2(scenario ~ column, scales = "free", independent = "x") +
  scale_colour_manual(
    name = NULL,
    values = c(setNames(RColorBrewer::brewer.pal(5, "Set1"), 1:5), truth = "grey30", baseline = "#D62728"),
    breaks = as.character(1:5), labels = paste("Replicate", 1:5)) +
  scale_linetype_manual(name = NULL, values = c(truth = "dotted", baseline = "dashed"),
                        breaks = c("truth", "baseline"), labels = line_labels) +
  scale_shape_manual(name = NULL, values = c(selected = 18),
                     labels = expression(paste("Selected ", n[uc]))) +
  scale_fill_manual(name = expression(paste("Selected ", n[uc], " across datasets:")), values = method_cols,
                    breaks = c("IV-TPS", "IV-GL")) +
  # Replicate colours only distinguish the five trajectories, so they get no legend.
  guides(colour = "none",
         shape = guide_legend(order = 2, override.aes = list(size = 6, colour = "black")),
         linetype = guide_legend(order = 3, override.aes = list(colour = c("grey30", "#D62728"))),
         fill = guide_legend(order = 4)) +
  labs(x = expression(paste("Number of basis functions ", n[uc],
                            " approximating the instrumental variation")),
       y = "Truncated exposure effect estimate") +
  theme_minimal(base_size = 22) +
  theme(strip.text = element_text(face = "bold"), panel.spacing.x = unit(3, "lines"),
        legend.position = "top", legend.box = "horizontal", legend.key.width = unit(2.5, "lines"))

ggsave(png_path, p, width = 20, height = 3.6 * length(plot_scenarios) + 4, dpi = 200, bg = "white",
       limitsize = FALSE)
message("Wrote ", png_path)
