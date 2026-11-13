# Map of the five confounding mechanisms: rows are Confounding mechanisms 1-5
# (LS, MS, twoU, coordnoise, stateJump), columns are the exposure
# A = A^c + A^uc, the instrumental variation A^uc, the non-instrumental
# variation A^c, the unmeasured confounder U, and the measured covariate X.
# Mechanism 3 shows both confounders U1 and U2 side by side.
#
# Evaluated on a 40 km grid over the US, with one draw of the non-spatial
# noise (eta, nu, eps) per cell, so X, A^c, A^uc, and A appear as they do in
# data. Each panel has its own symmetric colour scale (blue < 0 < red).
#
# Usage (from simulation/): Rscript make_decomp.R
# Requires sim.rds. Writes images/all-decomp.jpeg.

suppressMessages({
  library(sf)
  library(ggplot2)
})
source("funcs.R")

sim <- readRDS("sim.rds")
outline <- sim$outline

cell <- 40
grid_points <- st_coordinates(st_make_grid(outline, cellsize = cell, what = "centers"))
inside <- lengths(st_intersects(st_as_sf(as.data.frame(grid_points), coords = c("X", "Y"),
                                         crs = st_crs(outline)), outline)) > 0
grid_points <- grid_points[inside, , drop = FALSE]
colnames(grid_points) <- c("x", "y")

mechanisms <- c("LS", "MS", "twoU", "coordnoise", "stateJump")
components <- c("A", "Auc", "Ac", "U", "X")

# Components of one mechanism on the grid (list of named vectors; U may hold U1, U2).
components_of <- function(sc) {
  set.seed(1)  # same noise draw in every mechanism
  d <- draw_scenario_data(sim, sc, grid_points)
  U <- if (is.null(d$U1)) list(U = d$U) else list(U1 = d$U1, U2 = d$U2)
  list(A = d$A, Auc = d$Auc, Ac = d$Ac, U = U, X = d$X[, 1])
}

# Two confounders share one panel: each map is shrunk and placed side by side.
xr <- range(grid_points[, "x"]); yr <- range(grid_points[, "y"])
shrink <- 0.47
place <- function(k, of) {
  if (of == 1L) return(list(x = grid_points[, "x"], y = grid_points[, "y"], w = cell, h = cell))
  w <- diff(xr)
  list(x = xr[1] + (grid_points[, "x"] - xr[1]) * shrink + (k - 1) * w * (1 - shrink),
       y = mean(yr) + (grid_points[, "y"] - mean(yr)) * shrink,
       w = cell * shrink, h = cell * shrink)
}

rows <- list()
for (i in seq_along(mechanisms)) {
  comp <- components_of(mechanisms[i])
  for (cc in components) {
    maps <- if (cc == "U") comp$U else setNames(list(comp[[cc]]), cc)
    for (k in seq_along(maps)) {
      v <- maps[[k]]
      pl <- place(k, length(maps))
      rows[[length(rows) + 1]] <- data.frame(
        mechanism = i, component = cc, x = pl$x, y = pl$y, w = pl$w, h = pl$h,
        value = v / max(abs(v))  # per-map symmetric scale
      )
    }
  }
}
df <- do.call(rbind, rows)

row_labels <- setNames(sprintf('atop(bold("Confounding"), bold("mechanism %d"))', seq_along(mechanisms)),
                       seq_along(mechanisms))
col_labels <- c(A = '"Exposure"~A', Auc = '"Instrumental variation"~A[uc]',
                Ac = '"Non-instrumental variation"~A[c]', U = '"Confounder"~U',
                X = '"Measured covariate"~X')
df$mechanism <- factor(row_labels[as.character(df$mechanism)], levels = row_labels)
df$component <- factor(col_labels[df$component], levels = col_labels)

p <- ggplot(df, aes(x, y, fill = value, width = w, height = h)) +
  geom_tile() +
  facet_grid(mechanism ~ component, switch = "y", labeller = label_parsed) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1),
                       guide = "none") +
  coord_equal() +
  theme_void(base_size = 20) +
  theme(strip.text.x = element_text(size = 20, margin = margin(b = 8)),
        strip.text.y.left = element_text(size = 18, angle = 0, margin = margin(r = 10)),
        strip.placement = "outside", panel.spacing = unit(1, "lines"))

dir.create("images", showWarnings = FALSE)
ggsave("images/all-decomp.jpeg", p, width = 22.5, height = 3 * length(mechanisms), dpi = 250, bg = "white")
message("Wrote images/all-decomp.jpeg")
