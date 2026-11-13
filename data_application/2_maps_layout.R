# Rearranges images/PM_ivs_20pct_2010.png (written by 0_preprocessing.R; one
# column: exposure, then the TPS and GL estimates of A_uc and A_c) into two rows
# for the paper: exposure on the left, then A_uc and A_c for the thin plate
# spline basis (top row) and the graph Laplacian basis (bottom row).
# Crops the existing image, so it runs without the restricted data.
#
# Usage (from data_application/): Rscript 2_maps_layout.R

im <- png::readPNG("images/PM_ivs_20pct_2010.png")[, , 1:3]
exposure <- im[1:980, 420:2085, ]
tps_uc <- im[990:1985, 1:1250, ]
tps_c <- im[990:1985, 1251:2500, ]
gl_uc <- im[1990:2985, 1:1250, ]
gl_c <- im[1990:2985, 1251:2500, ]

gap <- 40
h <- 2 * nrow(tps_uc)
out <- array(1, c(h, ncol(exposure) + gap + 2 * ncol(tps_uc), 3))
put <- function(out, tile, row, col) {
  out[row + seq_len(nrow(tile)) - 1, col + seq_len(ncol(tile)) - 1, ] <- tile
  out
}
out <- put(out, exposure, (h - nrow(exposure)) %/% 2 + 1, 1)
col2 <- ncol(exposure) + gap + 1
col3 <- col2 + ncol(tps_uc)
out <- put(out, tps_uc, 1, col2)
out <- put(out, tps_c, 1, col3)
out <- put(out, gl_uc, nrow(tps_uc) + 1, col2)
out <- put(out, gl_c, nrow(tps_uc) + 1, col3)
png::writePNG(out, "images/PM_ivs_20pct_2010_2x3.png", dpi = 300)
