# Copies of the basis-selection figures with the x-axis symbol n_uc replaced
# by an italic k, matching the notation of the paper. The figures cannot be
# regenerated outside ReD, so the label is edited in the image: the trailing
# "n_uc" is erased, an italic k in Noto Sans (the font of the original) is drawn
# in its place, and the label is re-centred.
#
# Usage (from data_application/):
#   NOTO_SANS_DIR=<folder with NotoSans-{Regular,Italic}.ttf> Rscript 3_relabel_k.R
# Writes images/full_sample_basis_selection_cutoff{6,7,8}_k.png.

font_dir <- Sys.getenv("NOTO_SANS_DIR")
if (nzchar(font_dir)) {
  systemfonts::register_font("Noto Sans",
                             plain = file.path(font_dir, "NotoSans-Regular.ttf"),
                             italic = file.path(font_dir, "NotoSans-Italic.ttf"))
}
stopifnot("Noto Sans" %in% c(systemfonts::system_fonts()$family,
                             systemfonts::registry_fonts()$family))

# Italic k at the label's font size (64.5 px matches the label's "n").
render_k <- function(size = 64.5) {
  f <- tempfile(fileext = ".png")
  ragg::agg_png(f, 300, 300, res = 72, background = "white")
  grid::grid.text("k", x = 0.3, y = 0.3, just = c("left", "bottom"),
                  gp = grid::gpar(fontfamily = "Noto Sans", fontsize = size, fontface = 3))
  dev.off()
  m <- png::readPNG(f)[, , 1:3]
  d <- (m[, , 1] + m[, , 2] + m[, , 3]) / 3 < 0.99
  m[range(which(rowSums(d) > 0))[1]:range(which(rowSums(d) > 0))[2],
    range(which(colSums(d) > 0))[1]:range(which(colSums(d) > 0))[2], , drop = FALSE]
}
k_img <- render_k()

relabel <- function(in_path, out_path) {
  im <- png::readPNG(in_path)
  im <- im[, , 1:3]
  g <- (im[, , 1] + im[, , 2] + im[, , 3]) / 3
  # The x-axis title is the lowest band of text.
  rows <- which(rowSums(g[seq(nrow(g) - 250, nrow(g)), ] < 0.6) > 0) + nrow(g) - 251
  band <- split(rows, cumsum(c(1, diff(rows) > 1)))
  rows <- band[[length(band)]]
  rows <- seq(min(rows) - 5, max(rows) + 10)
  cols <- which(colSums(g[rows, ] < 0.6) > 0)
  clusters <- split(cols, cumsum(c(1, diff(cols) > 4)))
  n_cl <- length(clusters)
  # Last four clusters: comma, n, u, c.
  comma <- clusters[[n_cl - 3]]; n_glyph <- clusters[[n_cl - 2]]
  stopifnot(length(comma) < 15, length(n_glyph) > 20, length(n_glyph) < 35)
  baseline <- max(which(rowSums(g[rows, n_glyph] < 0.6) > 0)) + min(rows) - 1
  left <- min(cols); right_old <- max(cols)
  # New label: text up to and including the comma, then k.
  strip <- im[rows, left:max(comma), , drop = FALSE]
  k_left <- min(n_glyph) - left + 1
  new_w <- k_left + ncol(k_img) - 1
  new <- array(1, c(length(rows), new_w, 3))
  new[, seq_len(ncol(strip)), ] <- strip
  k_rows <- (baseline - min(rows) + 1) - nrow(k_img) + seq_len(nrow(k_img))
  new[k_rows, k_left - 1 + seq_len(ncol(k_img)), ] <- k_img
  # Erase the old label and paste the new one, centred where the old one was.
  im[rows, left:right_old, ] <- 1
  centre <- (left + right_old) / 2
  new_left <- round(centre - (new_w - 1) / 2)
  im[rows, new_left - 1 + seq_len(new_w), ] <- new
  png::writePNG(im, out_path)
  message("Wrote ", out_path)
}

for (cutoff in 6:8) {
  relabel(sprintf("images/full_sample_basis_selection_cutoff%d.png", cutoff),
          sprintf("images/full_sample_basis_selection_cutoff%d_k.png", cutoff))
}
