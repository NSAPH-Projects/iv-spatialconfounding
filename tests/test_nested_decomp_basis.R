source("funcs.R") # defines nested_decomp and spectral_decomp
library(tidyverse)

# Create a group matrix
L <- 5
n <- 2^(L + 3)
g <- matrix(0, nrow = n, ncol = L - 1)
for (l in 1:(L - 1)) {
    g[, l] <- 1 + ((0:(n - 1)) %/% (n %/% 2^l))
}

res <- nested_decomp_mats(g)

# example vector
x <- 1:n # just for example and visualization

# orthogonal components
x_averaged <- matrix(0, nrow = n, ncol = L)
x_orth <- matrix(0, nrow = n, ncol = L)
for (j in 1:L) {
    x_averaged[, j] <- res$proj_mat[[j]] %*% x
    x_orth[, j] <- res$decomp_mat[[j]] %*% x
}

# compare with nested decomp basis form
G <- nested_decomp(g)
print(all(abs(t(G) %*% G - diag(n)) < 1e-10))

projected <- list()
trunc_values <- c(0, 2, 4, 8, 16, 32, 64, 128, 256)

for (j in seq_along(trunc_values)) {
    xproj <- get_projection()(x, G, trunc_values[j])
    projected[[j]] <- data.frame(
        x = 1:n, y = xproj, truncation = sprintf("%03d", trunc_values[j])
    )
}

df <- bind_rows(projected)

ggplot(df, aes(x = x, y = y, color = truncation)) +
    geom_step()
ggsave("tests/nested_decomp_basis.png", width = 6, height = 3)
print("done")
