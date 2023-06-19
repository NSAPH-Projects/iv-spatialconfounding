source("funcs.R") # defines nested_decomp and spectral_decomp
library(tidyverse)

# Create a group matrix
L <- 4
n <- 128
g <- matrix(0, nrow = n, ncol = L - 1)
for (l in 1:(L - 1)) {
    g[, l] <- 1 + ((0:(n - 1)) %/% (n %/% 2^l))
}

res <- nested_decomp_mats(g)

# example vector
x <- 1:n # just for example and visualization

# orthogonal components
x_orth <- matrix(0, nrow = n, ncol = L)
x_averaged <- matrix(0, nrow = n, ncol = L)
for (j in 1:L) {
    x_averaged[, j] <- res$proj_mat[[j]] %*% x
    x_orth[, j] <- res$decomp_mat[[j]] %*% x
}

# print group averages
print(x_averaged)

# print orthogonal decomposition
print(x_orth)

# verify orthogonal decomposition
print(all(x == rowSums(x_orth)))

# compare with nested decomp fast
res2 <- nested_decomp(x, g)
print(all(x_averaged == res2$avs))
print(all(x_orth == res2$orth))

# eigen approach
# make group indicator matrix
tol <- 1e-4
G <- matrix(0, nrow = n, ncol = n)
num_vect <- 0

# starts with the coarses group
basis <- matrix(0, nrow = n, ncol = n)
gm <- g[, L - 1] # last column is the idenitty

# group size with zero at end
d <- numeric(L)
for (l in seq_len(L - 1)) {
    d[l] <- length(unique(g[, l]))
}
d[L] <- n


# complete the basis
for (l in (L - 1):1) {
    # get the group indicator matrix
    gl <- g[, l]

    # get the number of unique groups
    num_vect <- length(unique(gl))

    # get the projection matrix of this and the next level
    Al <- res$proj_mat[[l]]

    # get the eigen decomposition
    El <- eigen(Al)

    # get the basis
    ixs <- (d[l] + 1):n
    cols <- El$vectors[, ixs]
    for (i in seq_len(ncol(cols))) {
        cols[, i] <- cols[, i] / sqrt(sum(cols[, i]^2))
    }

    # make them orthogonal to the previous basis elements
    if (l < (L - 1)) {
        next_ixs <- (d[l + 1] + 1):n
        prev <- basis[, next_ixs]
        # compute all inner products to orthogonalize to the rank
        inner_prods <- t(prev) %*% cols
        cols <- cols - prev %*% inner_prods
        #
        El <- svd(cols)
        ixs <- (d[l] + 1):d[l + 1]
        k <- length(ixs)
        cols <- El$u[, 1:k]
        for (i in seq_len(ncol(cols))) {
            cols[, i] <- cols[, i] / sqrt(sum(cols[, i]^2))
        }
    }
    basis[, ixs] <- cols
}
# complete with the image of the first matrix
unique_vals <- unique(g[, 1])
for (i in seq_len(d[1])) {
    basis[, i] <- as.numeric(g[, 1] == unique_vals[i])
    basis[, i] <- basis[, i] / sqrt(sum(basis[, i]^2))
}


# test orthonormality
print(all(abs(t(basis) %*% basis - diag(n)) < 0.0001))


# get the orthogonal components
x_orth2 <- t(basis) %*% x

project <- function(x, basis, k) {
    xstar <- t(basis) %*% x
    n <- length(xstar)
    if (k < n) {
        xstar[(k + 1):n] <- 0
    }
    xproj <- basis %*% xstar
    return(xproj)
}

projected <- list(data.frame(x = 1:n, y = x, truncation = "original"))
for (j in seq_along(d)) {
    xproj <- project(x, basis, d[j])
    df <- data.frame(x = 1:n, y = xproj, truncation = as.character(d[j]))
    projected[[j + 1]] <- df
}
df <- bind_rows(projected)
ggplot(df, aes(x = x, y = y, color = truncation)) +
    geom_line() +
    theme_bw()
ggsave("tests/orthogonal_decomp.png", width = 6, height = 3)

print("done")
