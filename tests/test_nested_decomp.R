source("funcs.R") # defines nested_decomp and spectral_decomp

# Create a group matrix
L <- 3
n <- 128
g <- matrix(0, nrow = n, ncol = L)
for (l in 1:(L - 1)) {
    g[, l] <- 1 + ((0:(n - 1)) %/% (n %/% 2^l))
}
g[ ,L] <- 1:n

res <- nested_decomp_mats(g)

# example vector
x <- 1:n  # just for example and visualization

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
