Decompositions
================
Mauricio Tec
2023-06-08

``` r
source("funcs.R")
```

``` r
# Create a group matrix
L <- 4
n <- 128
g <- matrix(0, nrow = n, ncol = L + 1)
for (j in 1:(L - 1)) {
    g[, j] <- (0:(n - 1)) %/% (n %/% 2^j)
}
g[ ,L] <- 1:n
head(g)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]    0    0    0    1    0
    ## [2,]    0    0    0    2    0
    ## [3,]    0    0    0    3    0
    ## [4,]    0    0    0    4    0
    ## [5,]    0    0    0    5    0
    ## [6,]    0    0    0    6    0

Obtain two lists. The first one contains the projection matrices
$\{P_j\}_{j=1}^L$ and the second one contains the decomposition matrices
$\{A_j\}_{j=1}^L$.

``` r
res <- nested_decomp(g)
names(res)
```

    ## [1] "proj_mats"   "decomp_mats"

Apply decoposition to some vector $x$.

``` r
# example vector
x <- 1:n  # just for example and visualization

# orthogonal components
x_orth <- matrix(0, nrow = n, ncol = L)
x_averaged <- matrix(0, nrow = n, ncol = L)
for (j in 1:L) {
    x_averaged[, j] <- res$proj_mat[[j]] %*% x
    x_orth[, j] <- res$decomp_mat[[j]] %*% x
}
```

Verify that the averages are correct.

``` r
head(x_averaged, 16)
```

    ##       [,1] [,2] [,3] [,4]
    ##  [1,] 32.5 16.5  8.5    1
    ##  [2,] 32.5 16.5  8.5    2
    ##  [3,] 32.5 16.5  8.5    3
    ##  [4,] 32.5 16.5  8.5    4
    ##  [5,] 32.5 16.5  8.5    5
    ##  [6,] 32.5 16.5  8.5    6
    ##  [7,] 32.5 16.5  8.5    7
    ##  [8,] 32.5 16.5  8.5    8
    ##  [9,] 32.5 16.5  8.5    9
    ## [10,] 32.5 16.5  8.5   10
    ## [11,] 32.5 16.5  8.5   11
    ## [12,] 32.5 16.5  8.5   12
    ## [13,] 32.5 16.5  8.5   13
    ## [14,] 32.5 16.5  8.5   14
    ## [15,] 32.5 16.5  8.5   15
    ## [16,] 32.5 16.5  8.5   16

Verify the orthogonality.

``` r
cor(x_orth)
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    1    0    0    0
    ## [2,]    0    1    0    0
    ## [3,]    0    0    1    0
    ## [4,]    0    0    0    1

Verify the reconstruction.

``` r
rowSums(x_orth) == x
```

    ##   [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [46] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [76] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [91] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [106] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [121] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
