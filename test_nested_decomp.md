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
n <- 32
g <- matrix(0, nrow = n, ncol = L)
for (j in 1:(L - 1)) {
    g[, j] <- 1 + (0:(n - 1)) %/% (n %/% 2^j)
}
g[ ,L] <- 1:n
head(g, 16)
```

    ##       [,1] [,2] [,3] [,4]
    ##  [1,]    1    1    1    1
    ##  [2,]    1    1    1    2
    ##  [3,]    1    1    1    3
    ##  [4,]    1    1    1    4
    ##  [5,]    1    1    2    5
    ##  [6,]    1    1    2    6
    ##  [7,]    1    1    2    7
    ##  [8,]    1    1    2    8
    ##  [9,]    1    2    3    9
    ## [10,]    1    2    3   10
    ## [11,]    1    2    3   11
    ## [12,]    1    2    3   12
    ## [13,]    1    2    4   13
    ## [14,]    1    2    4   14
    ## [15,]    1    2    4   15
    ## [16,]    1    2    4   16

Obtain two lists. The first one contains the projection matrices
$(P_j)_{j=1}^L$ and the second one contains the decomposition matrices
$(A_j)_{j=1}^L$.

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
    ##  [1,]  8.5  4.5  2.5    1
    ##  [2,]  8.5  4.5  2.5    2
    ##  [3,]  8.5  4.5  2.5    3
    ##  [4,]  8.5  4.5  2.5    4
    ##  [5,]  8.5  4.5  6.5    5
    ##  [6,]  8.5  4.5  6.5    6
    ##  [7,]  8.5  4.5  6.5    7
    ##  [8,]  8.5  4.5  6.5    8
    ##  [9,]  8.5 12.5 10.5    9
    ## [10,]  8.5 12.5 10.5   10
    ## [11,]  8.5 12.5 10.5   11
    ## [12,]  8.5 12.5 10.5   12
    ## [13,]  8.5 12.5 14.5   13
    ## [14,]  8.5 12.5 14.5   14
    ## [15,]  8.5 12.5 14.5   15
    ## [16,]  8.5 12.5 14.5   16

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

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [31] TRUE TRUE
