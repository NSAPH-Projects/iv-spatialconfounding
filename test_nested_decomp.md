Decompositions
================
Mauricio Tec
2023-06-08

``` r
source("funcs.R")
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

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
