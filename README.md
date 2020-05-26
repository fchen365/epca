


<!-- README.md is generated from README.Rmd. Please edit that file -->

# Exploratory Principal Component Analysis

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

`epca` is an R package for comprehending any data matrix that contains
*low-rank* and *sparse* underlying signals of interest. The package
currently features two key tools:

  - `sca` for **s**parse principal **c**omponent **a**nalysis.
  - `sma` for **s**parse **m**atrix **a**pproximation, a two-way data
    analysis for simultaneously row and column dimensionality
    reductions.

## Installation

<!-- You can install the released version of epca from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("epca") -->

<!-- ``` -->

<!-- or the development version from [GitHub](https://github.com/) with: -->

`epca` is not yet on CRAN. You could install the development version
from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fchen365/epca")
```

## Example

Here is an example with simulated data. First, simulate a rank-5 data
matrix with some additive Gaussian noise.

``` r
n <- 300
p <- 50
k <- 5
Z <- shrinkage(polar(matrix(runif(n * k), n, k)), sqrt(n * k))
B <- diag(5) * 3
Y <- shrinkage(polar(matrix(runif(p * k), p, k)), sqrt(p * k))
E <- matrix(rnorm(n * p, sd = 0.01), n, p)
X <- scale(Z %*% B %*% t(Y) + E)
```

Now, perform the sparse PCA.

``` r
sca(X, k)
```

Similarly, we can do sparse matrix decomposition.

``` r
sma(X, k)
```

## Reference

Chen F and Rohe K, “A New Basis for Sparse PCA.”
