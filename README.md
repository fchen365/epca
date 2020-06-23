


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

The usage of `sca` and `sma` is straightforward. For example, to find
`k` sparse PCs of a data matrix `X`:

``` r
sca(X, k)
```

Similarly, we can find a rank-`k` sparse matrix decomposition by

``` r
sma(X, k)
```

For more examples, please see the vignette:

``` r
vignette("epca")
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/fchen365/epca/issues).

## Reference

Chen F and Rohe K, “A New Basis for Sparse PCA.”
