


<!-- README.md is generated from README.Rmd. Please edit that file -->

# Exploratory Principal Component Analysis

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

`epca` is an R package for **e**xploratory **p**rincipal **c**omponent
**a**nalysis. The goal of `epca` is to comprehend any data matrix that
contains *low-rank* and *sparse* underlying signals of interest.

`epca` features two key tools:

  - `sca` for **s**parse principal **c**omponent **a**nalysis.
  - `sma` for **s**parse **m**atrix **a**pproximation, a two-way data
    analysis for simultaneously row and column dimensionality
    reductions.

<!-- The SMA decomposes a data matrix $X$ into the form of $Z B Y^\top$, where $Z$ and $Y$ are sparse and orthogonal, and $B$ is low-rank. -->

<!-- This is accomplished by minimizing the reconstruction error under the [Frobenius norm](http://mathworld.wolfram.com/FrobeniusNorm.html). -->

## Installation

<!-- You can install the released version of epca from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("epca") -->

<!-- ``` -->

<!-- or the development version from [GitHub](https://github.com/) with: -->

`epca` is not yet on CRAN. You could install development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fchen365/epca")
```

<!-- ## Example -->

<!-- This is a basic example which shows you how to solve a common problem. -->

## Reference

Chen F and Rohe K, “A New Basis for Sparse PCA.”
