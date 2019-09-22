


<!-- README.md is generated from README.Rmd. Please edit that file -->
Exploratory Principle Component Analysis
========================================

<!-- badges: start -->
<!-- [![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) -->
<!-- badges: end -->
The goal of **e**xploratory **p**rinciple **c**omponent **a**nalysis (EPCA) is to comprehend any data with sparsity nature. The EPCA features sparse PCA via the **s**parse **m**ultivariate **d**ecomposition (SMD). The SMD factorizes a data matrix *X* into the form of *Z**B**Y*<sup>⊤</sup>, where *Z* and *Y* are sparse and orthogonal, and *B* is low-rank. This is accomplished by minimizing the restruction error under the [Frobenius norm](http://mathworld.wolfram.com/FrobeniusNorm.html).

Installation
------------

You can install the released version of epca from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("epca")
```

or the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fchen365/epca")
```

Example
-------

This is a basic example which shows you how to solve a common problem.
