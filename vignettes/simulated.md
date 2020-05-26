---
title: "Introduction to `epca`"
author: "Fan Chen (fan.chen@wisc.edu)"
date: "2020-05-25"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Introduction to `epca`}
  %\usepackage[UTF-8]{inputenc}
---



## Example

Let's simulate a rank-5 data matrix with some additive Gaussian noise.

```r
n <- 300
p <- 50
k <- 5
Z <- shrinkage(polar(matrix(runif(n * k), n, k)), sqrt(n))
B <- diag(5) * 3
Y <- shrinkage(polar(matrix(runif(p * k), p, k)), sqrt(p))
E <- matrix(rnorm(n * p, sd = 0.01), n, p)
X <- scale(Z %*% B %*% t(Y) + E)
```

Now, perform the sparse PCA.

```r
s.sca <- sca(X, k)
s.sca
```

```
## Call:sca(A = X, k = k)
## 
## 
## Num. non-zero loadings': 21 22 22 31 35 
## Abs. sum loadings': 3.467996 
## Cumulative proportion of variance explained (CPVE): 
##                      CPVE
## First component:    0.077
## First 2 components: 0.152
## First 3 components: 0.225
## First 4 components: 0.291
## First 5 components: 0.358
```

Similarly, we can do sparse matrix decomposition.

```r
s.sma <- sma(X, k)
s.sma
```

```
## Call: sma(A = X, k = k)
## 
## 
## Num. non-zero Z's:  201 194 201 208 196 
## Num. non-zero Y's:  23 24 30 36 20 
## Abs. sum Z's:  8.066438 
## Abs. sum Y's:  3.482323
```

