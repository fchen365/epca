## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      tidy = TRUE, 
                      tidy.opts = list(comment = FALSE))
library(epca)

## ----simu---------------------------------------------------------------------
## simulate a rank-5 data matrix with some additive Gaussian noise 
n <- 300
p <- 50
k <- 5 ## rank
Z <- shrinkage(polar(matrix(runif(n * k), n, k)), sqrt(n))
B <- diag(5) * 3
Y <- shrinkage(polar(matrix(runif(p * k), p, k)), sqrt(p))
E <- matrix(rnorm(n * p, sd = .01), n, p)
X <- scale(Z %*% B %*% t(Y) + E)

## -----------------------------------------------------------------------------
## perform sparse PCA
s.sca <- sca(X, k)
s.sca

## -----------------------------------------------------------------------------
## perform sparse matrix approximation
s.sma <- sma(X, k)
s.sma

