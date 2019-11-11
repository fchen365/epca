## ------------------------------------------------------------
## what: helper functions
## ------------------------------------------------------------

## ------------------------------------------------------------
## numerical
## ------------------------------------------------------------
#' Calculate fractional exponent/power
#'
#' Calculate fractional exponent/power, a^(num/den), where a could be negative.
#' @param a numeric, base (could be nagetive).
#' @param num positive integer, numerator of the exponent.
#' @param den positive integer, denominator of the exponent.
#' @return numeric, the evaluated a^(num/den)
#' @export
fracExp <- function(a, num, den) {
  if (num < 0 || den <= 0 || num %% 1 || den %% 1) {
    stop("num and den should be meaningful integers.")
  }
  if (den %% 2) {
    res = (abs(a) ^ (1/den)) * sign(a)
  } else {
    res = a ^ (1 / den)
  }
  return(res ^ num)
}

#' Calculate (1,1)-norm
#'
#' This is a helper function to \code{\link{soft}}.
norm11 = function(mat) {
  sum(abs(mat))
}


#' Soft-thresholding
#'
#' Evaluate the soft-thresholding function with either
#' (1) a cut-off value on elements, or
#' (2) a overall constraint on the (1,1)-norm.
#'
#' @param mat matrix or Matrix, to be soft-thresholded.
#' @param cutoff numeric, soft-thresholding cutoff. One and only one of cutoff and constraint should/can be specified.
#' @param constraint numeric, the overall (1,1)-norm constraint value (\eqn{\gamma}), i.e. \code{\link{norm11}}(mat) \eqn{\le \gamma} . One and only one of cutoff and constraint should/can be specified.
#' @param epsilon numeric, precision tolerance.
#' @details For the second option (constrainted soft-thresholding), this function uses a binary search to find the cut-off value.
#' @return a list that contains:
#' \item{matrix}{matrix, the matrix that results from soft-thresholding}
#' \item{norm}{numeric, the norm of the matrix after soft-thresholding. This value is close to constraint if using the second option.}
#' @export
soft <- function(mat, cutoff = NULL, constraint = NULL,
                 epsilon = 2*.Machine$double.eps) {
  if (is.null(cutoff) + is.null(constraint) != 1)
    stop("Either cutoff or constraint should/can be specified.")
  if (!is.null(cutoff)) {
    sign.mat = sign(mat)
    abs.mat = abs(mat)
    res = list()
    res$matrix = sign.mat * pmax(abs.mat - cutoff, 0)
    res$norm = norm11(res$matrix)
  }
  if (!is.null(constraint)) {
    ## no thresholding needed
    res = soft(mat, cutoff = 0)
    if (res$norm <= constraint) return(res)

    ## binary search
    lower = 0 ## whose res$norm > constraint
    upper = max(abs(mat)) ## whose res$norm = 0
    while (upper - lower > epsilon) {
      mid = (lower + upper) / 2
      res = soft(mat, cutoff = mid)
      if (res$norm > constraint) {
        lower = mid
      } else
        upper = mid
    }
  }
  return(res)
}


#' Polar decomposition
#'
#' Takes a matrix X, of dimensions n by p, and calculate polar components,
#' two matrices P and Z, such that X = UP,
#' assuming X has full col rank
#' @param X matrix or Matrix.
#' @return a matrix of the U part
#' @export
polar = function(X) {
  if (nrow(X) < ncol(X))
    stop('Input should be a tall (or at least square) matrix.')
  if (Matrix::rankMatrix(X) < ncol(X))
    stop('X is not full rank.')
  S = svd(X)
  U.part <- S$u %*% t(S$v)
  # P.part <- S$v %*% diag(S$d) %*% t(S$v)
  return(U.part)
}

#' Find root matrix
#'
#' Find X from the Gram matrix X^T X, i.e. sqrt(crossprod(X)).
#'
rootmatrix <- function(x) {
  x.eigen <- eigen(x)
  d <- x.eigen$values
  d <- (d + abs(d))/2
  v <- x.eigen$vectors
  root = v %*% diag(sqrt(d)) %*% t(v)
  dimnames(root) <- dimnames(x)
  return(root)
}

#' Orthogonal rotation for the SMD
#'
#' Perform varimax rotation on X from right.
#' Flip the signs of columns so that the resulting matrix is positive skewed.
#' @param X matrix or Matrix.
#' @param normalize logical, whether to (default to FALSE) first normalize rows before (and scale back afterward) the rotation.
#' @param eps pricision tolerance.
#' @return the rotated matrix
#' @export
rotation = function(X, normalize = F, eps = 1e-6) {
  if (is.null(dim(X)) || ncol(X) == 1) {
    warning('Input only one column; no rotation applied.')
    return(X)
  }
  R = varimax(X, normalize = normalize, eps = eps)$rotmat
  L = X %*% R
  ## (1) switch signs s.t. each column has positive third moment.
  signs = sign(colSums(L^3))
  L = t(t(L) * signs) ## L = sweep(L, 2, signs, '*')
  ## (2) order columns by Var(Z^2)
  dev = apply(R^2, 2, var)
  L = L[,order(dev, decreasing = T)]
  return(L)
}

#' Custom inner produce of matrices
#'
#' Calculate custom matrix inner product Z = <X, Y>, where Z[i,j] = FUN(X[,i], Y[,j]).
#'
#' @param X,Y matrix or Matrix.
#' @param FUN function, which takes in two vectors and returns a numeric value.
#' @param ... additional parameters for FUN.
#' @return a matrix of <X, Y>
#' @export
inner = function(X, Y, FUN = crossprod, ...) {
  t(apply(X, 2, function(x) {
    apply(Y, 2, function(y) FUN(x, y))
  }))
}

#' Proportion of variance explained
#'
#' Calculate the variance in a matrix explained by a set of linear transformation, (e.g. eigenvectors).
#' @param mat matrix or Matrix, the original data matrix X or cov(X) = crossprod(X) / (nrow(X) - 1)
#' @param V coefficients of linear transformation, e.g. a set of loadings in PCA.
#' @param type the input matrix for SMD (not used at the moment), either
#'           'predictor' (original data matrix) or 'Gram' (useful when n >> p).
#' @return a numeric value between [0,1], the proportion of variance in mat explained by Y.
#' @export
var.explain = function(mat, V, type = "predictor") {
  if (type == "predictor") {
    V <- as.matrix(V)
    XV = mat %*% V %*% solve(t(V) %*% V) %*% t(V)
    norm(XV, type = 'F')^2 / norm(mat, 'F')^2
  } else if (type == "Gram"){
    V <- as.matrix(V)
    TV <- sum(svd(mat)$d^2) ## total variance
    CL <- ifelse(sqrt(colSums(V^2)),sqrt(colSums(V^2)),0) ## column length, replace 0 with 1
    S <- mat %*% (V / CL) ## component scores
    EV <- diag(qr.R(qr(S))^2) ## explained variance
    sum(EV)/TV
    # CI = solve(crossprod(V))
    # sum(diag(V %*% CI %*% t(V) %*% mat %*% V %*% CI %*% t(V)))
  }
}

## ------------------------------------------------------------
## clustering
## ------------------------------------------------------------
#' Permute columns of a block membership matrix
#'
#' Perform column permutation of block membership matrix for aesthetical visualization.
#' That is, the k-th column gives k-th cluster. This is done by ranking the column sums of squares (by default).
#'
#' @param x a non-negative matrix, nNode x nBlock,
#' @param s integer, order of non-linear
#' @export
permColumn = function(x, s = 2) {
  if (!is.array(x)) stop("x is not array like.")
  x[,order(colMeans(row(x) * abs(x)^s))]
}


#' Label cluster
#'
#' Given the block membership matrix x, where x[i,j] is the likelihood of node i belongs to block j. Finds the cluster label.
#' @param x matrix, nNode x nBlock, non-negative, higher value means more likely
#' @param ties.method a character string specifying how ties are treated,
#'                  see \code{\link[base]{rank}}, only "random", "first", "last" are allowed.
#' @return integer vector %in% 1:ncol(x)
labelCluster <- function(x, ties.method = "random") {
  x1 <- t(apply(-abs(x), 1, rank, ties.method = ties.method))
  apply(x1, 1, which.min)
}


#' Match two matrices
#'
#' Match the columns of two matrices of the same size.
#' This is accomplished by solving an optimization over column permutation.
#' Given two matrices, X and X.new, find permutation p() that minimizes
#'      sum_i similarity(X.new[,p(i)], X[,i]),
#' where the similarity can be euclidean distance, cosine, maximum difference.
#'
#' @param X.new,X matrix or Matrix, with the same dimension.
#' @param method one of "maximum" (default), "cosine", or "euclidean".
#' @return X.new after permuting the columns.
#' @export
matchCols = function(X.new, X, method = 'maximum') {
  ## check
  if (!is.array(X.new) || dim(X.new)[2] == 1) return(X.new)

  ## (1) similarity (cosine)
  if (method == 'cosine') {
    norm.new = sqrt(colSums(X.new^2))
    norm = sqrt(colSums(X^2))
    cosine = t(X.new) %*% X / (norm.new %*% t(norm))
    x = 1 - cosine
  }
  ## (2) usual distance (euclidean)
  if (method == 'euclidean') {
    euclidean = function(x, y) sqrt(sum((x-y)^2))
    x = inner(X.new, X, euclidean)
  }
  ## (3) supremum norm (maximum)
  if (method == 'maximum') {
    maximum = function(x, y) max(abs(x-y))
    x = inner(X.new, X, maximum)
  }

  p = clue:: solve_LSAP(x, maximum = F)
  X.new[,as.vector(p)]
}


## FUN: calculate the empirical proportion of misclustered nodes
##      assuming #{cluster} = #{block},
##      allowing permutation on clusters
# INPUT:
##     clusters - vector of intergers or factors, cluster membership (estimated)
##     truth - true block labels
misClustRate = function(cluster, truth) {
  ## check
  if (length(unique(cluster)) != length(unique(truth))) {
    stop("#{cluster} != #{block}. \n",
         "Try send in cluster (and truch) as a factor vector, ",
         "with matching number of cluters.")
  }

  # count cluster labels in each block
  # clusterCounts[i,j] = #{cluster i ^ block j}
  clusterCounts <- table(cluster, truth)

  # maximum matches, linear solver
  p = clue:: solve_LSAP(clusterCounts, maximum = T)
  correct <- sum(clusterCounts[cbind(seq_along(p), p)])
  return(1 - correct / length(truth))
}


## ------------------------------------------------------------
## rotation variants via gradient projection
## ------------------------------------------------------------
#' The varimax criteria
#'
#' Calculate the varimax criteria.
#' @references \href{https://en.wikipedia.org/wiki/Varimax_rotation}{Varimax rotation (Wikipedia)}
#' @export
varimax.criteria = function(mat) {
  sum(apply(mat^2, 2, var))
}

#' Re-implementation of varimax
#' This is an extended version of \code{\link[stats]{varimax}}, allowing control of iteration.
#' @inheritParams stats::varimax
#' @export
varimax <- function(x, normalize = F, eps = 1e-05, maxit = 1000L) {
  nc <- ncol(x)
  if (nc < 2) return(x)
  if (normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  p <- nrow(x)
  TT <- diag(nc)
  d <- 0
  for (i in seq_len(maxit)) {
    z <- x %*% TT
    cm2 = colMeans(z^2)
    B <- t(x) %*% (z^3 - z %*% diag(cm2))
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    if (d < dpast * (1 + eps)) break
  }
  z <- x %*% TT
  if (normalize) z <- z * sc
  dimnames(z) <- dimnames(x)
  list(loadings = z, rotmat = TT, iter = i)
}


#' The absmax2 criteria
#'
#' Calculate the absmax2 criteria.
#' @references F. Chen, K Rohe. SPCA paper.
#' @export
varimax2.criteria = function(mat) {
  sum((apply(mat^2, 2, var) - 2 * colMeans(mat^2))^2)
}


## FUN: given x \in \R^{p \times k}, find the orthogonal matrix (rotation) T that
##      maximizes the varimax2 objective (z = x %*% T),
##      varimax2(z) = \sum_{j=1}^k \[\text{Var}(z_{\cdot j}^2) - 2\]^2
#' @export
varimax2 = function (x, normalize = F, eps = 1e-05, maxit = 1000L)
{
  nc <- ncol(x)
  if (nc < 2)
    return(x)
  if (normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  p <- nrow(x)
  TT <- diag(nc)
  d <- 0
  for (i in seq_len(maxit)) {
    z <- x %*% TT
    cs2 = colSums(z^2)
    cs4 = colSums(z^4)
    B <- t(x) %*% (z^3 - 3 * z %*% diag(cs2) / p) %*% diag(cs4 - 3 * cs2^2 / p)
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    if (d < dpast * (1 + eps))
      break
  }
  z <- x %*% TT
  if (normalize)
    z <- z * sc
  dimnames(z) <- dimnames(x)
  list(loadings = z, rotmat = TT, iter = i)
}

#' The absmin criteria
#'
#' Calculate the absmin criteria.
#' @references F. Chen, K Rohe. SPCA paper.
#' @export
absmin.criteria = function(mat) {
  sum(abs(mat))
}

## FUN: given x \in \R^{p \times k}, find the orthogonal matrix (rotation) T that
##      minimizes the absmin objective (z = x %*% T),
##      absmin_criteria(z) = \sum_{i=1}^p \sum_{j=1}^k |z_ij|
#' @export
absmin = function (L, Tmat = diag(ncol(L)), normalize = FALSE,
                   eps = 1e-5, maxit = 1000) {
  GPArotation::GPForth(L, Tmat = Tmat, method = "absmin",
                       normalize = normalize, ## al = 0.1,
                       eps = eps, maxit = maxit)
}

## FUN: gradient of absmin
vgQ.absmin = function (L) {
  list(Gq = sign(L), f = absmin.criteria(L), Method = "absmin")
}

