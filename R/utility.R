## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## what: helper functions
## who: fan chen (fan.chen@wisc.edu)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ clustering ------ 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Permute columns of a block membership matrix
#'
#' Perform column permutation of block membership matrix for aesthetic visualization.
#' That is, the k-th column gives k-th cluster. This is done by ranking the column sums of squares (by default).
#'
#' @param x a non-negative matrix, nNode x nBlock,
#' @param s integer, order of non-linear
permColumn = function(x, s = 2) {
  x[,order(colMeans(row(x) * abs(x)^s))]
} 


#' Label Cluster
#'
#' Assign cluster labels to each row from the membership matrix. 
#' @param x `matrix` with non-negative entries, where `x[i,j]` is the estimated likelihood (or any equivalent measure) of node i belongs to block j. The higher the more likely.
#' @param ties.method `character`, how should ties be handled, "random", "first", "last" are allowed. See [base::rank()] for details.
#' @return `integer` vector of the same length as `x`. Each entry is one of 1, 2, ..., `ncol(x)`.
labelCluster <- function(x, ties.method = "random") {
  rank <- apply(-abs(x), 1, rank, ties.method = ties.method)
  apply(t(rank), 1, which.min)
}


#' Mis-Classification Rate (MCR)
#' 
#' Compute the empirical MCR, assuming that #{cluster} = #{block}, 
#' This calculation allows a permutation on clusters.
#' 
#' @param cluster vector of `integer` or `factor`, estimated cluster membership.
#' @param truth a vector of the same length as `clusters`, the true cluster labels.
#' @return `numeric`, the MCR.
#' @examples 
#' truth = rep(1:3, each = 30) 
#' cluster = rep(3:1, times = c(25, 32, 33))
#' misClustRate(cluster, truth)
#' @export
misClustRate = function(cluster, truth) {
  ## check
  if (length(unique(cluster)) != length(unique(truth))) {
    stop("cluster and truth shoule be factor vectors with equal number of levels.")
  }
  
  # count cluster labels in each block
  # clusterCounts[i,j] = #{cluster i ^ block j}
  clusterCounts <- table(cluster, truth)
  
  # maximum matches, linear solver
  p = clue::solve_LSAP(clusterCounts, maximum = TRUE)
  correct <- sum(clusterCounts[cbind(seq_along(p), p)])
  return(1 - correct / length(truth))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ matrix ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Polar Decomposition
#' 
#' Perform the polar decomposition of an n x p (n > p) matrix X into U P, 
#' where U is an n x p matrix with orthogonal columns (i.e. `crossprod(U)` is the identity matrix), 
#' and P is a p x p positive-semidefinite Hermitian matrix.
#' The function returns the U matrix. 
#' This is a helper function of [prs()].
#' 
#' @param x a `matrix` or `Matrix`, which is presumed full-rank.
#' @return a `matrix` of the unitary part of the polar decomposition.
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse PCA." 
#' @examples 
#' x <- matrix(1:6, nrow = 3)
#' polar_x <- polar(x)
#' 
#' @export
polar = function(x) {
  stopifnot(nrow(x) >= ncol(x))
  stopifnot(Matrix::rankMatrix(x) == ncol(x))
  s = svd(x)
  U <- tcrossprod(s$u, s$v)
  # P <- s$v %*% diag(S$d) %*% t(s$v)
  return(U)
}

#' Find root matrix
#' 
#' Find `X` from the Gram matrix (i.e., `crossprod(X)`).
#' @param x a symmetric `matrix`.
rootmatrix <- function(x) {
  stopifnot(isSymmetric(x))
  x.eigen <- eigen(x)
  d <- x.eigen$values
  d <- (d + abs(d))/2
  v <- x.eigen$vectors
  root = v %*% diag(sqrt(d)) %*% t(v)
  dimnames(root) <- dimnames(x)
  return(root)
}

#' Matrix Inner Product
#'
#' Calculate the custom matrix inner product, 
#' `Z = crossprod(X, Y)`, 
#' where `Z[i,j] = FUN(X[,i], Y[,j]`).
#'
#' @param X,Y `matrix` or `Matrix`.
#' @param FUN `function` or a `character(1)` name of base function. 
#' The function should take in two vectors as input and ouput a `numeric(1)` result.
#' @param ... additional parameters for `FUN`.
#' @return `matrix`, inner product of `X` and `Y`.
#' @examples 
#' x <- matrix(1:6, 2, 3)
#' y <- matrix(7:12, 2, 3)
#' ## The default is equivalent to `crossprod(x, y)`
#' inner(x, y) 
#' ## We can compute the pair-wise Euclidean distance of columns.
#' EuclideanDistance = function(x, y) crossprod(x, y)^2
#' inner(x, y, EuclideanDistance)
#' 
#' @export
inner = function(X, Y, FUN = "crossprod", ...) {
  stopifnot(length(FUN) == 1) 
  if (is.character(FUN)) FUN <- match.fun(FUN)
  res <- apply(X, 2, function(x) {
    apply(Y, 2, function(y) FUN(x, y))
  })
  t(res)
}

#' Proportion of Variance Explained (PVE)
#'
#' Calculate the variance in a matrix explained by a set of linear transformation, (e.g. eigenvectors).
#' @param mat `matrix` or `Matrix`, the original data matrix `X` or `cov(X) = crossprod(X) / (nrow(X) - 1)`
#' @param V `matrix`, coefficients of linear transformation, e.g., loadings (in PCA).
#' @param is.cov `logical`, whether the input matrix is a covariance matrix or a Gram matrix.
#' @return a `numeric` value between 0 and 1, the proportion of variance in mat explained by Y.
#' @examples 
#' ## use the "swiss" data
#' ## find two sparse PCs 
#' s.sca <- sca(swiss, 2, gamma = sqrt(ncol(swiss))) 
#' ld <- loadings(s.sca)
#' pve(as.matrix(swiss), ld)
#' @export
pve = function(mat, V, is.cov = FALSE) {
  if (!is.cov) {
    # V <- as.matrix(V)
    XV = mat %*% V %*% solve(t(V) %*% V) %*% t(V)
    Matrix::norm(XV, type = 'F')^2 / 
      Matrix::norm(mat, 'F')^2
  } else {
    # V <- as.matrix(V)
    ## total variance
    TV <- sum(svd(mat)$d ^ 2) 
    ## column length, replace 0 with 1
    CL <- sqrt(colSums(V ^ 2)) 
    CL <- ifelse(CL, CL, 0) 
    ## component scores
    S <- mat %*% (V / CL) 
    ## explained variance
    EV <- diag(qr.R(qr(S)) ^ 2) 
    sum(EV) / TV
    # CI = solve(crossprod(V))
    # sum(diag(V %*% CI %*% t(V) %*% mat %*% V %*% CI %*% t(V)))
  }
}

#' Cumulative Proportion of Variance Explained (CPVE)
#' 
#' Calculate the CPVE. 
#' @inheritParams pve
#' @return a `numeric` vector of length `ncol(V)`, the i-th value is the CPVE of the first i columns in `V`.
#' @seealso [pve]
#' @examples 
#' ## use the "swiss" data
#' ## find two sparse PCs 
#' s.sca <- sca(swiss, 2, gamma = sqrt(ncol(swiss))) 
#' ld <- loadings(s.sca)
#' cpve(as.matrix(swiss), ld)
#' 
#' @export
cpve <- function(mat, V, is.cov = FALSE) {
  sdev = apply(mat %*% V, 2, stats::sd) 
  ord <- order(sdev, decreasing = TRUE)
  pve <- rep(0, ncol(V))
  for (i in 1:ncol(V)) {
    pve[i] <- pve(mat, V[,ord[1:i]], is.cov)
  }
  pve
}

#' Matrix Column Distance
#'
#' Compute the distance between two matrices.
#' The distance between two matrices is defined as the sum of distances between column pairs. 
#' This function matches the columns of two matrices, such that the matrix distance 
#' (i.e., the sum of paired column distances) is minimized.
#' This is accomplished by solving an optimization over column permutation.
#' Given two matrices, `x` and `y`, find permutation p() that minimizes
#'      sum_i similarity(`x[,p(i)], y[,i]`),
#' where the `similarity()` can be "euclidean" distance, 1 - "cosine", or "maximum" difference (manhattan distance).
#' The solution is computed by [clue::solve_LSAP()].
#'
#' @param x,y `matrix` or `Matrix`, of the same number of rows. The columns of `x` and `y` will be scaled to unit length.
#' @param method distance measure, "maximum", "cosine", or "euclidean" are implemented.
#' @return a `list` of four components:
#' \item{dist}{`dist`, the distance matrix.}
#' \item{match}{`solve_LSAP`, the column matches.}
#' \item{value}{`numeric` vector, the distance between pairs of columns.} 
#' \item{method}{`character`, the distance measure used.} 
#' \item{nrow}{`integer`, the dimension of the input matrices, i.e., `nrow(x)`.}
#' @seealso [clue::solve_LSAP]
#' @examples 
#' x <- diag(4) 
#' y <- x + rnorm(16, sd = 0.05) # add some noise
#' y = t(t(y) / sqrt(colSums(y ^ 2))) ## normalize the columns
#' ## euclidian distance between column pairs, with minimal matches
#' dist.matrix(x, y, "euclidean") 
#' 
#' @export
dist.matrix = function(x, y, method = 'euclidean') {
  ## check
  stopifnot(dim(x) == dim(y))
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    y <- as.matrix(y)
  }
  x = t(t(x) / sqrt(colSums(x ^ 2)))
  y = t(t(y) / sqrt(colSums(y ^ 2)))
  
  # s <- svd(crossprod(x, y))
  # sum(1 - s$d ^ 2)
  
  if (method == "sine") {
    ## sin ^ 2
    FUN = function(x, y) 1 - crossprod(x, y)^2
  } else if (method == "euclidean") {
    ## euclidean 
    FUN = function(x, y) sqrt(sum((x - y) ^ 2))
  } else if (method == "maximum") {
    ## supremum norm (manhattan)
    FUN = function(x, y) max(abs(x - y))
  }
  
  dist = inner(x, y, FUN)
  match = clue::solve_LSAP(dist, maximum = FALSE)
  value = dist[cbind(seq_along(match), match)]

  list(dist = dist,
       match = match,
       value = value,
       method = method,
       nrow = nrow(x))
}

#' Matrix Distance
#' @inheritParams dist.matrix
#' @return `numeric`, the distance between two matrices.
distance <- function(x, y, method = "euclidean") {
  dist <- dist.matrix(x, y, method)
  sum(dist$value)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ rotation ------ 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Varimax Rotation
#'
#' Perform varimax rotation.
#' Flip the signs of columns so that the resulting matrix is positive-skewed.
#' 
#' @param x a `matrix` or `Matrix`.
#' @param rotate `character(1)`, rotation method. Two options are currently available: "varimax" (default) or "absmin" (see details).
#' @param normalize `logical`, whether to rows normalization should be done before and undone afterward the rotation (see details). 
#' @param flip `logical`, whether to flip the signs of the columns of estimates such that all columns are positive-skewed (see details).
#' @param eps `numeric` precision tolerance.
#' @includeRmd man/rotate.md details
#' @includeRmd man/normalize.md details
#' @includeRmd man/flip.md details
#' @return the rotated matrix of the same dimension as `x`.
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse PCA." 
#' @seealso [prs], [varimax]
#' @examples 
#' ## use the "swiss" data
#' fa <- factanal( ~., 2, data = swiss, rotation = "none")
#' rotation(loadings(fa))
#' @export
rotation = function(x, 
                    rotate = c("varimax", "absmin"),
                    normalize = FALSE, 
                    flip = TRUE,
                    eps = 1e-6) {
  ## check input
  stopifnot(is.array(x)) 
  rotate <- match.arg(rotate)
  
  if (rotate == "varimax") {
    rotate.fun <- varimax
  } else {
    rotate.fun <- absmin
  }
  r = rotate.fun(x, normalize = normalize, eps = eps)
  y = x %*% r$rotmat
  
  if (flip) {
    skewness <- colSums(y ^ 3)
    y <- t(t(y) * sign(skewness))
  }
  
  return(y)
}

#' The varimax criterion
#'
#' Calculate the varimax criterion
#' @param mat a `matrix` or `Matrix`.
#' @return a `numeric` of evaluated varimax criterion.
#' @references \href{https://en.wikipedia.org/wiki/Varimax_rotation}{Varimax rotation (Wikipedia)}
#' @examples 
#' ## use the "swiss" data
#' fa <- factanal( ~., 2, data = swiss, rotation = "none")
#' lds <- loadings(fa)
#' 
#' ## compute varimax criterion:
#' varimax.criteria(lds) 
#' 
#' ## compute varimax criterion (after the varimax rotation):
#' rlds <- rotation(lds, rotate = "varimax")
#' varimax.criteria(rlds)
#' 
#' @export
varimax.criteria = function(mat) {
  sum(apply(mat^2, 2, stats::var))
}

#' Varimax Rotation 
#' 
#' This is a re-implementation of [stats::varimax], 
#' which (1) adds a parameter for the maximum number of iterations, 
#' (2) sets the default `normalize` parameter to `FALSE`,  
#' (3) outputs the number of iteration taken, and 
#' (4) returns regular `matrix` rather than in `loadings` class.
#' @inheritParams stats::varimax
#' @param maxit `integer`, maximum number of iteration (default to 1,000). 
#' @return A list with three elements: 
#' \item{rotated}{the rotated matrix.} 
#' \item{rotmat}{the (orthogonal) rotation matrix.}
#' \item{n.iter}{the number of iterations taken.}
#' @seealso [stats::varimax]
varimax = function(x,
                   normalize = FALSE,
                   eps = 1e-05,
                   maxit = 1000L) {
  nc <- ncol(x)
  if (nc < 2) {
    warnings("Rotation of a single-column matrix is invalid.")
    return(list(loadings = x, 
                rotmat = matrix(1), 
                n.iter = 1))
  }
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
    if (d < dpast * (1 + eps))
      break
  }
  z <- x %*% TT
  if (normalize)
    z <- z * sc
  dimnames(z) <- dimnames(x)
  list(loadings = z, rotmat = TT, n.iter = i)
}

#' Absmin Criteria
#'
#' Calculate the absmin criteria.
#' This is a helper function for [absmin].
#' @inheritParams absmin
absmin.criteria = function(L) {
  sum(abs(L))
}

#' Gradient of Absmin Criterion 
#' 
#' This is a helper function for [absmin] and is not to be used directly by users.
#' @inheritParams absmin
#' @return a list required by `GPArotation::GPForth` for the absmin rotation.
#' @examples 
#' \dontrun{
#' ## NOT RUN
#' ## NOT for users to call.
#' }
#' @export
vgQ.absmin = function (L) {
  list(Gq = sign(L), 
       f = absmin.criteria(L), 
       Method = "absmin")
}

#' Absmin Rotation
#' 
#' Given a p x k matrix `x`, 
#' finds the orthogonal matrix (rotation) that minimizes the [absmin.criteria].
#' 
#' @param L a `matrix` or `Matrix`. 
#' @param Tmat `matrix`, initial rotation matrix.
#' @inheritParams stats::varimax
#' @inheritParams varimax
#' @return A list with three elements: 
#' \item{rotated}{the rotated matrix.} 
#' \item{rotmat}{the (orthogonal) rotation matrix.}
#' \item{n.iter}{the number of iteration taken.}
#' @seealso `GPArotation::GPForth`
absmin <- function(L, 
                   Tmat = diag(ncol(L)), 
                   normalize = FALSE,
                   eps = 1e-5, 
                   maxit = 1000L) {
  nc <- ncol(L)
  if (nc < 2) {
    warnings("Rotation of a single-column matrix is invalid.")
    return(list(loadings = L, 
                rotmat = matrix(1), 
                n.iter = 1))
  }
  
  res <- GPArotation::GPForth(L, 
                              Tmat = Tmat, 
                              method = "absmin",
                              normalize = normalize, 
                              eps = eps, 
                              maxit = maxit)
  list(loadings = res$loadings, 
       rotmat = res$Th, 
       n.iter = nrow(res$Table))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ shrinkage ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Calculate fractional exponent/power
#'
#' Calculate fractional exponent/power, `a^(num/den)`, where a could be negative.
#' @param a `numeric(1)`, base (could be negative).
#' @param num a positive `integer`, numerator of the exponent.
#' @param den a positive `integer`, denominator of the exponent.
#' @return `numeric`, the evaluated a^(num/den)
exp.frac <- function(a, num, den) {
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

#' Element-wise Matrix Norm
#' 
#' Compute element-wise matrix Lp-norm. 
#' This is a helper function to [shrinkage()].
#' @param mat a `matrix` or `Matrix`.
#' @param p `numeric(1)`, the p for defining the Lp norm.
#' @return `numeric(1)`, the absolute sum of all elements.
norm.Lp = function(mat, p = 1) {
  stopifnot(p >= 0 && is.numeric(p))
  if (p == 0) {
    sum(!!mat)
  } else if (is.infinite(p)) {
    max(abs(mat))
  } else {
    sum(abs(mat) ^ p) ^ (1 / p)
  }
}

#' Soft-thresholding 
#' 
#' Perform soft-thresholding given the cut-off value.
#' @param x any numerical `matrix` or `vector`. 
#' @param t `numeric`, the amount to soft-threshold, i.e., \eqn{sgn(x_{ij}) (|x_{ij}-t|)_+}. 
soft <- function (x, t) {
  sign(x) * pmax(abs(x) - t, 0)
}

#' Hard-thresholding 
#' 
#' Perform hard-thresholding given the cut-off value.
#' @param x any numerical `matrix` or `vector`. 
#' @param t `numeric`, the amount to hard-threshold, i.e., \eqn{sgn(x_{ij}) (|x_{ij}-t|)_+}. 
hard <- function (x, t) {
  x * (abs(x) >= t)
}

#' Shrinkage
#'
#' Shrink a matrix using soft-thresholding or hard-thresholding.
#'
#' @param x `matrix` or `Matrix`, to be threshold.
#' @param gamma `numeric`, the constraint of Lp norm, i.e. \eqn{\|x\|\le \gamma}. 
#' @param shrink `character(1)`, shrinkage method, either "soft"- (default) or "hard"-thresholding (see details). 
#' @param epsilon `numeric`, precision tolerance. This should be greater than `.Machine$double.eps`.
#' @details A binary search to find the cut-off value.
#' @includeRmd man/shrink.md details
#' @return a `list` with two components:
#' \item{matrix}{matrix, the matrix that results from soft-thresholding}
#' \item{norm}{numeric, the norm of the matrix after soft-thresholding. This value is close to constraint if using the second option.}
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse PCA." 
#' @seealso [prs]
#' @examples 
#' x <- matrix(1:6, nrow = 3)
#' shrink_x <- shrinkage(x, 1)
#' 
#' @export
shrinkage = function(x,
                     gamma,
                     shrink = c("soft", "hard"),
                     epsilon = 1e-11) {
  ## check input
  stopifnot(gamma >= 0)
  shrink <- match.arg(shrink)
  if (shrink == "soft") {
    shrink.fun <- soft
  } else {
    shrink.fun <- hard
  } 
  p <- 1 * (shrink == "soft") 
  
  ## no shrinkage needed
  if (norm.Lp(x, p) <= gamma) return(x)
  
  ## binary search
  lower = 0 ## whose Lp norm > gamma
  upper = max(abs(x)) ## cut-off where L1 norm = 0
  while (upper - lower > epsilon) {
    mid = (lower + upper) / 2
    x.mid = shrink.fun(x, mid)
    if (norm.Lp(x.mid, p) > gamma) {
      lower = mid
    } else
      upper = mid
  }
  
  return(x.mid)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ image filters ------ 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## re-implement whitening procedure

#' ffshift
#' 
#' Rearrange the matrix so that zero frequency component is in the middle of the matrix. 
#' This is similar to the `fftshift` function in MATLAB 
#' 
#' @param x a `matrix` of image. 
#' @references Olshausen, Bruno A., and David J. Field. "Emergence of simple-cell receptive field properties by learning a sparse code for natural images." Nature 381.6583 (1996): 607.
fftshift <- function(x) {
  if (class(x) != 'matrix') {
    stop('sorry, this class of input x is not supported yet')
  }
  rd2 <- floor(nrow(x)/2)
  cd2 <- floor(ncol(x)/2)
  
  ## Identify the first, second, third, and fourth quadrants 
  q1 <- x[1:rd2,1:cd2]
  q2 <- x[1:rd2,(cd2+1):ncol(x)]
  q3 <- x[(rd2+1):nrow(x),(cd2+1):ncol(x)]
  q4 <- x[(rd2+1):nrow(x),1:cd2]
  
  ## rearrange the quadrants 
  centered.t <- rbind(q4,q1)
  centered.b <- rbind(q3,q2)
  centered <- cbind(centered.b,centered.t)
  
  return(Re(centered))
}

#' Pre-whiten Image
#' 
#' Mimic/simulate the processing in retina, the lateral geniculate nucleus (LGN), and V1
#'      (1) scale all pixels to constant variance (e.g., 1).
#'      (2) combined low-pass/whitening filter 
#' @param img a `matrix` of square image. 
#' @references Olshausen, Bruno A., and David J. Field. "Emergence of simple-cell receptive field properties by learning a sparse code for natural images." Nature 381.6583 (1996): 607.
prewhiten = function(img) {
  N = dim(img)[1]
  img = img / stats::sd(img) ## scale to variance 1
  
  px = matlab::meshgrid((-N / 2):(N / 2 - 1), (-N / 2):(N / 2 - 1))
  rho = sqrt(px$x^2 + px$y^2)
  f0 = 0.4 * N
  flt = rho * exp(-(rho / f0) ^ 4) ## frequency space
  ff = stats::fft(as.matrix(img))
  img.flt = stats::fft(ff * fftshift(flt), inverse = TRUE) 
  img.flt <- imager::as.cimg(Re(img.flt))
  img.flt = 0.1 * img.flt / stats::sd(img.flt) ## sd = 0.1
}

#' Trim Image
#' 
#' Trim the fringes of image
#' @param img a `matrix` of square image. 
#' @param hem `integer`, the width of trimmed fringes.
trim.fringe = function(img, hem = 6) {
  N = dim(img)[1]
  mat <- as.matrix(img)[(hem + 1):(N - hem), (hem + 1):(N - hem)]
  imager::as.cimg(mat)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ------ misc ------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' #' Varimin Rotation
#' varimin = function(x,
#'                    normalize = FALSE,
#'                    eps = 1e-05,
#'                    maxit = 1000L) {
#'   nc <- ncol(x)
#'   if (nc < 2)
#'     return(x)
#'   if (normalize) {
#'     sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
#'     x <- x/sc
#'   }
#'   p <- nrow(x)
#'   TT <- diag(nc)
#'   d <- 0
#'   for (i in seq_len(maxit)) {
#'     z <- x %*% TT
#'     cm2 = colMeans(z^2)
#'     B <- t(x) %*% (z^3 - z %*% diag(cm2)) 
#'     sB <- La.svd(TT - B)
#'     TT <- sB$u %*% sB$vt
#'     dpast <- d
#'     d <- sum(sB$d)
#'     if (d < dpast * (1 + eps))
#'       break
#'   }
#'   z <- x %*% TT
#'   if (normalize)
#'     z <- z * sc
#'   dimnames(z) <- dimnames(x)
#'   list(loadings = z, rotmat = TT, iter = i)
#' }
#' 
#' #' Varimax2 Criteria
#' #'
#' #' Calculate the absmax2 criteria.
#' #' @export
#' varimax2.criteria = function(mat) {
#'   sum((apply(mat^2, 2, stats::var) - 2 * colMeans(mat ^ 2)) ^ 2)
#' }
#' 
#' #' Varimax2 Rotation
#' #' 
#' #' Find the orthogonal matrix (rotation) T that
#' #' maximizes the varimax2 objective (z = x %*% T), 
#' #' varimax2(z) = \sum_{j=1}^k \[\text{Var}(z_{\cdot j}^2) - 2\]^2 
#' #' @inheritParams varimax
#' varimax2 = function (x, normalize = FALSE, eps = 1e-05, maxit = 1000L)
#' {
#'   nc <- ncol(x)
#'   if (nc < 2)
#'     return(x)
#'   if (normalize) {
#'     sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
#'     x <- x/sc
#'   }
#'   p <- nrow(x)
#'   TT <- diag(nc)
#'   d <- 0
#'   for (i in seq_len(maxit)) {
#'     z <- x %*% TT
#'     cs2 = colSums(z^2)
#'     cs4 = colSums(z^4)
#'     B <- t(x) %*% (z^3 - 3 * z %*% diag(cs2) / p) %*% diag(cs4 - 3 * cs2^2 / p)
#'     sB <- La.svd(B)
#'     TT <- sB$u %*% sB$vt
#'     dpast <- d
#'     d <- sum(sB$d)
#'     if (d < dpast * (1 + eps))
#'       break
#'   }
#'   z <- x %*% TT
#'   if (normalize)
#'     z <- z * sc
#'   dimnames(z) <- dimnames(x)
#'   list(loadings = z, rotmat = TT, iter = i)
#' }
#' 
#' #' Varimax2 Rotation
#' #' 
#' #' This function invokes [GPForth].
#' varimax2 = function (L, Tmat = diag(ncol(L)), normalize = FALSE,
#'                      eps = 1e-5, maxit = 1000) {
#'   GPForth(L, Tmat = Tmat, method = "varimax2",
#'           normalize = normalize, al = 1e3,
#'           eps = eps, maxit = maxit)
#' }
#' 
#' ## FUN: gradient of varimax2
#' vgQ.varimax2 = function (L) {
#'   p <- nrow(L)
#'   cs2 = colSums(L^2)
#'   cs4 = colSums(L^4)
#'   Gq <- (L^3 - 3 * L %*% diag(cs2) / p) %*% diag(cs4 - 3 * cs2^2 / p)
#'   list(Gq = -Gq, f = -varimax2.criteria(L),
#'        Method = "varimax2")
#' }
#' 
#' GPForth = function (A, Tmat = diag(ncol(A)), normalize = FALSE, al = 1,
#'                     eps = 1e-5, maxit = 1000,
#'                     method = "varimax", methodArgs = NULL) {
#'   if ((!is.logical(normalize)) || normalize) {
#'     W <- NormalizingWeight(A, normalize = normalize)
#'     normalize <- TRUE
#'     A <- A/W
#'   }
#'   if (1 >= ncol(A))
#'     stop("rotation does not make sense for single factor models.")
#'   alpha = al
#'   L <- A %*% Tmat
#'   Method <- paste("vgQ", method, sep = ".")
#'   VgQ <- do.call(Method, append(list(L), methodArgs))
#'   G <- crossprod(A, VgQ$Gq)
#'   f <- VgQ$f
#'   Table <- NULL
#'   VgQt <- do.call(Method, append(list(L), methodArgs))
#'   for (iter in 0:maxit) {
#'     M <- crossprod(Tmat, G)
#'     S <- (M + t(M))/2 ## symmetric part of M, the skew-symmetric part is skm(M) = (M-t(M))/2
#'     Gp <- G - Tmat %*% S ## gradient of the Lagrangian
#'     s <- sqrt(sum(diag(crossprod(Gp)))) ## Frobenius norm
#'     Table <- rbind(Table, c(iter, f, log10(s), al))
#'     if (s < eps || al < eps * 1e-3)
#'       break
#'     al <- 2 * al
#'     for (i in 0:10) {
#'       X <- Tmat - al * G ## dont know why the package uses Gp here
#'       UDV <- svd(X)
#'       Tmatt <- UDV$u %*% t(UDV$v)
#'       L <- A %*% Tmatt
#'       VgQt <- do.call(Method, append(list(L), methodArgs))
#'       if (VgQt$f < (f - 0.5 * s^2 * al))
#'         break
#'       al <- al/2
#'     }
#'     Tmat <- Tmatt
#'     f <- VgQt$f
#'     G <- crossprod(A, VgQt$Gq)
#'   }
#'   convergence <- (s < eps)
#'   if ((iter == maxit) & !convergence)
#'     warning("convergence not obtained in GPForth. ", maxit,
#'             " iterations used.")
#'   if (normalize)
#'     L <- L * W
#'   dimnames(L) <- dimnames(A)
#'   r <- list(loadings = L, Th = Tmat, Table = Table, method = VgQ$Method,
#'             orthogonal = TRUE, convergence = convergence, Gq = VgQt$Gq, iter = iter)
#'   class(r) <- "GPArotation"
#'   r
#' }
#' 
#' NormalizingWeight = function (A, normalize = FALSE)
#' {
#'   if ("function" == mode(normalize))
#'     normalize <- normalize(A)
#'   if (is.logical(normalize)) {
#'     if (normalize)
#'       normalize <- sqrt(rowSums(A^2))
#'     else return(array(1, dim(A)))
#'   }
#'   if (is.vector(normalize)) {
#'     if (nrow(A) != length(normalize))
#'       stop("normalize length wrong in NormalizingWeight")
#'     return(array(normalize, dim(A)))
#'   }
#'   stop("normalize argument not recognized in NormalizingWeight")
#' }
#' 
#' Random.Start = function (k) {
#'   qr.Q(qr(matrix(rnorm(k * k), k)))
#' }