## ------------------------------------------------------------
## what: helper functions
## ------------------------------------------------------------

## ------------------------------------------------------------
## numerical
## ------------------------------------------------------------
## FUN: pretty permutation of columns
##      s.t. k-th column gives k-th cluster (ad-hoc)
## INPUT:
##    x - matrix, nNode x nBlock, non-negative, higher value means more likely
##    s - integer, order of non-linear
permColumn = function(x, s = 2) {
  x[,order(colMeans(row(x) * abs(x)^s))]
}

## FUN: calculate a^(num/den)
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

## calculate (1,1)-norm
norm11 = function(mat) {
  sum(abs(mat))
}

#' evaluate soft-thresholding function given cut-off
#' @return $value and $matrix after soft-thresholding
evalThres = function(mat, cutoff, norm.fun = norm11) {
  ## mat - matrix or Matrix, to be soft-thresholded
  ## cutoff - soft cutoff
  ## norm.fun - matrix norm for $value, after soft-thresholding
  sign.mat = sign(mat)
  abs.mat = abs(mat)
  res = list()
  res$matrix = sign.mat * pmax(abs.mat - cutoff, 0)
  res$value = norm.fun(res$matrix)
  res
}

## use binary search to find sup(threshold) s.t.
## evalThres(mat, threshold) <= gamma
## return the matrix after thresholding
softThres = function(mat, gamma, epsilon = 2*.Machine$double.eps) {
  ## mat - matrix or Matix
  ## gamma - numeric, constraint parameter: sum_ij (mat_ij) <= gamma
  ## epsilon - tolerance for precision

  ## no thresholding
  res = evalThres(mat, 0)
  if (res$value <= gamma)
    return(res$matrix)

  lower = 0 ## whose res$value > gamma
  upper = max(abs(mat)) ## whose res$value = 0
  while (upper - lower > epsilon) {
    mid = (lower + upper) / 2
    res = evalThres(mat, mid)
    if (res$value > gamma) {
      lower = mid
    } else
      upper = mid
  }

  res$matrix
}

## FUN: takes a matrix X, of dimensions n by p, and
##      calculate polar components, two matrices P and Z, such that X = UP,
##      assuming X has full col rank
## OUTPUT: U part
polar = function(X) {
  if (nrow(X) < ncol(X))
    stop('Input should be a tall (or at least square) matrix.')
  if (rankMatrix(X) < ncol(X))
    stop('X is not full rank.')
  S = svd(X)
  U.part <- S$u %*% t(S$v)
  # P.part <- S$v %*% diag(S$d) %*% t(S$v)
  return(U.part)
}

## FUN: find X from the Gram matrix X^T X, i.e. sqrt(crossprod(X))
rootmatrix <- function(x) {
  x.eigen <- eigen(x)
  d <- x.eigen$values
  d <- (d + abs(d))/2
  v <- x.eigen$vectors
  root = v %*% diag(sqrt(d)) %*% t(v)
  dimnames(root) <- dimnames(x)
  return(root)
}

## FUN: rotate X from right via varimax
## OUTPUT: the rotated matrix
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

## FUN: calculate custom inner product Z = <X, Y>
##      Z[i,j] = FUN(X[,i], Y[,j])
inner = function(X, Y, FUN = crossprod, ...) {
  t(apply(X, 2, function(x) {
    apply(Y, 2, function(y) FUN(x, y))
  }))
}

## FUN: proportion of variance explained by V
## INPUT:
##    mat - the original data matrix X or cov(X) = crossprod(X) / (nrow(X) - 1)
##    V - loadings
##    type - the input matrix for SMD (not used at the moment)
##           'predictor' (original data matrix) or 'Gram' (useful when n >> p),
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
## FUN: permute columns of X.new, p(i), minimizes
##      sum_i similarity(X.new[,p(i)], X[,i])
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


## FUN: given theta matrix finds cluster label
## INPUT:
##    x - matrix, nNode x nBlock, non-negative, higher value means more likely
##    ties.method - a character string specifying how ties are treated,
##                  see base::rank, "random", "first", "last" are allowed
## OUTPUT: integer vector %in% 1:ncol(x)
labelCluster <- function(x, ties.method = "random") {
  x1 <- t(apply(-abs(x), 1, rank, ties.method = ties.method))
  apply(x1, 1, which.min)
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
## FUN: calculate varimax criteria
varimax.criteria = function(mat) {
  sum(apply(mat^2, 2, var))
}

# varimax = Varimax
varimax = function (x, normalize = F, eps = 1e-05, maxit = 1000L)
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
  list(loadings = z, rotmat = TT, iter = i)
}

varimin = function (x, normalize = F, eps = 1e-05, maxit = 1000L)
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
    cm2 = colMeans(z^2)
    B <- t(x) %*% (z^3 - z %*% diag(cm2))
    sB <- La.svd(TT - B)
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

## FUN: calculate varimax criteria
varimax2.criteria = function(mat) {
  sum((apply(mat^2, 2, var) - 2 * colMeans(mat^2))^2)
}


## FUN: given x \in \R^{p \times k}, find the orthogonal matrix (rotation) T that
##      maximizes the varimax2 objective (z = x %*% T),
##      varimax2(z) = \sum_{j=1}^k \[\text{Var}(z_{\cdot j}^2) - 2\]^2
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

## FUN: calculate absmin criteria
absmin.criteria = function(mat) {
  sum(abs(mat))
}

## FUN: given x \in \R^{p \times k}, find the orthogonal matrix (rotation) T that
##      minimizes the absmin objective (z = x %*% T),
##      absmin_criteria(z) = \sum_{i=1}^p \sum_{j=1}^k |z_ij|
absmin = function (L, Tmat = diag(ncol(L)), normalize = FALSE,
                   eps = 1e-5, maxit = 1000) {
  GPForth(L, Tmat = Tmat, method = "absmin",
          normalize = normalize, ## al = 0.1,
          eps = eps, maxit = maxit)
}

## FUN: gradient of absmin
vgQ.absmin = function (L) {
  list(Gq = sign(L), f = absmin.criteria(L), Method = "absmin")
}

# varimax2 = function (L, Tmat = diag(ncol(L)), normalize = FALSE,
#                      eps = 1e-5, maxit = 1000) {
#   GPForth(L, Tmat = Tmat, method = "varimax2",
#           normalize = normalize, al = 1e3,
#           eps = eps, maxit = maxit)
# }
#
# ## FUN: gradient of varimax2
# vgQ.varimax2 = function (L) {
#   p <- nrow(L)
#   cs2 = colSums(L^2)
#   cs4 = colSums(L^4)
#   Gq <- (L^3 - 3 * L %*% diag(cs2) / p) %*% diag(cs4 - 3 * cs2^2 / p)
#   list(Gq = -Gq, f = -varimax2.criteria(L),
#        Method = "varimax2")
# }

GPForth = function (A, Tmat = diag(ncol(A)), normalize = FALSE, al = 1,
                    eps = 1e-5, maxit = 1000,
                    method = "varimax", methodArgs = NULL) {
  if ((!is.logical(normalize)) || normalize) {
    W <- NormalizingWeight(A, normalize = normalize)
    normalize <- TRUE
    A <- A/W
  }
  if (1 >= ncol(A))
    stop("rotation does not make sense for single factor models.")
  alpha = al
  L <- A %*% Tmat
  Method <- paste("vgQ", method, sep = ".")
  VgQ <- do.call(Method, append(list(L), methodArgs))
  G <- crossprod(A, VgQ$Gq)
  f <- VgQ$f
  Table <- NULL
  VgQt <- do.call(Method, append(list(L), methodArgs))
  for (iter in 0:maxit) {
    M <- crossprod(Tmat, G)
    S <- (M + t(M))/2 ## symmetric part of M, the skew-symmetric part is skm(M) = (M-t(M))/2
    Gp <- G - Tmat %*% S ## gradient of the Lagrangian
    s <- sqrt(sum(diag(crossprod(Gp)))) ## Frobenius norm
    Table <- rbind(Table, c(iter, f, log10(s), al))
    if (s < eps || al < eps * 1e-3)
      break
    al <- 2 * al
    for (i in 0:10) {
      X <- Tmat - al * G ## dont know why the package uses Gp here
      UDV <- svd(X)
      Tmatt <- UDV$u %*% t(UDV$v)
      L <- A %*% Tmatt
      VgQt <- do.call(Method, append(list(L), methodArgs))
      if (VgQt$f < (f - 0.5 * s^2 * al))
        break
      al <- al/2
    }
    Tmat <- Tmatt
    f <- VgQt$f
    G <- crossprod(A, VgQt$Gq)
  }
  convergence <- (s < eps)
  if ((iter == maxit) & !convergence)
    warning("convergence not obtained in GPForth. ", maxit,
            " iterations used.")
  if (normalize)
    L <- L * W
  dimnames(L) <- dimnames(A)
  r <- list(loadings = L, Th = Tmat, Table = Table, method = VgQ$Method,
            orthogonal = TRUE, convergence = convergence, Gq = VgQt$Gq, iter = iter)
  class(r) <- "GPArotation"
  r
}

NormalizingWeight = function (A, normalize = FALSE)
{
  if ("function" == mode(normalize))
    normalize <- normalize(A)
  if (is.logical(normalize)) {
    if (normalize)
      normalize <- sqrt(rowSums(A^2))
    else return(array(1, dim(A)))
  }
  if (is.vector(normalize)) {
    if (nrow(A) != length(normalize))
      stop("normalize length wrong in NormalizingWeight")
    return(array(normalize, dim(A)))
  }
  stop("normalize argument not recognized in NormalizingWeight")
}

Random.Start = function (k) {
  qr.Q(qr(matrix(rnorm(k * k), k)))
}


