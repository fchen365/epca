
#' Polar-Rotation-Shrinkage
#' 
#' This function is a helper function of [sma()]. 
#' It performs polar docomposition, orthgonal rotation, and soft-thresholding shrinkage in order. 
#' The three steps together enable sparse estimates of the SMA and SCA. 
#' 
#' @param X,Z.hat the matrix product `crossprod(X, Z.hat)` is the input.
#' @param gamma `numeric`, the sparsity parameter.
#' @param rotate `character(1)`, rotation method, "varimax" (default) or "absmin".
#' @param shrink `character(1)`, shrinkage method, "soft" (default) or "hard". 
#' @param normalize `logical`, whether to normalize rows before rotation, then scale back after.
#' @param order `logical`, whether to re-order decreasingly the columns of estimates by PVE.
#' @param flip `logical`, whether to flip the signs of the columns of estimates such that all columns are positive-skewed, i.e., the sums of cubic elements are positive.
#' @param epsilon `numeric`, tolerance of convergence precision (dafault to 0.00001).
#' 
#' @return a `matrix` of the same dimension as `A`, sparse estimate.
prs <- function(X, Z.hat,
                gamma,
                rotate, 
                shrink,
                normalize,
                order,
                flip,
                epsilon = 1e-5) {
  ## input
  A <- crossprod(X, Z.hat)

  ## three steps  
  Y.tilde <- polar(A)
  Y.star = rotation(Y.tilde, 
                    rotate = rotate,
                    normalize = normalize, 
                    flip = flip, 
                    eps = epsilon / sqrt(nrow(A)))  
  Y.hat = shrinkage(Y.star, gamma, shrink = shrink)

  ## re-order by std dev of scores
  if (order) {
    sdev <- apply(X %*% Y.hat, 2, sd) 
    ord <- order(sdev, decreasing = T)
    Y.hat <- Y.hat[,ord,drop=F]
  }

  Y.hat
}

#' Sparse matrix approximation 
#'
#' Perform the sparse matrix approximation (SMA) of a data matrix X as three components: Z B Y', 
#' where Z and Y are sparse, and B is low-rank but not necessarily diagonal.
#'
#' @param A `matrix` or `Matrix` to be analyzed. 
#' @param k `integer`, rank of approximation.
#' @param gamma `numeric(2)`, sparsity parameters. If `gamma` is `numeric(1)`, it is used for both left and right sparsity component (i.e, Z and Y). If absent, the two parameters are set as (default): √(nk) and √(pk) for Z and Y respoectively, where n x p is the dimension of `A`.
#' @param center `logical`, whether to center columns of `A` (default to `TRUE`).
#' @param scale `logical`, whether to scale columns of `A` (default to `TRUE`).
#' @param max.iter `integer`, maximum number of iteration (default to 1,000). 
#' @param quiet `logical`, whether to mute the process report (default to `TRUE`)
#' @inheritParams prs
#'
#' @return a `list` of four elements:
#' \item{Z}{an n x k `matrix` , sparse component Z of `A`'srow space.}
#' \item{B}{a k x k `matrix`, the middle B matrix. This is like the scores of the SMA, because the objective of SMA is to maximize the Frobenius norm of B.}
#' \item{Y}{an n x k `matrix`, sparse component Y of `A`'s column space.}
#' \item{score}{optimal objective value.}
#' \item{n.iter}{an `integer`, the number of iteration taken.}

#' @export
sma = function(A,
               k = min(5, dim(A)),
               gamma = NULL,
               rotate = "varimax", 
               shrink = "soft",
               center = F,
               scale = F,
               normalize = F,
               order = F,
               flip = F,
               max.iter = 1e3,
               epsilon = 1e-5,
               quiet = T) {
  
  ## check gamma
  if (length(gamma) == 1) {
    warning('Same sparsity parameters for both sides.')
    gamma = rep(gamma, 2)
  } 
  if (length(gamma) == 0) 
    gamma = sqrt(dim(A)) * sqrt(k)
  if (length(gamma) > 2) 
    stop('Too many sparsity parameters.')
  if (gamma[1] < k || gamma[1] > k * sqrt(nrow(A))) 
    message("Improper sparsity parameter (gamma) for Z.") 
  if (gamma[2] < k || gamma[2] > k * sqrt(ncol(A))) 
    message("Improper sparsity parameter (gamma) for Y.") 
  names(gamma) = c('Z', 'Y')
  
  ## center and scale
  A = scale(x = A, center = center, scale = scale)
  
  ## initialize
  S = RSpectra::svds(A, k)
  # S = irlba::irlba(A, k, tol = 1e-10) 
  Z = S$u
  B = diag(S$d)
  Y = S$v
  score = sqrt(sum(S$d^2))
  diff = c(Z = Inf, Y = Inf)
  
  for (n.iter in seq_len(max.iter)) {
    ## update Z
    Z.new <- prs(t(A), Y, 
                 gamma = gamma["Z"],
                 rotate = rotate, 
                 shrink = shrink,
                 normalize = normalize, 
                 order = order, 
                 flip = flip,
                 epsilon = epsilon)
    diff['Z'] = distance(Z, Z.new, "maximum")
    Z = Z.new
    
    ## update Y
    Y.new = prs(A, Z, 
                gamma = gamma["Y"],
                rotate = rotate, 
                shrink = shrink,
                normalize = normalize, 
                order = order, 
                flip = flip,
                epsilon = epsilon)
    diff['Y'] = distance(Y, Y.new, "maximum")
    Y = Y.new
    
    ## update B and score
    B <- t(Z.new) %*% A %*% Y.new
    # score.new = norm(B, 'F')
    # diff['score'] = score.new - score
    score = norm(B, 'F')
    
    ## report progress
    if (!quiet) {
      cat("Iter", n.iter, "score =", format(score, digit = 7), "and the differences:\n") 
      print(diff, digit = 4, quote = F)
    }
    
    ## convergence?
    if (max(diff) < epsilon) break
  }
  
  list(Z = Z, B = B, Y = Y, score = score, n.iter = n.iter)
}


#' Sparse component analysis
#'
#' Perform sparse component analysis (SCA). 
#' SCA is a method of sparse PCA where the loadings are column sparse. 
#'
#' @param gamma `numeric(1)`, sparsity parameter, default to √(pk), where n x p is the dimension of `A`.
#' @param is.cov  `logical`, whether the `A` is a covariance matrix or Gram matrix (i.e., `crossprod(X)`). This function presumes that `A` is *not* covariance matrix.
#' @inheritParams prs
#' @inheritParams sma
#'
#' @return a list of
#' \item{loadings}{`matrix`, sparse loadings of PCs.}
#' \item{scores}{an n x k `matrix`, the component scores.}
#' \item{sdev}{a `numeric` vector of length `k`, standard deviation of each columns of scores. These may not sum to exactly 1 because of a slight loss of orthogonality.}
#' \item{pve}{a `numeric` vector of length `k`, cumulative proportion of variance in `A` explained by the top PCs.}
#' \item{center}{`logical`, this records the `center` parameter.}
#' \item{scale}{`logical`, this records the `scale` parameter.}
#' \item{n.iter}{`integer`, number of iteration tabke.}
#' \item{n.obs}{`integer`, sample size, that is, `nrow(A)`.}
#' 
#' @export
sca = function(A, 
               k = min(5, dim(A)), 
               gamma = NULL, 
               is.cov = F,
               rotate = "varimax", 
               shrink = "soft",
               center = T, 
               scale = T, 
               normalize = F, 
               order = T,
               flip = T,
               max.iter = 1e3, 
               epsilon = 1e-5, 
               quiet = T) {
  ## check gamma
  if (!length(gamma)) 
    gamma = sqrt(ncol(A)) * sqrt(k)
  stopifnot(length(gamma) == 1) 
  if (gamma < k || gamma > k * sqrt(ncol(A))) 
    message("Improper sparsity parameter (gamma) for Z.") 
  
  ## covariance or Gram matrix 
  if (is.cov) {
    A = rootmatrix(A)
    gamma2 <- rep(gamma, 2) 
  } else {
    gamma2 <- c(k * sqrt(nrow(A)), gamma) ## for SMA
  }
  
  ## center and scale
  A = scale(x = A, center = center, scale = scale)
  
  S = sma(A, k = k, 
          gamma = gamma2, 
          rotate = rotate, 
          shrink = shrink,
          center = F, 
          scale = F, 
          normalize = normalize,
          order = order,
          flip = flip,
          max.iter = max.iter, 
          epsilon = epsilon, 
          quiet = quiet)
  loadings = S$Y
  rownames(loadings) <- colnames(A)
  colnames(loadings) <- paste0("PC", 1:k)
  scores = A %*% S$Y
  rownames(scores) <- rownames(A)
  colnames(scores) <- paste0("PC", 1:k)
  
  n.iter = S$n.iter ## number of iterations
  sdev = apply(scores, 2, sd) ## standard deviation of scores (column)
  pve = cpve(A, loadings)
  
  list(loadings = loadings, 
       scores = scores, 
       pve = pve, 
       sdev = sdev, 
       center = center, 
       scale = scale, 
       n.iter = n.iter, 
       n.obs = ifelse(is.cov, NA, nrow(A)))
}

# ## ------ test ------
# n = 1500
# p = 50
# k = 15 ## number of PCs to calculate
# S = matrix(runif(n * k), n, k)
# Y = shrinkage(polar(matrix(runif(p * k), p, k)), 
#               sqrt(p * k))$matrix
# A = tcrossprod(S, Y) + matrix(rnorm(n * p, sd = .02), n, p)
# 
# ## SMA
# s.sma = sma(A, k = k, quiet = F)
# 
# ## SCA 
# s.sca = sca(A, k = k, quiet = F)

