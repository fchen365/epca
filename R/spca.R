#  sparse.pca
#' Sparse PCA via SMD
#'
#' Perform sparse PCA via the sparse multivariate decomposition (SMD)
#'
#' @param A  matrix (or Matrix), data matrix
#' @param k  number of PCs to compute
#' @param lambda  penalized parameter, default is sqrt(nrow(A)) and sqrt(ncol(A))
#' @param input.type  'predictor' (original) or 'Gram' (covariance, when n >> p), the input matrix for SMD
#' @param centering  logical, whether to (TRUE) center columns of A
#' @param scaling logical, whether to (TRUE) scale columns of A
#' @param normalizing logical, whether to (TRUE) normalize rows before rotation, then scale back after
#' @param ordering logical or "sdev" or "pve", whether to (TRUE) reorder columns of loadings by decreasing P.V.E.
#' @param epsilon tolerance
#'
#' @return a list of
#'    $sdev - numeric vector, standard deviation of PC,
#'    $loadings - (sparse) matrix, loadings of PC,
#'    $center - logical (=centering),
#'    $scale - logical (=scaling),
#'    $n.obs - integer, sample size, = nrow(A),
#'    $n.iter - integer, number of iterations taken,
#'    $scores - matrix of n*k, the scores of the supplied data on the principal components.
#' @export
sparse.pca = function(A, k = 5,
                      lambda = NULL,
                      input.type = 'predictor',
                      centering = T, scaling = F, normalizing = F, ordering = T,
                      max.iter = 1e3, epsilon = 1e-5,
                      quiet = T) {

  ## input original data matrix
  if (input.type == 'predictor') {
    ## centering and scaling
    A1 = scale(x = A, center = centering, scale = scaling)
    S = smd(A1, k = k,
            lambda = lambda,
            side = 'right',
            centering = F, scaling = F, normalizing = normalizing,
            max.iter = max.iter, epsilon = epsilon,
            quiet = quiet)
    loadings = as.matrix(S$Y)
    scores = S$Z %*% S$B
  }

  ## input Gram matrix, i.e. Gram(A) = crossprod(A)
  if (input.type == 'Gram') {
    A1 = rootmatrix(A)
    S = smd(A1, k = k,
            lambda = lambda,
            side = 'both',
            centering = F, scaling = F,
            max.iter = max.iter, epsilon = epsilon,
            quiet = quiet)
    loadings = (S$Z + S$Y) / 2
    scores = A1 %*% S$Y
  }

  n.iter = S$n.iter ## number of iterations
  sdev = apply(scores, 2, sd) ## standard deviation of scores (column)
  pve = apply(loadings, 2, var.explain, mat = A1, type = input.type) ## percent variance explained

  ## permute U and V s.t. PVE's are decreasing
  if (ordering == T | ordering == 'pve') {
    ord = order(pve, decreasing = T)
  } else if (ordering == 'sdev') {
    ord = order(sdev, decreasing = T)
  } else {
    ord = seq_len(k)
  }
  loadings = loadings[,ord,drop = F]
  dimnames(loadings) <- list(colnames(A1), paste("PC", 1:k, sep = ""))
  scores = scores[,ord]
  sdev = sdev[ord]
  ## re-calculate PVE by groups to avoid inflation (since loadings are not exactly orthogonal)
  for (i in seq_len(k)) {
    pve[i] = var.explain(A1, loadings[,1:i], type = input.type)
  }
  pve = c(pve[1], diff(pve))

  list(sdev = sdev, loadings = loadings, pve = pve, scores = scores,
       center = centering, scale = scaling, n.iter = n.iter,
       n.obs = ifelse(input.type == "predictor", nrow(A1), NA))
}

## test ----
# A = matrix(runif(3000), 60, 50)
# A = scale(A, center = T, scale = F) ## WLOG assume A is centered by col
# S1 = sparse.pca(A, k = 7, input.type = 'ori')
# S2 = sparse.pca(A, k = 7, input.type = 'cov')
