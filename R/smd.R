#  smd
#' Sparse multivariate decomposition (SMD)
#'
#' Perform sparse multivariate decomposition (SMD): A = Z B Y^T
#'
#' @param A matrix or Matrix to be decomposed, of dimensions n x p
#' @param k integer, rank of approximation
#' @param lambda numeric, penalized parameter, default is √n√k (for Z) and/or √p√k (for Y)
#' @param centering logical, whether to (TRUE) center columns of A
#' @param scaling logical, whether to (TRUE) scale columns of A
#' @param normalizing logical, whether to (TRUE) normalize rows before rotation, then scale back after
#' @param epsilon numeric, tolerance of precision
#'
#' @return a list() containing
#' $Z - matrix n x k, sparse left singular vector (in column)
#' $B - matrix k x k, Z^T A Y
#' $Y - matrix n x k, sparse right singluar vector (in column)
#' $n.iter - integer, number of iteration taken
#' @export
smd = function(A, k = 5,
                lambda = NULL,
                side = c('both', 'left', 'right'),
                centering = F, scaling = F, normalizing = F,
                max.iter = 1e3, epsilon = 1e-5,
                quiet = T) {

  ## arguments check
  if (side == 'both') {
    if (length(lambda) == 1) {
      warning('Same sparsity parameters for both sides.')
      lambda = rep(lambda, 2)
    }
    if (length(lambda) == 0)
      lambda = sqrt(dim(A)) * sqrt(k)
    if (length(lambda) > 2)
      stop('Too many sparsity parameters.')
    lambda = setNames(lambda, c('left', 'right'))
  }
  if (side %in% c('left', 'right')) {
    if (length(lambda) > 1)
      stop('Too many sparsity parameters.')
    if (length(lambda) == 0)
      lambda = ifelse(side == 'left', sqrt(nrow(A)), sqrt(ncol(A))) * sqrt(k)
    lambda = setNames(lambda, side)
  }

  ## centering and scaling
  A = scale(x = A, center = centering, scale = scaling)

  ## No rotation, get result from PMD or SPC
  if (k == 1 && FALSE) {
    if (side == 'both') {
      S = PMD(A, type = 'standard',
              sumabsu = min(lambda['left'], sqrt(nrow(A))),
              sumabsv = min(lambda['right'], sqrt(nrow(A))),
              niter = max.iter, K = 1, trace = F)
    } else if (side == 'right') {
      S = SPC(A, sumabsv = min(lambda['right'], sqrt(nrow(A))),
              niter = max.iter, K = 1, trace = F)
    } else {
      S = SPC(t(A), sumabsv = min(lambda['left'], sqrt(ncol(A))),
              niter = max.iter, K = 1, trace = F)
    }
    return(list(Z = S$u,
                B = S$d,
                Y = S$v,
                n.iter = NULL))
  }

  ## initialize
  if (k < 0.5 * min(dim(A))) {
    S = irlba::irlba(A, nu = k, nv = k)
  } else {
    S = svd(A, k, k)
  }
  U = S$u
  Z = rotation(U, normalize = normalizing, eps = 0.1 * epsilon)
  if (side != 'right') Z = softThres(Z, lambda['left'])
  V = S$v
  Y = rotation(V, normalize = normalizing, eps = 0.1 * epsilon)
  if (side != 'left') Y = softThres(Y, lambda['right'])
  obj = norm(t(Z) %*% A %*% Y, 'F')
  obj.best = -Inf
  n.iter.best = 0

  converge = c(Z = Inf, Y = Inf, Obj = Inf)
  n.iter = 0
  while (max(converge[c('Z', 'Y')]) > epsilon &&
         n.iter < n.iter.best + max(1e2, max.iter/10) &&
         n.iter < max.iter) {
    n.iter = n.iter + 1

    ## given Y, update Z
    U = polar(A %*% Y)
    Z.new = rotation(U, normalize = normalizing, eps = 0.1 * epsilon)
    if (side != 'right')
      Z.new = softThres(Z.new, lambda['left'])
    obj.new = norm(t(Z.new) %*% A %*% Y, 'F')
    if (converge['Z'] < epsilon || obj.new < obj) {
      Z.new = matchCols(Z.new, Z)
    }
    converge['Z'] = norm(Z - Z.new, type = 'M')
    Z = Z.new

    ## given Z, update Y
    V = polar(t(A) %*% Z)
    Y.new = rotation(V, normalize = normalizing, eps = 0.1 * epsilon)
    if (side != 'left')
      Y.new = softThres(Y.new, lambda['right'])
    obj.new = norm(t(Z) %*% A %*% Y.new, 'F')
    if (converge['Y'] < epsilon || obj.new < obj) {
      Y.new = matchCols(Y.new, Y)
    }
    converge['Y'] = norm(Y - Y.new, type = 'M')
    Y = Y.new

    ## update obj value
    obj.new = norm(t(Z) %*% A %*% Y, 'F')
    converge['Obj'] = obj.new - obj
    obj = obj.new
    ## obj may decrease when close to converging
    if (obj > obj.best) {
      Y.best = Y
      Z.best = Z
      B.best = t(Z) %*% A %*% Y
      obj.best = obj
      n.iter.best = n.iter
    }

    if ((n.iter < 6 || !n.iter %% round(max.iter/50)) && !quiet) {
      message('Iterate ', n.iter, '/', max.iter, ', obj value = ', format(obj, digit = 7),
              ',\ndiff{Z, Y, Obj} = ', paste0(format(converge, digit = 3), collapse = ', '))
    }
  }

  # if (max(converge[c('Z', 'Y')]) <= epsilon) {
  #   list(Z = Z, B = B, Y = Y, n.iter = n.iter, obj = obj)
  # } else
  #   list(Z = Z.best, B = B.best, Y = Y.best, n.iter = n.iter.best, obj = obj.best)
  list(Z = Z.best, B = B.best, Y = Y.best, n.iter = n.iter, obj = obj.best)
}

