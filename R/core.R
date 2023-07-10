

#' Polar-Rotate-Shrink
#'
#' This function is a helper function of [sma()].
#' It performs polar docomposition, orthogonal rotation, and soft-thresholding shrinkage in order.
#' The three steps together enable sparse estimates of the SMA and SCA.
#'
#' @param x,z.hat the matrix product `crossprod(x, z.hat)` is the actual Polar-Rotate-Shrink object. `x` and `z.hat` are input separatedly because the former is additionally used to compute the proportion of variance explained, in the case when `order = TRUE`.
#' @param gamma `numeric`, the sparsity parameter.
#' @param order `logical`, whether to re-order the columns of the estimates (see Details below).
#' @param epsilon `numeric`, tolerance of convergence precision (default to 0.00001).
#' @inheritParams rotation
#' @inheritParams shrinkage
#' @includeRmd man/rotate.md details
#' @includeRmd man/shrink.md details
#' @includeRmd man/normalize.md details
#' @includeRmd man/order.md details
#' @includeRmd man/flip.md details
#' @return a `matrix` of the sparse estimate, of the same dimension as `crossprod(x, z.hat)`.
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse Principal Component Analysis."
#' @seealso [sma], [sca], [polar], [rotation], [shrinkage]
prs <- function(x,
                z.hat,
                gamma,
                rotate,
                shrink,
                normalize,
                order,
                flip,
                epsilon) {
  ## input
  a <- crossprod(x, z.hat)

  ## three steps
  y.tilde <- polar(a)
  y.star = rotation(
    y.tilde,
    rotate = rotate,
    normalize = normalize,
    flip = flip,
    eps = epsilon / sqrt(nrow(a))
  )
  y.hat = shrinkage(y.star, gamma, shrink = shrink)
  # y.hat <- t(t(y.hat) / sqrt(colSums(y.hat ^ 2)))

  ## re-order by std dev of scores
  if (order) {
    sdev <- colSums((x %*% y.hat) ^ 2)
    ord <- order(sdev, decreasing = TRUE)
    y.hat <- y.hat[, ord, drop = FALSE]
  }

  y.hat
}

#' Sparse Matrix Approximation
#'
#' Perform the sparse matrix approximation (SMA) of a data matrix `x` as three multiplicative components: `z`, `b`, and `t(y)`,
#' where `z` and `y` are sparse, and `b` is low-rank but not necessarily diagonal.
#'
#' @param x `matrix` or `Matrix` to be analyzed.
#' @param k `integer`, rank of approximation.
#' @param gamma `numeric(2)`, sparsity parameters. If `gamma` is `numeric(1)`, it is used for both left and right sparsity component (i.e, `z` and `y`). If absent, the two parameters are set as (default): `sqrt(nk)` and `sqrt(pk)` for `z` and `y` respectively, where n x p is the dimension of `x`.
#' @param center `logical`, whether to center columns of `x` (see [scale()]).
#' @param scale `logical`, whether to scale columns of `x` (see [scale()]).
#' @param max.iter `integer`, maximum number of iteration (default to 1,000).
#' @param quiet `logical`, whether to mute the process report (default to `TRUE`)
#' @inheritParams prs
#' @includeRmd man/rotate.md details
#' @includeRmd man/shrink.md details
#' @includeRmd man/normalize.md details
#' @includeRmd man/order.md details
#' @includeRmd man/flip.md details
#'
#' @return an `sma` object that contains:
#' \item{`z`, `b`, `t(y)`}{the three parts in the SMA.
#' `z` is a sparse n x k `matrix` that contains the row components (loadings).
#' The row names of `z` inherit the row names of `x`.
#' `b` is a k x k `matrix` that contains the scores of SMA;
#' the Frobenius norm of `b` equals to the total variance explained by the SMA.
#' `y` is a sparse n x k `matrix`that contains the column components (loadings).}
#' The row names of `y` inherit the column names of `x`.
#' \item{score}{the total variance explained by the SMA.
#' This is the optimal objective value obtained.}
#' \item{n.iter}{`integer`, the number of iteration taken.}
#'
#' @seealso [sca], [prs]
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse Principal Component Analysis."
#' @examples
#' ## simulate a rank-5 data matrix with some additive Gaussian noise
#' n <- 300
#' p <- 50
#' k <- 5 ## rank
#' z <- shrinkage(polar(matrix(runif(n * k), n, k)), sqrt(n))
#' b <- diag(5) * 3
#' y <- shrinkage(polar(matrix(runif(p * k), p, k)), sqrt(p))
#' e <- matrix(rnorm(n * p, sd = .01), n, p)
#' x <- scale(z %*% b %*% t(y) + e)
#'
#' ## perform sparse matrix approximation
#' s.sma <- sma(x, k)
#' s.sma
#'
#' @export
sma = function(x,
               k = min(5, dim(x)),
               gamma = NULL,
               rotate = c("varimax", "absmin"),
               shrink = c("soft", "hard"),
               center = FALSE,
               scale = FALSE,
               normalize = FALSE,
               order = FALSE,
               flip = FALSE,
               max.iter = 1e3,
               epsilon = 1e-5,
               quiet = TRUE) {
  ## check gamma
  if (length(gamma) == 1) {
    warning('Same sparsity parameters for both sides.')
    gamma = rep(gamma, 2)
  }
  if (length(gamma) == 0)
    gamma = sqrt(dim(x)) * sqrt(k)
  if (length(gamma) > 2)
    stop('Too many sparsity parameters.')
  if (gamma[1] < k || gamma[1] > k * sqrt(nrow(x)))
    message("Improper sparsity parameter (gamma) for 'z'.")
  if (gamma[2] < k || gamma[2] > k * sqrt(ncol(x)))
    message("Improper sparsity parameter (gamma) for 'y'.")
  names(gamma) = c('z', 'y')

  ## center and scale
  x = scale(x = x,
            center = center,
            scale = scale)

  ## initialize
  # s = RSpectra::svds(x, k)
  s = irlba::irlba(x, k, tol = 1e-10)
  z = s$u
  b = diag(s$d)
  y = s$v
  score = sqrt(sum(s$d ^ 2))
  diff = c(z = Inf, y = Inf)

  for (n.iter in seq_len(max.iter)) {
    ## update z
    z.new <- prs(
      t(x),
      y,
      gamma = gamma["z"],
      rotate = rotate,
      shrink = shrink,
      normalize = normalize,
      order = order,
      flip = flip,
      epsilon = epsilon
    )
    diff['z'] = distance(z, z.new, "maximum")
    z = z.new

    ## update y
    y.new = prs(
      x,
      z,
      gamma = gamma["y"],
      rotate = rotate,
      shrink = shrink,
      normalize = normalize,
      order = order,
      flip = flip,
      epsilon = epsilon
    )
    diff['y'] = distance(y, y.new, "maximum")
    y = y.new

    ## update b and score
    b <- t(z.new) %*% x %*% y.new
    # score.new = norm(b, 'F')
    # diff['score'] = score.new - score
    score = norm(b, 'F')

    ## report progress
    if (!quiet) {
      cat("Iter",
          n.iter,
          "score =",
          format(score, digit = 7),
          "and the differences:\n")
      print(diff, digit = 4, quote = F)
    }

    ## convergence?
    if (max(diff) < epsilon)
      break
  }

  rownames(z) <- rownames(x)
  rownames(y) <- colnames(x)
  res <- list(
    z = z,
    b = b,
    y = y,
    score = score,
    n.iter = n.iter,
    call = match.call()
  )
  class(res) <- "sma"
  return(res)
}

#' Print SMA
#'
#' @method print sma
#' @param x an `sma` object.
#' @param verbose `logical(1)`, whether to print out loadings.
#' @param ... additional input to generic [print].
#' @return Print an `sma` object interactively.
#' @export
print.sma <- function(x, verbose = FALSE, ...) {
  cat("Call: ")
  dput(x$call)
  cat("\n\n")
  cat("Num. non-zero z's: ", colSums(!!x$z), "\n")
  cat("Num. non-zero y's: ", colSums(!!x$y), "\n")
  cat("Abs. sum z's (L1-norm): ", norm(x$z, "1"), "\n")
  cat("Abs. sum y's (L1-norm): ", norm(x$y, "1"), "\n")
  if (verbose) {
    rn <- rownames(x$z)
    cn <- rownames(x$y)
    if (is.null(rn))
      rn <- 1:nrow(x$z)
    if (is.null(cn))
      cn <- 1:nrow(x$y)
    for (k in 1:ncol(x$y)) {
      cat("\n Component ", k, ":\n")
      u <- x$z[, k]
      v <- x$y[, k]
      cat(fill = TRUE)
      us <- cbind(rn[!!u], round(u[!!u], 3))
      dimnames(us) <-
        list(1:sum(!!u), c("row feature", "row weight"))
      vs <- cbind(cn[!!v], round(v[!!v], 3))
      dimnames(vs) <-
        list(1:sum(!!v), c("column feature", "column weight"))
      print(us, quote = FALSE, sep = "\t")
      cat(fill = TRUE)
      print(vs, quote = FALSE, sep = "\t")
    }
  }
}

#' Sparse Component Analysis
#'
#' `sca` performs sparse principal components analysis on the given numeric data matrix.
#' Choices of rotation techniques and shrinkage operators are available.
#'
#' @param gamma `numeric(1)`, sparsity parameter, default to `sqrt(pk)`, where n x p is the dimension of `x`.
#' @param is.cov  `logical`, default to `FALSE`, whether the `x` is a covariance matrix (or Gram matrix, i.e., `crossprod()` of some design matrix). If `TRUE`, both `center` and `scale` will be ignored/skipped.
#' @inheritParams prs
#' @inheritParams sma
#' @includeRmd man/rotate.md details
#' @includeRmd man/shrink.md details
#' @includeRmd man/normalize.md details
#' @includeRmd man/order.md details
#' @includeRmd man/flip.md details
#'
#' @return an `sca` object that contains:
#' \item{loadings}{`matrix`, sparse loadings of PCs.}
#' \item{scores}{an n x k `matrix`, the component scores, calculated using centered (and/or scaled) `x`. This will only be available when `is.cov = FALSE`.}
# \item{sss}{a `numeric` vector of length `k`, sum of squared scores (by column) corresponding to each sparse PC. These numbers were used to re-order sparse PCs, if `order = TRUE`.}
#' \item{cpve}{a `numeric` vector of length `k`, cumulative proportion of variance in `x` explained by the top PCs (after center and/or scale).}
#' \item{center}{`logical`, this records the `center` parameter.}
#' \item{scale}{`logical`, this records the `scale` parameter.}
#' \item{n.iter}{`integer`, number of iteration taken.}
#' \item{n.obs}{`integer`, sample size, that is, `nrow(x)`.}
#'
#' @seealso [sma], [prs]
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse Principal Component Analysis."
#' @examples
#' ## ------ example 1 ------
#' ## simulate a low-rank data matrix with some additive Gaussian noise
#' n <- 300
#' p <- 50
#' k <- 5 ## rank
#' z <- shrinkage(polar(matrix(runif(n * k), n, k)), sqrt(n))
#' b <- diag(5) * 3
#' y <- shrinkage(polar(matrix(runif(p * k), p, k)), sqrt(p))
#' e <- matrix(rnorm(n * p, sd = .01), n, p)
#' x <- scale(z %*% b %*% t(y) + e)
#'
#' ## perform sparse PCA
#' s.sca <- sca(x, k)
#' s.sca
#'
#' ## ------ example 2 ------
#' ## use the `pitprops` data from the `elasticnet` package
#' data(pitprops)
#'
#' ## find 6 sparse PCs
#' s.sca <- sca(pitprops, 6, gamma = 6, is.cov = TRUE)
#' print(s.sca, verbose = TRUE)
#'
#' @export
sca = function(x,
               k = min(5, dim(x)),
               gamma = NULL,
               is.cov = FALSE,
               rotate = c("varimax", "absmin"),
               shrink = c("soft", "hard"),
               center = TRUE,
               scale = FALSE,
               normalize = FALSE,
               order = TRUE,
               flip = TRUE,
               max.iter = 1e3,
               epsilon = 1e-5,
               quiet = TRUE) {
  ## check gamma
  if (!length(gamma))
    gamma = sqrt(ncol(x)) * sqrt(k)
  stopifnot(length(gamma) == 1)
  if (gamma < k || gamma > k * sqrt(ncol(x)))
    message("Improper sparsity parameter (gamma).")

  ## covariance or Gram matrix
  if (is.cov) {
    stopifnot("'x' must be symmetric when is.cov = TRUE." = isSymmetric(x))
    # x = rootmatrix(x)
    gamma2 <- rep(gamma, 2)
    center <- FALSE
    scale <- FALSE
  } else {
    gamma2 <- c(k * sqrt(nrow(x)), gamma) ## for SMA
  }

  ## center and scale
  x = scale(x = x,
            center = center,
            scale = scale)

  s = sma(
    x = x,
    k = k,
    gamma = gamma2,
    rotate = rotate,
    shrink = shrink,
    center = FALSE,
    scale = FALSE,
    normalize = normalize,
    order = order,
    flip = flip,
    max.iter = max.iter,
    epsilon = epsilon,
    quiet = quiet
  )

  ## construction result object
  res <- list()
  loadings = s$y
  rownames(loadings) <- colnames(x)
  colnames(loadings) <- paste0("PC", 1:k)
  res$loadings = loadings
  if (!is.cov) {
    scores = x %*% loadings
    rownames(scores) <- rownames(x)
    colnames(scores) <- paste0("PC", 1:k)
    res$scores = scores
    res$sss = colSums(scores ^ 2) ## Sum of squared scores
  }
  res$cpve = cpve(x, loadings)
  res$center = center
  res$scale = scale
  res$n.iter = s$n.iter ## number of iterations
  res$n.obs = ifelse(is.cov, NA, nrow(x))
  res$z = s$z
  res$y = s$y
  res$call = match.call()

  class(res) <- "sca"
  return(res)
}


#' Print SCA
#'
#' @method print sca
#'
#' @param x an `sca` object.
#' @param verbose `logical(1)`, whether to print out loadings.
#' @param ... additional input to generic [print].
#' @return Print an `sca` object interactively.
#' @export
print.sca <- function(x, verbose = FALSE, ...) {
  cat("Call:")
  dput(x$call)
  cat("\n\n")
  cat("Num. non-zero loadings':", colSums(!!x$loadings), "\n")
  cat("Abs. sum loadings' (L1-norm):", norm(x$loadings, "1"), "\n")
  cat("Cumulative proportion of variance explained (CPVE): \n")
  tab <- matrix(round(x$cpve, 3), dimnames =
                  list(paste(
                    "First", seq_along(x$cpve), "components:"
                  ), "CPVE"))
  rownames(tab)[1] = "First component:"
  print(tab,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE)
  if (verbose) {
    nm <- rownames(x$loadings)
    if (is.null(nm))
      nm <- 1:nrow(x$loadings)
    for (k in 1:ncol(x$loadings)) {
      cat("\n Component ", k, ":\n")
      v <- x$loadings[, k]
      cat(fill = TRUE)
      vs <- cbind(nm[!!v], round(v[!!v], 3))
      dimnames(vs) <- list(1:sum(!!v), c("feature", "loadings"))
      print(vs, quote = FALSE, sep = "\t")
    }
  }
}
