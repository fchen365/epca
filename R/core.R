
#' Polar-Rotate-Shrink
#'
#' This function is a helper function of [sma()].
#' It performs polar decomposition, orthogonal rotation, and soft-thresholding shrinkage in order.
#' The three steps together enable sparse estimates of the SMA and SCA.
#'
#' @param X,Z.hat the matrix product `A <- crossprod(X, Z.hat)` is the input. `X` and `Z.hat` are separated because the former is additionally used to compute the proportion of variance explained, in the case when `order = TRUE`.
#' @param gamma `numeric`, the sparsity parameter.
#' @param order `logical`, whether to re-order the columns of the estimates (see details).
#' @param epsilon `numeric`, tolerance of convergence precision (default to 0.00001).
#' @inheritParams rotation
#' @inheritParams shrinkage
#' @includeRmd man/rotate.md details
#' @includeRmd man/shrink.md details
#' @includeRmd man/normalize.md details
#' @includeRmd man/order.md details
#' @includeRmd man/flip.md details
#' @return a `matrix` of the sparse estimate, of the same dimension as `crossprod(X, Z.hat)`.
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse PCA."
#' @seealso [sma], [sca], [polar], [rotation], [shrinkage]
prs <- function(X, Z.hat,
                gamma,
                rotate,
                shrink,
                normalize,
                order,
                flip,
                epsilon) {
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
    sdev <- apply(X %*% Y.hat, 2, stats::sd)
    ord <- order(sdev, decreasing = TRUE)
    Y.hat <- Y.hat[,ord,drop = FALSE]
  }

  Y.hat
}

#' Sparse Matrix Approximation
#'
#' Perform the sparse matrix approximation (SMA) of a data matrix X as three components: Z B Y',
#' where Z and Y are sparse, and B is low-rank but not necessarily diagonal.
#'
#' @param A `matrix` or `Matrix` to be analyzed.
#' @param k `integer`, rank of approximation.
#' @param gamma `numeric(2)`, sparsity parameters. If `gamma` is `numeric(1)`, it is used for both left and right sparsity component (i.e, Z and Y). If absent, the two parameters are set as (default): `sqrt(nk)` and `sqrt(pk)` for Z and Y respectively, where n x p is the dimension of `A`.
#' @param center `logical`, whether to center columns of `A` (see [scale()]).
#' @param scale `logical`, whether to scale columns of `A` (see [scale()]).
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
#' \item{Z, B, Y}{the three parts in the SMA (i.e., *ZBY'*).
#' Z is a sparse n x k `matrix` that contains the row components (loadings).
#' The row names of Z inherit the row names of `A`.
#' B is a k x k `matrix` that contains the scores of SMA;
#' the Frobenius norm of B equals to the total variance explained by the SMA.
#' Y is a sparse n x k `matrix`that contains the column components (loadings).}
#' The row names of Y inherit the column names of `A`.
#' \item{score}{the total variance explained by the SMA.
#' This is the optimal objective value obtained.}
#' \item{n.iter}{`integer`, the number of iteration taken.}
#'
#' @seealso [sca], [prs]
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse PCA."
#' @examples
#' ## simulate a rank-5 data matrix with some additive Gaussian noise
#' n <- 300
#' p <- 50
#' k <- 5 ## rank
#' Z <- shrinkage(polar(matrix(runif(n * k), n, k)), sqrt(n))
#' B <- diag(5) * 3
#' Y <- shrinkage(polar(matrix(runif(p * k), p, k)), sqrt(p))
#' E <- matrix(rnorm(n * p, sd = .01), n, p)
#' X <- scale(Z %*% B %*% t(Y) + E)
#'
#' ## perform sparse matrix approximation
#' s.sma <- sma(X, k)
#' s.sma
#'
#' @export
sma = function(A,
               k = min(5, dim(A)),
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
  # S = RSpectra::svds(A, k)
  S = irlba::irlba(A, k, tol = 1e-10)
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

  rownames(Z) <- rownames(A)
  rownames(Y) <- colnames(A)
  res <- list(
    Z = Z,
    B = B,
    Y = Y,
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
  cat("Num. non-zero Z's: ", colSums(!!x$Z), "\n")
  cat("Num. non-zero Y's: ", colSums(!!x$Y), "\n")
  cat("Abs. sum Z's: ", norm(x$Z, "1"), "\n")
  cat("Abs. sum Y's: ", norm(x$Y, "1"), "\n")
  if(verbose){
    rn <- rownames(x$Z)
    cn <- rownames(x$Y)
    if(is.null(rn)) rn <- 1:nrow(x$Z)
    if(is.null(cn)) cn <- 1:nrow(x$Y)
    for(k in 1:ncol(x$Y)){
      cat("\n Component ", k, ":\n")
      u <- x$Z[,k]
      v <- x$Y[,k]
      cat(fill = TRUE)
      us <- cbind(rn[!!u], round(u[!!u], 3))
      dimnames(us) <- list(1:sum(!!u), c("row feature", "row weight"))
      vs <- cbind(cn[!!v], round(v[!!v], 3))
      dimnames(vs) <- list(1:sum(!!v), c("column feature", "column weight"))
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
#' @param gamma `numeric(1)`, sparsity parameter, default to `sqrt(pk)`, where n x p is the dimension of `A`.
#' @param is.cov  `logical`, whether the `A` is a covariance matrix or Gram matrix (i.e., `crossprod(X)`). This function presumes that `A` is *not* covariance matrix.
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
#' \item{scores}{an n x k `matrix`, the component scores.}
#' \item{sdev}{a `numeric` vector of length `k`, standard deviation of each columns of scores. These may not sum to exactly 1 because of a slight loss of orthogonality.}
#' \item{pve}{a `numeric` vector of length `k`, cumulative proportion of variance in `A` explained by the top PCs.}
#' \item{center}{`logical`, this records the `center` parameter.}
#' \item{scale}{`logical`, this records the `scale` parameter.}
#' \item{n.iter}{`integer`, number of iteration taken.}
#' \item{n.obs}{`integer`, sample size, that is, `nrow(A)`.}
#'
#' @seealso [sma], [prs]
#' @references Chen, F. and Rohe, K. (2020) "A New Basis for Sparse PCA."
#' @examples
#' ## ------ example 1 ------
#' ## simulate a low-rank data matrix with some additive Gaussian noise
#' n <- 300
#' p <- 50
#' k <- 5 ## rank
#' Z <- shrinkage(polar(matrix(runif(n * k), n, k)), sqrt(n))
#' B <- diag(5) * 3
#' Y <- shrinkage(polar(matrix(runif(p * k), p, k)), sqrt(p))
#' E <- matrix(rnorm(n * p, sd = .01), n, p)
#' X <- scale(Z %*% B %*% t(Y) + E)
#'
#' ## perform sparse PCA
#' s.sca <- sca(X, k)
#' s.sca
#'
#' ## ------ example 2 ------
#' ## use the `pitprops` data from the `elasticnet` package
#' data(pitprops)
#'
#' ## find 3 sparse PCs
#' s.sca <- sca(pitprops, 3, gamma = 4.5)
#' print(s.sca, verbose = TRUE)
#'
#' ## find 6 sparse PCs
#' s.sca <- sca(pitprops, 6, gamma = 6)
#' print(s.sca, verbose = TRUE)
#'
#' @export
sca = function(A,
               k = min(5, dim(A)),
               gamma = NULL,
               is.cov = FALSE,
               rotate = c("varimax", "absmin"),
               shrink = c("soft", "hard"),
               center = TRUE,
               scale = TRUE,
               normalize = FALSE,
               order = TRUE,
               flip = TRUE,
               max.iter = 1e3,
               epsilon = 1e-5,
               quiet = TRUE) {
  ## check gamma
  if (!length(gamma))
    gamma = sqrt(ncol(A)) * sqrt(k)
  stopifnot(length(gamma) == 1)
  if (gamma < k || gamma > k * sqrt(ncol(A)))
    message("Improper sparsity parameter (gamma).")

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
          center = FALSE,
          scale = FALSE,
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
  sdev = apply(scores, 2, stats::sd) ## standard deviation of scores (column)
  pve = cpve(A, loadings)

  res <- list(loadings = loadings,
              scores = scores,
              pve = pve,
              sdev = sdev,
              center = center,
              scale = scale,
              n.iter = n.iter,
              n.obs = ifelse(is.cov, NA, nrow(A)),
              call = match.call())

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
#' @return Print an `sma` object interactively.
#' @export
print.sca <- function(x, verbose = FALSE, ...) {
  cat("Call:")
  dput(x$call)
  cat("\n\n")
  cat("Num. non-zero loadings':", colSums(!!x$loadings), "\n")
  cat("Abs. sum loadings':", norm(x$loadings, "1"), "\n")
  cat("Cumulative proportion of variance explained (CPVE): \n")
  tab <- matrix(round(x$pve, 3), dimnames =
                  list(paste("First", seq_along(x$pve), "components:"), "CPVE"))
  rownames(tab)[1] = "First component:"
  print(tab, quote = FALSE, sep = "\t", row.names = FALSE)
  if(verbose){
    nm <- rownames(x$loadings)
    if (is.null(nm)) nm <- 1:nrow(x$loadings)
    for(k in 1:ncol(x$loadings)){
      cat("\n Component ", k, ":\n")
      v <- x$loadings[, k]
      cat(fill = TRUE)
      vs <- cbind(nm[!!v], round(v[!!v], 3))
      dimnames(vs) <- list(1:sum(!!v), c("feature", "loadings"))
      print(vs, quote = FALSE, sep = "\t")
    }
  }
}



