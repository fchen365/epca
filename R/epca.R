#' Exploratory Principal Component Analysis
#'
#' `epca` is for comprehending any data matrix that contains *low-rank* and *sparse* underlying signals of interest.
#' The package currently features two key tools: (1) `sca` for **s**parse principal **c**omponent **a**nalysis and
#' (2) `sma` for **s**parse **m**atrix **a**pproximation, a two-way data analysis for simultaneously row and column dimensionality reductions.
#'
#' @name epca-package
#' @docType package
#' @references Chen, F. and Rohe K. (2020) "A New Basis for Sparse PCA".
#' @keywords package
NULL


#' Pitprops correlation data
#'
#' The `pitprops` data is a correlation matrix that was calculated from 180 observations.
#' There are 13 explanatory variables.
#' Jeffers (1967) tried to interpret the first six PCs.
#' This is a classical example showing the difficulty of interpreting principal components.
#'
#' @name pitprops
#' @docType data
#'
#' @references
#' Jeffers, J. (1967) "Two case studies in the application of principal component", *Applied Statistics*, 16, 225-236.
#'
#' @keywords datasets
#' @examples
#' \donttest{
#' ## NOT TEST
#' data(pitprops)
#' ggcorrplot::ggcorrplot(pitprops)
#' }
#'
NULL


