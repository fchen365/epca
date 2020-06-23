#' `sca` result for scRNA-seq data
#'
#' The sparse PCA result (by `sca`) for the pancreas islet cell gene expression (types) from Baron et al. (2017).
#' We removed the genes that do not have any variation across samples (i.e., zero standard deviation)
#' and the cell types that contain fewer than 100 cells.
#' This resulted in the expression (counts) of 17499 genes in 8451 cells.
#'
#' @name scar
#' @docType data
#'
#' @references
#'
#' Baron M et al. (2017). Single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell population structure. *Cell Syst.* 3(4), 346-360.
#'
#' @keywords datasets
#' @examples
#' data(scar, package = "epca")
#' print(scar)
NULL

#' Cell type labels in scRNA-seq data
#'
#' The cell type labels in the pancreas islet cell gene expression (types) from Baron et al. (2017).
#' We removed the genes that do not have any variation across samples (i.e., zero standard deviation)
#' and the cell types that contain fewer than 100 cells.
#' This resulted in 8451 cells across 9 cell types.
#'
#' @name label
#' @docType data
#'
#' @references
#'
#' Baron M et al. (2017). Single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell population structure. *Cell Syst.* 3(4), 346-360.
#'
#' @keywords datasets
#' @examples
#' data(label, package = "epca")
#' table("cell type" = label)
#'
NULL

