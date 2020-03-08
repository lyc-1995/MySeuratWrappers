
# Preprocess ---------------------------

#' Calculate cut-off value
#'
#' Compute cut-off for outliars of an given feature.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a vector storing cut-off(s) for outliars
#'
#' @export
#'
#' @rdname CalculateCutOff
#' @export CalculateCutOff
#'
CalculateCutOff <- function (object, ...) {
  UseMethod(generic = 'CalculateCutOff', object = object)
}

#
# #' Regress out technical effects and cell cycle using regularized Negative Binomial regression (old Seurat function)
# #'
# #' Remove unwanted effects from umi data and set scale.data to Pearson residuals
# #' Uses future_lapply; you can set the number of cores it will use to n with plan(strategy = 'multicore', workers = n).
# #' If n.features is set, only a (somewhat-random) subset of features is used for estimating the initial model parameters.
# #'
# #' @param object An object
# #' @param ... Arguments passed to other methods
# #'
# #' @return Returns a Seurat object with pearson residuals. Intermediate results are saved in a misc slot
# #'
# #' @references Mayer, C., Hafemeister, C., Bandler, R. C., Machold, R., Batista Brito, R., Jaglin, X., … Satija, R. (2018). Developmental diversification of cortical inhibitory interneurons. Nature, 555(7697), 457–462. https://doi.org/10.1038/nature25999
# #'
# #' @rdname RegressOutNBreg
# #' @export RegressOutNBreg
# #'
# RegressOutNBreg <- function (object, ...) {
#   UseMethod(generic = 'RegressOutNBreg', object = object)
# }

# Dimensional reduction ---------------------------

#' Run diffusion map
#'
#' Run a diffusion map dimensionality reduction. For details about stored difussion map calculation
#' parameters, see \code{PrintDiffusionParams}.
#'
#' @param object An object
#' @param ... Arguments passed to other methods and destiny
#'
#' @return Returns Seurat object with the diffusion map calculation stored in the reductions slot
#'
#' @export
#'
#' @rdname RunDiffusion
#' @export RunDiffusion
#'
RunDiffusion <- function (object, ...) {
  UseMethod(generic = 'RunDiffusion', object = object)
}

# Clustering ---------------------------

#' Perform spectral density clustering on single cells
#'
#' Find point clounds single cells in a two-dimensional space using density clustering (DBSCAN).
#'
#' @param object An object
#' @param ... Arguments passed to other methods and destiny
#'
#' @return Returns a Seurat object where the idents have been updated with new cluster info;
#' latest clustering results will be stored in object metadata under 'DB_clusters'.
#' Note that 'DB_clusters' will be overwritten everytime FindClusters is run
#'
#' @export
#'
#' @rdname DBclustering
#' @export DBclustering
#'
#' @examples
#' pbmc_small
#' # Density based clustering on the first two tSNE dimensions
#' pbmc_small <- DBclustering(pbmc_small)
#'
DBclustering <- function(object, ...)
{
  UseMethod(generic = 'DBclustering', object = object)
}
