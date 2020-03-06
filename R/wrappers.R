#' @include internal.R
#'
NULL

#' Wrapper for SCTransform
#'
#' Wrapper for SCTransform and CellCycleEnrichScore
#'
#' @param object A seurat object
#' @param assay Name of assay to pull the count data from; default is 'RNA'
#' @param get.cellcycle.score Whether or not to calculate cell-cycle enrich score
#' @param min.cells Include features detected in at least this many cells
#' @param min.cells.cc Include features detected in at least this many cells when assign cell-cycle phase. If null, will be as same as min.cells
#' @param s.features A vector of features associated with S phase
#' @param g2m.features A vector of features associated with G2M phase
#' @param bg.N Number of random background feature sets, default is 100
#' @param seed Set a random seed
#' @param vars.to.regress Variables to regress out in a second non-regularized linear
#' regression. For example, percent.mito. Default is NULL
#' @param variable.features.to.remove A vector of features to be removed from variable features
#' @param future.multicore Whether or not to use multicore computation
#' @param workers Number of workers used for future.multicore
#' @param ... Additional parameters passed to \code{\link{SCTransform}}
#'
#' @import Seurat
#' @import future
#'
#' @export
#'
SCTransformWrapper <- function(
  object,
  assay = NULL,
  get.cellcycle.score = FALSE,
  min.cells = 5,
  min.cells.cc = NULL,
  s.features = NULL,
  g2m.features = NULL,
  bg.N = 100,
  seed = 42,
  vars.to.regress = NULL,
  variable.features.to.remove = NULL,
  future.multicore = TRUE,
  workers = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (get.cellcycle.score) {
    min.cells.cc <- min.cells.cc %||% min.cells
    object <- CellCycleEnrichScore(object = object, min.cells = min.cells.cc, s.features = s.features, g2m.features = g2m.features, N = bg.N, seed = seed)
  }
  if (future.multicore) {
    workers <- workers %||% 4
    future::plan(strategy = 'multiprocess', workers = workers)
  }
  object <- SCTransform(object = object, min_cells = min.cells, vars.to.regress = vars.to.regress, assay = assay, ...)
  plan("sequential");gc(reset = TRUE)

  if (!is.null(variable.features.to.remove)) {
    message('Removing ', length(intersect(variable.features.to.remove, VariableFeatures(object))), ' features from variable features')
    VariableFeatures(object) <- VariableFeatures(object)[!VariableFeatures(object) %in% variable.features.to.remove]
  }
  return(object)
}
