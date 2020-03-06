#' @include internal.R
#'
NULL

#' @param assay Assay to use in DBSCAN clustering
#' @param reduction Which dimensional reduction to use, default is tsne.
#'
#' @import Seurat
#'
#' @rdname DBclustering
#' @export
#' @method DBclustering Seurat
#'
DBclustering.Seurat <- function(
  object,
  assay = NULL,
  reduction = "tsne",
  key = NULL,
  dim.1 = 1,
  dim.2 = 2,
  eps,
  minPts = 5,
  seed.use = 42,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  DimReduction <- Reductions(
    object = object,
    slot = reduction
  )
  if (assay != DefaultAssay(object = DimReduction)) {
    stop("The provided assay dosen't match the assay.used of DimReduc object.")
  }
  res <- DBclustering(object = DimReduction, dim.1 = dim.1, dim.2 = dim.2, key = key, eps = eps, minPts = minPts, seed.use = seed.use, ...)
  colnames(x = res) <- paste(assay, colnames(x = res), sep = '_')
  object <- AddMetaData(object = object, metadata = res)
  Idents(object = object) <- colnames(x = res)[ncol(x = res)]
  levels <- levels(x = object)
  levels <- tryCatch(
    expr = as.numeric(x = levels),
    warning = function(...) {
      return(levels)
    },
    error = function(...) {
      return(levels)
    }
  )
  Idents(object = object) <- factor(x = Idents(object = object), levels = sort(x = levels))
  object[['DB_clusters']] <- Idents(object = object)
  object <- LogSeuratCommand(object)
  return(object)
}

#' @param key Key used in provided DimReduc object
#' @param dim.1 First dimension to use
#' @param dim.2 second dimension to use
#'
#' @import Seurat
#'
#' @rdname DBclustering
#' @export
#' @method DBclustering DimReduc
#'
DBclustering.DimReduc <- function(
  object,
  key = NULL,
  dim.1 = 1,
  dim.2 = 2,
  eps,
  minPts = 5,
  seed.use = 42,
  ...
) {
  key <- key %||% Key(object = object)
  cell.embeddings <- Embeddings(object = object)
  x1 <- paste0(key, dim.1)
  x2 <- paste0(key, dim.2)
  key.check <- table(c(x1, x2) %in% colnames(cell.embeddings))
  if (all(key.check) != TRUE) {
    stop("Provided key dosen't match column names of cell embeddings.")
  }
  cell.embeddings <- cell.embeddings[, c(x1, x2)]
  res <- DBclustering(object = cell.embeddings, eps = eps, minPts = minPts, seed.use = seed.use, ...)
  colnames(x = res) <- paste(key, colnames(x = res), sep = '_')
  return(res)
}

#' @param eps size of the epsilon neighborhood. Lower value to get more fine-scale clustering.
#' @param minPts number of minimum points in the eps region (for core points). Default is 5 points.
#' @param seed.use Random seed for the dbscan function
#'
#' @import dbscan
#'
#' @rdname DBclustering
#' @export
#'
DBclustering.default <- function(
  object,
  eps,
  minPts = 5,
  seed.use = 42,
  ...
) {
  set.seed(seed = seed.use)
  db.res <- dbscan::dbscan(x = object, eps = eps, minPts = minPts, ...)
  res <- data.frame(cluster = as.factor(x = db.res$cluster),
                    row.names = rownames(x = object))
  colnames(res) <- paste(eps, minPts, sep = '_')
  return(res)
}
