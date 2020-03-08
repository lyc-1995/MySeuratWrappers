#' @param assay Name of Assay diffusion map is being run on
#' @param dim.embed Total Number of diffusion map components (DMCs) to compute and store (2 by default)
#' @param q Quantile to clip diffusion map components at. This addresses an
#' issue where 1-2 cells will have extreme values that obscure all other points.
#' 0.01 by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. DMC by default
#'
#' @importFrom utils installed.packages
#' @importFrom stats dist quantile
#'
#' @rdname RunDiffusion
#' @export
#' @method RunDiffusion matrix
#'
RunDiffusion.matrix <- function (
  object,
  assay = NULL,
  dim.embed = 2,
  q = 0.01,
  reduction.key = "DMC_",
  ...
) {
  if (!"destiny" %in% rownames(x = installed.packages())) {
    stop("Please install destiny - learn more at https://bioconductor.org/packages/release/bioc/html/destiny.html")
  }
  dist.data <- dist(object)
  diffusion <- destiny::DiffusionMap(data = as.matrix(dist.data), n_eigs = dim.embed, ...)
  diffusion.data <- data.frame(diffusion@eigenvectors)
  colnames(x = diffusion.data) <- paste0(reduction.key, 1:ncol(x = diffusion.data))
  rownames(x = diffusion.data) <- rownames(x = object)
  for (i in 1:dim.embed) {
    x <- diffusion.data[, i]
    x <- MinMax(data = x, min = quantile(x = x, probs = q),
                quantile(x = x, probs = 1 - q))
    diffusion.data[, i] <- x
  }
  diffusion.reduction <- CreateDimReducObject(
    embeddings = as.matrix(diffusion.data),
    stdev = diffusion@eigenvalues,
    key = reduction.key,
    assay = assay
  )
  rm(dist.data, diffusion, diffusion.data)
  return(diffusion.reduction)
}

#' @param cells Which cells to analyze (default, all cells)
#' @param dims Which dimensions to use as input features.
#'
#' @import Seurat
#' @rdname RunDiffusion
#' @export
#' @method RunDiffusion DimReduc
#'
RunDiffusion.DimReduc <- function(
  object,
  cells = NULL,
  dims = 1:5,
  q = 0.01,
  dim.embed = 2,
  reduction.key = "DMC_",
  ...
) {
  args <- as.list(x = sys.frame(which = sys.nframe()))
  args <- c(args, list(...))
  args$object <- args$object[[cells, args$dims]]
  args$dims <- NULL
  args$cells <- NULL
  args$assay <- DefaultAssay(object = object)
  return(do.call(what = 'RunDiffusion', args = args))
}

#' @param features If set, run the diffusion map on this subset of features
#' (instead of running on a set of reduced dimensions). Not set (NULL) by default;
#' \code{dims} will automaticly be NULL to run on features
#' @param reduction Which dimensional reduction (e.g. PCA, ICA) to use for
#' the diffusion map. Default is PCA
#' @param normalize.method 'LogNormalize' or 'SCTransform'. If NULL, set 'LogNormalize'
#' by default
#' @param scale.clip Max/min value for scaled data. Default is 10
#' @param reduction.name dimensional reduction name, specifies the position in the object$dr list. diffusionmap by default
#'
#' @import Seurat
#' @rdname RunDiffusion
#' @export
#' @method RunDiffusion Seurat
#'
RunDiffusion.Seurat <- function (
  object,
  assay = NULL,
  cells = NULL,
  dims = 1:5,
  reduction = "diffusionmap",
  features = NULL,
  normalize.method = c("LogNormalize", "SCTransform"),
  q = 0.01,
  dim.embed = 2,
  scale.clip = 10,
  reduction.name = "diffusionmap",
  reduction.key = "DMC_",
  ...
) {
  cells <- cells %||% Cells(x = object)
  if (!normalize.method == "SCTransform") {
    normalize.method <- "LogNormalize"
  }
  if (!is.null(features)) {
    dims <- NULL
  }
  diffusion.reduction <- if (!is.null(x = dims)) {
    RunDiffusion(
      object = object[[reduction]],
      cells = cells,
      dims = dims,
      q = q,
      dim.embed = dim.embed,
      reduction.key = reduction.key,
      ...
    )} else if (!is.null(x = features)) {
      RunDiffusion(
        object = MinMax(data = t(as.matrix(x = if (normalize.method == "LogNormalize") {
          GetAssayData(object = object, assay = assay)[features, cells]
        } else if (normalize.method == "SCTransform") {
          GetAssayData(object = object, assay = assay, slot = "scale.data")[features, cells]
        })), min = -1 * scale.clip, max = scale.clip),
        assay = DefaultAssay(object = object),
        q = q,
        dim.embed = dim.embed,
        reduction.key = reduction.key,
        ...
      )} else {
        stop("Unknown way of running DiffusionMap")
      }
  object[[reduction.name]] <- diffusion.reduction
  object <- LogSeuratCommand(object = object)
  return(object)
}
