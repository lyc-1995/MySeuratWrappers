#' @include internal.R
#'
NULL

#' Calculate some QC metrics
#'
#' Compute quality control (QC) metrics for each control feature and cell in a Seurat object, accounting for specified control sets.
#'
#' @param object Seurat object
#' @param assay Name of assay to use
#' @param calculate.mapped.reads Whether or not to calculate mapped reads metrics. Default is FALSE. If TRUE, will need to set molecule information
#' @param mol.info Molecule information for mapped reads calculation (molecule_info.h5 files of 10x)
#' @param total.reads a data.frame with column named 'total.reads', collecting total count of reads per barcode.This can be extracted from .BAM file that CellRanger output and is used to calculate mapping rate per cell.
#' @param control.feature.list A list containing one or more vectors (a character vector of ensembl IDs named by feature symbols), used to identify feature controls
#' @param use.names Whether or not to use feature names instead of ID. If FALSE, shall specify symbol.to.id data.frame for transferation
#' @param symbol.to.id Data.frame for convertion from gene symbols to such as ENSEMBL ID, must contain at least two columns (use 'ID' and 'Symbol' as colnames)
#' @param do.log10 Whether or not to calculate log10 transformed QC metrics
#'
#' @importFrom Seurat DefaultAssay AddMetaData
#'
#' @return Returns Seurat object with new meta.data storing QC metrics
#'
#' @export
#'
CalculateSeuratQC <- function(
  object,
  assay = NULL,
  calculate.mapped.reads = FALSE,
  mol.info = NULL,
  total.reads = NULL,
  control.feature.list = NULL,
  use.names = FALSE,
  symbol.to.id = NULL,
  do.log10 = TRUE
) {
  assay <- assay %||% DefaultAssay(object)
  new.meta.data <- data.frame(row.names = colnames(object))
  new.meta.data[['nCount']] <- object@meta.data[, paste0('nCount_', assay)]
  new.meta.data[['nFeature']] <- object@meta.data[, paste0('nFeature_', assay)]
  object@meta.data <- object@meta.data[, which(!colnames(object@meta.data) %in% c(paste0('nFeature_', assay), paste0('nCount_', assay))), drop = FALSE]
  if (calculate.mapped.reads) {
    if (!is.null(total.reads)) {
      new.meta.data[['total.reads']] <- total.reads[['total.reads']][match(rownames(new.meta.data), rownames(total.reads))]
    } else {
      message('No available information about total reads per cell, will skip calculating mapping rate.')
    }
    if (!is.null(mol.info)) {
      mol.info$data <- mol.info$data[which(mol.info$data$cell %in% rownames(object@meta.data)), ]
      mapped.reads <- aggregate(mol.info$data[, 'reads'], by = list(mol.info$data$cell), FUN=sum)
      new.meta.data[['mapped.reads']] <- mapped.reads$x[match(rownames(new.meta.data), mapped.reads$Group.1)]
      if (!is.null(total.reads)) {
        new.meta.data[['mapping.rate']] <- new.meta.data[['mapped.reads']]*100 / new.meta.data[['total.reads']]
      }
      new.meta.data[['reads.per.umi']] <- new.meta.data[['mapped.reads']] / new.meta.data[['nCount']]
      new.meta.data[['reads.per.gene']] <- new.meta.data[['mapped.reads']] / new.meta.data[['nFeature']]
    } else {
      message('No available molecular information, will skip calculating mapped reads')
    }
  }
  new.meta.data[['umi.per.gene']] <- new.meta.data[['nCount']] / new.meta.data[['nFeature']]
  if (do.log10) {
    log10.new.meta.data <- log10(new.meta.data)
    colnames(log10.new.meta.data) <- paste0('log10.', colnames(new.meta.data))
    new.meta.data <- cbind(new.meta.data, log10.new.meta.data)
  }
  colnames(new.meta.data) <- sapply(colnames(new.meta.data), function(x) {
    paste0(x, '_', assay)
  })
  object <- AddMetaData(object, metadata = new.meta.data)
  if (!is.null(control.feature.list)) {
    if (use.names) {
      all.features <- rownames(object[[assay]]@counts)
      control.feature.list <- lapply(control.feature.list, names)
    } else {
      if (!is.null(symbol.to.id)) {
        all.features <- symbol.to.id$ID[match(rownames(object[[assay]]@counts), symbol.to.id$Symbol)]
        if (any(is.na(all.features))) {
          stop("The provided symbol.to.id data.frame dosen't contain right IDs for some features")
        }
      } else {
        stop('Cannot convert gene symbol to ensembl id without available symbol.to.id')
      }
    }
    total <- lapply(control.feature.list, function(.ele){
      Matrix::colSums(object[[assay]]@counts[which(all.features %in% .ele), , drop = FALSE])
    })
    percent <- lapply(total, function(.ele){
      .ele*100 / Matrix::colSums(object[[assay]]@counts)
    })
    for (j in 1:length(control.feature.list)) {
      object <- AddMetaData(object, metadata = total[[j]], col.name = paste0('total.', names(control.feature.list[j]), '_', assay))
      object <- AddMetaData(object, metadata = percent[[j]], col.name = paste0('percent.', names(control.feature.list[j]), '_', assay))
    }
    if (do.log10) {
      log10.total <- lapply(total, log10)
      for (j in 1:length(control.feature.list)) {
        object <- AddMetaData(object, metadata = log10.total[[j]], col.name = paste0('log10.total.', names(control.feature.list[j]), '_', assay))
      }
    }
  }
  rm(new.meta.data, log10.new.meta.data, total, percent, log10.total, all.features);gc(reset = TRUE)
  object <- LogSeuratCommand(object = object)
  return(object)
}

#' @param metrics Default is 'mean', and standard deviation will be used as discrete parameter. Median absolute deviation will be used for 'median' as metrics
#' @param n.dev Number for discrete parameter
#' @param cut.off.end Default is paired-end. Can set with 'lower' or 'higher' end
#' @param cut.off.given Add extra cut-off values
#' @param cut.off.given.only Whether or not to use the given cut-offs only, default is FALSE
#'
#' @rdname CalculateCutOff
#' @export
#'
CalculateCutOff.default <- function(
  object,
  metrics = c('mean', 'median'),
  n.dev = 3,
  cut.off.end = c('paired', 'lower', 'higher'),
  cut.off.given = NULL,
  cut.off.given.only = FALSE,
  ...
) {
  metrics <- match.arg(arg = metrics)
  cut.off.end <- match.arg(arg = cut.off.end)
  outliars <- c(lower = mean(object) - n.dev*sd(object),
                higher = mean(object) + n.dev*sd(object))
  if (metrics == 'median') {
    outliars <- c(lower = median(object) - n.dev*mad(object),
                  higher = median(object) + n.dev*mad(object))
  }
  if (cut.off.end != 'paired') {
    outliars <- outliars[cut.off.end]
  }
  if (!is.null(cut.off.given)) {
    names(cut.off.given) <- sapply(1:length(cut.off.given), function(x) {
      paste0('given.', as.character(x))
    })
    outliars <- c(outliars, cut.off.given)
    if (cut.off.given.only) {
      outliars <- cut.off.given
    }
  }
  return(outliars)
}

#' @param assay Which assay to use
#' @param feature Which feature to calculate. Must be from either rownames of Seurat data or colnames of meta.data
#' @param cells Which cells to use. Default is all cells
#' @param slot Which slot to use. Default is counts slot
#'
#' @import Seurat
#' @rdname CalculateCutOff
#' @export
#' @method CalculateCutOff Seurat
#'
CalculateCutOff.Seurat <- function(
  object,
  assay = NULL,
  feature = NULL,
  cells = NULL,
  slot = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  slot <- slot %||% 'counts'
  cells <- cells %||% colnames(GetAssayData(object = object, assay = assay, slot = slot))
  if (!is.null(feature)) {
    if (feature %in% rownames(GetAssayData(object = object, assay = assay, slot = slot))) {
      x <- t(GetAssayData(object = object, assay = assay, slot = slot)[feature, cells])
    }
    if (feature %in% colnames(object@meta.data)) {
      x <- object@meta.data[cells, feature]
    }
  } else {
    stop('No specified feature')
  }
  outliars <- CalculateCutOff(object = x, ...)
  return(outliars)
}
