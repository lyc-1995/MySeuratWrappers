#' @include internal.R
#'
NULL

# Read 10X data -----------------------------------------------------------------

#' Load in data from 10X Cellranger aggr
#'
#' Enables easy loading of sparse data matrices provided by 10X Cellranger aggr.
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
#' files provided by by 10X Cellranger aggr
#' @param sample.names A vector of characters indicating the sample names. The length of sample.names must be equal
#' to the number of libraries input to aggr pipeline
#' @param split Whether or not to split the matrix by each sample
#' @param gene.column Specify which column of genes.tsv or features.tsv to use for gene names; default is 2
#' @param unique.features Make feature names unique (default TRUE)
#'
#' @return If split=FALSE, will return a sparse matrix with cell barcodes attched to each sample name. If split=TRUE,
#' will return a list of splited matrices
#'
#' @importFrom Seurat Read10X
#' @importFrom stringr str_extract
#'
#' @export
#'
Read10X_aggr <- function(
  data.dir,
  sample.names,
  split = FALSE,
  gene.column = 2,
  unique.features = TRUE
) {
  raw.counts <- Read10X(data.dir = data.dir, gene.column = gene.column, unique.features = unique.features)
  rownames(x = raw.counts) <- gsub(pattern = '_', replacement = '-', x = rownames(x = raw.counts))
  raw.bc <- colnames(x = raw.counts)
  lib.id <- str_extract(string = raw.bc, pattern = '[0-9]*$')
  if (length(x = sample.names) != dim(x = table(lib.id))) {
    stop("The number of sample names dosen't match the number of libraries")
  }
  trimed.bc <- gsub(pattern = '[0-9]*$', replacement = '', x = raw.bc)
  trimed.bc <- gsub(pattern = '-', replacement = '', x = trimed.bc)
  bc <- data.frame(raw.bc = raw.bc,
                   lib.id = as.numeric(lib.id),
                   trimed.bc = trimed.bc)
  bc[['sample']] <- sample.names[bc[['lib.id']]]
  bc[['new.bc']] <- paste0(bc[['sample']], '_', bc[['trimed.bc']])
  colnames(x = raw.counts) <- bc[['new.bc']]
  if (split) {
    bc[['sample']] <- factor(x = bc[['sample']], levels = sample.names)
    bc <- split(x = bc, f = bc[, 'sample'], drop = FALSE)
    new.counts <- lapply(X = bc, function(x) {
      raw.counts[, x[['new.bc']], drop = FALSE]
    })
    return(new.counts)
  } else {
    return(raw.counts)
  }
}

# Cell cycle scoring -----------------------------------------------------------------

#' Calculate cell-cycle z-score
#'
#' Another method to assign cell-cycle phase by computing z-score
#'
#' @param object A Seurat object
#' @param assay Which assay to use
#' @param s.features A vector of features associated with S phase
#' @param g2m.features A vector of features associated with G2M phase
#' @param min.cells Include features detected in at least this many cells
#' @param N Number of random background feature sets, default is 100
#' @param seed Set a random seed
#'
#' @return A Seurat object with the following columns added to object meta data: s.score, g2m.score, and phase
#'
#' @references Mayer, C., Hafemeister, C., Bandler, R. C., Machold, R., Batista Brito, R., Jaglin, X., … Satija, R. (2018). Developmental diversification of cortical inhibitory interneurons. Nature, 555(7697), 457–462. https://doi.org/10.1038/nature25999
#'
#' @import Seurat
#'
#' @export
#'
CellCycleEnrichScore <- function(
  object,
  assay = NULL,
  s.features = NULL,
  g2m.features = NULL,
  min.cells = 3,
  N = 100,
  seed = 42
) {
  set.seed(seed)
  if (is.null(s.features) | is.null(g2m.features)) {
    stop('No available cell-cycle features.')
  }
  assay = assay %||% DefaultAssay(object = object)
  message('Getting cell-cycle score, number of random background feature sets set to ', N, '\n')
  counts <- GetAssayData(object = object, assay = assay, slot = 'counts')
  nCells <- apply(counts > 0, 1, sum)
  counts <- counts[nCells >= min.cells, ]

  feature.mean <- apply(counts, 1, mean)
  breaks <- unique(quantile(log10(feature.mean), probs = seq(0,1, length.out = 50)))
  feature.bin <- cut(log10(feature.mean), breaks = breaks, labels = FALSE)
  names(feature.bin) <- rownames(counts)
  feature.bin[is.na(feature.bin)] <- 0

  feature.list <- list('S' = rownames(counts)[!is.na(match(toupper(rownames(counts)), s.features))],
                       'G2M' = rownames(counts)[!is.na(match(toupper(rownames(counts)), g2m.features))])

  n <- min(40, min(sapply(feature.list, length)))
  feature.list <- lapply(feature.list, function(x) x[order(feature.mean[x], decreasing = TRUE)[1:n]])
  bg.list <- list('S'=GetBackgroundList(feature.list[['S']], N, feature.bin),
                  'G2M'=GetBackgroundList(feature.list[['G2M']], N, feature.bin))
  all.features <- sort(unique(c(unlist(feature.list, use.names=FALSE), unlist(bg.list, use.names=FALSE))))

  expr <- log10(counts[all.features, ]+1)
  rm(counts);gc(reset = TRUE)

  s.score <- EnrichScore(expr, feature.list[['S']], bg.list[['S']])
  g2m.score <- EnrichScore(expr, feature.list[['G2M']], bg.list[['G2M']])
  rm(expr);gc(reset = TRUE)

  phase <- rep('G1', times = length(s.score))
  phase[g2m.score > 2 & s.score <= 2] <- 'G2/M'
  phase[g2m.score <= 2 & s.score > 2] <- 'S'
  phase <- factor(phase, levels = c('G1', 'S', 'G2/M'))
  cc <- data.frame(score = s.score - g2m.score, s.score, g2m.score, phase)

  object <- AddMetaData(object, metadata = cc)
  rm(cc);gc(reset = TRUE)
  object <- LogSeuratCommand(object = object)
  return(object)
}

