#' @include internal.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Dimensional reduction plot
#'
#' Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a
#' cell and it's positioned based on the cell embeddings determined by the reduction technique. By
#' default, cells are colored by their identity class (can be changed with the group.by parameter).
#'
#' @param object Seurat object
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
#' or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
#' By default, ggplot2 assigns colors. We also include a number of palettes from the pals package.
#' See \code{\link{DiscretePalette}} for details.
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of one or more metadata columns to group (color) cells by
#' (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by;
#' see \code{\link{FetchData}} for more details
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param order Specify the order of plotting for the idents. This can be
#' useful for crowded plots if points of interest are being buried. Provide
#' either a full list of valid idents or a subset to be plotted last (on top)
#' @param label Whether to label the clusters
#' @param label.size Sets size of labels
#' @param repel Repel labels
#' @param cells.highlight A list of character or numeric vectors of cells to
#' highlight. If only one group of cells desired, can simply
#' pass a vector instead of a list. If set, colors selected cells to the color(s)
#' in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#' will also resize to the size(s) passed to \code{sizes.highlight}
#' @param cols.highlight A vector of colors to highlight the cells as; will
#' repeat to the length groups in cells.highlight
#' @param sizes.highlight Size of highlighted cells; will repeat to the length
#' groups in cells.highlight
#' @param na.value Color value for NA points when using custom scale
#' @param axis.type Select one of four axis types
#' @param line.size Size of axis line
#' @param axis.title.size,axis.ticks.size,axis.text.size Parameters about axis to be passed to \code{\link{theme}}
#' @param show.legend Whether or not to show legend
#' @param legend.title A vector of user-defined legend titles.
#' @param legend.title.size,legend.text.size,legend.key.size Parameters about legend to be passed to \code{\link{theme}}
#' @param ncol Number of columns for display when combining plots
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @importFrom rlang !!
#' @importFrom ggplot2 facet_wrap vars sym
#' @importFrom patchwork wrap_plots
#'
#' @export
#'
#' @note For the old \code{do.hover} and \code{do.identify} functionality, please see
#' \code{HoverLocator} and \code{CellSelector}, respectively.
#'
#' @aliases TSNEPlot PCAPlot ICAPlot
#' @seealso \code{\link{FeaturePlot}} \code{\link{HoverLocator}}
#' \code{\link{CellSelector}} \code{\link{FetchData}}
#'
#' @examples
#' DimPlot(object = pbmc_small)
#' DimPlot(object = pbmc_small, split.by = 'ident')
#'
DimPlot <- function(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  axis.type = c('default', 'keep line', 'keep title', 'no axis'),
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  legend.position = c('right', 'bottom', 'top', 'none'),
  legend.title = NULL,
  legend.title.size = NULL,
  legend.text.size = NULL,
  legend.key.size = 0.5,
  ncol = NULL,
  combine = TRUE
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[['ident']] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% 'ident'
  if (!is.null(legend.title)) {
    if (length(x = legend.title) != length(x = group.by)) {
      stop("If specifying the legend titles, please specify one title per grouping variable")
    }
    names(legend.title) <- group.by
  }
  data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        axis.type = axis.type,
        line.size = line.size,
        axis.title.size = axis.title.size,
        axis.ticks.size = axis.ticks.size,
        axis.text.size = axis.text.size,
        legend.position = legend.position,
        legend.title = legend.title[[x]],
        legend.title.size = legend.title.size,
        legend.text.size = legend.text.size,
        legend.key.size = legend.key.size
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + Seurat:::FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      return(plot)
    }
  )
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}



#' Visualize 'features' on a dimensional reduction plot
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#' @inheritParams DimPlot
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if
#' cells expressing given feature are getting buried.
#' @param features Vector of features to plot. Features can come from:
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#'     \item A column name from a \code{DimReduc} object corresponding to the cell embedding values
#'     (e.g. the PC 1 scores - "PC_1")
#' }
#' @param cols The two colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high. Also accepts a Brewer
#' color scale or vector of colors. Note: this will bin the data into number of colors provided.
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param split.by A factor in object metadata to split the feature plot by, pass 'ident'
#'  to split by cell identity'; similar to the old \code{FeatureHeatmap}
#' @param slot Which slot to pull expression data from?
#' @param sort.cell If \code{TRUE}, the positive cells will overlap the negative cells
#' @param features.position Where the feature label should be placed
#' @param title.size,title.face Size and face of feature labels
#'
#' @import Seurat
#'
#' @importFrom grDevices rgb
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 labs scale_x_continuous scale_y_continuous theme element_rect facet_grid facet_wrap theme_bw
#' dup_axis guides element_blank element_text margin scale_color_brewer scale_color_gradientn
#' scale_color_manual coord_fixed ggtitle
#'
#' @export
#'
MultiFeaturePlot <- function(
  object,
  features,
  dims = c(1, 2),
  cells = NULL,
  cols = c("lightgrey", "blue"),
  pt.size = NULL,
  order = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  shape.by = NULL,
  slot = "data",
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  ncol = NULL,
  sort.cell = FALSE,
  features.position = c('row', 'col'),
  title.size = 12,
  title.face = c('bold', 'plain', 'italic', 'bold.italic'),
  axis.type = c('default', 'keep line', 'keep title', 'no axis'),
  line.size = 0.5,
  legend.title = NULL
) {
  features.position <- match.arg(arg = features.position)
  title.face <- match.arg(arg = title.face)
  # Set a theme to remove right-hand Y axis lines
  # Also sets right-hand Y axis text label formatting
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank()
  )
  # Get the DimReduc to use
  reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  # Name the reductions
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  # Get plotting data
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident', features),
    cells = cells,
    slot = slot
  )
  # Check presence of features/dimensions
  if (ncol(x = data) < 4) {
    stop(
      "None of the requested features were found: ",
      paste(features, collapse = ', '),
      " in slot ",
      slot,
      call. = FALSE
    )
  } else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols, ]$maxcolors,
    no = length(x = cols)
  )
  # Apply cutoffs
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
      max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      }
      else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    Seurat:::RandomName()
  } else {
    switch(
      EXPR = split.by,
      ident = Idents(object = object)[cells],
      object[[split.by, drop = TRUE]][cells]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  # Set shaping variable
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  # Apply common limits
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  # Combine data
  data.plot <- data
  data.plot <- lapply(
    X = features,
    FUN = function(x) {
      df <- data.plot[, c(dims, 'ident', 'split')]
      df$features <- x
      df$value <- (data.plot[, x] - min(data.plot[, x])) / (max(data.plot[, x]) - min(data.plot[, x]))
      if (sort.cell) {
        df <- df[order(df$value), ]
      }
      return(df)
    }
  )
  data.plot <- do.call(rbind, data.plot)
  data.plot$features <- factor(data.plot$features, levels = features)
  breaks <- c(min(data.plot$value), max(data.plot$value))

  # Set expression fomula for facet
  if (length(x = unique(x = data.plot$split)) != 1) {
    split.by <- 'split'
  } else {
    split.by <- NULL
  }
  row <- switch(
    EXPR = features.position,
    'row' = 'features',
    'col' = split.by
  )
  col <- switch(
    EXPR = features.position,
    'row' = split.by,
    'col' = 'features'
  )
  if (!is.null(x = split.by)) {
    expr <- paste0(col, ' ~ ', row)
    facet <- facet_grid(expr, switch = 'y')
    y <- scale_y_continuous(position = 'right')
  } else {
    facet <- facet_wrap('~features', strip.position = 'top', ncol = ncol)
    y <- NULL
  }
  # Make plot
  plot <- suppressMessages(
    SingleDimPlot(
      data = data.plot,
      dims = dims,
      col.by = 'value',
      order = order,
      pt.size = pt.size,
      cols = cols,
      shape.by = shape.by,
      label = FALSE
    ) + y + facet +
      theme_bw(base_line_size = line.size, base_rect_size = line.size) +
      guides(color = NULL) +
      labs(colour = legend.title) +
      scale_color_gradientn(colors = cols, breaks = breaks, labels = c('Min', 'Max')) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )
  )
  if (label) {
    plot <- LabelClusters(
      plot = plot,
      id = 'ident',
      repel = repel,
      size = label.size
    )
  }
  facet.theme <- switch(
    EXPR = features.position,
    'row' = theme(strip.text.x = element_text(size = title.size, face = title.face),
                  strip.text.y = element_text(size = title.size, face = 'plain')),
    'col' = theme(strip.text.y = element_text(size = title.size, face = title.face),
                  strip.text.x = element_text(size = title.size, face = 'plain'))
  )
  axis.type <- match.arg(arg = axis.type)
  axis.theme <- switch(
    EXPR = axis.type,
    'default' = theme(),
    'keep line' = theme(axis.text = element_blank(), axis.ticks = element_blank()),
    'keep title' = theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()),
    'no axis' = NoAxes()
  )
  plot <- plot + facet.theme + axis.theme
  return(plot)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Expression by identity plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Single cell ridge plot
#'
#' Draws a ridge plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData)
#' @param cols Colors to use for plotting
#' @param idents Which classes to include in the plot (default is all)
#' @param sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted, can also pass 'increasing' or 'decreasing' to change sort direction
#' @param assay Name of assay to use, defaults to the active assay
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param log plot the feature axis on log scale
#' @param ncol Number of columns if multiple plots are displayed
#' @param x.lab Title of x axis
#' @param y.lab Title of y axis
#' @param line.size Size of axis line
#' @param features.face,features.size,axis.title.size,axis.ticks.size,axis.text.size Parameters to be passed to \code{\link{theme}}
#' @param slot Use non-normalized counts data for plotting
#' @param stacked Make stacked plot
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @export
#'
#' @examples
#' RidgePlot(object = pbmc_small, features = 'PC_1')
#'
RidgePlot <- function(
  object,
  features,
  cols = NULL,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  group.by = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  x.lab = NULL,
  y.lab = NULL,
  line.size = 0.5,
  features.face = c('bold', 'plain', 'italic', 'bold.italic'),
  features.size = NULL,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  slot = 'data',
  stacked = FALSE,
  combine = TRUE
) {
  features.face <- match.arg(arg = features.face)
  return(ExIPlot(
    object = object,
    type = 'ridge',
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    cols = cols,
    group.by = group.by,
    log = log,
    x.lab = x.lab,
    y.lab = y.lab,
    line.size = line.size,
    features.face = features.face,
    features.size = features.size,
    axis.title.size = axis.title.size,
    axis.ticks.size = axis.ticks.size,
    axis.text.size = axis.text.size,
    slot = slot,
    stacked = stacked,
    combine = combine
  ))
}



#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @inheritParams RidgePlot
#' @param pt.size Point size for geom_violin
#' @param split.by A variable to split the violin plots by,
#' @param multi.group  plot each group of the split violin plots by multiple or single violin shapes
#' see \code{\link{FetchData}} for more details
#' @param adjust Adjust parameter for geom_violin
#'
#' @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
#' \code{combine = TRUE}; otherwise, a list of ggplot objects
#'
#' @export
#'
#' @seealso \code{\link{FetchData}}
#'
#' @examples
#' VlnPlot(object = pbmc_small, features = 'PC_1')
#' VlnPlot(object = pbmc_small, features = 'LYZ', split.by = 'groups')
#'
VlnPlot <- function(
  object,
  features,
  cols = NULL,
  pt.size = 1,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  group.by = NULL,
  split.by = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  slot = 'data',
  x.lab = NULL,
  y.lab = NULL,
  features.face = c('bold', 'plain', 'italic', 'bold.italic'),
  features.size = NULL,
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  direction = c('vertical', 'horizontal'),
  stacked = FALSE,
  multi.group = FALSE,
  combine = TRUE
) {
  direction <- match.arg(arg = direction)
  features.face <- match.arg(arg = features.face)
  return(ExIPlot(
    object = object,
    type = ifelse(test = multi.group, yes = 'multiViolin', no = 'violin'),
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    adjust = adjust,
    pt.size = pt.size,
    cols = cols,
    group.by = group.by,
    split.by = split.by,
    log = log,
    x.lab = x.lab,
    y.lab = y.lab,
    features.face = features.face,
    features.size = features.size,
    line.size = line.size,
    axis.title.size = axis.title.size,
    axis.ticks.size = axis.ticks.size,
    axis.text.size = axis.text.size,
    slot = slot,
    direction = direction,
    stacked = stacked,
    combine = combine
  ))
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Heatmaps
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Feature expression heatmap
#'
#' Draws a heatmap of single cell feature expression. (Modified from Seurat)
#'
#' @param object Seurat object
#' @param features A vector of features to plot, defaults to \code{VariableFeatures(object = object)}
#' @param cells A vector of cells to plot
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped); defaults to 2.5
#' if \code{slot} is 'scale.data', 6 otherwise
#' @param group.by A vector of variables to group cells by; pass 'ident' to group by cell identity classes
#' @param group.bar Add a color bar showing group status for cells
#' @param group.colors Colors to use for the color bar
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param assay Assay to pull from
# @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param label Label the cell identies above the color bar
#' @param size Size of text above color bar
#' @param hjust Horizontal justification of text above color bar
#' @param angle Angle of text above color bar
#' @param raster If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on
#' some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE
#' if you are encountering that issue (note that plots may take longer to produce/render).
#' @param draw.lines Include white lines to separate the groups
#' @param lines.width Integer number to adjust the width of the separating white lines.
#' Corresponds to the number of "cells" between each group.
#' @param group.bar.height Scale the height of the color bar
#' @param direction Direction of heatmap, option can be "vertical" (default) or "horizontal".
#' @param expr.title Title of scaled expression legend bar
#' @param identity.legend Whether or not to show legend of identity
#' @param ident.title Customed legend title of identity
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
#' ggplot object. If \code{FALSE}, return a list of ggplot objects
#'
#' @return A ggplot object
#'
#' @import Seurat
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian scale_color_manual
#' ggplot_build aes_string
#' @export
#'
#' @examples
#' DoHeatmap(object = pbmc_small)
#'
DoHeatmap <- function(
  object,
  features = NULL,
  cells = NULL,
  group.by = 'ident',
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = 'scale.data',
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  direction = c('vertical', 'horizontal'),
  expr.title = 'Expression',
  identity.legend = TRUE,
  ident.title = 'Identity',
  combine = TRUE
) {
  direction <- match.arg(arg = direction)
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  # make sure features are present
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot,
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(
    object = object,
    slot = slot)[features, cells, drop = FALSE])))

  object <- suppressMessages(expr = StashIdent(object = object, save.name = 'ident'))
  group.by <- group.by %||% 'ident'
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  # group.use <- switch(
  #   EXPR = group.by,
  #   'ident' = Idents(object = object),
  #   object[[group.by, drop = TRUE]]
  # )
  # group.use <- factor(x = group.use[cells])
  plots <- vector(mode = 'list', length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      # create fake cells to serve as the white lines, fill with NAs
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 0.0025)
      placeholder.cells <- sapply(
        X = 1:(length(x = levels(x = group.use)) * lines.width),
        FUN = function(x) {
          return(Seurat:::RandomName(length = 20))
        }
      )
      placeholder.groups <- rep(x = levels(x = group.use), times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups), levels = group.levels)
      na.data.group <- matrix(
        data = NA,
        nrow = length(x = placeholder.cells),
        ncol = ncol(x = data.group),
        dimnames = list(placeholder.cells, colnames(x = data.group))
      )
      data.group <- rbind(data.group, na.data.group)
    }
    lgroup <- length(levels(group.use))
    plot <- SingleRasterMap(
      data = data.group,
      raster = raster,
      disp.min = disp.min,
      disp.max = disp.max,
      feature.order = features,
      cell.order = names(x = sort(x = group.use)),
      group.by = group.use,
      direction = direction,
      identity.legend = identity.legend,
      ident.title = ident.title,
      expr.title = expr.title
    )
    if (group.bar) {
      # TODO: Change group.bar to annotation.bar
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      cols <- group.colors[1:length(x = levels(x = group.use))] %||% default.colors
      if (any(is.na(x = cols))) {
        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
        cols <- Col2Hex(cols)
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
          x = cols,
          start = 1,
          stop = 7
        )))))
        through <- length(x = default.colors)
        while (length(x = col.dups) > 0) {
          pal.max <- length(x = col.dups) + through
          cols.extra <- hue_pal()(pal.max)[(through + 1):pal.max]
          cols[col.dups] <- cols.extra
          col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
            x = cols,
            start = 1,
            stop = 7
          )))))
        }
      }
      group.use2 <- sort(x = group.use)
      if (draw.lines) {
        na.group <- Seurat:::RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2), na.group)
        group.use2[placeholder.cells] <- na.group
        cols <- c(cols, "#FFFFFF")
      }
      pbuild <- ggplot_build(plot = plot)
      names(x = cols) <- levels(x = group.use2)
      # scale the height of the bar
      switch(
        EXPR = direction,
        'vertical' = {
          y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
          y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
          y.max <- y.pos + group.bar.height * y.range
          plot <- plot + annotation_raster(
            raster = t(x = cols[group.use2]),
            xmin = -Inf,
            xmax = Inf,
            ymin = y.pos,
            ymax = y.max
          ) +
            coord_cartesian(
              ylim = c(0, y.max),
              clip = 'off'
            ) +
            scale_color_manual(values = cols)
          if (label) {
            x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
            # Attempt to pull xdivs from x.major in ggplot2 < 3.3.0; if NULL, pull from the >= 3.3.0 slot
            x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
            x <- data.frame(group = sort(x = group.use), x = x.divs)
            label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = median) * x.max
            label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
            plot <- plot + geom_text(
              stat = "identity",
              data = label.x.pos,
              aes_string(label = 'group', x = 'label.x.pos'),
              y = y.max + y.max * 0.03 * 0.5,
              angle = angle,
              hjust = hjust,
              size = size
            )
            plot <- suppressMessages(plot + coord_cartesian(
              ylim = c(0, y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * size),
              clip = 'off')
            )
          }
        },
        'horizontal' = {
          y.range <- diff(x = pbuild$layout$panel_params[[1]]$x.range)
          y.pos <- min(pbuild$layout$panel_params[[1]]$x.range) - y.range * 0.015
          y.max <- y.pos - group.bar.height * y.range
          plot <- plot + annotation_raster(
            raster = cols[group.use2],
            ymin = -Inf,
            ymax = Inf,
            xmin = y.pos,
            xmax = y.max
          ) +
            coord_cartesian(
              xlim = c(y.max, max(pbuild$layout$panel_params[[1]]$x.range)),
              clip = 'off'
            ) +
            scale_color_manual(values = cols)
          if (label) {
            x.max <- max(pbuild$layout$panel_params[[1]]$y.range)
            x.divs <- rev(x = pbuild$layout$panel_params[[1]]$y.major) %||% pbuild$layout$panel_params[[1]]$y$break_positions()
            x <- data.frame(group = sort(x = rev(group.use)), x = rev(x.divs))
            blank <- (diff(pbuild$layout$panel_params[[1]]$y.range) - length(x = x.divs)) / 2
            label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = median) * diff(pbuild$layout$panel_params[[1]]$y.range) +
              min(pbuild$layout$panel_params[[1]]$y.range) + blank
            label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
            max.char <- max(nchar(x = as.character(x = label.x.pos$group)))
            min.char <- min(nchar(x = as.character(x = label.x.pos$group)))
            plot <- plot + geom_text(
              stat = "identity",
              data = label.x.pos,
              aes_string(label = 'group', y = 'label.x.pos'),
              x = y.max - (max(pbuild$layout$panel_params[[1]]$x.range) - y.max) * 0.5 * 0.03,
              angle = 0,
              hjust = 'right',
              vjust = -hjust,
              size = size
            ) + theme(plot.margin = unit(c(0, 0, 0, max.char - min.char), 'lines'))
            plot <- suppressMessages(plot + coord_cartesian(
              xlim = c(y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * size, max(pbuild$layout$panel_params[[1]]$x.range)),
              clip = 'off')
            )
          }
        }
      )
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Plot a single dimension
#
# @param data Data to plot
# @param dims A two-length numeric vector with dimensions to use
# @param pt.size Adjust point size for plotting
# @param col.by ...
# @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character
# or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
# By default, ggplot2 assigns colors
# @param shape.by If NULL, all points are circles (default). You can specify any cell attribute
# (that can be pulled with FetchData) allowing for both different colors and different shapes on
# cells.
# @param order Specify the order of plotting for the idents. This can be useful for crowded plots if
# points of interest are being buried. Provide either a full list of valid idents or a subset to be
# plotted last (on top).
# @param label Whether to label the clusters
# @param repel Repel labels
# @param label.size Sets size of labels
# @param cells.highlight A list of character or numeric vectors of cells to
# highlight. If only one group of cells desired, can simply
# pass a vector instead of a list. If set, colors selected cells to the color(s)
# in \code{cols.highlight} and other cells black (white if dark.theme = TRUE);
#  will also resize to the size(s) passed to \code{sizes.highlight}
# @param cols.highlight A vector of colors to highlight the cells as; will
# repeat to the length groups in cells.highlight
# @param sizes.highlight Size of highlighted cells; will repeat to the length
# groups in cells.highlight
# @param na.value Color value for NA points when using custom scale.
#
#' @import Seurat
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 ggplot aes_string labs geom_text guides
#'  scale_color_brewer scale_color_manual element_rect guide_legend discrete_scale
#'
SingleDimPlot <- function(
  data,
  dims,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  shape.by = NULL,
  order = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  cells.highlight = NULL,
  cols.highlight = '#DE2D26',
  sizes.highlight = 1,
  na.value = 'grey50',
  axis.type = c('default', 'keep line', 'keep title', 'no axis'),
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  legend.position = c('right', 'bottom', 'top', 'none'),
  legend.title = NULL,
  legend.title.size = NULL,
  legend.text.size = NULL,
  legend.key.size = 0.5
) {
  pt.size <- pt.size %||% Seurat:::AutoPointSize(data = data)
  axis.type <- match.arg(arg = axis.type)
  legend.position <- match.arg(arg = legend.position)
  axis.theme1 <- theme(
    axis.title = element_text(size = axis.title.size),
    axis.ticks = element_line(size = axis.ticks.size),
    axis.text = element_text(size = axis.text.size)
  )
  axis.theme2 <- switch(
    EXPR = axis.type,
    'default' = theme(),
    'keep line' = theme(axis.text = element_blank(), axis.ticks = element_blank()),
    'keep title' = theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()),
    'no axis' = Seurat::NoAxes()
  )
  legend.theme <- theme(
    legend.position = legend.position,
    legend.text = element_text(size = legend.text.size),
    legend.title = element_text(size = legend.title.size),
    legend.key.size = unit(legend.key.size, 'lines')
  )
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  } else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- Seurat:::SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = sizes.highlight %||% pt.size,
      cols.highlight = cols.highlight,
      col.base = cols[1] %||% '#C3C3C3',
      pt.size = pt.size
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    } else {
      order <- rev(x = c(
        order,
        setdiff(x = unique(x = data[, col.by]), y = order)
      ))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  } else {
    # col.index <- grep(pattern = col.by, x = colnames(x = data), fixed = TRUE)
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = '^\\d', x = col.by)) {
      # Do something for numbers
      col.by <- paste0('x', col.by)
    } else if (grepl(pattern = '-', x = col.by)) {
      # Do something for dashes
      col.by <- gsub(pattern = '-', replacement = '.', x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  plot <- ggplot(data = data) +
    geom_point(
      mapping = aes_string(
        x = dims[1],
        y = dims[2],
        color = paste0("`", col.by, "`"),
        shape = shape.by
      ),
      size = pt.size
    ) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(colour = legend.title)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(
      plot = plot,
      id = col.by,
      repel = repel,
      size = label.size
    )
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    } else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot(line_size = line.size) + axis.theme1 + axis.theme2 + legend.theme
  return(plot)
}



# Plot feature expression by identity
#
# Basically combines the codebase for VlnPlot and RidgePlot
#
# @param object Seurat object
# @param type Plot type, choose from 'ridge', 'violin', or 'multiViolin'
# @param features Features to plot (gene expression, metrics, PC scores,
# anything that can be retreived by FetchData)
# @param idents Which classes to include in the plot (default is all)
# @param ncol Number of columns if multiple plots are displayed
# @param sort Sort identity classes (on the x-axis) by the average expression of the attribute being potted
# @param y.max Maximum y axis value
# @param same.y.lims Set all the y-axis limits to the same values
# @param adjust Adjust parameter for geom_violin
# @param pt.size Point size for geom_violin
# @param cols Colors to use for plotting
# @param group.by Group (color) cells in different ways (for example, orig.ident)
# @param split.by A variable to split the plot by
# @param log plot Y axis on log scale
# @param slot Use non-normalized counts data for plotting
# @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed}
# ggplot object. If \code{FALSE}, return a list of ggplot objects
#
# @return A \code{\link[patchwork]{patchwork}ed} ggplot object if
# \code{combine = TRUE}; otherwise, a list of ggplot objects
#
#' @importFrom scales hue_pal
#' @importFrom ggplot2 xlab ylab
#' @importFrom patchwork wrap_plots
#
ExIPlot <- function(
  object,
  features,
  type = 'violin',
  idents = NULL,
  ncol = NULL,
  sort = FALSE,
  assay = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  adjust = 1,
  cols = NULL,
  pt.size = 0,
  group.by = NULL,
  split.by = NULL,
  log = FALSE,
  x.lab = NULL,
  y.lab = NULL,
  features.face = c('bold', 'plain', 'italic', 'bold.italic'),
  features.size = NULL,
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  slot = 'data',
  direction = c('vertical', 'horizontal'),
  stacked = FALSE,
  combine = TRUE
) {
  direction = match.arg(arg = direction)
  features.face <- match.arg(arg = features.face)
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  data <- FetchData(object = object, vars = features, slot = slot)
  features <- colnames(x = data)
  ncol <- ncol %||% ifelse(
    test = length(x = features) > 9,
    yes = 4,
    no = min(length(x = features), 3)
  )
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (type == 'ridge') {
    direction <- 'horizontal'
    split.by <- NULL
  }
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Seurat:::Interleave(cols, Seurat:::InvertHex(hexadecimal = cols))
    } else if (length(x = cols) == 1 && cols == 'interaction') {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    } else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Seurat:::Interleave(cols, Seurat:::InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- sort(x = levels(x = split))
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  label.fxn <- switch(
    EXPR = type,
    'violin' = ylab,
    "multiViolin" = ylab,
    'ridge' = xlab,
    stop("Unknown ExIPlot type ", type, call. = FALSE)
  )
  if (length(x = features) == 1){
    stacked <- FALSE
  }
  if (stacked) {
    data <- lapply(
      X = features,
      FUN = function(x) {
        df <- data.frame(cells = as.character(cells), ident = idents)
        df$features <- x
        df$value <- data[, x]
        return(df)
      }
    )
    data <- do.call(rbind, data)
    data$features <- factor(x = data$features, levels = features)
    plots <- SingleExIPlot(
      type = type,
      data = data,
      idents = idents,
      split = split,
      sort = sort,
      y.max = y.max,
      adjust = adjust,
      cols = cols,
      pt.size = pt.size,
      log = log,
      x.lab = x.lab,
      y.lab = y.lab,
      features.face = features.face,
      features.size = features.size,
      line.size = line.size,
      axis.title.size = axis.title.size,
      axis.ticks.size = axis.ticks.size,
      axis.text.size = axis.text.size,
      stacked = stacked,
      direction = direction
    )
    obj <- sapply(
      X = features,
      FUN = function(f) {
        key <- paste0(unlist(x = strsplit(x = f, split = '_'))[1], '_')
        obj <- names(x = which(x = Key(object = object) == key))
        return(obj)
      }
    )
    type.chk <- lapply(
      X = obj,
      FUN = function(x) {
        if (length(x = x) == 1) {
          if (inherits(x = object[[x]], what = 'DimReduc')) {
            return('DimReduc')
          } else if (inherits(x = object[[x]], what = 'Assay')) {
            return('Assay')
          }
        } else {
          return('Unknown')
        }
      }
    )
    if (all(type.chk == 'DimReduc')) {# From DimReduc
      plots <- plots + label.fxn(label = 'Embeddings Value')
    } else if (all(type.chk == 'Unknown') && any(!names(x = obj)[which(type.chk != 'Assay')] %in% rownames(x = object))) {# From meta.data
      warning("Unknown object type ", class(x = object), immediate. = TRUE, call. = FALSE)
      plots <- plots + label.fxn(label = NULL)
    } else if (all(names(x = obj)[which(type.chk != 'Assay')] %in% rownames(x = object))) {
      plots <- plots
    } else {
      plots <- plots + label.fxn(label = NULL)
    }
    return(plots)
  } else {
    plots <- lapply(
      X = features,
      FUN = function(x) {
        return(
          SingleExIPlot(
            type = type,
            data = data[, x, drop = FALSE],
            idents = idents,
            split = split,
            sort = sort,
            y.max = y.max,
            adjust = adjust,
            cols = cols,
            pt.size = pt.size,
            log = log,
            x.lab = x.lab,
            y.lab = y.lab,
            features.face = features.face,
            features.size = features.size,
            line.size = line.size,
            axis.title.size = axis.title.size,
            axis.ticks.size = axis.ticks.size,
            axis.text.size = axis.text.size,
            stacked = stacked,
            direction = direction
          )
        )
      }
    )
    for (i in 1:length(x = plots)) {
      key <- paste0(unlist(x = strsplit(x = features[i], split = '_'))[1], '_')
      obj <- names(x = which(x = Key(object = object) == key))
      if (length(x = obj) == 1) {
        if (inherits(x = object[[obj]], what = 'DimReduc')) {
          plots[[i]] <- plots[[i]] + label.fxn(label = 'Embeddings Value')
        } else if (inherits(x = object[[obj]], what = 'Assay')) {
          next
        } else {
          warning("Unknown object type ", class(x = object), immediate. = TRUE, call. = FALSE)
          plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
        }
      } else if (!features[i] %in% rownames(x = object)) {
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    }
    if (combine) {
      plots <- wrap_plots(plots, ncol = ncol)
      if (length(x = features) > 1) {
        plots <- plots & NoLegend()
      }
    }
    return(plots)
  }
}



# Plot a single expression by identity on a plot
#
# @param type Make either a 'ridge' or 'violin' plot
# @param data Data to plot
# @param idents Idents to use
# @param sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param adjust Adjust parameter for geom_violin
# @param cols Colors to use for plotting
# @param log plot Y axis on log scale
# @param seed.use Random seed to use. If NULL, don't set a seed
#
# @return A ggplot-based Expression-by-Identity plot
#
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes_string theme labs geom_violin geom_jitter ylim theme_bw facet_grid
#' scale_fill_manual scale_y_log10 scale_x_log10 scale_x_discrete scale_y_discrete scale_x_continuous waiver
#' @importFrom cowplot theme_cowplot
#'
SingleExIPlot <- function(
  data,
  idents,
  split = NULL,
  type = 'violin',
  sort = FALSE,
  y.max = NULL,
  adjust = 1,
  pt.size = 0,
  cols = NULL,
  seed.use = 42,
  log = FALSE,
  x.lab = NULL,
  y.lab = NULL,
  features.face = c('bold', 'plain', 'italic', 'bold.italic'),
  features.size = NULL,
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  stacked = FALSE,
  direction = c('vertical', 'horizontal')
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  features.face <- match.arg(arg = features.face)
  direction <- match.arg(arg = direction)
  x.lab <- x.lab %||% 'identity'
  y.lab <- y.lab %||% 'Expression Level'
  if (stacked) {
    if (!is.data.frame(x = data) || any(!c('features', 'ident', 'value') %in% colnames(x = data))) {
      stop('Stacked plot requires a data frame with "features", "idents" and "value" as column names')
    }
    sort <- FALSE
    for (f in levels(x = data$features)) {
      if (log) {
        noise <- rnorm(n = length(x = data[data$features == f, 'value'])) / 200
        data[data$features == f, 'value'] <- data[data$features == f, 'value'] + 1
      } else {
        noise <- rnorm(n = length(x = data[data$features == f, 'value'])) / 100000
      }
      if (all(data[data$features == f, 'value'] == data[data$features == f, 'value'][1])) {
        warning(paste0("All cells have the same value of ", f, "."))
      } else{
        data[data$features == f, 'value'] <- data[data$features == f, 'value'] + noise
      }
    }
    if (!is.null(x = split)) {
      data$split <- NA
      for (f in levels(x = data$features)) {
        data[data$features == f, 'split'] <- as.character(split)
      }
      data$split <- factor(x = data$split, levels = levels(x = split))
    }
    y.min <- min(data[, 'value'])
    y.max <- y.max %||% max(data[, 'value'])
    feature <- NULL
  } else {
    if (!is.data.frame(x = data) || ncol(x = data) != 1) {
      stop("SingleExIPlot requires a data frame with 1 column")
    }
    feature <- colnames(x = data)
    data$ident <- idents
    if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
      data$ident <- factor(
        x = data$ident,
        levels = names(x = rev(x = sort(
          x = tapply(
            X = data[, feature],
            INDEX = data$ident,
            FUN = mean
          ),
          decreasing = grepl(pattern = paste0('^', tolower(x = sort)), x = 'decreasing')
        )))
      )
    }
    if (log) {
      noise <- rnorm(n = length(x = data[, feature])) / 200
      data[, feature] <- data[, feature] + 1
    } else {
      noise <- rnorm(n = length(x = data[, feature])) / 100000
    }
    if (all(data[, feature] == data[, feature][1])) {
      warning(paste0("All cells have the same value of ", feature, "."))
    } else{
      data[, feature] <- data[, feature] + noise
    }
    y.min <- min(data[, feature])
    y.max <- y.max %||% max(data[, feature])
    if (!is.null(x = split)) {
      data$split <- split
    }
  }
  title <- feature
  if (type == 'violin' && !is.null(x = split)) {
    vln.geom <- Seurat:::geom_split_violin
    fill <- 'split'
  } else if (type == 'multiViolin' && !is.null(x = split )) {
    vln.geom <- geom_violin
    fill <- 'split'
    type <- 'violin'
  } else {
    vln.geom <- geom_violin
    fill <- 'ident'
  }
  switch(
    EXPR = type,
    'violin' = {
      x <- 'ident'
      if (stacked) {
        y <- 'value'
      } else {
        y <- paste0("`", feature, "`")
      }
      xlab <- x.lab
      ylab <- y.lab
      geom <- list(
        vln.geom(scale = 'width', adjust = adjust, trim = TRUE),
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      jitter <- geom_jitter(height = 0, size = pt.size)
      log.scale <- scale_y_log10()
      axis.scale <- ylim
    },
    'ridge' = {
      if (stacked) {
        x <- 'value'
      } else {
        x <- paste0("`", feature, "`")
      }
      y <- 'ident'
      xlab <- y.lab
      ylab <- x.lab
      geom <- list(
        geom_density_ridges(scale = 4),
        theme_bw(base_line_size = line.size),
        scale_y_discrete(expand = c(0.01, 0), limits = rev(levels(x = idents))),
        scale_x_continuous(expand = c(0, 0))
      )
      jitter <- geom_jitter(width = 0, size = pt.size)
      log.scale <- scale_x_log10()
      axis.scale <- function(...) {
        invisible(x = NULL)
      }
    },
    stop("Unknown plot type: ", type)
  )
  if (stacked) {
    base.theme <- theme_bw(base_line_size = line.size)
  } else {
    base.theme <- theme_cowplot(line_size = line.size)
  }
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]
  ) +
    labs(x = xlab, y = ylab, title = title, fill = NULL) +
    base.theme
  if (type == 'violin' && direction == 'horizontal') {
    plot <- plot + coord_flip() + scale_x_discrete(limits = rev(levels(x = idents)))
  }
  plot <- do.call(what = '+', args = list(plot, geom))
  plot <- plot +
    theme(
      plot.title = element_text(hjust = 0.5, face = features.face, size = features.size),
      axis.title = element_text(size = axis.title.size),
      axis.ticks = element_line(size = axis.ticks.size),
      axis.text = element_text(size = axis.text.size)
    )
  plot <- plot + if (log) {
    log.scale
  } else {
    if (!stacked) {
      axis.scale(y.min, y.max)
    }
  }
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      } else {
        unlist(x = lapply(
          X = idents,
          FUN = function(pattern, x) {
            x.mod <- gsub(
              pattern = paste0(pattern, '.'),
              replacement = paste0(pattern, ': '),
              x = x,
              fixed = TRUE
            )
            x.keep <- grep(pattern = ': ', x = x.mod, fixed = TRUE)
            x.return <- x.mod[x.keep]
            names(x = x.return) <- x[x.keep]
            return(x.return)
          },
          x = unique(x = as.vector(x = data$split))
        ))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    } else {
      labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  if (stacked) {
    plot <- switch(
      EXPR = direction,
      'vertical' = plot + facet_grid(features ~ ., scales = 'free_y'),
      'horizontal' = plot + facet_grid(~ features, scales = 'free_x')
    )
    angle <- switch(
      EXPR = direction,
      'vertical' = 0,
      'horizontal' = 90,
    )
    if (type == 'ridge') {
      angle <- 0
    }
    stacked.theme <- theme(
      strip.text = element_text(angle = angle, face = features.face, size = features.size),
      strip.placement = 'outside',
      strip.background.x = element_rect(colour = "red", fill = "#FFFFFF"),
      strip.background.y = element_rect(colour = "red", fill = "#FFFFFF"),
      panel.spacing = unit(0.01, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    )
    plot <- plot + stacked.theme
  }
  return(plot)
}



# A single heatmap from ggplot2 using geom_raster
#
# @param data A matrix or data frame with data to plot
# @param raster switch between geom_raster and geom_tile
# @param cell.order ...
# @param feature.order ...
# @param cols A vector of colors to use
# @param disp.min Minimum display value (all values below are clipped)
# @param disp.max Maximum display value (all values above are clipped)
# @param limits A two-length numeric vector with the limits for colors on the plot
# @param group.by A vector to group cells by, should be one grouping identity per cell
#
#' @import Seurat
#'
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient
#' scale_fill_gradientn theme element_blank labs geom_point guides guide_legend geom_tile
#
SingleRasterMap <- function(
  data,
  raster = TRUE,
  cell.order = NULL,
  feature.order = NULL,
  colors = PurpleAndYellow(),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  group.by = NULL,
  direction = c('vertical', 'horizontal'),
  expr.title = 'Expression',
  identity.legend = TRUE,
  ident.title = 'Identity'
) {
  direction <- match.arg(arg = direction)
  data <- Seurat::MinMax(data = data, min = disp.min, max = disp.max)
  data <- Seurat:::Melt(x = t(x = data))
  colnames(x = data) <- c('Feature', 'Cell', expr.title)
  if (!is.null(x = feature.order)) {
    data$Feature <- factor(x = data$Feature, levels = unique(x = feature.order))
  }
  if (!is.null(x = cell.order)) {
    data$Cell <- factor(x = data$Cell, levels = unique(x = cell.order))
  }
  if (!is.null(x = group.by)) {
    data$Identity <- group.by[data$Cell]
    colnames(x = data) <- c('Feature', 'Cell', expr.title, ident.title)
  }
  limits <- limits %||% c(min(data$Expression), max(data$Expression))
  if (length(x = limits) != 2 || !is.numeric(x = limits)) {
    stop("limits' must be a two-length numeric vector")
  }
  my_geom <- ifelse(test = raster, yes = geom_raster, no = geom_tile)
  if (direction == 'vertical') {
    plot <- ggplot(data = data) +
      my_geom(mapping = aes_string(x = 'Cell', y = 'Feature', fill = expr.title)) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_fill_gradientn(limits = limits, colors = colors, na.value = "white") +
      labs(x = NULL, y = NULL, fill = group.by %iff% expr.title) +
      WhiteBackground() + NoAxes(keep.text = TRUE)
    if (!is.null(x = group.by)) {
      if (identity.legend) {
        plot <- plot + geom_point(
          mapping = aes_string(x = 'Cell', y = 'Feature', color = ident.title),
          alpha = 0
        ) +
          guides(color = guide_legend(override.aes = list(alpha = 1)))
      }
    }
  } else if (direction == 'horizontal') {
    plot <- ggplot(data = data) +
      my_geom(mapping = aes_string(x = 'Feature', y = 'Cell', fill = expr.title)) +
      scale_x_discrete(limits = rev(levels(data$Feature)), position = 'top') +  scale_y_discrete(limits = rev(levels(data$Cell))) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
            legend.position = 'bottom') +
      scale_fill_gradientn(limits = limits, colors = colors, na.value = "white") +
      labs(x = NULL, y = NULL, fill = group.by %iff% expr.title) +
      WhiteBackground() + NoAxes(keep.text = TRUE)
  }
  return(plot)
}
