#' @include internal.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Expression by identity plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Single cell ridge plot (Modified version)
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
#' @param x.lab Title of x axis
#' @param y.lab Title of y axis
#' @param check.italic Whether or not to use italic text on gene names
#' @param do.bold Whether or not to use bold font on plot title
#' @param line.size Size of axis line
#' @param plot.title.size,axis.title.size,axis.ticks.size,axis.text.size Parameters about axis to be passed to \code{\link{theme}}
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param log plot the feature axis on log scale
#' @param ncol Number of columns if multiple plots are displayed
#' @param combine Combine plots using cowplot::plot_grid or patchwork package. If set 'None', will return a list of ggplot objects.
#' @param slot Use non-normalized counts data for plotting
#' @param ... Extra parameters passed on to \code{\link{CombinePlots}}
#'
#' @export
#'
RidgePlot <- function(
  object,
  features,
  cols = NULL,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  group.by = NULL,
  x.lab = NULL,
  y.lab = NULL,
  check.italic = FALSE,
  do.bold = TRUE,
  plot.title.size = NULL,
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  combine = c('CombinePlots', 'patchwork', 'None'),
  slot = 'data',
  ...
) {
  combine <- match.arg(arg = combine)
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
    x.lab = x.lab,
    y.lab = y.lab,
    check.italic = check.italic,
    do.bold = do.bold,
    plot.title.size = plot.title.size,
    line.size = line.size,
    axis.title.size = axis.title.size,
    axis.ticks.size = axis.ticks.size,
    axis.text.size = axis.text.size,
    log = log,
    slot = slot,
    combine = combine,
    ...
  ))
}

#' Single cell violin plot (Modified version)
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.).
#'
#' @inheritParams RidgePlot
#' @param direction Direction of violin plot, option can be "vertical" (default) or "horizontal".
#' @param pt.size Point size for geom_violin
#' @param split.by A variable to split the violin plots by,
#' see \code{\link{FetchData}} for more details
#' @param adjust Adjust parameter for geom_violin
#'
#' @export
#'
#' @seealso \code{\link{FetchData}}
#'
VlnPlot <- function (
  object,
  features,
  cols = NULL,
  pt.size = 1,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  group.by = NULL,
  split.by = NULL,
  x.lab = NULL,
  y.lab = NULL,
  check.italic = FALSE,
  do.bold = TRUE,
  plot.title.size = NULL,
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  combine = c('CombinePlots', 'patchwork', 'None'),
  direction = c('vertical', 'horizontal'),
  slot = "data",
  ...
) {
  direction <- match.arg(arg = direction)
  combine <- match.arg(arg = combine)
  return(ExIPlot(
    object = object,
    type = 'violin',
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
    x.lab = x.lab,
    y.lab = y.lab,
    check.italic = check.italic,
    do.bold = do.bold,
    plot.title.size = plot.title.size,
    line.size = line.size,
    axis.title.size = axis.title.size,
    axis.ticks.size = axis.ticks.size,
    axis.text.size = axis.text.size,
    log = log,
    slot = slot,
    combine = combine,
    direction = direction,
    ...
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Visualize 'features' on a dimensional reduction plot, with a modified type of legend
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#' @param object A Seurat object
#' @param features Vector of features to plot. Features can come from:
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#'     \item A column name from a \code{DimReduc} object corresponding to the cell embedding values
#'     (e.g. the PC 1 scores - "PC_1")
#' }
#' @param ncol Number of columns to combine multiple feature plots to, ignored if split.by is not NULL
#' @param check.italic Whether or not to use italic text on gene names
#' @param do.bold Whether or not to use bold font on plot title
#' @param plot.title.size Plot title size
#' @param axis.type Select one of four axis types
#' @param line.size Size of axis line
#' @param axis.title.size,axis.ticks.size,axis.text.size Parameters about axis to be passed to \code{\link{theme}}
#' @param show.legend Whether or not to show legend
#' @param legend.title A user-defined title of legend
#' @param legend.title.size,legend.text.size,legend.key.size Parameters about legend to be passed to \code{\link{theme}}
#' @param ... Extra parameters passed on to \code{\link{FeaturePlot}}
#'
#' @return A ggplot object
#'
#' @import Seurat
#' @import ggplot2
#' @importFrom cowplot get_legend plot_grid theme_cowplot
#'
#' @export
#'
MultiFeaturePlot <- function(
  object,
  features,
  ncol = NULL,
  check.italic = FALSE,
  do.bold = TRUE,
  plot.title.size = NULL,
  axis.type = c('default', 'keep line', 'keep title', 'no axis'),
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  show.legend = TRUE,
  legend.title = NULL,
  legend.title.size = NULL,
  legend.text.size = NULL,
  legend.key.size = 0.5,
  ...
) {
  axis.type <- match.arg(arg = axis.type)
  ncol <- ncol %||% min(length(x = features), 4)

  plot.tmp <- FeaturePlot(object = object, features = features[1], ...)
  plot.legend <- get_legend(plot = plot.tmp + theme(legend.position = 'right'))
  cols <- as.character(plot.legend$grobs[[1]]$grobs[[2]]$raster)

  data.tmp <- plot.tmp[['data']][, features[1]]
  breaks <- c(min(data.tmp), max(data.tmp))
  plot.tmp <- suppressMessages(plot.tmp + labs(colour = legend.title) + scale_color_gradientn(colors = rev(cols), breaks = breaks, labels = c('Min', 'Max')) +
                                 theme(legend.text = element_text(size = legend.text.size),
                                       legend.title = element_text(size = legend.title.size),
                                       legend.key.size = unit(legend.key.size, 'lines')))
  plot.legend <- get_legend(plot = plot.tmp + theme(legend.position = 'right'))
  rm(plot.tmp, data.tmp, cols, breaks);gc(reset = TRUE)

  plots <- lapply(features, function(f) {
    face <- ifelse(
      test = do.bold,
      yes = 'bold',
      no = 'plain'
    )
    if (check.italic) {
      if (f %in% rownames(object)) {
        face <- ifelse(
          test = do.bold,
          yes = 'bold.italic',
          no = 'italic'
        )
      }
    }
    p <- FeaturePlot(object = object, features = f, ...) + theme_cowplot(line_size = line.size) + NoLegend() +
      theme(plot.title = element_text(face = face, size = plot.title.size, hjust = 0.5),
            axis.title = element_text(size = axis.title.size),
            axis.ticks = element_line(size = axis.ticks.size),
            axis.text = element_text(size = axis.text.size))
    axis.theme <- switch(
      EXPR = axis.type,
      'default' = theme(),
      'keep line' = theme(axis.text = element_blank(), axis.ticks = element_blank()),
      'keep title' = theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank()),
      'no axis' = NoAxes()
    )
    p <- p + axis.theme
    return(p)
  })

  plots.combined <- plot_grid(
    plotlist = plots,
    ncol = ncol,
    align = 'hv'
  )
  if (show.legend) {
    plots.combined <- plot_grid(
        plots.combined,
        plot.legend,
        rel_widths = c(3, 0.3)
    )
    return(plots.combined)
  } else {
    return(plots.combined)
  }
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
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work
#' when plotting multiple dimensions
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
            x.divs <- pbuild$layout$panel_params[[1]]$x.major
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
            x.divs <- rev(x = pbuild$layout$panel_params[[1]]$y.major)
            x <- data.frame(group = sort(x = group.use), x = x.divs)
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
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Plot feature expression by identity
#
# Basically combines the codebase for VlnPlot and RidgePlot
#
# @param object Seurat object
# @param features Features to plot (gene expression, metrics, PC scores,
# anything that can be retreived by FetchData)
# @param type Plot type, choose from 'ridge' or 'violin'
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
# @param x.lab Title of x axis
# @param y.lab Title of y axis
# @param check.italic Whether or not to use italic text on gene names
# @param do.bold Whether or not to use bold font on plot title
# @param line.size Size of axis line
# @param plot.title.size,axis.title.size,axis.ticks.size,axis.text.size Parameters about axis to be passed to \code{\link{theme}}
# @param log plot Y axis on log scale
# @param combine Combine plots using cowplot::plot_grid or patchwork package. If set 'None', will return a list of ggplot objects.
# @param direction Direction of violin plot, option can be "vertical" (default) or "horizontal".
# @param slot Use non-normalized counts data for plotting
# @param ... Extra parameters passed to \code{\link{CombinePlots}}
#
#' @importFrom scales hue_pal
#' @importFrom ggplot2 xlab ylab
#' @importFrom cowplot get_plot_component
#'
#' @import patchwork
#'
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
  x.lab = NULL,
  y.lab = NULL,
  check.italic = FALSE,
  do.bold = TRUE,
  plot.title.size = NULL,
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  log = FALSE,
  combine = c('CombinePlots', 'patchwork', 'None'),
  direction = c('vertical', 'horizontal'),
  slot = 'data',
  ...
) {
  direction <- match.arg(arg = direction)
  combine <- match.arg(arg = combine)
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  ncol <- ncol %||% ifelse(
    test = length(x = features) > 9,
    yes = 4,
    no = min(length(x = features), 3)
  )
  x.lab <- x.lab %||% group.by
  data <- FetchData(object = object, vars = features, slot = slot)
  features <- colnames(x = data)
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
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    } else if (length(x = cols) == 1 && cols == 'interaction') {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    } else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- sort(x = levels(x = split))
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  patchwork <- FALSE
  do.italic <- FALSE
  if (combine == 'patchwork') {
    patchwork <- TRUE
    pt.size <- 0
    if (check.italic) {
      if (features[1] %in% rownames(x = object)) {
        do.italic <- TRUE
      }
    }
    plot1 <- SingleExIPlot(
      type = type,
      data = data[, features[1], drop = FALSE],
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
      do.italic = do.italic,
      do.bold = do.bold,
      plot.title.size = plot.title.size,
      line.size = line.size,
      axis.title.size = axis.title.size,
      axis.ticks.size = axis.ticks.size,
      axis.text.size = axis.text.size,
      patchwork = patchwork,
      is.plot1 = TRUE,
      direction = direction
    )
    features <- features[-1]
  }
  plots <- lapply(
    X = features,
    FUN = function(x) {
      if (check.italic) {
        if (x %in% rownames(x = object)) {
          do.italic <- TRUE
        }
      }
      return(SingleExIPlot(
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
        do.italic = do.italic,
        do.bold = do.bold,
        plot.title.size = plot.title.size,
        line.size = line.size,
        axis.title.size = axis.title.size,
        axis.ticks.size = axis.ticks.size,
        axis.text.size = axis.text.size,
        patchwork = patchwork,
        is.plot1 = FALSE,
        direction = direction
      ))
    }
  )
  if (combine != 'patchwork') {
    label.fxn <- switch(
      EXPR = type,
      'violin' = ylab,
      'ridge' = xlab,
      stop("Unknown ExIPlot type ", type, call. = FALSE)
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
  }
  if (combine == 'patchwork') {
    all.plots <- switch(
      EXPR = type,
      'violin' = if (direction == 'vertical') {
        rev(c(list(plot1), plots))
      } else if (direction == 'horizontal') {
        c(list(plot1), plots)
      },
      'ridge' = c(list(plot1), plots)
    )
    plots <- all.plots[[1]]
    for (i in 2:length(all.plots)) {
      plots <- plots + all.plots[[i]]
    }
    plot.tmp <- suppressMessages(SingleExIPlot(
      type = type,
      data = data[, features[1], drop = FALSE],
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
      do.italic = do.italic,
      do.bold = do.bold,
      plot.title.size = plot.title.size,
      line.size = line.size,
      axis.title.size = axis.title.size,
      axis.ticks.size = axis.ticks.size,
      axis.text.size = axis.text.size,
      patchwork = FALSE,
      is.plot1 = FALSE,
      direction = direction
    ) +
      theme(axis.title.x.top = element_text(angle = 0, hjust = 0.5)))
    if (type == 'violin') {
      if (direction == 'horizontal') {
        y <- get_plot_component(plot = plot.tmp, pattern = "xlab-b")
        plots <- plots + plot_layout(nrow = 1, widths = rep(1, length(x = all.plots)))
        plots <- (plots) - y + plot_layout(nrow = 2, heights = c(1, 0.05))
      } else if (direction == 'vertical') {
        y <- get_plot_component(plot = plot.tmp, pattern = "ylab-l")
        plots <- plots + plot_layout(ncol = 1, heights = rep(1, times = length(x = all.plots)))
        plots <- (plots) - y + plot_layout(nrow = 1, widths = c(1, 0.05))
      }
    } else if (type == 'ridge') {
      y <- get_plot_component(plot = plot.tmp, pattern = "xlab-b")
      plots <- {plots + plot_layout(nrow = 1, widths = rep(1, times = length(x = all.plots)))} -
        y + plot_layout(nrow = 2, heights = c(1, 0.05), widths = 1)
    }
  }
  if (combine == 'CombinePlots') {
    combine.args <- list(
      'plots' = plots,
      'ncol' = ncol
    )
    combine.args <- c(combine.args, list(...))
    if (!'legend' %in% names(x = combine.args)) {
      combine.args$legend <- 'none'
    }
    plots <- do.call(what = 'CombinePlots', args = combine.args)
  }
  return(plots)
}

# Plot a single expression by identity on a plot
#
# @param data Data to plot
# @param idents Idents to use
# @param type Make either a 'ridge' or 'violin' plot
# @param sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param adjust Adjust parameter for geom_violin
# @param cols Colors to use for plotting
# @param log plot Y axis on log scale
# @param x.lab Title of x axis
# @param y.lab Title of y axis
# @param do.italic Whether or not to use italic text on gene names
# @param do.bold Whether or not to use bold font on plot title
# @param line.size Size of axis line
# @param plot.title.size,axis.title.size,axis.ticks.size,axis.text.size Parameters about axis to be passed to \code{\link{theme}}
# @param patchwork Whether or not to generate plot used for patchwork combination
# @param is.plot1 Whether or not to generate the initial plot used for patchwork combination
# @param direction Direction of violin plot, option can be "vertical" (default) or "horizontal".
#
# @return A ggplot-based Expression-by-Identity plot
#
# @import ggplot2
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggplot2 ggplot aes_string theme labs geom_violin geom_jitter ylim
#' scale_fill_manual scale_y_log10 scale_x_log10 scale_y_discrete scale_x_continuous waiver
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
  log = FALSE,
  x.lab = NULL,
  y.lab = NULL,
  do.italic = FALSE,
  do.bold = TRUE,
  plot.title.size = NULL,
  line.size = 0.5,
  axis.title.size = NULL,
  axis.ticks.size = NULL,
  axis.text.size = NULL,
  patchwork = FALSE,
  is.plot1 = FALSE,
  direction = c('vertical', 'horizontal')
) {
  direction <- match.arg(arg = direction)
  set.seed(seed = 42)
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
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
  x.lab <- x.lab %||% 'identity'
  y.lab <- y.lab %||% ifelse(test = log, yes = 'Log Expression Level', no = 'Expression Level')
  y.max <- y.max %||% max(data[, feature])
  if (is.null(x = split) || type != 'violin') {
    vln.geom <- geom_violin
    fill <- 'ident'
  } else {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- 'split'
  }
  face <- ifelse(
    test = do.bold,
    yes = 'bold',
    no = 'plain'
  )
  if (do.italic) {
    face <- ifelse(
      test = do.bold,
      yes = 'bold.italic',
      no = 'italic'
    )
  }
  theme <- theme_cowplot(line_size = line.size)
  if (patchwork) {
    theme <- theme +
      theme(plot.title = switch (
        EXPR = type,
        'violin' = element_blank(),
        'ridge' = element_text(hjust = 0.5, face = face, size = plot.title.size)
        ),
        legend.position = 'none')
    pt.size <- 0
    y.lab <- feature
  } else {
    theme <- theme + theme(plot.title = element_text(hjust = 0.5, face = face, size = plot.title.size))
  }
  switch(
    EXPR = type,
    'violin' = {
      x <- 'ident'
      y <- paste0("`", feature, "`")
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
      x <- paste0("`", feature, "`")
      y <- 'ident'
      xlab <- y.lab
      ylab <- x.lab
      geom <- list(
        geom_density_ridges(scale = 4),
        theme_ridges(line_size = line.size),
        theme(axis.title.x = element_text(hjust = 0.5),
              plot.title = element_text(face = face, size = plot.title.size)),
        scale_y_discrete(expand = c(0.01, 0)),
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
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]
  ) +
    labs(x = xlab, y = ylab, title = feature, fill = NULL) +
    theme
  plot <- do.call(what = '+', args = list(plot, geom))
  plot <- plot + if (log) {
    log.scale
  } else {
    axis.scale(min(data[, feature]), y.max)
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
  plot <- plot + theme(
    axis.title = element_text(size = axis.title.size),
    axis.ticks = element_line(size = axis.ticks.size),
    axis.text = element_text(size = axis.text.size)
  )
  if (patchwork) {
    if (type == 'violin') {
      face <- NULL
      if (do.italic) {
        face <- 'italic'
      }
      plot <- plot + theme(
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = axis.title.size, face = face),
        plot.margin = unit(c(0.5, 0, 0, 0), units = 'line')
      )
      if (!is.plot1) {
        plot <- plot + theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0.5, 0, -0.75, 0), units = 'line')
        )
      }
    } else if (type == 'ridge') {
      plot <- plot + theme(
        axis.title.x = element_blank(),
        legend.position = 'none'
      )
      if (!is.plot1) {
        plot <- plot + theme(
          axis.title.y = element_blank(),
          axis.line.y = element_line(),
          axis.text.y = element_blank()
        )
      }
    }
  }
  if (direction == 'horizontal' && type == 'violin') {
    idents.level <- levels(x = idents)[which(table(idents) > 0)]
    plot <- plot + coord_flip() + scale_x_discrete(limits = rev(idents.level)) +
      theme(plot.margin = unit(c(0, -0.1, 0, 0), units = 'lines'),
            axis.title = element_text(size = axis.title.size),
            axis.ticks = element_line(size = axis.ticks.size),
            axis.text = element_text(size = axis.text.size)
            )
    if (patchwork) {
      face <- NULL
      if (do.italic) {
        face <- 'italic'
      }
      plot <- suppressMessages(
        plot + scale_y_discrete(position = 'right') +
        theme(axis.title.x.top = element_text(angle = 30, vjust = 0, hjust = 0.1, face = face, size = axis.title.size),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = unit(c(0, 0.5, 0, -0.75), 'line'))
      )
      if (is.plot1) {
        plot <- plot + theme(
          axis.title.y = element_text(angle = 90, vjust = 2.5, size = axis.title.size),
          axis.text.y = element_text(hjust = 1, size = axis.text.size),
          plot.margin = unit(c(0, 0.5, 0, 0), 'line')
        )
      }
    }
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
