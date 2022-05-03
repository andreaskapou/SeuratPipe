#' @title Tailored dimensional reduction plot
#'
#' @description This function extends the DimPlot Seurat function by
#'   providing additional plotting options.
#'
#' @param reduction Which dimensionality reduction to use (required).
#' @param group.by Name of metadata column to group (color) cells by (required).
#' @param split.by Name of a metadata column to split plot by.
#' @param col_pal Discrete colour palette to use, default is Hue palette
#' (hue_pal) from 'scales' package. Should be of equal length to number of
#' groups in 'group.by'.
#' @param ncol Number of columns for display when combining plots.
#' @param dims_plot Dimensions to plot, must be a two-length numeric vector
#' specifying x- and y-dimensions.
#' @param combine Combine plots into a single patchworked ggplot object.
#'    If FALSE, return a list of ggplot objects
#' @param ... Additional parameters passed to Seurat's DimPlot.
#' @inheritParams cluster_analysis
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
dim_plot <- function(seu, reduction = "umap", group.by = "active.ident",
                     split.by = NULL, ncol = NULL, legend.position = "right",
                     col_pal = NULL, dims_plot = c(1, 2), pt.size = 1.4,
                     label = FALSE, label.size = 7, combine = TRUE, ...) {
  # So CMD passes
  D1 = D2 = ident = x = y <- NULL
  assertthat::assert_that(!is.null(reduction))
  assertthat::assert_that(!is.null(group.by))

  if (reduction %in% c("umap", "tsne")) {
    dims_plot <- c(1,2)
  } else {
    dims_plot <- dims_plot
  }

  # Extract dimensionally reduced data
  dim_dt <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  # Keep only specified dimensions and provide generic column names
  dim_dt <- dim_dt[, dims_plot]
  colnames(dim_dt) <- c("D1", "D2")
  # Extract x and y plotting limits
  xlim <- c(min(dim_dt[,1]) - 0.5, max(dim_dt[,1]) + 0.5)
  ylim <- c(min(dim_dt[,2]) - 0.5, max(dim_dt[,2]) + 0.5)

  if (!(is.character(seu@meta.data[[group.by]]) ||
        is.factor(seu@meta.data[[group.by]]) ||
        is.logical(seu@meta.data[[group.by]]) ) ) {
    stop("Error: 'group.by' should not be continuous in metadata slot")
  }
  if (!is.null(split.by)) {
    label <- FALSE
    if (!(is.character(seu@meta.data[[split.by]]) ||
          is.factor(seu@meta.data[[split.by]]) ||
          is.logical(seu@meta.data[[split.by]]) ) ) {
      stop("Error: split.by should not be continuous in metadata slot")
    }
  } else {
    ncol <- 1
  }
  group <- as.factor(seu@meta.data[[group.by]])
  if (is.null(col_pal)) {
    col_pal <- scales::hue_pal()(nlevels(group))
    names(col_pal) <- levels(group)
  } else {
    col_pal <- col_pal[1:nlevels(group)]
    names(col_pal) <- levels(group)
  }
  t <- Seurat::DimPlot(seu, reduction = reduction, dims = dims_plot,
        group.by = group.by, split.by = split.by, ncol = ncol, pt.size = 0,
        combine = combine, label.size = label.size, label = FALSE, ...)
  t <- t &
    ggplot2::theme(legend.position = legend.position,
          legend.key.size = ggplot2::unit(1.3,"line"),
          legend.text = ggplot2::element_text(size = 11),
          axis.line = ggplot2::element_line(colour = "black",
                                            size = 1.25, linetype = "solid"),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(size = 16),
          axis.title.y = ggplot2::element_text(size=16),
          axis.ticks.length = ggplot2::unit(0.2,"cm")) &
    ggplot2::geom_point(ggplot2::aes(t$data[, 1], t$data[, 2], fill = group),
               shape = 21, size = pt.size, stroke = 0.1) &
    ggplot2::guides(size = "none", colour = "none",
            fill = ggplot2::guide_legend(override.aes = list(size = 2))) &
    ggplot2::scale_fill_manual(name = NULL, values = col_pal) &
    ggplot2::scale_x_continuous(name = NULL, minor_breaks = NULL,
                                limits = xlim) &
    ggplot2::scale_y_continuous(name = NULL, minor_breaks = NULL,
                                limits = ylim)

  if (label) {
    # Compute cluster centres
    dim_dt$ident <- group
    centres <- dim_dt |> dplyr::group_by(ident) |>
      dplyr::summarize(x = median(D1), y = median(D2))
    t <- t + ggplot2::geom_text(centres,
              mapping = ggplot2::aes(x = x, y = y, label = ident),
              colour = "black", size = label.size)
  }
  return(t)
}


#' @title Tailored feature plot
#'
#' @description This function adapts the FeaturePlot Seurat function by
#'   providing additional plotting options.
#'
#' @param reduction Which dimensionality reduction to use (required).
#' @param features Vector of features to plot.
#' @param ncol Number of columns for display when having multiple features.
#' @param col_pal Continuous colour palette to use, default "RdYlBu".
#' @param dims_plot Dimensions to plot, must be a two-length numeric vector
#' specifying x- and y-dimensions.
#' @param combine Combine plots into a single patchworked ggplot object.
#'    If FALSE, return a list of ggplot objects
#' @param ... Additional parameters passed to Seurat's FeaturePlot.
#' @inheritParams cluster_analysis
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
feature_plot <- function(
    seu, reduction = "umap", features = "nFeature_RNA", max.cutoff = "q98",
    min.cutoff = NA, ncol = NULL, legend.position = "right", col_pal = NULL,
    dims_plot = c(1, 2), pt.size = 1.4, combine = TRUE, ...) {
  assertthat::assert_that(!is.null(reduction))
  assertthat::assert_that(!is.null(features))

  # Extract features present in the Seurat object
  features <- features[(features %in% rownames(seu)) |
                         (features %in% colnames(seu@meta.data))]
  if (length(features) == 0) {
    message("No features present to plot.")
    return(ggplot2::ggplot())
  }

  if (reduction %in% c("umap", "tsne")) {
    dims_plot <- c(1,2)
  } else {
    dims_plot <- dims_plot
  }

  # Extract dimensionally reduced data
  dim_dt <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  # Keep only specified dimensions and provide generic column names
  dim_dt <- dim_dt[, dims_plot]
  # Extract x and y plotting limits
  xlim <- c(min(dim_dt[,1]) - 0.5, max(dim_dt[,1]) + 0.5)
  ylim <- c(min(dim_dt[,2]) - 0.5, max(dim_dt[,2]) + 0.5)

  if (is.null(col_pal)) { col_pal = "RdYlBu" }

  key_height = 0.6
  key_width = 0.2
  if (legend.position == "top") {
    key_height = 0.2
    key_width = 0.7
  }

  t <- Seurat::FeaturePlot(
    seu, features = features, reduction = reduction, max.cutoff = max.cutoff,
    min.cutoff = min.cutoff, ncol = ncol, combine = combine, pt.size = NULL, ...)
  t <- t &
    ggplot2::theme(legend.position = legend.position,
          legend.key.size = ggplot2::unit(0.7,"line"),
          legend.key.height = ggplot2::unit(key_height, 'cm'),
          legend.key.width = ggplot2::unit(key_width, "cm"),
          legend.text = ggplot2::element_text(size = 5),
          axis.line = ggplot2::element_line(colour = "black", size = 1.25,
                                            linetype = "solid"),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(size = 16),
          axis.title.y = ggplot2::element_text(size = 16),
          axis.ticks.length = ggplot2::unit(0.2,"cm")) &
    #ggplot2::geom_point(ggplot2::aes(t$data[, 1], t$data[, 2]), shape = 3,
    #           size = pt.size, stroke = 0) &
    ggplot2::scale_colour_distiller(palette = col_pal) &
    ggplot2::scale_fill_distiller(palette = col_pal) &
    ggplot2::scale_x_continuous(name = NULL, minor_breaks = NULL,
                                limits = xlim) &
    ggplot2::scale_y_continuous(name = NULL, minor_breaks = NULL,
                                limits = ylim)
  return(t)
}


#' @title Subset dim plot
#'
#' @description This function splits the dim plots according to a grouping
#' variable defined in `subset.by`. The main difference is that the whole
#' dataset is also shown in the background as smaller points.
#'
#' @param subset.by Metadata column whose unique values will generate split
#' dim plots.
#' @param col_pal Discrete colour palette to use.
#' @inheritParams subset_feature_plot
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
subset_dim_plot <- function(
    seu, subset.by, reduction = "umap", ncol = NULL, col_pal = NULL,
    pt.size = 2, stroke = 0.05, back.pt.size = 0.5, back.alpha = 0.1,
    back.color = "grey", combine = TRUE) {

  # Extract dimensionally reduced data from full object
  dim_dt <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  # Keep only specified dimensions and provide generic column names
  dim_dt <- dim_dt[, c(1,2)]
  # Extract x and y plotting limits
  xlim <- c(min(dim_dt[,1]) - 0.5, max(dim_dt[,1]) + 0.5)
  ylim <- c(min(dim_dt[,2]) - 0.5, max(dim_dt[,2]) + 0.5)
  # Add cell ID as column
  dim_dt[["subset"]] <- seu[[subset.by]][[1]]

  # Split to list
  dim_dt_list <- dim_dt |> dplyr::group_by(subset) |> dplyr::group_split()
  names(dim_dt_list) <- unique(seu[[subset.by]][[1]])

  if (is.null(ncol)) {
    if (length(dim_dt_list) <= 3) {
      ncol <- length(dim_dt_list)
    } else {
      ncol <- 4
    }
  }

  if (is.null(col_pal)) {
    col_pal <- scales::hue_pal()(length(dim_dt_list))
    names(col_pal) <- names(dim_dt_list)
  } else {
    col_pal <- col_pal[1:names(dim_dt_list)]
    names(col_pal) <- names(dim_dt_list)
  }

  gg_list <- list()
  for (s in names(dim_dt_list)) {
    gg_list[[s]] <- local({
      s <- s
      dim_subset <- as.data.frame(dim_dt_list[[s]])
      gg <- ggplot2::ggplot(dim_dt, ggplot2::aes(x = dim_dt[, 1],
                                                 y = dim_dt[, 2],
                                                 fill = dim_dt[, 3])) +
        ggplot2::geom_point(size = back.pt.size, alpha = back.alpha,
                            color = back.color) +
        ggplot2::geom_point(data = dim_subset,
                            mapping = aes(x = dim_subset[, 1],
                                          y = dim_subset[, 2],
                                          fill = dim_subset[, 3]),
                            shape = 21, size = pt.size, stroke = stroke) +
        ggplot2::scale_fill_manual(name = NULL, values = col_pal) +
        ggplot2::scale_x_continuous(name = NULL, minor_breaks = NULL,
                                    limits = xlim) +
        ggplot2::scale_y_continuous(name = NULL, minor_breaks = NULL,
                                    limits = ylim) +
        ggplot2::labs(title = s) +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none",
                       plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
      return(gg)
    })
  }
  if (combine) {
    gg_list <- patchwork::wrap_plots(gg_list, ncol = ncol)
  }
  return(gg_list)
}


#' @title Subset feature plot
#'
#' @description This function splits the feature plots according to a grouping
#' variable defined in `subset.by`. The main difference is that the whole
#' dataset is also shown in the background as smaller points.
#'
#' @param seu Seurat object
#' @param subset.by Metadata column whose unique values will generate split
#' feature plots.
#' @param feature Feature to plot.
#' @param min.cutoff Minimum cutoff value for feature, may specify quantile
#' in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10').
#' @param max.cutoff Maximum cutoff value for feature, may specify quantile
#' in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10').
#' @param reduction Dimensionality reduction to use.
#' @param slot Slot to extract data from.
#' @param ncol Number of columns for display when having multiple features.
#' @param col_pal Continuous colour palette to use, default "RdYlBu".
#' @param pt.size Adjust point size for plotting.
#' @param stroke Stroke value for each point.
#' @param legend.position Position of legend, default "right" (set to "none"
#' for clean plot).
#' @param back.pt.size Adjust background point size for plotting.
#' @param back.alpha Adjust opacity for background points.
#' @param back.color Colour for background points.
#' @param combine Combine plots into a single patchworked ggplot object.
#'    If FALSE, return a list of ggplot objects
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
subset_feature_plot <- function(
    seu, subset.by, feature, max.cutoff = "q98", min.cutoff = NA,
    reduction = "umap", slot = "data", ncol = NULL, col_pal = NULL,
    pt.size = 2, stroke = 0.05, legend.position = "right",
    back.pt.size = 0.5, back.alpha = 0.1, back.color = "grey", combine = TRUE) {

  # Extract dimensionally reduced data from full object
  dim_dt <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  # Keep only specified dimensions and provide generic column names
  dim_dt <- dim_dt[, c(1,2)]
  # Extract x and y plotting limits
  xlim <- c(min(dim_dt[,1]) - 0.5, max(dim_dt[,1]) + 0.5)
  ylim <- c(min(dim_dt[,2]) - 0.5, max(dim_dt[,2]) + 0.5)
  # Add feature to plot as 3rd column
  dim_dt[, 3] <- Seurat::FetchData(object = seu, vars = feature, slot = slot)
  # Add cell ID as column
  dim_dt[["subset"]] <- seu[[subset.by]][[1]]

  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(dim_dt[, 3]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = feature
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(dim_dt[, 3]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = feature
  )

  # Apply cutoffs
  data.feature <- as.vector(x = dim_dt[, 3])
  min.use <- Seurat::SetQuantile(cutoff = min.cutoff, data.feature)
  max.use <- Seurat::SetQuantile(cutoff = max.cutoff, data.feature)
  data.feature[data.feature < min.use] <- min.use
  data.feature[data.feature > max.use] <- max.use
  dim_dt[, 3] <- data.feature

  # Split to list
  dim_dt_list <- dim_dt |> dplyr::group_by(subset) |> dplyr::group_split()
  names(dim_dt_list) <- unique(seu[[subset.by]][[1]])

  if (is.null(ncol)) {
    if (length(dim_dt_list) <= 3) {
      ncol <- length(dim_dt_list)
    } else {
      ncol <- 4
    }
  }

  if (is.null(col_pal)) { col_pal = "RdYlBu" }

  key_height = 0.6
  key_width = 0.2
  if (legend.position == "top") {
    key_height = 0.2
    key_width = 0.7
  }

  gg_list <- list()
  for (s in names(dim_dt_list)) {
    gg_list[[s]] <- local({
      s <- s
      dim_subset <- as.data.frame(dim_dt_list[[s]])
      gg <- ggplot2::ggplot(dim_dt, ggplot2::aes(x = dim_dt[, 1], y = dim_dt[, 2])) +
        ggplot2::geom_point(size = back.pt.size, alpha = back.alpha, color = back.color) +
        ggplot2::geom_point(data = dim_subset,
                            mapping = aes(x = dim_subset[, 1],
                                          y = dim_subset[, 2],
                                          fill = dim_subset[, 3]),
                            shape = 21, size = pt.size, stroke = stroke) +
        ggplot2::scale_fill_distiller(palette = col_pal, limits = c(min.use, max.use)) +
        ggplot2::scale_x_continuous(name = NULL, minor_breaks = NULL,
                                    limits = xlim) +
        ggplot2::scale_y_continuous(name = NULL, minor_breaks = NULL,
                                    limits = ylim) +
        ggplot2::labs(fill = feature, title = s) +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = legend.position,
                       legend.key.size = ggplot2::unit(0.7, "line"),
                       legend.key.height = ggplot2::unit(key_height, "cm"),
                       legend.key.width = ggplot2::unit(key_width, "cm"),
                       legend.text = ggplot2::element_text(size = 3),
                       plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
      return(gg)
    })
  }
  if (combine) {
    gg_list <- patchwork::wrap_plots(gg_list, ncol = ncol)
  }
  return(gg_list)
}


#' @title Tailored feature plot
#'
#' @description This function generates the same plot as `feature_plot`,
#' although it focuses on a single feature and generates slightly better
#' looking plot.
#'
#' @param seu Seurat object
#' @param feature Feature to plot.
#' @param min.cutoff Minimum cutoff value for feature, may specify quantile
#' in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10').
#' @param max.cutoff Maximum cutoff value for feature, may specify quantile
#' in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10').
#' @param reduction Dimensionality reduction to use.
#' @param slot Slot to extract data from.
#' @param col_pal Continuous colour palette to use, default "RdYlBu".
#' @param pt.size Adjust point size for plotting.
#' @param stroke Stroke value for each point.
#' @param legend.position Position of legend, default "right" (set to "none"
#' for clean plot).
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
feature_plot_tailored <- function(seu, feature, max.cutoff = "q98", min.cutoff = NA,
    reduction = "umap", slot = "data", col_pal = NULL, pt.size = 2,
    stroke = 0.05, legend.position = "right") {

  # Extract dimensionally reduced data from full object
  dim_dt <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  # Keep only specified dimensions and provide generic column names
  dim_dt <- dim_dt[, c(1,2)]
  # Extract x and y plotting limits
  xlim <- c(min(dim_dt[,1]) - 0.5, max(dim_dt[,1]) + 0.5)
  ylim <- c(min(dim_dt[,2]) - 0.5, max(dim_dt[,2]) + 0.5)
  # Add feature to plot as 3rd column
  dim_dt[, 3] <- Seurat::FetchData(object = seu, vars = feature, slot = slot)

  # Determine cutoffs
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(dim_dt[, 3]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = feature
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(dim_dt[, 3]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = feature
  )

  # Apply cutoffs
  data.feature <- as.vector(x = dim_dt[, 3])
  min.use <- Seurat::SetQuantile(cutoff = min.cutoff, data.feature)
  max.use <- Seurat::SetQuantile(cutoff = max.cutoff, data.feature)
  data.feature[data.feature < min.use] <- min.use
  data.feature[data.feature > max.use] <- max.use
  dim_dt[, 3] <- data.feature

  if (is.null(col_pal)) { col_pal = "RdYlBu" }

  key_height = 0.6
  key_width = 0.2
  if (legend.position == "top") {
    key_height = 0.2
    key_width = 0.7
  }

  gg <- ggplot2::ggplot(dim_dt, ggplot2::aes(x = dim_dt[, 1], y = dim_dt[, 2],
                                             fill = dim_dt[, 3])) +
    ggplot2::geom_point(shape = 21, size = pt.size, stroke = stroke) +
    ggplot2::scale_fill_distiller(palette = col_pal, limits = c(min.use, max.use)) +
    ggplot2::scale_x_continuous(name = NULL, minor_breaks = NULL,
                                limits = xlim) +
    ggplot2::scale_y_continuous(name = NULL, minor_breaks = NULL,
                                limits = ylim) +
    ggplot2::labs(fill = feature, title = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = legend.position,
                   legend.key.size = ggplot2::unit(0.7, "line"),
                   legend.key.height = ggplot2::unit(key_height, "cm"),
                   legend.key.width = ggplot2::unit(key_width, "cm"),
                   legend.text = ggplot2::element_text(size = 3),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  return(gg)
}



#' @title Tailored spatial dim plot
#'
#' @description This function adapts the SpatialDimPlot Seurat function by
#'   providing additional plotting options.
#'
#' @param seu Seurat object (required).
#' @param group.by Name of meta.data column to group the data by.
#' @param title Plot title
#' @param alpha Controls opacity of spots. Provide a single alpha value for
#' each plot.
#' @param pt.size.factor Scale the size of the spots.
#' @param crop Crop the plot in to focus on points plotted. Set to FALSE to
#' show entire background image.
#' @param col_pal Discrete colour palette to use, default is Hue palette
#' (hue_pal) from 'scales' package. Should be of equal length to number of
#' groups in 'group.by'.
#' @param legend.position Position of legend, default "right" (set to "none"
#' for clean plot).
#' @param combine Combine plots into a single patchworked ggplot object.
#'    If FALSE, return a list of ggplot objects
#' @param ... Additional parameters passed to Seurat's SpatialDimPlot.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
spatial_dim_plot <- function(
    seu, group.by = "active.ident", title = NULL, alpha = 0.6,
    pt.size.factor = 1.4, crop = TRUE, col_pal = NULL,
    legend.position = "top", combine = TRUE, ...) {

  assertthat::assert_that(!is.null(group.by))

  if (!(is.character(seu@meta.data[[group.by]]) ||
        is.factor(seu@meta.data[[group.by]]) ||
        is.logical(seu@meta.data[[group.by]]) ) ) {
    stop("Error: 'group.by' should not be continuous in metadata slot")
  }

  group <- as.factor(seu@meta.data[[group.by]])
  if (is.null(col_pal)) {
    col_pal <- scales::hue_pal()(nlevels(group))
    names(col_pal) <- levels(group)
  } else {
    col_pal <- col_pal[1:nlevels(group)]
    names(col_pal) <- levels(group)
  }

  t <- Seurat::SpatialDimPlot(
    seu, group.by = group.by, crop = crop, combine = combine,
    pt.size.factor = pt.size.factor, alpha = alpha, ...)
  t <- t &
    ggplot2::theme(legend.position = legend.position,
                   legend.key.size = ggplot2::unit(1.3, "line"),
                   legend.text = ggplot2::element_text(size = 11)) &
    ggplot2::scale_color_manual(name = NULL, values = col_pal) &
    ggplot2::scale_fill_manual(name = NULL, values = col_pal) &
    ggplot2::ggtitle(label = title)
    return(t)
}


#' @title Tailored spatial feature plot
#'
#' @description This function adapts the SpatialFeaturePlot Seurat function by
#'   providing additional plotting options.
#'
#' @param features Vector of features to plot.
#' @param alpha Controls opacity of spots. Provide as a vector specifying the
#' min and max range of values (between 0 and 1).
#' @param ncol Number of columns for display when having multiple features.
#' @param max.cutoff Vector of maximum cutoff values for each feature, may
#' specify quantile in the form of 'q##' where '##' is the quantile
#' (eg, 'q1', 'q10').
#' @param min.cutoff Vector of minimum cutoff values for each feature, may
#' specify quantile in the form of 'q##' where '##' is the quantile
#' (eg, 'q1', 'q10').
#' @param col_pal Continuous colour palette to use from viridis package,
#' default "inferno".
#' @param ... Additional parameters passed to Seurat's SpatialFeaturePlot.
#'
#' @inheritParams spatial_dim_plot
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
spatial_feature_plot <- function(
    seu, features = "nFeature_Spatial", title = NULL, alpha = c(0.1, 0.9),
    pt.size.factor = 1.4, ncol = NULL, max.cutoff = "q98", min.cutoff = NA,
    crop = TRUE, col_pal = "inferno", legend.position = "top", combine = TRUE, ...) {

  assertthat::assert_that(!is.null(features))
  # Extract features present in the Seurat object
  features <- features[(features %in% rownames(seu)) |
                         (features %in% colnames(seu@meta.data))]
  if (length(features) == 0) {
    message("No features present to plot.")
    return(ggplot2::ggplot())
  }

  if (is.null(col_pal)) { col_pal = "inferno" }

  key_height = 0.6
  key_width = 0.2
  if (legend.position == "top") {
    key_height = 0.2
    key_width = 0.8
  }

  t <- Seurat::SpatialFeaturePlot(
    seu, features = features, crop = crop, max.cutoff = max.cutoff, min.cutoff = min.cutoff,
    ncol = ncol, combine = combine, pt.size.factor = pt.size.factor, alpha = alpha, ...)
  t <- t &
    ggplot2::theme(legend.position = legend.position,
                   legend.key.size = ggplot2::unit(0.7,"line"),
                   legend.key.height = ggplot2::unit(key_height, 'cm'),
                   legend.key.width = ggplot2::unit(key_width, "cm"),
                   legend.text = ggplot2::element_text(size = 5)) &
    viridis::scale_fill_viridis(option = col_pal) &
    ggplot2::ggtitle(label = title)
  return(t)
}


#' @title Tailored scatter plot of metadata
#'
#' @description This function creates a scatter plot to assess feature
#' to feature relationships
#'
#' @param seu Seurat object (required).
#' @param features Vector of features in metadata to plot.
#'
#' @return A ggplot2 list with all combinations of elements in features
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
scatter_meta_plot <- function(
    seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mito")) {
  x = y <- NULL
  # Get all combinations of QCs to plot
  combs <- combn(features, 2)
  # Create scatter plots
  tmp <- seu@meta.data
  gg_list <- list()
  for (comb in 1:NCOL(combs)) {
    feats <- combs[, comb]
    to_plot <- data.frame(x = tmp[[feats[1]]], y = tmp[[feats[2]]])

    dens <- try(.get_density(x = tmp[[feats[1]]],
                             y = tmp[[feats[2]]], n = 100))
    if (!inherits(dens, "try-error")) {
      to_plot$dens <- dens
      gg_list[[comb]] <- ggplot2::ggplot(to_plot) +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = dens),
                            pointsize = 1.2) +
        viridis::scale_color_viridis() + ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab(feats[1]) + ggplot2::ylab(feats[2])
    }
  }
  return(gg_list)
}


#' @title PCA and feature metadata correlation heatmap plot
#'
#' @description This function computes correlation of principal components
#' with specific feature metadata (e.g. nFeature) to identify potential
#' technical/batch effects in the data.
#'
#' @param seu Seurat object (required).
#' @param features Vector of features in metadata to plot.
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
pca_feature_cor_plot <- function(seu, features) {
  PC = QC = cor <- NULL
  # Heatmap showing correlation of metadata with Principal Components
  if ("pca" %in% names(seu@reductions)) {
    df_corr <- abs(stats::cor(seu@reductions[["pca"]]@cell.embeddings,
                              seu[[features]])) |>
      dplyr::as_tibble(rownames = "PC") |>
      tidyr::pivot_longer(cols = -c(PC), names_to = "QC", values_to = "cor")
    df_corr$PC <- factor(df_corr$PC,
                         levels = paste0("PC_", seq(1, NCOL(seu@reductions[["pca"]]))))
    # Create heatmap
    gg <- ggplot2::ggplot(df_corr, ggplot2::aes(x = PC, y = QC, fill = cor)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_distiller(name = "abs(cor)", type="div", palette="RdYlBu") +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                         hjust = 1, size = 6))
    return(gg)
  } else {
    return(ggplot2::ggplot())
  }
}


#' @title Tailored dot plot
#'
#' @description This function adapts the DotPlot Seurat function by
#'   providing additional plotting options.
#'
#' @param features Vector of features to plot.
#' @param group.by Name of metadata column to group cells by (required).
#' @param labels If we want to have different labels plotted for each feature.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param col_pal Continuous colour palette to use, default "RdYlBu".
#' @param ... Additional parameters passed to Seurat's DotPlot.
#' @inheritParams cluster_analysis
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
dot_plot <-  function(seu, features, group.by = NULL, labels = NULL,
                      xlab = "Signature", ylab = "Cluster",
                      legend.position = "right", col_pal = NULL, ...) {
  # Keep only features that are present in the metadata
  idx <- features %in% colnames(seu@meta.data)
  if (sum(idx) == 0) { return(ggplot2::ggplot()) }
  features <- features[idx]

  # If we want to have different labels plotted
  if (!is.null(labels)) {
    labels <- labels[idx]
    if (length(labels) != length(features)) { labels <- features }
  } else {
    labels <- features
  }
  if (is.null(col_pal)) { col_pal = "RdYlBu" }

  p <- Seurat::DotPlot(seu, features = features, group.by = group.by,  ...) +
    ggplot2::coord_flip() +
    ggplot2::scale_colour_distiller(name = NULL, type = "div",
                                    palette = col_pal, limits = c(-2.5,2.5)) +
    ggplot2::scale_y_discrete(name = ylab) +
    ggplot2::scale_x_discrete(name = xlab, labels = labels) +
    ggplot2::theme(legend.position = legend.position)
  return(p)
}

#' @title Tailored heatmap plot
#'
#' @description Heatmap plot showing marker genes after the clustering process.
#' Expression levels are taken from the `scaled.data` slot. Current
#' implementation is 'hard' coded: that is for heatmap column annotation,
#' it assumes the Seurat object contains 'condition', 'sample' and
#' 'seurat_clusters' as columns in meta.data slot.
#'
#' @param seu Seurat object (required).
#' @param markers Data frame with marker genes for each cluster.
#' Expects a format as the output of Seurat's FindAllMarkers.
#' @param topn_genes Top N marker genes to plot for each cluster.
#' @param filename Filename for saving the heatmap plot. If null, the heatmap
#' is just plotted in device.
#' @param col_pal Discrete colour palette to use, default is Hue palette
#' (hue_pal) from 'scales' package.
#' @param ... Additional parameters passed to 'pheatmap' function
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
heatmap_plot <- function(seu, markers, topn_genes = 10,
                         filename = NULL, col_pal = NULL, ...) {
  # TODO: Make the function more generic so the user defines what
  # column annotation is show in the heatmap.

  # So CMD passes without NOTES
  p_val_adj = cluster = avg_log2FC <- NULL

  if (NROW(markers) == 0) { return(-1) }

  # Colours for each group (by default cluster, sample, condition)
  if (is.null(col_pal)) {
    colours_cluster <- scales::hue_pal()(nlevels(seu$seurat_clusters))
    names(colours_cluster) <- levels(seu$seurat_clusters)
    colours_sample <- scales::hue_pal()(nlevels(seu$sample))
    names(colours_sample) <- levels(seu$sample)
    colours_condition <- scales::hue_pal()(nlevels(seu$condition))
    names(colours_condition) <- levels(seu$condition)
  } else {
    colours_cluster <- col_pal[1:nlevels(seu$seurat_clusters)]
    names(colours_cluster) <- levels(seu$seurat_clusters)
    colours_sample <- col_pal[1:nlevels(seu$sample)]
    names(colours_sample) <- levels(seu$sample)
    colours_condition <- col_pal[1:nlevels(seu$condition)]
    names(colours_condition) <- levels(seu$condition)
  }

  annotation_colours <- list(Cluster = colours_cluster,
                             Sample = colours_sample,
                             Condition = colours_condition)

  # Extract scale data
  scale_data <- Seurat::GetAssayData(object = seu, slot = "scale.data")
  # Make cluster a factor
  markers$cluster <- factor(markers$cluster)
  #markers <- markers[!(grepl("^Rpl|^Rps|^mt-", markers$gene)),]
  #markers <- markers[markers$pct.2 < 0.2,]
  markers <- markers[markers$gene %in% rownames(scale_data), ]
  markers <- markers |> dplyr::filter(p_val_adj < 0.05) |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(n = topn_genes, order_by = avg_log2FC)
  markers <- markers[!duplicated(markers$gene), ]

  # Shuffle cells (columns) so we don't get sample
  # specific block signal in heatmap
  cell_order <- unlist(lapply(split(Seurat::Cells(seu), seu$seurat_clusters),
                              FUN = sample))
  annotation_col <- data.frame(Cluster = seu$seurat_clusters[cell_order],
                               Sample = seu$sample[cell_order],
                               Condition = seu$condition[cell_order])
  rownames(annotation_col) <- cell_order
  # Annotate genes (rows)
  annotation_row <- data.frame(Cluster = factor(markers$cluster))
  rownames(annotation_row) <- markers$gene

  # Subset matrix for heatmap plot
  seu_filt <- scale_data[as.character(markers$gene[order(annotation_row[,1])]),
                         rownames(annotation_col)[order(annotation_col[,1])]]
  seu_filt[seu_filt > 2.5] <- 2.5
  seu_filt[seu_filt < -2.5] <- -2.5
  gaps_col <- cumsum(summary(seu$seurat_clusters))
  gaps_row <- cumsum(summary(markers$cluster))

  if (length(dev.list()) != 0) { dev.off() }

  if (!is.null(filename)) {
    png(filename, width = 18, height = nrow(markers)*3/13,
        res = 200, units = "in")
    print(pheatmap::pheatmap(seu_filt,
                             color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11,name = "RdBu")))(1000),
                             breaks = seq(-2.5, 2.5, 0.005), gaps_row = gaps_row, gaps_col = gaps_col,
                             tree_height_row = NA, tree_height_col = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                             annotation_col = annotation_col, annotation_row = annotation_row,
                             annotation_colors = annotation_colours,
                             annotation_legend = TRUE, annotation_names_row = FALSE, annotation_names_col = TRUE,
                             show_rownames = TRUE, show_colnames = FALSE, legend = TRUE, ...))
    dev.off()
  } else {
    print(pheatmap::pheatmap(seu_filt,
                             color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11,name = "RdBu")))(1000),
                             breaks = seq(-2.5, 2.5, 0.005), gaps_row = gaps_row, gaps_col = gaps_col,
                             tree_height_row = NA, tree_height_col = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                             annotation_col = annotation_col, annotation_row = annotation_row,
                             annotation_colors = annotation_colours,
                             annotation_legend = TRUE, annotation_names_row = FALSE, annotation_names_col = TRUE,
                             show_rownames = TRUE, show_colnames = FALSE, legend = TRUE, ...))
  }
  if (length(dev.list()) != 0) { dev.off() }
  return(0)
}
