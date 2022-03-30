#' @title Tailored dimensional reduction plot
#'
#' @description This function extends the DimPlot Seurat function by
#'   providing additional plotting options.
#'
#' @param reduction Which dimensionality reduction to use (required).
#' @param group.by Name of metadata column to group (color) cells by (required).
#' @param split.by Name of a metadata column to split plot by
#' @param ncol Number of columns for display when combining plots
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
                     split.by = NULL, ncol = 1, legend.position = "top",
                     col_pal = NULL, dims_plot = c(1, 2), pt.size = 1.4,
                     label = TRUE, label.size = 7, combine = TRUE, ...) {
  # So CMD passes
  D1 = D2 = ident = x = y <- NULL
  assertthat::assert_that(!is.null(reduction))
  assertthat::assert_that(!is.null(group.by))

  if (reduction != c("umap", "tsne")) { dims_plot <- dims_plot }
  else { dims_plot <- c(1,2) }

  # Extract dimensionally reduced data
  dim_dt <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  # Keep only specified dimensions and provide generic column names
  dim_dt <- dim_dt[, dims_plot]
  colnames(dim_dt) <- c("D1", "D2")
  # Extract x and y plotting limits
  xlim <- c(min(dim_dt@cell.embeddings[,1]) - 0.5,
            max(dim_dt@cell.embeddings[,1]) + 0.5)
  ylim <- c(min(dim_dt@cell.embeddings[,2]) - 0.5,
            max(dim_dt@cell.embeddings[,2]) + 0.5)

  if (!(is.character(seu@meta.data[[group.by]]) ||
        is.factor(seu@meta.data[[group.by]]) ) ) {
    stop("Error: 'group.by' should be a factor in Seurat metadata slot")
  }
  if (!is.null(split.by)) {
    label <- FALSE
    if (!(is.character(seu@meta.data[[split.by]]) ||
          is.factor(seu@meta.data[[split.by]]) ) ) {
      stop("Error: group.by should be a factor in Seurat metadata slot")
    }
  } else {
    ncol <- 1
  }
  group <- as.factor(seu@meta.data[[group.by]])
  if (is.null(col_pal)) {
    col_pal <- scales::hue_pal()(nlevels(group))
    names(col_pal) <- levels(group)
  }
  t <- Seurat::DimPlot(seu, dims = dims_plot, group.by = group.by,
        split.by = split.by, ncol = ncol, pt.size = pt.size, combine = combine,
        label.size = label.size, label = FALSE, ...)
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
    ggplot2::guides(size = FALSE, colour = FALSE,
            fill = ggplot2::guide_legend(override.aes = list(size = 2))) &
    ggplot2::scale_fill_manual(name = NULL, values = col_pal) &
    ggplot2::scale_x_continuous(name = NULL, minor_breaks = NULL,
                                limits = xlim) &
    ggplot2::scale_y_continuous(name = NULL, minor_breaks = NULL,
                                limits = ylim)

  if (label) {
    # Compute cluster centres
    dim_dt$ident <- group
    centres <- dim_dt %>% dplyr::group_by(ident) %>%
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
#' @param ncol Number of columns for display when combining plots
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
feature_plot <- function(seu, reduction = "umap", features = "nFeature_RNA",
                         max.cutoff = "q98", ncol = 1, legend.position = "top",
                         col_pal = NULL, dims_plot = c(1, 2), pt.size = 1.4,
                         combine = TRUE, ...) {
  # Extract features present in the Seurat object
  features <- features[(features %in% rownames(seu)) |
                         (features %in% colnames(seu@meta.data))]
  if (length(features) == 0) {
    message("No features present to plot.")
    return(ggplot2::ggplot())
  }
  assertthat::assert_that(!is.null(reduction))
  assertthat::assert_that(!is.null(features))

  if (reduction != c("umap", "tsne")) { dims_plot <- dims_plot }
  else { dims_plot <- c(1,2) }

  # Extract dimensionally reduced data
  dim_dt <- as.data.frame(seu@reductions[[reduction]]@cell.embeddings)
  # Keep only specified dimensions and provide generic column names
  dim_dt <- dim_dt[, dims_plot]
  # Extract x and y plotting limits
  xlim <- c(min(dim_dt@cell.embeddings[,1]) - 0.5,
            max(dim_dt@cell.embeddings[,1]) + 0.5)
  ylim <- c(min(dim_dt@cell.embeddings[,2]) - 0.5,
            max(dim_dt@cell.embeddings[,2]) + 0.5)

  if (is.null(col_pal)) { col_pal = "RdYlBu" }

  t <- Seurat::FeaturePlot(seu, features = features, reduction = reduction,
                           max.cutoff = max.cutoff, ncol = ncol,
                           combine = combine, pt.size = pt.size, ...)
  t <- t &
    ggplot2::theme(legend.position = legend.position,
          legend.key.size = ggplot2::unit(0.3,"line"),
          legend.key.width = ggplot2::unit(0.6, "cm"),
          legend.text = ggplot2::element_text(size = 5),
          axis.line = ggplot2::element_line(colour = "black", size = 1.25,
                                   linetype = "solid"),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(size = 16),
          axis.title.y = ggplot2::element_text(size = 16),
          axis.ticks.length = ggplot2::unit(0.2,"cm")) &
    ggplot2::geom_point(ggplot2::aes(t$data[, 1], t$data[, 2]), shape = 21,
               size = pt.size, stroke = 0) &
    ggplot2::scale_colour_distiller(palette = col_pal) &
    ggplot2::scale_x_continuous(name = NULL, minor_breaks = NULL,
                                limits = xlim) &
    ggplot2::scale_y_continuous(name = NULL, minor_breaks = NULL,
                                limits = ylim)
  return(t)
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
                                    palette = col_pal, limits = c(-2.5, 2.5)) +
    ggplot2::scale_y_discrete(name = ylab) +
    ggplot2::scale_x_discrete(name = xlab, labels = labels) +
    ggplot2::theme(legend.position = legend.position)
  return(p)
}

#' @title Tailored heatmap plot
#'
#' @description Hetmap plot showing marker genes after the clustering process.
#' Expression levels are taken from the `scaled.data` slot.
#'
#' @param seu Seurat object (required).
#' @param markers Data frame with marker genes for each cluster.
#' Expects a format as the output of Seurat's FindAllMarkers.
#' @param topn Top N marker genes to plot for each cluster.
#' @param assay Assay to extract gene expression data from (default 'RNA').
#' @param filename Filename for saving the heatmap plot. If null, the
#' @param ... Additional parameters passed to 'pheatmap' function
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
heatmap_plot <- function(seu, markers, topn = 10, assay = "RNA",
                         filename = NULL, ...){
  # So CMD passes without NOTES
  p_val_adj = cluster = avg_log2FC <- NULL
  # TODO: Make the function more generic so the user defines what
  # column annotation is show in the heatmap.
  if (NROW(markers) == 0) { return(-1) }

  # Colours for each group (by default cluster, sample, condition)
  colours_cluster <- scales::hue_pal()(nlevels(Seurat::Idents(seu)))
  names(colours_cluster) <- levels(Seurat::Idents(seu))
  colours_sample <- scales::hue_pal()(nlevels(seu$sample))
  names(colours_sample) <- levels(seu$sample)
  colours_condition <- scales::hue_pal()(nlevels(seu$condition))
  names(colours_condition) <- levels(seu$condition)

  annotation_colours <- list(Cluster = colours_cluster,
                             Sample = colours_sample,
                             Condition = colours_condition)

  # Make cluster a factor
  markers$cluster <- factor(markers$cluster)
  #markers <- markers[!(grepl("^Rpl|^Rps|^mt-", markers$gene)),]
  #markers <- markers[markers$pct.2 < 0.2,]
  markers <- markers[markers$gene %in% rownames(seu[[assay]]@scale.data),]
  markers <- markers %>% dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>% dplyr::top_n(topn, avg_log2FC)
  markers <- markers[!duplicated(markers$gene), ]

  # Shuffle cells (columns) so we don't get sample
  # specific block signal in heatmap
  cell_order <- unlist(lapply(split(Seurat::Cells(seu), Seurat::Idents(seu)),
                              FUN = sample))
  annotation_col <- data.frame(Cluster = Seurat::Idents(seu)[cell_order],
                               Sample = seu$sample[cell_order],
                               Condition = seu$condition[cell_order])
  rownames(annotation_col) <- cell_order
  # Annotate genes (rows)
  annotation_row <- data.frame(Cluster = factor(markers$cluster))
  rownames(annotation_row) <- markers$gene

  # Subset matrix for heatmap plot
  seu_filt <- seu[[assay]]@scale.data[
    as.character(markers$gene[order(annotation_row[,1])]),
    rownames(annotation_col)[order(annotation_col[,1])]]
  seu_filt[seu_filt > 2.5] <- 2.5
  seu_filt[seu_filt < -2.5] <- -2.5
  gaps_col <- cumsum(summary(Seurat::Idents(seu)))
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
}
