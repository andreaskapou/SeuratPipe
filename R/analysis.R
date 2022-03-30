#' @rdname dimred_qc_plots
#'
#' @title QC and general metadata plots visualised on dimensiona
#'
#' @description This function generates QC and general metadata plots that
#' are visualised in dimensional reduced space (e.g. PCA and UMAP). The aim is
#' to visualise weather there are any technical factors driving the variability
#' in the data, e.g. integration driven by sample, correlation of specific
#' PCs with number of genes expressed in a cell, etc.
#'
#' @param reductions Vector with dimensional reductions objects to use for
#' plotting the metadata and QC features.
#' @param metadata_to_plot Vector with metadata names to plot, they should be
#' present in the  meta.data slot.
#' @param qc_to_plot Vector with QC names to plot, they should be
#' present in the  meta.data slot.
#' @param dims_plot Dimensions to plot. For UMAP and TSNE this is set to c(1,2).
#' @param ... Additional named parameters passed to Seurat's DimPlot
#' and FeaturePlot.
#' @inheritParams cluster_analysis
#'
#' @return Nothing, plots are saved in the specified folder if not NULL.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
dimred_qc_plots <- function(seu, reductions = c("umap", "pca"),
                            metadata_to_plot = NULL, qc_to_plot = NULL,
                            plot_dir, max.cutoff = "q98", legend.position = "top",
                            col_pal = NULL, dims_plot = c(1,2), pt.size = 1.4,
                            fig.res = 200, ...) {
  # TODO: Compute variance of first PCs and show them in the plot

  # Plots based on metadata columns
  for (meta in metadata_to_plot) {
    for (red in reductions) {
      if (red != c("umap", "tsne")) { dims_plot <- dims_plot }
      else { dims_plot <- c(1,2) }
      png(paste0(plot_dir,"01_", red, "_", meta, ".png"), width = 7,
          height = 7, res = fig.res, units = "in")
      plot(dim_plot(seu = seu, reduction = red, group.by = meta,
                    split.by = NULL, ncol = 1, legend.position = legend.position,
                    col_pal = col_pal, dims_plot = dims_plot, pt.size = pt.size,
                    label = FALSE, combine = TRUE, ...) & Seurat::NoAxes())
      dev.off()

      # Split by metadata plots
      # Plot dimensions
      plot_dim <- .plot_dims(feat_len = length(unique(seu@meta.data[[meta]])))
      png(paste0(plot_dir, "01_", red, "_", meta, "_split.png"),
          width = plot_dim$width, height = plot_dim$height, res = fig.res,
          units = "in")
      plot(dim_plot(seu = seu, reduction = red, group.by = meta,
                    split.by = meta, ncol = plot_dim$ncols,
                    legend.position = legend.position,
                    col_pal = col_pal, dims_plot = dims_plot, pt.size = pt.size,
                    label = FALSE, combine = TRUE, ...) & Seurat::NoAxes())
      dev.off()
    }
  }

  # Major QCs plots
  if (!is.null(qc_to_plot)) {
    for (red in reductions) {
      if (red != c("umap", "tsne")) { dims_plot <- dims_plot }
      else { dims_plot <- c(1,2) }

      plot_dim <- .plot_dims(feat_len = length(qc_to_plot))
      png(paste0(plot_dir, "02_qc_", red, ".png"), width = plot_dim$width,
          height = plot_dim$height, res = fig.res, units = "in")
      plot(feature_plot(seu = seu, reduction = red, features = qc_to_plot,
                        max.cutoff = max.cutoff, ncol = plot_dim$ncols,
                        legend.position = legend.position,
                        col_pal = col_pal, dims_plot = dims_plot,
                        pt.size = pt.size, ...) &
             Seurat::NoAxes())
      dev.off()
    }
  }
}


#' @rdname harmony_analysis
#'
#' @title Analysis steps for Harmony integration
#'
#' @description This function implements all the analysis steps for performing
#' Harmony integration on a Seurat object. These include, 1. merge all samples
#' in a single Seurat object (if a list of Seurat objects is provided) 2.
#' data normalisation, 2. identification of HVGs, 3. scaling of expression
#' values, 4. Computing PCA, 5. running Harmony integration, 6. Generating
#' plots.
#'
#' @param npcs Number of Principal Components (PCs) to compute.
#' @param dims.use Vector with PCs to use for Harmony integration, e.g. 1:50.
#' @param n_hvgs Number of highly variable genes (HVGs) to compute, which will
#' be used as inpute to PCA.
#' @param max.iter.harmony Maximum number of iterations for Harmony integration.
#' @param ... Additional named parameters passed to RunHarmony and other Seurat
#' processing functions, such as RunPCA and ScaleData.
#' @inheritParams cluster_analysis
#'
#' @return Updated Seurat object with integrated and processed data.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
harmony_analysis <- function(seu, npcs, dims.use, plot_dir = NULL, n_hvgs = 3000,
                             max.iter.harmony = 50, assay = "RNA", seed = 1,
                             fig.res = 200, ...) {
  # If list, then we have un-merged independent samples
  if (is.list(seu)) {
    # Normalise data
    seu <- lapply(X = seu, FUN = function(x) {
      x <- Seurat::NormalizeData(x, normalization.method = "LogNormalize",
                                 assay = assay,
                                 scale.factor = 10000)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst",
                                        assay = assay,
                                        nfeatures = n_hvgs)
      # Plot HVGs
      png(paste0(plot_dir,"hvgs_", x$sample[1], ".png"), width = 10,
          height = 5, res = fig.res, units = "in")
      print(Seurat::LabelPoints(plot = Seurat::VariableFeaturePlot(x),
                                points = head(Seurat::VariableFeatures(x), 10),
                                repel = TRUE))
      dev.off()
      return(x)
    })
    seu <- merge(x = seu[[1]], y = seu[2:length(seu)], project = "Liver")
    # Seurat merge bug which converts factors to character.
    seu <- .as_factor_metadata_seurat(seu)
  }
  seu <- seu %>% Seurat::NormalizeData(normalization.method = "LogNormalize",
                                       assay = assay,
                                       scale.factor = 10000) %>%
    Seurat::FindVariableFeatures(selection.method = "vst",
                                 assay = assay,
                                 nfeatures = n_hvgs) %>%
    Seurat::ScaleData(assay = assay, ...) %>%
    # PCA reduction -> Harmony integration
    Seurat::RunPCA(features = Seurat::VariableFeatures(.),
                   npcs = npcs, assay = assay, ...)

  # Run Harmony
  set.seed(seed) # Set seed due to Harmony being stochastic
  seu <- run_harmony(seu, group.by.vars = c("sample"), dims.use = dims.use,
                     max.iter.harmony = max.iter.harmony, assay.use = assay, ...)

  # Plots
  png(paste0(plot_dir, "pca_heatmap.png"), width = 15, height = 15,
      res = fig.res, units = "in")
  Seurat::DimHeatmap(seu, dims = 1:9, nfeatures = 30, cells = 300,
                     reduction = "pca", balanced = TRUE)
  dev.off()
  png(paste0(plot_dir, "harmony_heatmap.png"), width = 15, height = 15,
      res = fig.res, units = "in")
  Seurat::DimHeatmap(seu, dims = 1:9, nfeatures = 30, cells = 300,
                     reduction = "harmony", balanced = TRUE)
  dev.off()
  png(paste0(plot_dir, "pca_elbow.png"), width = 15, height = 15,
      res = fig.res, units = "in")
  print(Seurat::ElbowPlot(seu, ndims = npcs, reduction = "pca"))
  dev.off()
  png(paste0(plot_dir, "harmony_elbow.png"), width = 15, height = 15,
      res = fig.res, units = "in")
  print(Seurat::ElbowPlot(seu, ndims = npcs, reduction = "harmony"))
  dev.off()

  # Return processed and integrated object
  return(seu)
}


#' @rdname module_score_analysis
#'
#' @title Module score analysis
#'
#' @description This function takes a group of modules and calculates module
#' scores for the expression programs. In addition, if `plot_dir` is not NULL
#' creates feature plots with marker genes and the module score and saved them
#' in `plot_dir` directory.
#'
#' @param reduction Which dimensionality reduction to use (required).
#' @param dims_plot Dimensions to plot, must be a two-length numeric vector
#' specifying x- and y-dimensions.
#' @param ... Additional named parameters passed to Seurat's AddModuleScore
#' and FeaturePlot
#' @inheritParams cluster_analysis
#'
#' @return Updated Seurat object with calculated module scores in metadata.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
module_score_analysis <- function(seu, modules_group, plot_dir = NULL,
                                  reduction = "umap", max.cutoff = "q98",
                                  ctrl = 100, seed = 1, fig.res = 200,
                                  legend.position = "top", col_pal = NULL,
                                  dims_plot = c(1, 2), ...) {
  # Required for the really low quality samples so they pass without errors
  ctrl <- ifelse(NCOL(seu) <= ctrl, NCOL(seu) - 10, ctrl)
  # Iterate over the group of modules
  for (mg in names(modules_group)) {
    # Extract list of modules within the group
    modules <- modules_group[[mg]]
    # For each module compute module scores and plot
    for (m in names(modules)) {
      # Extract genes that are present in data
      features <- modules[[m]][modules[[m]] %in% rownames(seu)]
      if (length(features) == 0) { next }
      seu <- compute_module_score(seu = seu, features = features, name = m,
                                  ctrl = ctrl, seed = seed, ...)

      if (!is.null(plot_dir)) {
        plot_dim <- .plot_dims(feat_len = length(features))
        png(paste0(plot_dir, "01_markers_", m, ".png"), width = plot_dim$width,
            height = plot_dim$height, res = fig.res, units = "in")
        plot(feature_plot(seu = seu, reduction = reduction, features = modules[[m]],
                          max.cutoff = max.cutoff, ncol = plot_dim$ncols,
                          legend.position = legend.position, col_pal = col_pal,
                          dims_plot = dims_plot, combine = TRUE, ...) & Seurat::NoAxes())
        dev.off()
      }
    }
    if (!is.null(plot_dir)) {
      # Module scores in one plot
      plot_dim <- .plot_dims(feat_len = length(modules))
      png(paste0(plot_dir, "02_score_", mg, ".png"), width = plot_dim$width,
          height = plot_dim$height, res = fig.res, units = "in")
      plot(feature_plot(seu = seu, reduction = reduction, features = names(modules),
                        max.cutoff = max.cutoff, ncol = plot_dim$ncols,
                        legend.position = legend.position, col_pal = col_pal,
                        dims_plot = dims_plot, combine = TRUE, ...) & Seurat::NoAxes())
      dev.off()
    }
  }
  return(seu)
}


#' @rdname cluster_analysis
#'
#' @title Common clustering analysis steps
#'
#' @description This function implements all the analysis steps for performing
#' clustering on a Seurat object. These include, 1. finding neighbours in lower
#' dimensional space (defined in 'cluster_reduction' parameter) 2. obtaining
#' clusters, 3. identifying marker genes (NOTE: to speed up re-analysis it
#' first checks if file with marker genes is already present, if yes reads the
#' file instead of calling FinaAllMarkers) and 4. generating plots, which
#' include heatmap with (scaled) expression of marker genes in each cluster,
#' marker gene expression on feature plots (e.g. UMAP space, defined in
#' plot_reduction'  parameter), dot / feature plots with
#' pre-computed module scores on each cluster (assumes we have first run
#' 'module_score_analysis' function). This step could be useful for
#' lineage annotation.
#'
#' @param seu Seurat object (required).
#' @param ndim Number of dimensions to use for clustering (from
#' cluster_reduction' object).
#' @param res Vector with clustering resolutions (e.g. seq(0.1, 0.6, by = 0.1)).
#' @param logfc.threshol Limit testing to genes which show, on average, at
#' least X-fold difference (log-scale) between the two groups of cells.
#' @param min.pct Only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations.
#' @param only.pos Only return positive markers (TRUE by default).
#' @param topn Top cluster marker genes to use for plot (in heatmap and
#' feature plots), default is 10.
#' @param plot_dir Directory to save generated plots. If NULL, plots are
#' not saved.
#' @param plot_cluster_markers Logical, wheather to create feature plots with
#' 'topn' cluster markers. Added mostly to reduce number of files (and size)
#' in analysis folders. Default is TRUE.
#' @param modules_group Group of modules (named list of lists) storing features
#' (e.g. genes) to compute module score for each identified cluster. This step
#' can be useful for annotating the different clusters by saving dot
#' plots for each group. Assumes that we already have computed the modules e.g.
#' by calling the 'module_score_analysis' function. If 'plot_dir' is NULL,
#' no plots will be generated.
#' @param cluster_reduction Dimensionality reduction to use for performing
#' clustering. Default is 'harmony', should be set to 'pca' if we do not
#' perform data integration.
#' @param plot_reduction Dimensionality reduction to use for plotting
#' functions. Default is 'umap'.
#' @param max.cutoff Vector of maximum cutoff values for each feature,
#' may specify quantile in the form of 'q##' where '##' is the quantile
#' (eg, 'q1', 'q10').
#' @param legend.position Position of legend, default "top" (set to "none"
#' for clean plot).
#' @param cl.legend.position Position of legend for clusters, default "right"
#' (set to "none" for clean plot)
#' @param col_pal Colour palette to use.
#' @param pt.size Adjust point size for plotting.
#' @param label Whether to label the clusters.
#' @param label.size Sets size of labels.
#' @param assay Assay to perform analysis, default "RNA".
#' @param ctrl Number of control features selected from the same bin per
#' analyzed feature.
#' @param seed Set a random seed.
#' @param fig.res Figure resolution in ppi (see 'png' function).
#' @param ... Additional named parameters passed to Seurat
#' analysis and plotting functions, such as FindClusters, FindAllMarkers,
#' DimPlot and FeaturePlot.
#'
#' @return Updated Seurat object clustered cells
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
cluster_analysis <- function(seu, ndim = 20, res = seq(0.1, 0.1, by = 0.1),
                             logfc.threshol = 0.5, min.pct = 0.25,
                             only.pos = TRUE, topn = 10, plot_dir = NULL,
                             plot_cluster_markers = TRUE, modules_group = NULL,
                             cluster_reduction = "harmony",
                             plot_reduction = "umap", max.cutoff = "q98",
                             cl.legend.position = "right",
                             legend.position = "top", col_pal = NULL,
                             pt.size = 1.4, label = TRUE, label.size = 8,
                             assay = "RNA", ctrl = 100, seed = 1,
                             fig.res = 200, ...) {
  # So CMD passes without NOTES
  cluster = avg_log2FC <- NULL
  # Required for the really low quality samples so they pass without errors
  ctrl <- ifelse(NCOL(seu) <= ctrl, NCOL(seu) - 10, ctrl)

  # Iterate over each resolution
  for (r in res) {
    cat("Res", r, "\n")
    # Identify clusters
    seu <- Seurat::FindNeighbors(seu, dims = 1:ndim, assay = assay,
                                 reduction = cluster_reduction, ...)
    seu <- Seurat::FindClusters(seu, resolution = r, random.seed = seed, ...)

    if (!is.null(plot_dir)) {
      png(paste0(plot_dir,"z_", plot_reduction, "_res", r, ".png"), width = 7,
          height = 7, res = fig.res, units = "in")
      plot(dim_plot(seu = seu, reduction = plot_reduction, split.by = NULL,
                    group.by = "seurat_clusters", ncol = 1, col_pal = col_pal,
                    legend.position = cl.legend.position, dims_plot = c(1,2),
                    pt.size = pt.size, label = label, label.size = label.size,
                    combine = TRUE, ...) & Seurat::NoAxes())
      dev.off()
    }

    # Create dot plot for each module group, assumes we already have computed
    # the modules e.g. by calling the 'module_score_analysis' function.
    if (!is.null(plot_dir) & !is.null(modules_group)) {
      # Iterate over the group of modules
      for (mg in names(modules_group)) {
        # Extract module names to use for plotting
        m_names <- names(modules_group[[mg]])
        # Check that they indeed are computed in Seurat object
        if (sum(m_names %in% colnames(seu@meta.data)) > 0) {
          png(paste0(plot_dir,"dotplot_module_", mg, "_res", r, ".png"),
              width = 6 + 0.3*nlevels(seu@meta.data[["seurat_clusters"]]),
              height = 3 + 0.15*length(m_names), res = fig.res, units = "in")
          plot(dot_plot(seu = seu, features = m_names, labels = NULL,
                        group.by = "seurat_clusters", xlab = "Signature",
                        ylab = "Cluster", legend.position = "right",
                        col_pal = col_pal, ...))
          dev.off()
        }
      }
    }

    # Identify cluster markers
    if (!file.exists(paste0(plot_dir, "seu_markers_res", r, ".csv"))) {
      mark <- Seurat::FindAllMarkers(seu, only.pos = only.pos,
                                     min.pct = min.pct,
                                     logfc.threshold = logfc.threshol, ...)
      write.csv(mark, file = paste0(plot_dir, "seu_markers_res", r, ".csv"))
      write.csv(seu@meta.data,
                file = paste0(plot_dir, "seu_meta_res", r, ".csv"))
    } else {
      mark <- read.csv(file = paste0(plot_dir, "seu_markers_res", r, ".csv"))
    }

    # Heatmap of marker genes
    heatmap_plot(seu = seu, markers = mark, topn = topn, assay = assay,
                 filename = paste0(plot_dir, "z_heatmap_res", r, ".png"), ...)

    ## Feature and violin plots
    if (plot_cluster_markers & !is.null(plot_dir)) {
      # Extract top marker genes
      top_mark <- mark[mark$gene %in% rownames(seu[[assay]]@data), ] %>%
        dplyr::group_by(cluster) %>% dplyr::top_n(topn, avg_log2FC)
      # For each cluster plot marker genes
      for (cl in levels(top_mark$cluster)) {
        if (nrow(top_mark[top_mark$cluster == cl, ]) > 0) {
          # get genes for the specific cluster
          genes <- top_mark$gene[top_mark$cluster == cl]
          plot_dim <- .plot_dims(feat_len = length(genes))
          png(paste0(plot_dir, "01_feature_seu_res", r, "_cl", cl, ".png"),
              width = plot_dim$width, height = plot_dim$height, res = 150,
              units = "in")
          plot(feature_plot(seu = seu, reduction = plot_reduction,
                            features = genes, max.cutoff = max.cutoff,
                            ncol = plot_dim$ncols, col_pal = col_pal,
                            legend.position = legend.position, dims_plot = c(1,2),
                            pt.size = pt.size, ...) & Seurat::NoAxes())
          dev.off()

          png(paste0(plot_dir, "02_violin_seu_res", r, "_cl", cl,".png"),
              width = plot_dim$width, height = plot_dim$height, res = 100,
              units = "in")
          plot(Seurat::VlnPlot(seu, features = genes, ncol = plot_dim$ncols,
                               pt.size = 0.03, ...))
          dev.off()

          # Each cluster as a module and compute score
          seu <- compute_module_score(seu = seu, features = genes, ctrl = ctrl,
                                      name = paste0("Cluster", cl), ...)
        }
      }
      # Module scores in one plot
      plot_dim <- .plot_dims(feat_len = length(unique(top_mark$cluster)))
      png(paste0(plot_dir, "03_score_seu_res", r, ".png"),
          width = plot_dim$width, height = plot_dim$height, res = fig.res,
          units = "in")
      plot(feature_plot(seu = seu, reduction = plot_reduction,
                        features = paste0("Cluster", unique(top_mark$cluster)),
                        max.cutoff = max.cutoff,
                        ncol = plot_dim$ncols, col_pal = col_pal,
                        legend.position = legend.position, dims_plot = c(1,2),
                        pt.size = pt.size, ...) & Seurat::NoAxes())
      dev.off()
    }
  }
  return(seu)
}
