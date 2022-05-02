#' @rdname dimred_qc_plots
#'
#' @title QC and general metadata plots visualised on dimensional reduced space
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
dimred_qc_plots <- function(seu, reductions = c("pca"),
                            metadata_to_plot = NULL, qc_to_plot = NULL,
                            plot_dir, max.cutoff = "q98",
                            legend.position = "right",
                            cont_col_pal = NULL, discrete_col_pal = NULL,
                            dims_plot = c(1,2), pt.size = 1.4,
                            fig.res = 200, ...) {
  PC = QC = cor = NULL

  # Plots based on metadata columns
  for (meta in metadata_to_plot) {
    for (red in reductions) {
      if (red %in% c("umap", "tsne")) {
        dims_plot <- c(1,2)
      } else {
        dims_plot <- dims_plot
      }
      plot_dim <- .plot_dims(feat_len = 1)
      png(paste0(plot_dir,"01_", red, "_", meta, ".png"), width = plot_dim$width,
          height = plot_dim$height, res = fig.res, units = "in")
      plot(dim_plot(seu = seu, reduction = red, group.by = meta,
                    split.by = NULL, ncol = plot_dim$ncols,
                    legend.position = legend.position,
                    col_pal = discrete_col_pal, dims_plot = dims_plot,
                    pt.size = pt.size, label = FALSE,
                    combine = TRUE, ...) & Seurat::NoAxes())
      dev.off()

      # Plot dimensions
      plot_dim <- .plot_dims(feat_len = length(unique(seu@meta.data[[meta]])))
      # Split by metadata plots
      png(paste0(plot_dir, "01_", red, "_", meta, "_split.png"),
          width = plot_dim$width, height = plot_dim$height, res = fig.res,
          units = "in")
      plot(dim_plot(seu = seu, reduction = red, group.by = meta,
                    split.by = meta, ncol = plot_dim$ncols,
                    legend.position = legend.position,
                    col_pal = discrete_col_pal, dims_plot = dims_plot,
                    pt.size = pt.size, label = FALSE,
                    combine = TRUE, ...) & Seurat::NoAxes())
      dev.off()
    }
  }

  # Major QCs plots
  if (!is.null(qc_to_plot)) {
    for (red in reductions) {
      if (red %in% c("umap", "tsne")) {
        dims_plot <- c(1,2)
      } else {
        dims_plot <- dims_plot
      }

      plot_dim <- .plot_dims(feat_len = length(qc_to_plot))
      png(paste0(plot_dir, "02_qc_", red, ".png"), width = plot_dim$width,
          height = plot_dim$height, res = fig.res, units = "in")
      plot(feature_plot(seu = seu, reduction = red, features = qc_to_plot,
                        max.cutoff = max.cutoff, ncol = plot_dim$ncols,
                        legend.position = legend.position,
                        col_pal = cont_col_pal, dims_plot = dims_plot,
                        pt.size = pt.size, ...) &
             Seurat::NoAxes())
      dev.off()
    }

    # Heatmap showing correlation of QCs with Principal Components
    if ("pca" %in% reductions) {
      png(paste0(plot_dir, "pca_qc_cor.png"), width = 12,
          height = 5, res = fig.res, units = "in")
      plot(pca_feature_cor_plot(seu = seu, features = qc_to_plot))
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
#' data normalisation, 3. identification of HVGs, 4. Scaling of expression
#' values and computing PCA, 5. running Harmony integration, 8. Generating
#' plots.
#'
#' @param seu Seurat object (required).
#' @param npcs Number of principal components.
#' @param dims.use Vector with PCs to use for Harmony integration, e.g. 1:50.
#' @param plot_dir Directory to save generated plots. If NULL, plots are
#' not saved.
#' @param ... Additional named parameters passed to RunHarmony and other Seurat
#' processing functions, such as RunPCA and ScaleData.
#' @inheritParams run_harmony_pipeline
#'
#' @return Updated Seurat object with integrated and processed data.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
harmony_analysis <- function(
    seu, batch_id = "sample", npcs, dims.use = NULL, plot_dir = NULL,
    n_hvgs = 3000, max.iter.harmony = 50, seed = 1, fig.res = 200, ...) {
  # If list, then we have un-merged independent samples
  if (is.list(seu)) {
    # Normalise, obtain HVGs and run PCA
    seu <- lognormalize_and_pca(seu, npcs = NULL, n_hvgs = n_hvgs, ...)
    # Plot HVGs
    for (s in names(seu)) {
      png(paste0(plot_dir, "hvgs_", s, ".png"), width = 10, height = 5,
          res = fig.res, units = "in")
      print(Seurat::LabelPoints(plot = Seurat::VariableFeaturePlot(seu[[s]]),
               points = head(Seurat::VariableFeatures(seu[[s]]), 10), repel = TRUE))
      dev.off()
    }
    # Merge all samples
    seu <- merge(x = seu[[1]], y = seu[2:length(seu)])
    # Seurat merge bug which converts factors to character.
    seu <- .as_factor_metadata(seu)
  }
  # Drop unused factor levels - mostly done when a Seurat object is filtered
  seu <- .drop_factors(seu)

  # Process merged data and run PCA
  seu <- lognormalize_and_pca(seu, npcs = npcs, n_hvgs = n_hvgs, ...)

  # Run Harmony
  set.seed(seed) # Set seed due to Harmony being stochastic
  seu <- run_harmony(object = seu, group.by.vars = batch_id,
                     reduction = "pca", dims.use = dims.use,
                     max.iter.harmony = max.iter.harmony, ...)

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
  png(paste0(plot_dir, "pca_elbow.png"), width = 10, height = 6,
      res = fig.res, units = "in")
  print(Seurat::ElbowPlot(seu, ndims = npcs, reduction = "pca"))
  dev.off()
  png(paste0(plot_dir, "harmony_elbow.png"), width = 10, height = 6,
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
#' @param col_pal Continuous colour palette to use, default "RdYlBu".
#' @param dims_plot Dimensions to plot, must be a two-length numeric vector
#' specifying x- and y-dimensions.
#' @param alpha Controls opacity of spots. Provide as a vector specifying the
#' min and max range of values (between 0 and 1).
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
                                  legend.position = "top", col_pal = NULL,
                                  dims_plot = c(1, 2), seed = 1, ctrl = 100,
                                  pt.size = 1.4, fig.res = 200,
                                  alpha = c(0.1, 0.9), pt.size.factor = 1.6,
                                  spatial_col_pal = "inferno", crop = TRUE, ...) {
  # If no modules are given, return Seurat object
  if (is.null(modules_group)) {
    return(seu)
  }

  # Extract default assay
  assay <- Seurat::DefaultAssay(object = seu)

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
        plot(feature_plot(seu = seu, reduction = reduction,
                          features = features,
                          max.cutoff = max.cutoff, ncol = plot_dim$ncols,
                          legend.position = legend.position, col_pal = col_pal,
                          dims_plot = dims_plot, pt.size = pt.size,
                          combine = TRUE, ...) & Seurat::NoAxes())
        dev.off()

        # If assay is "Spatial" plot expression on tissue as well
        if (assay == "Spatial") {
          spat_plot_dim <- .spatial_plot_dims(feat_len = length(features))
          pdf(paste0(plot_dir, "01_markers_spatial_", m, ".pdf"), width = spat_plot_dim$width,
              height = spat_plot_dim$height, useDingbats = FALSE)
          print(spatial_feature_plot(
            seu, features = features, alpha = alpha, pt.size.factor = pt.size.factor,
            ncol = spat_plot_dim$ncols, max.cutoff = max.cutoff,
            crop = TRUE, col_pal = "inferno", legend.position = "top", ...))
          dev.off()
        }
      }
    }
    if (!is.null(plot_dir)) {
      # Module scores in one plot
      plot_dim <- .plot_dims(feat_len = length(modules))
      png(paste0(plot_dir, "02_score_", mg, ".png"), width = plot_dim$width,
          height = plot_dim$height, res = fig.res, units = "in")
      plot(feature_plot(seu = seu, reduction = reduction,
                        features = names(modules),
                        max.cutoff = max.cutoff, ncol = plot_dim$ncols,
                        legend.position = legend.position, col_pal = col_pal,
                        dims_plot = dims_plot, pt.size = pt.size,
                        combine = TRUE, ...) & Seurat::NoAxes())
      dev.off()

      # If assay is "Spatial" plot expression on tissue as well
      if (assay == "Spatial") {
        spat_plot_dim <- .spatial_plot_dims(feat_len = length(modules))
        pdf(paste0(plot_dir, "02_score_spatial_", mg, ".png"), width = spat_plot_dim$width,
            height = spat_plot_dim$height, useDingbats = FALSE)
        print(spatial_feature_plot(
          seu, features = features, alpha = alpha, pt.size.factor = pt.size.factor,
          ncol = spat_plot_dim$ncols, max.cutoff = max.cutoff,
          crop = crop, col_pal = spatial_col_pal, legend.position = "top", ...))
        dev.off()
      }
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
#' @param dims Vector denoting dimensions to use for nearest neighnors and
#' clustering (from 'cluster_reduction' parameter below).
#' @param res Vector with clustering resolutions (e.g. seq(0.1, 0.6, by = 0.1)).
#' @param logfc.threshold Limit testing to genes which show, on average, at
#' least X-fold difference (log-scale) between the two groups of cells.
#' @param min.pct Only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations.
#' @param only.pos Only return positive markers (TRUE by default).
#' @param topn_genes Top cluster marker genes to use for plot (in heatmap and
#' feature plots), default is 10.
#' @param plot_dir Directory to save generated plots. If NULL, plots are
#' not saved.
#' @param plot_cluster_markers Logical, wheather to create feature plots with
#' 'topn_genes' cluster markers. Added mostly to reduce number of files (and size)
#' in analysis folders. Default is TRUE.
#' @param modules_group Group of modules (named list of lists) storing features
#' (e.g. genes) to compute module score for each identified cluster. This step
#' can be useful for annotating the different clusters by saving dot
#' plots for each group. Assumes that we already have computed the modules e.g.
#' by calling the 'module_score_analysis' function. If 'plot_dir' is NULL,
#' no plots will be generated.
#' @param cluster_reduction Dimensionality reduction to use for performing
#' clustering. Default is 'pca', should be set to 'harmony' if we perform data
#' integration.
#' @param plot_reduction Dimensionality reduction to use for plotting
#' functions. Default is 'umap'.
#' @param max.cutoff Vector of maximum cutoff values for each feature,
#' may specify quantile in the form of 'q##' where '##' is the quantile
#' (eg, 'q1', 'q10').
#' @param seed Set a random seed, for reproducibility.
#' @param ctrl Number of control features selected from the same bin per
#' analyzed feature.
#' @param force_reanalysis Logical, if cluster marker genes file
#' exists and force_reanalysis = FALSE, run identification of cluster markers.
#' Otherwise, read cluster markers from file. Added for computing time
#' efficiency purposes.
#' @param label Whether to label the clusters in 'plot_reduction' space.
#' @param label.size Sets size of labels.
#' @param legend.position Position of legend, default "right" (set to "none"
#' for clean plot).
#' @param pt.size Adjust point size for plotting.
#' @param cont_col_pal Continuous colour palette to use, default "RdYlBu".
#' @param discrete_col_pal Discrete colour palette to use, default is Hue palette
#' (hue_pal) from 'scales' package.
#' @param fig.res Figure resolution in ppi (see 'png' function).
#' @param cont_alpha Controls opacity of spots. Provide as a vector specifying the
#' min and max range of values (between 0 and 1).
#' @param discrete_alpha Controls opacity of spots. Provide a single alpha value.
#' @param pt.size.factor Scale the size of the spots.
#' @param spatial_col_pal Continuous colour palette to use from viridis package to
#' colour spots on tissue, default "inferno".
#' @param crop Crop the plot in to focus on points plotted. Set to FALSE to
#' show entire background image.
#' @param ... Additional named parameters passed to Seurat
#' analysis and plotting functions, such as FindClusters, FindAllMarkers,
#' DimPlot and FeaturePlot.
#'
#' @return Updated Seurat object clustered cells
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
cluster_analysis <- function(seu, dims = 1:20, res = seq(0.1, 0.1, by = 0.1),
                             logfc.threshold = 0.5, min.pct = 0.25,
                             only.pos = TRUE, topn_genes = 10, plot_dir = NULL,
                             plot_cluster_markers = TRUE, modules_group = NULL,
                             cluster_reduction = "pca", plot_reduction = "umap",
                             max.cutoff = "q98", seed = 1,
                             ctrl = 100, force_reanalysis = TRUE,
                             label = TRUE, label.size = 8,
                             legend.position = "right", pt.size = 1.4,
                             cont_col_pal = NULL, discrete_col_pal = NULL,
                             fig.res = 200, cont_alpha = c(0.1, 0.9),
                             discrete_alpha = 0.6, pt.size.factor = 1.6,
                             spatial_col_pal = "inferno", crop = TRUE, ...) {
  # So CMD passes without NOTES
  cluster = avg_log2FC = seurat_clusters = condition = sample = freq = n <- NULL
  assertthat::assert_that(methods::is(seu, "Seurat"))
  # Drop unused factor levels - mostly done when a Seurat object is filtered
  seu <- .drop_factors(seu)
  # Extract default assay
  assay <- Seurat::DefaultAssay(object = seu)

  # Iterate over each clustering resolution
  for (r in res) {
    cat("Res", r, "\n")
    # Identify nearest neighbors and perform clustering
    seu <- find_neighbors(seu = seu, dims = dims,
                          reduction = cluster_reduction, ...)
    seu <- find_clusters(seu = seu, resolution = r, random.seed = seed, ...)
    # Store metadata with obtained clusters locally
    write.csv(seu@meta.data, file = paste0(plot_dir, "seu_meta_res",r,".csv"))

    if (!is.null(plot_dir)) {
      plot_dim <- .plot_dims(feat_len = 1)
      png(paste0(plot_dir,"z_", plot_reduction, "_res", r, ".png"),
          width = plot_dim$width, height = plot_dim$height, res = fig.res,
          units = "in")
      plot(dim_plot(seu = seu, reduction = plot_reduction, split.by = NULL,
                    group.by = "seurat_clusters", ncol = 1,
                    col_pal = discrete_col_pal, legend.position=legend.position,
                    dims_plot = c(1,2), pt.size = pt.size, label = label,
                    label.size=label.size, combine=TRUE,...) & Seurat::NoAxes())
      dev.off()

      # If assay is "Spatial" plot expression on tissue as well
      if (assay == "Spatial") {
        spat_plot_dim <- .spatial_plot_dims(feat_len = 1)
        pdf(paste0(plot_dir, "z_cluster_res", r, ".pdf"), width = spat_plot_dim$width,
            height = spat_plot_dim$height, useDingbats = FALSE)
        print(spatial_dim_plot(
          seu, group.by = "seurat_clusters", alpha = discrete_alpha,
          pt.size.factor = pt.size.factor, ncol = spat_plot_dim$ncols,
          crop = crop, col_pal = discrete_col_pal, legend.position = "top", ...))
        dev.off()
      }
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
          png(paste0(plot_dir, "dotplot_module_", mg, "_res", r, ".png"),
              width = 6 + 0.3*nlevels(seu@meta.data[["seurat_clusters"]]),
              height = 3 + 0.15*length(m_names), res = fig.res, units = "in")
          plot(dot_plot(seu = seu, features = m_names, labels = NULL,
                        group.by = "seurat_clusters", xlab = "Signature",
                        ylab = "Cluster", legend.position = "right",
                        col_pal = cont_col_pal, ...))
          dev.off()
        }
      }
    }

    # Identify cluster markers
    if (!file.exists(paste0(plot_dir, "seu_markers_res", r, ".csv")) ||
        force_reanalysis == TRUE) {
      mark <- find_all_markers(seu = seu, random.seed = seed,
                               only.pos = only.pos, min.pct = min.pct,
                               logfc.threshold = logfc.threshold, ...)
      write.csv(mark, file = paste0(plot_dir, "seu_markers_res", r, ".csv"))
    } else {
      mark <- read.csv(file = paste0(plot_dir, "seu_markers_res", r, ".csv"))
      # Make cluster a factor from character
      mark$cluster <- factor(mark$cluster)
    }

    # Heatmap of marker genes
    heatmap_plot(seu = seu, markers = mark, topn_genes = topn_genes,
                 filename = paste0(plot_dir, "z_heatmap_res", r, ".png"), ...)

    ## Feature and violin plots
    if (plot_cluster_markers & !is.null(plot_dir)) {
      # Extract top marker genes
      top_mark <- mark[mark$gene %in% rownames(x = seu), ] |>
        dplyr::group_by(cluster) |>
        dplyr::slice_max(n = topn_genes, order_by = avg_log2FC)
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
                            ncol = plot_dim$ncols, col_pal = cont_col_pal,
                            legend.position = legend.position, dims_plot = c(1,2),
                            pt.size = pt.size, ...) & Seurat::NoAxes())
          dev.off()

          png(paste0(plot_dir, "02_violin_seu_res", r, "_cl", cl,".png"),
              width = plot_dim$width, height = plot_dim$height, res = 100,
              units = "in")
          plot(Seurat::VlnPlot(seu, features = genes, ncol = plot_dim$ncols,
                               pt.size = 0.03, ...))
          dev.off()

          # If assay is "Spatial" plot expression on tissue as well
          if (assay == "Spatial") {
            spat_plot_dim <- .spatial_plot_dims(feat_len = length(genes))
            pdf(paste0(plot_dir, "01_feature_spatial_seu_res", r, "_cl", cl, ".pdf"),
                width = spat_plot_dim$width, height = spat_plot_dim$height, useDingbats = FALSE)
            plot(spatial_feature_plot(
              seu, features = genes, alpha = cont_alpha, pt.size.factor = pt.size.factor,
              ncol = spat_plot_dim$ncols, max.cutoff = max.cutoff,
              crop = TRUE, col_pal = spatial_col_pal, legend.position = "top", ...))
            dev.off()
          }

          # Each cluster as a module and compute score
          seu <- compute_module_score(seu = seu, features = genes, ctrl = ctrl,
                                      name = paste0("Cluster", cl), ...)
        }
      }
      # Module scores in one plot
      if (NROW(top_mark) >0) {
        plot_dim <- .plot_dims(feat_len = length(unique(top_mark$cluster)))
        png(paste0(plot_dir, "03_score_seu_res", r, ".png"),
            width = plot_dim$width, height = plot_dim$height, res = fig.res,
            units = "in")
        plot(feature_plot(seu = seu, reduction = plot_reduction,
                          features = paste0("Cluster", unique(top_mark$cluster)),
                          max.cutoff = max.cutoff, ncol = plot_dim$ncols,
                          col_pal = cont_col_pal,
                          legend.position = legend.position, dims_plot = c(1,2),
                          pt.size = pt.size, ...) & Seurat::NoAxes())
        dev.off()

        # If assay is "Spatial" plot expression on tissue as well
        if (assay == "Spatial") {
          spat_plot_dim <- .spatial_plot_dims(feat_len = length(unique(top_mark$cluster)))
          pdf(paste0(plot_dir, "03_score_spatial_seu_res", r, ".pdf"),
              width = spat_plot_dim$width, height = spat_plot_dim$height, useDingbats = FALSE)
          plot(spatial_feature_plot(
            seu, features = paste0("Cluster", unique(top_mark$cluster)),
            alpha = cont_alpha, pt.size.factor = pt.size.factor,
            ncol = spat_plot_dim$ncols, max.cutoff = max.cutoff,
            crop = TRUE, col_pal = spatial_col_pal, legend.position = "top", ...))
          dev.off()
        }
      }
    }

    # Condition specific histogram plots
    if (length(unique(seu$condition)) > 1) {
      # Update discrete colour palette if NULL
      if (is.null(discrete_col_pal)) {
        if (nlevels(seu$condition) > 35) {
          col_pal <- scales::hue_pal()(nlevels(seu$condition))
        } else {
          col_pal <- .internal_col_pal()[1:nlevels(seu$condition)]
        }
        names(col_pal) <- levels(seu$condition)
      } else {
        col_pal <- discrete_col_pal
      }
      cl_condition <- seu@meta.data |> dplyr::group_by(seurat_clusters, condition) |>
        dplyr::summarise(n = n()) |> dplyr::mutate(freq = n / sum(n))

      png(paste0(plot_dir, "hist_condition_contr_per_cluster_res", r, ".png"),
          width = 7 + 0.3*nlevels(seu@meta.data[["seurat_clusters"]]),
          height = 5, res = 150, units = "in")
      plot(ggplot2::ggplot(data = cl_condition,
                           ggplot2::aes(x = seurat_clusters, y = freq, fill = condition)) +
        ggplot2::geom_bar(stat = "identity", color = "black") +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
              plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)) +
        ggplot2::xlab(NULL) + ggplot2::ylab("Cluster contribution") + ggplot2::ggtitle(NULL) +
        ggplot2::scale_fill_manual(values = col_pal))
      dev.off()

      cl_condition <- seu@meta.data |> dplyr::group_by(condition, seurat_clusters) |>
        dplyr::summarise(n = n()) |> dplyr::mutate(freq = n / sum(n))

      png(paste0(plot_dir, "hist_condition_contr_res", r, ".png"),
          width = 7 + 0.3*nlevels(seu@meta.data[["seurat_clusters"]]),
          height = 5, res = 150, units = "in")
      plot(ggplot2::ggplot(data = cl_condition,
                           ggplot2::aes(x = seurat_clusters, y = freq, fill = condition)) +
             ggplot2::geom_bar(stat = "identity", color = "black", position = ggplot2::position_dodge()) +
             ggplot2::theme_classic() +
             ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)) +
             ggplot2::xlab(NULL) + ggplot2::ylab("Cluster contribution") + ggplot2::ggtitle(NULL) +
          ggplot2::scale_fill_manual(values = col_pal))
      dev.off()
    }

    # Sample specific histogram plots
    if (length(unique(seu$sample)) > 1) {
      # Update discrete colour palette if NULL
      if (is.null(discrete_col_pal)) {
        if (nlevels(seu$sample) > 35) {
          col_pal <- scales::hue_pal()(nlevels(seu$sample))
        } else {
          col_pal <- .internal_col_pal()[1:nlevels(seu$sample)]
        }
        names(col_pal) <- levels(seu$sample)
      } else {
        col_pal <- discrete_col_pal
      }

      cl_condition <- seu@meta.data |> dplyr::group_by(seurat_clusters, sample) |>
        dplyr::summarise(n = n()) |> dplyr::mutate(freq = n / sum(n))

      png(paste0(plot_dir, "hist_sample_contr_per_cluster_res", r, ".png"),
          width = 7 + 0.3*nlevels(seu@meta.data[["seurat_clusters"]]),
          height = 5, res = 150, units = "in")
      plot(ggplot2::ggplot(data = cl_condition,
                           ggplot2::aes(x = seurat_clusters, y = freq, fill = sample)) +
             ggplot2::geom_bar(stat = "identity", color="black") +
             ggplot2::theme_classic() +
             ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)) +
             ggplot2::xlab(NULL) + ggplot2::ylab("Cluster contribution") + ggplot2::ggtitle(NULL) +
             ggplot2::scale_fill_manual(values = col_pal))
      dev.off()

      cl_condition <- seu@meta.data |> dplyr::group_by(sample, seurat_clusters) |>
        dplyr::summarise(n = n()) |> dplyr::mutate(freq = n / sum(n))

      png(paste0(plot_dir, "hist_sample_contr_res", r, ".png"),
          width = 7 + 0.3*nlevels(seu@meta.data[["seurat_clusters"]]),
          height = 5, res = 150, units = "in")
      plot(ggplot2::ggplot(data = cl_condition,
                           ggplot2::aes(x = seurat_clusters, y = freq, fill = sample)) +
             ggplot2::geom_bar(stat = "identity", color="black", position = ggplot2::position_dodge()) +
             ggplot2::theme_classic() +
             ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)) +
             ggplot2::xlab(NULL) + ggplot2::ylab("Cluster contribution") + ggplot2::ggtitle(NULL) +
             ggplot2::scale_fill_manual(values = col_pal))
      dev.off()
    }
  }
  return(seu)
}
