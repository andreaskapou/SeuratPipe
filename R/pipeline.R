#' @rdname run_qc_pipeline
#'
#' @title QC pipeline
#'
#' @description This function implements all the analysis steps for perfoming
#' QC. These include: 1. reading all sample information from metadata object/file
#' and generating one Seurat object per sample. 2. Performs SoupX (ambient
#' RNA removal) and Scrublet (doublet detection) if user defines the
#' corresponding parameters. 3. Filter Seurat object according to QC
#' criteria 4. Generate correspond QC plots.
#'
#' @param data_dir Parent directory where all sample 10x files are stored.
#' Think of it as project directory.
#' @param sample_meta Sample metadata information in a Data.frame like object.
#' Columns should at least contain 'sample', 'donor', 'condition' and 'pass_qc'.
#' @param sample_meta_filename Filename of sample metadata information, same as 'meta'
#' parameter above. User should provide one of 'meta' or 'meta_filename'.
#' @param nfeat_thresh Filter cells that have less than 'nfeat_thresh' counts
#' expressed.
#' @param mito_thresh Filter cells with more than 'mito_thresh'% counts.
#' @param meta_colnames Sample metadata column names to store in Seurat metadata.
#' @param out_dir Output directory for storing analysis results.
#' @param qc_to_plot Vector of features in metadata to plot.
#' @param use_scrublet Logical, wether to use Scrublet for doublet detection.
#' @param use_soupx Logical, wether to use SoupX for ambient RNA removal.
#' @param tenx_dir Name of 10x base directory, e.g. with outputs after running
#' cellranger. Default 'premrna_outs', i.e. assumes single-nuclei RNA-seq.
#' @param tenx_counts_dir Name of 10x directory where count matrices are stored.
#' Default 'filtered_feature_bc_matrix'
#' @param obj_filename Filename of the stored Seurat object, default 'seu_qc'.
#' @param expected_doublet_rate The expected fraction of transcriptomes that
#' are doublets, typically 0.05 - 0.1
#' @param force_reanalysis Logical, if intermediate file 'seu_preqc.rds' (with
#' created Seurat object) exists and force_reanalysis = FALSE, read object
#' instead of re-running whole analysis with soupX. Added for computing time
#' efficiency purposes and intermediate object will be created only when
#' 'use_soupx = TRUE'.
#' @param min.cells Include features/genes detected in at least this many cells.
#' @param min.features Include cells where at least this many features/genes
#' are detected.
#' @param ... Additional named parameters passed to Seurat, Scrublet or SoupX.
#'
#' @return List of Seurat objects as the length of the number of samples in
#' the sample metadata file. If a single sample, return a Seurat object
#' instead of a list.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
run_qc_pipeline <- function(
    data_dir, sample_meta, sample_meta_filename = NULL, nfeat_thresh = 500,
    mito_thresh = 5, meta_colnames = c("donor", "condition", "pass_qc"),
    out_dir = NULL, qc_to_plot = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
    use_scrublet = TRUE, use_soupx = FALSE, tenx_dir = "premrna_outs",
    tenx_counts_dir = "filtered_feature_bc_matrix", obj_filename = "seu_qc",
    expected_doublet_rate = 0.06, force_reanalysis = TRUE, min.cells = 10,
    min.features = 100, ...) {

  # Store all parameters for reproducibility
  opts <- c(as.list(environment()), list(...))
  # f  <- rlang::fn_fmls()
  # mc <- as.list(match.call())[-1]
  # opts <- append(mc, f[!names(f) %in% names(mc)])[names(f)]

  # SO CMD passes without NOTES
  x = y = dens = plot_dir = seu <- NULL

  # Check and return if metadata object or filename is given as input
  opts$sample_meta <- .get_metadata(sample_meta = sample_meta,
                                    sample_meta_filename = sample_meta_filename)
  # Create plot directory if it doesn't exist
  if (!is.null(out_dir)) {
    plot_dir <- paste0(out_dir, "/plots/")
    if (!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  }

  if (!file.exists(paste0(out_dir, "/", "seu_preqc.rds")) ||
      force_reanalysis == TRUE) {
    # Create Seurat object
    seu <- create_seurat_object(
      data_dir = data_dir, sample_meta = opts$sample_meta, sample_meta_filename = NULL,
      meta_colnames = meta_colnames, plot_dir = plot_dir, use_scrublet = use_scrublet,
      use_soupx = use_soupx, tenx_dir = tenx_dir, tenx_counts_dir = tenx_counts_dir,
      expected_doublet_rate = expected_doublet_rate,
      min.cells = min.cells, min.features = min.features, ...)
    # Save pre-QC Seurat object and opts if we use SoupX that takes a long time
    if (use_soupx) {
      saveRDS(object = list(seu = seu, opts = opts),
              file = paste0(out_dir, "/", "seu_preqc.rds"))
    }
  } else {
    obj <- readRDS(file = paste0(out_dir, "/", "seu_preqc.rds"))
    seu <- obj$seu
    # Keep only samples that are present in current metadata file
    if (is.list(seu)) { seu <- seu[opts$sample_meta$sample] }
    base::rm(obj)
    base::gc(verbose = FALSE)
  }

  # In case we have one sample, we make the Seurat object a list
  if (!is.list(seu)) {
    # Extract sample name
    sample <- opts$sample_meta$sample
    seu <- list(seu)
    names(seu) <- sample
  }
  # Instantiate default assay by looking at DefaultAssay of first sample
  assay <- Seurat::DefaultAssay(object = seu[[1]])

  # Plot QC prior to filtering
  if (!is.null(plot_dir) && !is.null(qc_to_plot)) {
    plot_dim <- .plot_dims(feat_len = length(qc_to_plot))
    pdf(paste0(plot_dir, "vln_preqc.pdf"), width = 10, height = 6,
        useDingbats = FALSE)
    for (s in names(seu)) {
      print(Seurat::VlnPlot(seu[[s]], features = qc_to_plot,
                            ncol = plot_dim$ncols, pt.size = 0.05, ...))
    }
    dev.off()

    # Get all combinations of QCs to plot
    combs <- utils::combn(qc_to_plot, 2)
    plot_dim <- .plot_dims(feat_len = NCOL(combs))

    # Create scatter plots for feature-to-feature relationships
    pdf(paste0(plot_dir, "scatter_preqc.pdf"), width = plot_dim$width,
        height = plot_dim$height, useDingbats = FALSE)
    for (s in names(seu)) {
      gg_list <- scatter_meta_plot(seu = seu[[s]], features = qc_to_plot)
      plot(patchwork::wrap_plots(gg_list, ncol = plot_dim$ncols) +
             patchwork::plot_annotation(title = s, theme = ggplot2::theme(
               plot.title = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                  size = 20))))
    }
    dev.off()
  }

  # Pre-QC summary
  if (!is.null(plot_dir)) {
    preqc_summary <- data.frame(sample = character(), cells = numeric(),
                                median_nGenes = numeric(),
                                median_nCounts = numeric(),
                                median_mito = numeric(),
                                prop_filt_nGenes = numeric(),
                                prop_filt_mito = numeric())
    for (s in names(seu)) {
      tmp <- seu[[s]][[]]
      preqc_summary <- rbind(preqc_summary, data.frame(
        sample = s, cells = NROW(tmp),
        median_nGenes = round(median(tmp[[paste0("nFeature_", assay)]])),
        median_nCounts = round(median(tmp[[paste0("nCount_", assay)]])),
        median_mito = round(median(tmp$percent.mito), 2),
        prop_filt_nGenes = round(
          sum(tmp[[paste0("nFeature_", assay)]] <= nfeat_thresh) / NROW(tmp), 2),
        prop_filt_mito = round(
          sum(tmp$percent.mito >= mito_thresh &
                tmp[[paste0("nFeature_", assay)]] > nfeat_thresh) / NROW(tmp), 2)
      ))
    }
    write.csv(preqc_summary, file = paste0(plot_dir, "preqc_sample_summary.csv"))
  }

  # Perform actual QC filtering
  seu <- qc_filter_seurat_object(seu = seu, nfeat_thresh = nfeat_thresh,
                                 mito_thresh = mito_thresh)

  if (!is.null(plot_dir)) {
    # Plot QC after filtering
    if (!is.null(plot_dir) && !is.null(qc_to_plot)) {
      plot_dim <- .plot_dims(feat_len = length(qc_to_plot))
      pdf(paste0(plot_dir, "vln_qc.pdf"), width = 10, height = 6,
          useDingbats = FALSE)
      for (s in names(seu)) {
        print(Seurat::VlnPlot(seu[[s]], features = qc_to_plot,
                              ncol = plot_dim$ncols, pt.size = 0.05, ...))
      }
      dev.off()
    }

    # Post-QC summary
    qc_summary <- data.frame(sample = character(), cells = numeric(),
                             median_nGenes = numeric(),
                             median_nCounts = numeric(),
                             median_mito = numeric())
    for (s in names(seu)) {
      tmp <- seu[[s]][[]]
      qc_summary <- rbind(qc_summary, data.frame(
        sample = s, cells = NROW(tmp),
        median_nGenes = round(median(tmp[[paste0("nFeature_", assay)]])),
        median_nCounts = round(median(tmp[[paste0("nCount_", assay)]])),
        median_mito = round(median(tmp$percent.mito), 2)
      ))
    }
    write.csv(qc_summary, file = paste0(plot_dir, "qc_sample_summary.csv"))
  }

  # If we have single sample, remove it from list so we return Seurat object
  if (length(seu) == 1) { seu <- seu[[1]] }

  # Save Seurat object and opts
  if (!is.null(out_dir)) {
    if (obj_filename == "") { obj_filename <- "seu_qc.rds" }
    saveRDS(object = list(seu = seu, opts = opts),
                          file = paste0(out_dir,"/",obj_filename,".rds"))
  }
  # Return QCed Seurat object
  return(seu)
}


#' @rdname run_harmony_pipeline
#'
#' @title Pipeline for Harmony integration
#'
#' @description This function implements all the analysis steps for perfoming
#' data integration using Harmony. These include: 1. data processing, e.g.
#' normalisation, PCA. 2. Running Harmony 3. Perform UMAP and clustering after
#' data integration. Analysis outputs are stored in corresponding directories.
#'
#' @param seu_obj Seurat object or list of Seurat objects(required).
#' @param out_dir Output directory for storing analysis results.
#' @param batch_id Name of batch to try and remove with data integration
#' (required). Can also be a vector if multiple batch information are present.
#' Should be a column name in Seurat 'meta.data'. Default is "sample".
#' This parameter is called 'group.by.vars' in Harmony.
#' @param npcs Number of principal components, can be a vector e.g. c(50, 70).
#' @param ndims Top Harmony dimensions to perform UMAP and clustering,
#' can be a vector e.g. c(50, 70).
#' @param res Vector with clustering resolutions (e.g. seq(0.1, 0.6, by = 0.1)).
#' @param modules_group Group of modules (named list of lists) storing features
#' (e.g. genes) to compute module score for each identified cluster. This step
#' can be useful for annotating the different clusters by saving dot/feature
#' plots for each group.
#' @param metadata_to_plot Vector with metadata names to plot, they should be
#' present in the meta.data slot of the Seurat object.
#' @param qc_to_plot Vector with QC names to plot, they should be
#' present in the meta.data slot of the Seurat object.
#' @param logfc.threshold Limit testing to genes which show, on average, at
#' least X-fold difference (log-scale) between the two groups of cells.
#' @param min.pct Only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations.
#' @param only.pos Only return positive markers (TRUE by default).
#' @param topn_genes Top cluster marker genes to use for plotting (in heatmap and
#' feature plots), default is 10.
#' @param pcs_to_remove Which PCs should be removed prior to running Harmony.
#' Possibly due to being correlated with technical/batch effects. If NULL,
#' all PCs are used.
#' @param obj_filename Filename of the stored Seurat object, default
#' 'seu_harmony'. Number of PCs will be added to the filename automatically.
#' @param force_reanalysis Logical, if intermediate object 'seu_harmony_<>.rds'
#' exists and force_reanalysis = FALSE, read object instead of re-running
#' Harmony integration. Added for computing time efficiency purposes.
#' @param plot_cluster_markers Logical, wheather to create feature plots with
#' 'topn_genes' cluster markers. Added mostly to reduce number of files
#' (and size) in analysis folders. Default is TRUE.
#' @param max.cutoff Maximum cutoff values for plotting each continuous
#' feature, e.g. gene expression levels. May specify quantile in the form of
#' 'q##' where '##' is the quantile (eg, 'q1', 'q10').
#' @param min.cutoff Maximum cutoff values for plotting each continuous
#' feature, e.g. gene expression levels. May specify quantile in the form of
#' 'q##' where '##' is the quantile (eg, 'q1', 'q10').
#' @param n_hvgs Number of highly variable genes (HVGs) to compute, which will
#' be used as input to PCA.
#' @param max.iter.harmony Maximum number of iterations for Harmony integration.
#' @param seed Set a random seed, for reproducibility.
#' @param label Whether to label the clusters in 'plot_reduction' space.
#' @param label.size Sets size of labels.
#' @param pt.size Adjust point size for plotting.
#' @param fig.res Figure resolution in ppi (see 'png' function).
#' @param cont_col_pal Continuous colour palette to use, default "RdYlBu".
#' @param discrete_col_pal Discrete colour palette to use, default is Hue palette
#' (hue_pal) from 'scales' package.
#' @param cont_alpha Controls opacity of spots. Provide as a vector specifying the
#' min and max range of values (between 0 and 1).
#' @param discrete_alpha Controls opacity of spots. Provide a single alpha value.
#' @param pt.size.factor Scale the size of the spots.
#' @param spatial_col_pal Continuous colour palette to use from viridis package to
#' colour spots on tissue, default "inferno".
#' @param crop Crop the plot in to focus on points plotted. Set to FALSE to
#' show entire background image.
#' @param plot_spatial_markers Logical, whether to create spatial feature plots
#' with expression of individual genes.
#' @param ... Additional named parameters passed to Seurat's or Harmony
#' functions.
#'
#' @return An updated Seurat object. Note that if multiple `npcs` and `ndims`
#' are given, only the last setting will be returned. All analysis results are
#' also stored on disk.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
run_harmony_pipeline <- function(
    seu_obj, out_dir, batch_id = "sample", npcs = c(50), ndims = c(30),
    res = seq(0.1, 0.3, by = 0.1), modules_group = NULL,
    metadata_to_plot = c("sample", "condition"), qc_to_plot = NULL,
    logfc.threshold = 0.5, min.pct = 0.25, only.pos = TRUE, topn_genes = 10,
    pcs_to_remove = NULL, obj_filename = "seu_harmony", force_reanalysis = TRUE,
    plot_cluster_markers = TRUE, max.cutoff = "q98", min.cutoff = NA, n_hvgs = 3000,
    max.iter.harmony = 50, seed = 1, label = TRUE, label.size = 8,
    pt.size = 1.4, fig.res = 200, cont_col_pal = NULL, discrete_col_pal = NULL,
    cont_alpha = c(0.1, 0.9), discrete_alpha = 0.9, pt.size.factor = 1.6,
    spatial_col_pal = "inferno", crop = TRUE,
    plot_spatial_markers = FALSE, ...) {

  # Store all parameters for reproducibility
  opts <- c(as.list(environment()), list(...))
  # Do not store the Seurat object in opts
  opts$seu_obj <- NULL
  # So CMD passes without NOTES
  seu <- NULL

  if (obj_filename == "") { obj_filename <- "seu_harmony.rds" }

  pcs_remove_name <- ""
  # Iterate over the number of principal components
  for (npc in npcs) {
    # Define number of dimensions to use
    dims.use <- seq(from = 1, to = npc, by = 1)
    if (!is.null(pcs_to_remove)) {
      pcs_remove_name <- paste0("_r", paste(pcs_to_remove, collapse = ""))
      dims.use <- dims.use[-pcs_to_remove]
    }

    obj_name <- paste0("_npcs", npc, pcs_remove_name)
    npc_dir <- paste0(out_dir, "/", obj_filename, obj_name, "/")
    if (!dir.exists(npc_dir)) {dir.create(npc_dir, recursive = TRUE)}

    qc_dir <- paste0(npc_dir, "/qc/")
    if (!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}

    # Run Harmony integration analysis
    if (!file.exists(paste0(out_dir, "/", obj_filename, obj_name, ".rds")) ||
        force_reanalysis == TRUE) {
      seu <- harmony_analysis(seu = seu_obj, batch_id = batch_id,
                              npcs = npc, dims.use = dims.use,
                              plot_dir = qc_dir, n_hvgs = n_hvgs,
                              max.iter.harmony = max.iter.harmony,
                              seed = seed, fig.res = fig.res, ...)
      opts$npcs <- npc
      # Store final object prior to exploring the data
      saveRDS(list(seu = seu, opts = opts),
              file = paste0(out_dir, "/", obj_filename, obj_name, ".rds"))
    } else{
      obj <- readRDS(paste0(out_dir, "/", obj_filename, obj_name, ".rds"))
      seu <- obj$seu
      # Number of principal components
      opts$npcs <- npc
      base::rm(obj)
      base::gc(verbose = FALSE)
    }

    # Explore the data
    for (ndim in ndims) {
      if (ndim > length(dims.use)) {
        message("Skipping analysis: ndim larger than harmony dimensions.")
        next
      }
      # Create directories
      dim_dir <- paste0(npc_dir, "/dim", ndim, "/")
      qc_dir <- paste0(dim_dir, "/qc/")
      if (!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}

      # Run UMAP
      seu <- run_umap(seu, dims = 1:ndim, reduction="harmony", seed=seed, ...)
      # Store UMAP embedding for reproducibility
      write.csv(seu@reductions[["umap"]]@cell.embeddings,
                file = paste0(dim_dir, "umap_embedding.csv"))

      # QC plots
      dimred_qc_plots(seu = seu, reductions = c("umap", "pca", "harmony"),
                      metadata_to_plot = metadata_to_plot,
                      qc_to_plot = qc_to_plot, plot_dir = qc_dir,
                      max.cutoff = max.cutoff, min.cutoff = min.cutoff,
                      legend.position = "right",
                      cont_col_pal = cont_col_pal,
                      discrete_col_pal = discrete_col_pal,
                      dims_plot = c(1,2), pt.size = pt.size,
                      fig.res = fig.res, ...)

      # Compute module scores for modules of interest
      module_dir <- paste0(dim_dir, "/modules/")
      if (!dir.exists(module_dir)) {dir.create(module_dir, recursive = TRUE)}
      seu <- module_score_analysis(seu = seu, modules_group = modules_group,
                                   plot_dir = module_dir, reduction = "umap",
                                   max.cutoff = max.cutoff, min.cutoff = min.cutoff,
                                   legend.position = "right", col_pal = cont_col_pal,
                                   dims_plot = c(1,2), seed = seed, ctrl = 100,
                                   fig.res = fig.res, alpha = cont_alpha,
                                   pt.size.factor = pt.size.factor,
                                   spatial_col_pal = spatial_col_pal, crop = crop,
                                   plot_spatial_markers = plot_spatial_markers,
                                   spatial_legend_position = "top", ...)
      # Different clustering resolutions and DGE
      cl_dir <- paste0(dim_dir, "/clusters/")
      if (!dir.exists(cl_dir)) {dir.create(cl_dir, recursive = TRUE)}
      seu <- cluster_analysis(seu = seu, dims = 1:ndim, res = res,
                              logfc.threshold = logfc.threshold,
                              min.pct = min.pct, only.pos = only.pos,
                              topn_genes = topn_genes, plot_dir = cl_dir,
                              plot_cluster_markers = plot_cluster_markers,
                              modules_group = modules_group,
                              cluster_reduction = "harmony",
                              plot_reduction = "umap", max.cutoff = max.cutoff,
                              min.cutoff = min.cutoff,
                              force_reanalysis = force_reanalysis,
                              seed = seed, ctrl = 100,
                              label = label, label.size = label.size,
                              legend.position = "right", pt.size = pt.size,
                              cont_col_pal = cont_col_pal,
                              discrete_col_pal = discrete_col_pal,
                              fig.res = fig.res, cont_alpha = cont_alpha,
                              discrete_alpha = discrete_alpha,
                              pt.size.factor = pt.size.factor,
                              spatial_col_pal = spatial_col_pal, crop = crop,
                              plot_spatial_markers = plot_spatial_markers,
                              spatial_legend_position = "top", ...)
    }
  }
  return(seu)
}


#' @rdname run_cluster_pipeline
#'
#' @title Pipeline for clustering analysis
#'
#' @description This function wraps the most common Seurat analysis pipeline
#' for cell type identification. These include: 1. data processing, e.g.
#' normalisation. 2. Running PCA. 3. Perform UMAP and clustering. Analysis
#' outputs are stored in corresponding directories. Note this pipeline
#' requires a single Seurat object/sample.
#'
#' @param seu_obj Seurat object (required).
#' @param ndims Top PCA dimensions to perform UMAP and clustering,
#' can be a vector e.g. c(50, 70).
#' @param pcs_to_remove Which PCs should be removed prior to performing clustering.
#' Possibly due to being correlated with technical/batch effects. If NULL,
#' all PCs are used.
#' @param ... Additional named parameters passed to Seurat functions.
#'
#' @inheritParams run_harmony_pipeline
#'
#' @return An updated Seurat object. Note that if multiple `npcs` and `ndims`
#' are given, only the last setting will be returned. All analysis results are
#' also stored on disk.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
run_cluster_pipeline <- function(
    seu_obj, out_dir, npcs = c(50), ndims = c(30), res = seq(0.1, 0.3, by = 0.1),
    modules_group = NULL, metadata_to_plot = c("sample", "condition"),
    qc_to_plot = NULL, logfc.threshold = 0.5, min.pct = 0.25, only.pos = TRUE,
    topn_genes = 10, pcs_to_remove = NULL, plot_cluster_markers = TRUE,
    max.cutoff = "q98", min.cutoff = NA, n_hvgs = 3000, seed = 1, label = TRUE, label.size = 8,
    pt.size = 1.4, fig.res = 200, cont_col_pal = NULL, discrete_col_pal = NULL,
    cont_alpha = c(0.1, 0.9), discrete_alpha = 0.9,
    pt.size.factor = 1.6, spatial_col_pal = "inferno", crop = TRUE,
    plot_spatial_markers = FALSE, ...) {
  # Store all parameters for reproducibility
  opts <- c(as.list(environment()), list(...))
  # Do not store the Seurat object in opts
  opts$seu_obj <- NULL

  assertthat::assert_that(methods::is(seu_obj, "Seurat"))

  # So CMD passes without NOTES
  seu <- NULL

  pcs_remove_name <- ""
  # Iterate over the number of principal components
  for (npc in npcs) {
    # Define number of dimensions to use
    dims.use <- seq(from = 1, to = npc, by = 1)
    if (!is.null(pcs_to_remove)) {
      pcs_remove_name <- paste0("_r", paste(pcs_to_remove, collapse = ""))
      dims.use <- dims.use[-pcs_to_remove]
    }
    obj_name <- paste0("npcs", npc, pcs_remove_name)
    npc_dir <- paste0(out_dir, "/", obj_name, "/")
    if (!dir.exists(npc_dir)) {dir.create(npc_dir, recursive = TRUE)}

    qc_dir <- paste0(npc_dir, "/qc/")
    if (!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}

    # Process and run PCA
    seu <- lognormalize_and_pca(seu = seu_obj, npcs = npc, n_hvgs=n_hvgs, ...)
    opts$npcs <- npc

    # Plot HVGs and PCs information
    png(paste0(qc_dir, "hvgs.png"), width = 10, height = 5,
        res = fig.res, units = "in")
      print(Seurat::LabelPoints(
        plot = Seurat::VariableFeaturePlot(seu),
        points = head(Seurat::VariableFeatures(seu), 10), repel = TRUE))
    dev.off()

    png(paste0(qc_dir, "pca_heatmap.png"), width = 15, height = 15,
        res = fig.res, units = "in")
    Seurat::DimHeatmap(seu, dims = 1:9, nfeatures = 30, cells = 300,
                       reduction = "pca", balanced = TRUE)
    dev.off()
    png(paste0(qc_dir, "pca_elbow.png"), width = 12, height = 8,
        res = fig.res, units = "in")
    print(Seurat::ElbowPlot(seu, ndims = npcs, reduction = "pca"))
    dev.off()

    # Explore the data
    for (ndim in ndims) {
      if (ndim > length(dims.use)) {
        message("Skipping analysis: ndim larger than pca dimensions.")
        next
      }
      # Create directories
      dim_dir <- paste0(npc_dir, "/dim", ndim, "/")
      qc_dir <- paste0(dim_dir, "/qc/")
      if (!dir.exists(qc_dir)) {dir.create(qc_dir, recursive = TRUE)}

      # Run UMAP
      seu <- run_umap(seu, dims = dims.use[1:ndim], reduction="pca", seed=seed, ...)
      # Store UMAP embedding for reproducibility
      write.csv(seu@reductions[["umap"]]@cell.embeddings,
                file = paste0(dim_dir, "umap_embedding.csv"))

      # QC plots
      dimred_qc_plots(seu = seu, reductions = c("umap", "pca"),
                      metadata_to_plot = metadata_to_plot,
                      qc_to_plot = qc_to_plot, plot_dir = qc_dir,
                      max.cutoff = max.cutoff, min.cutoff = min.cutoff,
                      cont_col_pal = cont_col_pal, discrete_col_pal = discrete_col_pal,
                      legend.position = "right",
                      dims_plot = c(1,2), pt.size = pt.size,
                      fig.res = fig.res, ...)

      # Compute module scores for modules of interest
      module_dir <- paste0(dim_dir, "/modules/")
      if (!dir.exists(module_dir)) {dir.create(module_dir, recursive = TRUE)}
      seu <- module_score_analysis(seu = seu, modules_group = modules_group,
                                   plot_dir = module_dir, reduction = "umap",
                                   max.cutoff = max.cutoff, min.cutoff = min.cutoff,
                                   legend.position = "right",
                                   col_pal = cont_col_pal,
                                   dims_plot = c(1,2), seed = seed, ctrl = 100,
                                   fig.res = fig.res, alpha = cont_alpha,
                                   pt.size.factor = pt.size.factor,
                                   spatial_col_pal = spatial_col_pal, crop = crop,
                                   plot_spatial_markers = plot_spatial_markers,
                                   spatial_legend_position = "top", ...)
      # Different clustering resolutions and DGE
      cl_dir <- paste0(dim_dir, "/clusters/")
      if (!dir.exists(cl_dir)) {dir.create(cl_dir, recursive = TRUE)}
      seu <- cluster_analysis(seu = seu, dims = dims.use[1:ndim], res = res,
                              logfc.threshold = logfc.threshold,
                              min.pct = min.pct, only.pos = only.pos,
                              topn_genes = topn_genes, plot_dir = cl_dir,
                              plot_cluster_markers = plot_cluster_markers,
                              modules_group = modules_group,
                              cluster_reduction = "pca",
                              plot_reduction = "umap", max.cutoff = max.cutoff,
                              min.cutoff = min.cutoff,
                              force_reanalysis = TRUE,
                              seed = seed, ctrl = 100,
                              label = label, label.size = label.size,
                              legend.position = "right", pt.size = pt.size,
                              cont_col_pal = cont_col_pal,
                              discrete_col_pal = discrete_col_pal,
                              fig.res = fig.res, cont_alpha = cont_alpha,
                              discrete_alpha = discrete_alpha,
                              pt.size.factor = pt.size.factor,
                              spatial_col_pal = spatial_col_pal, crop = crop,
                              plot_spatial_markers = plot_spatial_markers,
                              spatial_legend_position = "top", ...)
    }
  }
  return(seu)
}
