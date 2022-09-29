#' @rdname create_seurat_object
#'
#' @title Create Seurat object based on sample metadata.
#'
#' @description This function creates Seurat object for each sample
#' in the metadata file. It also allows the user to perform SoupX (ambient
#' RNA removal) and Scrublet (doublet detection) analysis.
#' @param plot_dir Directory for storing QC plots. Used if use_soupx = TRUE.
#'
#' @inheritParams run_qc_pipeline
#'
#' @return List of Seurat objects equal the number of samples in
#' the sample metadata file. If a single sample, returns a Seurat object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
create_seurat_object <- function(data_dir, sample_meta = NULL, sample_meta_filename = NULL,
                                 meta_colnames = c("donor", "condition", "pass_qc"),
                                 plot_dir = NULL, use_scrublet = TRUE,
                                 use_soupx = FALSE, tenx_dir = "premrna_outs",
                                 tenx_counts_dir = "filtered_feature_bc_matrix",
                                 expected_doublet_rate = 0.06,
                                 min.cells = 10, min.features = 200, ...) {
  # All parameters passed to ...
  dot_params <- rlang::list2(...)

  pass_qc = sample <- NULL
  # Check if metadata object or filename is given as input
  sample_meta <- .get_metadata(sample_meta = sample_meta,
                               sample_meta_filename = sample_meta_filename)

  # Create plot directoy if it doesn't exist
  if (!is.null(plot_dir) && use_soupx) {
    if (!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  }

  # Create Seurat object for each sample
  seu <- list()
  for (i in 1:NROW(sample_meta)) {
    s <- sample_meta$sample[i]
    print(paste("Sample ", s, "... \n"))
    if (is.null(tenx_dir)) {
      if (is.null(sample_meta$technology[i])) {
        stop("Stopping: 'technology' column is not defined in metadata
             and 'tenx_dir' is NULL.")
      }
      if (sample_meta$technology[i] == "whole-cell") {
        tenx_dir <- "/outs/"
      } else if (sample_meta$technology[i] == "single-nuclei") {
        tenx_dir <- "/premrna_outs/"
      } else {
        stop("Stopping: 'technology' should be 'whole-cell' or 'single-nuclei'")
      }
    }

    ##
    # Load data and create intermediate Seurat object
    params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat::Read10X))]
    params <- params[which(!(names(params) %in% c("data.dir") ) )]
    tmp <-  eval(rlang::expr(Seurat::Read10X(
      data.dir = paste0(data_dir, sample_meta$path[i], "/", tenx_dir, "/", tenx_counts_dir), !!!params)))
    tmp <- Seurat::CreateSeuratObject(counts = tmp, project = s, min.cells = 0, min.features = 0)
    if (use_soupx) {
      # Run SoupX for ambient RNA correction
      params <- dot_params[which(names(dot_params) %in% methods::formalArgs(SoupX::load10X))]
      params <- params[which(!(names(params) %in% c("dataDir", "keepDroplets") ) )]
      soup <- eval(rlang::expr(SoupX::load10X(
        dataDir = paste0(data_dir, sample_meta$path[i], "/", tenx_dir),
        keepDroplets = TRUE, !!!params)))
      # Estimate contamination
      params <- dot_params[which(names(dot_params) %in% methods::formalArgs(SoupX::autoEstCont))]
      params <- params[which(!(names(params) %in% c("sc", "forceAccept", "doPlot") ) )]
      if (!is.null(plot_dir)) {
        png(paste0(plot_dir, "soupX_", s, ".png"), width = 7, height = 5,
            res = 150, units = "in")
        soup <- eval(rlang::expr(SoupX::autoEstCont(sc = soup, forceAccept = TRUE, doPlot = TRUE, !!!params)))
        dev.off()
      } else {
        soup <- eval(rlang::expr(SoupX::autoEstCont(sc = soup, forceAccept = TRUE, doPlot = FALSE, !!!params)))
      }
      # Adjust counts
      params <- dot_params[which(names(dot_params) %in% methods::formalArgs(SoupX::autoEstCont))]
      params <- params[which(!(names(params) %in% c("sc") ) )]
      gex <- eval(rlang::expr(SoupX::adjustCounts(sc = soup, !!!params)))
    } else {
      # No ambient RNA correction, reduction counts matrix
      gex <- eval(rlang::expr(Seurat::GetAssayData(object = tmp, slot = "counts")))
    }

    ##
    # Create final Seurat object and add metadata information
    params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat::CreateSeuratObject))]
    params <- params[which(!(names(params) %in% c("counts", "project", "min.cells", "min.features") ) )]
    seu[[s]] <- eval(rlang::expr(Seurat::CreateSeuratObject(
      counts = gex, project = s, min.cells = min.cells, min.features = min.features, !!!params)))
    # Rename cells to distinguish across datasets
    seu[[s]] <- Seurat::RenameCells(object = seu[[s]], add.cell.id = s)
    # Compute Mitochondrial percentage
    # TODO: Hack assuming either human OR mouse, and since for one of them we
    # will have all 0s, OR (|) will give us the same information we need.
    seu[[s]]$percent.mito <- Seurat::PercentageFeatureSet(seu[[s]], pattern = "^MT-|^mt-")

    # Add required sample metadata information
    seu[[s]]$sample <- as.factor(s)
    seu[[s]]$donor <- as.factor(sample_meta$donor[i])
    seu[[s]]$condition <- as.factor(sample_meta$condition[i])

    # Add user defined sample metadata information
    for (m in meta_colnames) {
      seu[[s]][[m]] <- as.factor(sample_meta[[m]][i])
    }

    if (use_scrublet) {
      # Run scrublet to quantify doublets
      params <- dot_params[which(names(dot_params) %in% c("sim_doublet_ratio", "n_neighbors",
                                                          "stdev_doublet_rate", "random_state"))]
      scr <- eval(rlang::expr(scrub$Scrublet(
        counts_matrix = reticulate::r_to_py(Seurat::GetAssayData(object = seu[[s]], slot = "counts"))$T$tocsc(),
        expected_doublet_rate = expected_doublet_rate, !!!params)))

      params <- dot_params[which(names(dot_params) %in% c(
        "synthetic_doublet_umi_subsampling", "use_approx_neighbors", "distance_metric",
        "get_doublet_neighbor_parents", "min_counts", "min_cells", "min_gene_variability_pctl",
        "log_transform", "mean_center", "normalize_variance", "n_prin_comps"))]
      doublet_results <- eval(rlang::expr(reticulate::py_to_r(scr$scrub_doublets(!!!params))))
      doublet_score <- doublet_results[[1]]
      names(doublet_score) <- colnames(seu[[s]])
      doublet_prediction <- doublet_results[[2]]
      names(doublet_prediction) <- colnames(seu[[s]])
      # Add scrublet output
      seu[[s]] <- Seurat::AddMetaData(seu[[s]], metadata = doublet_score, col.name = "doublet_score")
      seu[[s]] <- Seurat::AddMetaData(seu[[s]], metadata = doublet_prediction, col.name = "doublet_prediction")
    }
  }
  # We have a single sample, so removing from list
  if (length(seu) == 1) { seu <- seu[[1]] }
  # Return generated Seurat object or list of Seurat objects
  return(seu)
}


#' @name qc_filter_seurat_object
#' @rdname qc_filter_seurat_object
#'
#' @title Filter Seurat object based on QC metrics.
#'
#' @description This function filters a Seurat object (i.e. removes cells)
#' based on QC metric thresholds. For now only based on nFeatures and
#' mitochondrial percentage.
#'
#' @param seu Seurat object or list of Seurat objects(required).
#'
#' @inheritParams run_qc_pipeline
#'
#' @return List of QC filtered Seurat objects. If a single sample, returns
#' a QC filtered Seurat object instead of a list.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
NULL

#' @rdname qc_filter_seurat_object
#'
#' @export
qc_filter_seurat_object <- function(seu, nfeat_thresh, mito_thresh){
  UseMethod("qc_filter_seurat_object")
}

# Default function for the generic function 'qc_filter_seurat_object'
qc_filter_seurat_object.default <- function(seu, nfeat_thresh,
                                            mito_thresh){
  stop("Object 'seu' should be either list or Seurat object!")
}


#' @rdname qc_filter_seurat_object
#'
#' @export
qc_filter_seurat_object.list <- function(seu, nfeat_thresh, mito_thresh) {
  # Perform actual QC filtering
  for (s in names(seu)) {
    seu[[s]] <- qc_filter_seurat_object.Seurat(
      seu = seu[[s]], nfeat_thresh=nfeat_thresh, mito_thresh=mito_thresh)
  }
  return(seu)
}

#' @rdname qc_filter_seurat_object
#'
#' @export
qc_filter_seurat_object.Seurat <- function(seu, nfeat_thresh, mito_thresh) {
  assay <- Seurat::DefaultAssay(object = seu)
  # Perform actual QC filtering
  cells <- seu@meta.data[paste0("nFeature_", assay)] > nfeat_thresh &
    seu@meta.data["percent.mito"] < mito_thresh
  cells <- rownames(cells[which(cells == TRUE), , drop = FALSE])
  seu <- subset(seu, cells = cells)
}



#' @name lognormalize_and_pca
#' @rdname lognormalize_and_pca
#'
#' @title Log normalisation and PCA computation
#'
#' @description This function wrapps the default Seurat analysis of 1.
#' data normalisation, 2. identification of HVGs, 3. scaling expression
#' and computing PCA. If a list of Seurat objects is provided
#' this is performed for each element in the list.
#'
#' @param seu Seurat object or list of Seurat objects(required).
#' @param npcs Number of Principal Components (PCs) to compute.
#' If NULL, PCA analysis is skipped.
#' @param n_hvgs Number of highly variable genes (HVGs) to compute. Only HVGs
#' will be used as input to PCA.
#' @param ... Additional named parameters passed to Seurat functions.
#'
#' @return List of processed Seurat objects. If a single sample, returns
#' a processed Seurat object instead of a list.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
NULL

#' @rdname lognormalize_and_pca
#'
#' @export
lognormalize_and_pca <- function(seu, npcs = NULL, n_hvgs = 3000, ...){
  UseMethod("lognormalize_and_pca")
}

# Default function for the generic function 'lognormalize_and_pca'
lognormalize_and_pca.default <- function(seu, npcs = NULL, n_hvgs = 3000, ...){
  stop("Object 'seu' should be either list or Seurat object!")
}

#' @rdname lognormalize_and_pca
#'
#' @export
lognormalize_and_pca.list <- function(seu, npcs = NULL, n_hvgs = 3000, ...) {
  for (s in names(seu)) {
    seu[[s]] <- lognormalize_and_pca.Seurat(seu = seu[[s]], npcs = npcs,
                                       n_hvgs = n_hvgs, ...)
  }
  return(seu)
}


#' @rdname lognormalize_and_pca
#'
#' @export
lognormalize_and_pca.Seurat <- function(seu, npcs = NULL, n_hvgs = 3000, ...) {
  # All parameters passed to ...
  dot_params <- rlang::list2(...)

  params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat:::NormalizeData.Seurat))]
  seu <- eval(rlang::expr(Seurat::NormalizeData(seu, !!!params)))

  params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat:::FindVariableFeatures.Seurat))]
  params <- params[which(!(names(params) %in% c("nfeatures") ) )]
  seu <- eval(rlang::expr(Seurat::FindVariableFeatures(object = seu, nfeatures = n_hvgs, !!!params)))
  # PCA reduction
  if (!is.null(npcs)) {
    params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat:::ScaleData.Seurat))]
    seu <- eval(rlang::expr(Seurat::ScaleData(object = seu, !!!params)))

    params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat:::RunPCA.Seurat))]
    params <- params[which(!(names(params) %in% c("npcs") ) )]
    seu <- eval(rlang::expr(Seurat::RunPCA(object = seu, npcs = npcs, !!!params)))
  }
  return(seu)
}


#' @name run_umap
#' @rdname run_umap
#'
#' @title Tailored UMAP function
#'
#' @description This function adapts the Seurat RunUMAP function by also
#' checking if 'dims' is larger than dimensions in 'reduction' object.
#' If yes, it uses the maximum columns of 'reduction' object.
#'
#' @param seu Seurat object (required).
#' @param dims Vector denoting which dimensions to use as input
#' features, e.g. dims = 1:30
#' @param reduction Reduction object to run UMAP. E.g. can be 'pca'
#' or 'harmony' for integrated data (required).
#' @param seed Set a random seed, for reproducibility.
#' @param ... Additional named parameters passed to Seurat RunUMAP function.
#'
#' @return Seurat object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
run_umap <- function(seu, dims, reduction, seed = 1, ...) {
  # All parameters passed to ...
  dot_params <- rlang::list2(...)

  # Number of dimensions to perform UMAP
  dims <- dims[dims <= NCOL(seu@reductions[[reduction]])]
  # Run UMAP
  params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat:::RunUMAP.Seurat))]
  params <- params[which(!(names(params) %in% c("dims", "reduction", "seed.use") ) )]
  seu <- eval(rlang::expr(Seurat::RunUMAP(seu, dims = dims, seed.use = seed,
                                          reduction = reduction, !!!params)))
  return(seu)
}


#' @name find_neighbors
#' @rdname find_neighbors
#' @title Tailored function for finding neighbors in lower dimensional space.
#'
#' @description This function adapts the Seurat FindNeighbors function by also
#' checking if 'dims' is larger than dimensions in 'reduction' object.
#' If yes, it takes maximum columns of 'reduction' object.
#'
#' @param seu Seurat object (required).
#' @param dims Vector denoting which dimensions to use as input
#' features, e.g. dims = 1:30
#' @param reduction Reduction object to run UMAP. E.g. can be 'pca'
#' or 'harmony' for integrated data (required).
#' @param ... Additional named parameters to Seurat FindNeighbors function.
#'
#' @return Seurat object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
find_neighbors <- function(seu, dims, reduction, ...) {
  assertthat::assert_that(methods::is(reduction, "character"))
  # All parameters passed to ...
  dot_params <- rlang::list2(...)
  # Number of dimensions to find neighbors
  dims <- dims[dims <= NCOL(seu@reductions[[reduction]])]
  # Identify neighbors
  params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat:::FindNeighbors.Seurat))]
  params <- params[which(!(names(params) %in% c("dims", "reduction")) )]
  seu <- eval(rlang::expr(Seurat::FindNeighbors(seu, dims = dims, reduction = reduction, !!!params)))
  return(seu)
}


#' @name find_clusters
#' @rdname find_clusters
#'
#' @title Wrapper FindClusters function
#'
#' @description This function just calls the Seurat FindClusters function.
#' Used only for naming consistency in this package.
#'
#' @param seu Seurat object (required).
#' @param resolution Value of the resolution parameter, use a value
#' above (below) 1.0 if you want to obtain a larger (smaller) number of
#' communities.
#' @param random.seed Seed to use for reproducibility purposes.
#' @param ... Additional named parameters to Seurat FindClusters function.
#'
#' @return Seurat object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
find_clusters <- function(seu, resolution, random.seed = 1, ...) {
  # All parameters passed to ...
  dot_params <- rlang::list2(...)
  # Identify neighbors
  params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat:::FindClusters.Seurat))]
  params <- params[which(!(names(params) %in% c("resolution", "random.seed")) )]
  seu <- eval(rlang::expr(Seurat::FindClusters(seu, resolution = resolution,
                                               random.seed = random.seed, !!!params)))
  return(seu)
}


#' @name find_all_markers
#' @rdname find_all_markers
#'
#' @title Wrapper FindAllMarkers function
#'
#' @description This function just calls the Seurat FindAllMarkers function.
#' Used only for naming consistency in this package.
#'
#' @param seu Seurat object (required).
#' @param random.seed Seed to use for reproducibility purposes.
#' @param ... Additional named parameters to Seurat FindAllMarkers function.
#'
#' @return Data frame with marker genes for each cluster and associated
#' statistics.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
find_all_markers <- function(seu, random.seed = 1, ...) {
  # All parameters passed to ...
  dot_params <- rlang::list2(...)
  # Identify neighbors
  params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat::FindAllMarkers))]
  params <- params[which(!(names(params) %in% c("random.seed")) )]
  markers <- eval(rlang::expr(Seurat::FindAllMarkers(seu, random.seed = 1, !!!params)))
  return(markers)
}


#' @title Tailored module score calculation
#'
#' @description This function adapts the AddModuleScore Seurat function to
#' compute module scores based on a group of genes of interest.
#'
#' @param features Vector of features (e.g. genes) to compute module score.
#' @param name Name of module score to be stored in metadata slot.
#' @param ... Additional parameters passed to Seurat's AddModuleScore.
#' @inheritParams cluster_analysis
#'
#' @return A ggplot2 object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
compute_module_score <- function(seu, features, name, seed = 1, ...) {
  # All parameters passed to ...
  dot_params <- rlang::list2(...)
  params <- dot_params[which(names(dot_params) %in% methods::formalArgs(Seurat::AddModuleScore))]
  params <- params[which(!(names(params) %in% c("features", "name", "seed") ) )]

  # Required for the really low quality samples so they pass without errors
  if (any(names(dot_params) %in% "ctrl")) {
    if (NROW(seu) <= dot_params[[which(names(dot_params) %in% "ctrl")]]) {
      dot_params[[which(names(dot_params) %in% "ctrl")]] <- 20
    }
  }

  features <- list(markers = features[features %in% rownames(seu)])[1]
  if (length(features$markers) == 0) {
    message("No features to compute module score.")
    return(seu)
  }
  seu <- eval(rlang::expr(Seurat::AddModuleScore(object = seu, features = features,
                                                 name = name, seed = seed, !!!params)))
  seu[[name]] <- seu[[paste0(name, "1")]]
  seu[[paste0(name, "1")]] <- NULL
  return(seu)
}

#' @title Local implementation of RunHarmony function
#'
#' @description This is a local implementation of RunHarmony function
#' (from `harmony` package) which solves the bug where the 'dims.use'
#' parameter is ignored.
#'
#' @param object Seurat object.
#' @param group.by.vars Which variable(s) to remove (character vector).
#' @param reduction Dimensionality reduction to use.
#' @param dims.use Which PCA dimensions to use for Harmony. By default, use all
#' @param theta Diversity clustering penalty parameter. Specify for each
#' variable in group.by.vars. Default theta=2. theta=0 does not encourage any
#'  diversity. Larger values of theta result in more diverse clusters.
#' @param lambda Ridge regression penalty parameter. Specify for each variable
#' in group.by.vars. Default lambda=1. Lambda must be strictly positive.
#' Smaller values result in more aggressive correction.
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
#' the distance from a cell to cluster centroids. Larger values of sigma result
#'  in cells assigned to more clusters. Smaller values of sigma make soft
#'  kmeans cluster approach hard clustering.
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple
#'  linear regression.
#' @param tau Protection against overclustering small datasets with large ones.
#'  tau is the expected number of cells per cluster.
#' @param block.size What proportion of cells to update during clustering.
#'  Between 0 to 1, default 0.05. Larger values may be faster but less accurate
#' @param max.iter.cluster Maximum number of rounds to run clustering at each
#' round of Harmony.
#' @param epsilon.cluster Convergence tolerance for clustering round of Harmony
#'  Set to -Inf to never stop early.
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round
#'  of Harmony involves one clustering and one correction step.
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to
#'  never stop early.
#' @param plot_convergence Whether to print the convergence plot of the
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging.
#' @param verbose Whether to print progress messages. TRUE to print, FALSE to
#'  suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s). Cells
#' that have batch variables values matching reference_values will not be moved
#' @param reduction.save Keyword to save Harmony reduction. Useful if you want
#' to try Harmony with multiple parameters and save them as e.g.
#' 'harmony_theta0', 'harmony_theta1', 'harmony_theta2'
#' @param assay.use (Seurat V3 only) Which assay to Harmonize with
#' (Default Assay as default).
#' @param project.dim Project dimension reduction loadings. Default TRUE.
#' @param ... other parameters
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
run_harmony <- function(object, group.by.vars, reduction = 'pca',
                        dims.use = NULL, theta = NULL, lambda = NULL,
                        sigma = 0.1, nclust = NULL, tau = 0, block.size = 0.05,
                        max.iter.harmony = 10, max.iter.cluster = 20,
                        epsilon.cluster = 1e-5, epsilon.harmony = 1e-4,
                        plot_convergence = FALSE, verbose = TRUE,
                        reference_values = NULL, reduction.save = "harmony",
                        assay.use = NULL, project.dim = TRUE, ...) {
  assay.use <- assay.use %||% Seurat::DefaultAssay(object)
  if (reduction == "pca") {
    tryCatch(
      embedding <- Seurat::Embeddings(object, reduction = "pca"),
      error = function(e) {
        if (verbose) { message("Harmony needs PCA. Trying to run PCA now.") }
        tryCatch(
          object <- Seurat::RunPCA(object, assay = assay.use, verbose = verbose),
          error = function(e) {
            stop("Harmony needs PCA. Tried to run PCA and failed.")
          })
      })
  } else {
    available.dimreduc <- names(methods::slot(object = object,
                                              name = "reductions"))
    if (!(reduction %in% available.dimreduc)) {
      stop("Requested dimension reduction is not present in the Seurat object")
    }
    embedding <- Seurat::Embeddings(object, reduction = reduction)
  }
  if (is.null(dims.use)) { dims.use <- seq_len(ncol(embedding)) }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rerun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  # Fix Andreas: Subset by dims.use
  embedding <- embedding[, dims.use, drop = FALSE]
  # Extract grouping variable
  metavars_df <- Seurat::FetchData(object, group.by.vars)

  harmonyEmbed <- harmony::HarmonyMatrix(embedding, metavars_df, group.by.vars,
                                         FALSE, 0, theta, lambda, sigma, nclust,
                                         tau, block.size, max.iter.harmony,
                                         max.iter.cluster, epsilon.cluster,
                                         epsilon.harmony, plot_convergence,
                                         FALSE, verbose, reference_values
  )
  rownames(harmonyEmbed) <- row.names(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))

  suppressWarnings({
    harmonydata <- Seurat::CreateDimReducObject(
      embeddings = harmonyEmbed,
      stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
      assay = assay.use, key = reduction.save)
  })

  object[[reduction.save]] <- harmonydata
  if (project.dim) {
    object <- Seurat::ProjectDim(object, reduction = reduction.save,
                                 overwrite = TRUE, verbose = FALSE)
  }
  return(object)
}


#' Add UMAP embedding in Seurat object
#'
#' Add UMAP embedding in existing Seurat object. This is the case when running
#' the pipeline and then you want to attach the stored UMAP embedding, stored
#' as CSV file, instead of recomputing the UMAP.
#'
#' @param seu Seurat object
#' @param embedding UMAP embedding. If a character, the function will assume a
#' filename is given and will read the corresponding UMAP embedding. Otherwise,
#' it assumes a UMAP embedding matrix is provided as input containing only the
#' embeddings as columns.
#'
#' @export
add_umap_embedding <- function(seu, embedding) {
  # If character assume we are given a file name
  if (is.character(embedding)) {
    embedding <- read.csv(file = embedding) |>
      tibble::column_to_rownames(var = "X") |> as.matrix()
  }
  if (NROW(embedding) != NCOL(seu)) {
    stop("Embedding matrix not equal to number of cells in Seurat object.")
  }
  # Create DimReducObject
  seu@reductions$umap <- Seurat::CreateDimReducObject(
    embeddings = embedding, assay = Seurat::DefaultAssay(object = seu), key = "UMAP_")
  return(seu)
}


#' Install Scrublet Python Package
#'
#' Install Scrublet Python package into a virtualenv or conda env.
#'
#' On Linux and OS X the "virtualenv" method will be used by default
#' ("conda" will be used if virtualenv isn't available). On Windows,
#' the "conda" method is always used.
#'
#' @param envname Name of environment to install packages into
#' @param method Installation method. By default, "auto" automatically finds
#' a method that will work in the local environment. Change the default to
#' force a specific installation method. Note that the "virtualenv" method
#' is not available on Windows.
#' @param conda Path to conda executable (or "auto" to find conda using the PATH
#'  and other conventional install locations).
#' @param pip Install from pip, if possible.
#' @param ... Additional arguments passed to conda_install() or
#' virtualenv_install().
#'
#' @export
install_scrublet <- function(envname = "r-reticulate", method = "auto",
                             conda = "auto", pip = TRUE, ...) {
  # Check: https://github.com/KrishnaswamyLab/phateR/blob/master/R/utils.R
  tryCatch({
    message("Attempting to install Scrublet Python package with reticulate")
    reticulate::py_install("scrublet",
                           envname = envname, method = method,
                           conda = conda, pip = pip, ...
    )
    message("Install complete. Please restart R and try again.")
  },
  error = function(e) {
    stop(paste0(
      "Cannot locate Scrublet Python package, please install through pip ",
      "(e.g. ", reticulate::py_config()$python,
      " -m pip install --user scrublet) and then restart R."
    ))
  }
  )
}

# Helper function for checking metadata data.frame or filename and returning
# the object
.get_metadata <- function(sample_meta = NULL, sample_meta_filename = NULL) {
  pass_qc <- NULL
  # Check if metadata object or filename is given as input
  if (!is.null(sample_meta) && !is.null(sample_meta_filename)) {
    message("Both 'sample_meta' and 'sample_meta_filename' provided.
            Using 'sample_meta'.")
  } else if (is.null(sample_meta) && !is.null(sample_meta_filename)) {
    sample_meta <- read.csv(file = paste0(sample_meta_filename)) |>
      dplyr::filter(pass_qc == TRUE)
  }
  if (is.null(sample_meta) && is.null(sample_meta_filename)) {
    stop("One of 'sample_meta' or 'sample_meta_filename' should be provided.")
  }

  # Check that required colnames exist in metadata file
  if (!all(c("sample", "donor", "path", "condition", "pass_qc")
           %in% colnames(sample_meta))) {
    stop("Required columns 'sample', 'donor', 'path', 'condition', 'pass_qc'
         are not present in metadata file. Stopping.")
  }
  if (NROW(sample_meta) != length(unique(sample_meta$sample))) {
    stop("Stopping 'sample' ids are not unique.")
  }

  # Order metadata by sample name
  idx <- gtools::mixedorder(sample_meta$sample)
  sample_meta <- sample_meta[idx, ]

  # Return sample metadata
  return(sample_meta)
}

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
.get_density <- function(x, y, ...) {
  # All parameters passed to ...
  dot_params <- rlang::list2(...)

  params <- dot_params[which(names(dot_params) %in% methods::formalArgs(MASS::kde2d))]
  dens <- eval(rlang::expr(MASS::kde2d(x, y, !!!params)))

  ix <- base::findInterval(x, dens$x)
  iy <- base::findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# General function for converting all columns that are string to Factors, due
# to a bug when merging Seurat objects that they get converted to characters.
.as_factor_metadata <- function(seu) {
  for (meta_ids in colnames(seu@meta.data)) {
    meta_vals <- seu@meta.data[[meta_ids]]
    if (is.character(meta_vals) || is.logical(meta_vals)) {
      meta_levels <- gtools::mixedsort(unique(meta_vals))
      seu@meta.data[[meta_ids]] <- factor(meta_vals, levels = meta_levels)
    }
  }
  return(seu)
}

# General function for dropping unwanted factors in a (filterred) Seurat dataset.
.drop_factors <- function(seu) {
  for (meta_ids in colnames(seu@meta.data)) {
    meta_vals <- seu@meta.data[[meta_ids]]
    if (is.factor(meta_vals)) {
      seu@meta.data[[meta_ids]] <- droplevels(seu@meta.data[[meta_ids]])
    }
  }
  return(seu)
}

# General function for obtaining consistent plotting dimensions
.plot_dims <- function(feat_len) {
  if (!is.numeric(feat_len)) {
    message("Invalid argument:  Setting 'feat_len' for plotting to 2")
    feat_len <- 2
  } else if (feat_len == 0){
    feat_len <- 2
  }
  if (feat_len == 1) {
    width = 9.5
    height = 7
    ncols = feat_len
  } else if (feat_len <= 3 & feat_len > 1) {
    width = 7.5 * feat_len
    height = 6
    ncols = feat_len
  } else {
    width = 30
    height = 6 * ceiling(feat_len / 4)
    ncols = 4
  }
  return(list(width = width, height = height, ncols = ncols))
}


# General function for obtaining consistent spatial plotting dimensions
.spatial_plot_dims <- function(feat_len, sample_len = 1, technology = "VisiumV1") {
  if (!is.numeric(feat_len)) {
    message("Invalid argument:  Setting 'feat_len' for plotting to 2")
    feat_len <- 2
  } else if (feat_len == 0){
    feat_len <- 2
  }
  if (sample_len == 1) {
    if (feat_len == 1) {
      width <- 6
      height <- 7
      ncols <- 1
    } else if (feat_len <= 3 & feat_len > 1) {
      width <- 6 * feat_len
      height <- 7
      ncols <- feat_len
    } else {
      width <- 24
      height <- 7 * ceiling(feat_len / 4)
      ncols <- 4
    }
  } else {
    if (technology == "Visium") {
      width <- 6 * sample_len
      height <- 7 * feat_len
      ncols <- 1
    } else if (technology == "SlideSeq") {
      width <- 6 * sample_len
      height <- 7 * feat_len
      ncols <- sample_len
    } else {
      width <- 6 * sample_len
      height <- 7 * feat_len
      ncols <- 1
    }
  }
  return(list(width = width, height = height, ncols = ncols))
}


.internal_col_pal <- function() {
  # https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
  return(c(
    "indianred", # red
    "#6699CB",
    "#FDBF6F", # lt orange
    "#CAB2D6", # lt purple
    "#FB9A99", # lt pink
    "tan3",
    "darkolivegreen4", # darkgreen
    "darkgrey", # darkgrey
    "skyblue2", # lightblue
    "mediumpurple1",
    "darkseagreen3",
    "khaki2", "ivory3",
    "steelblue4",
    "#6A3D9A", # purple
    "seagreen", "orchid1", "blue1", "deeppink1", "gold1",
    "darkturquoise", "darkorange4", "#FF7F00",
    "dodgerblue", "yellow3", "mediumorchid1", "firebrick4",
    "wheat4", "maroon", "grey30", "red2", "burlywood2", "cyan",
    "darkolivegreen2", "yellowgreen"
  ))
}
