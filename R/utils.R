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
    tmp <-  Seurat::Read10X(
      data.dir = paste0(data_dir, sample_meta$path[i], "/", tenx_dir, "/",
                        tenx_counts_dir), ...)
    tmp <- Seurat::CreateSeuratObject(counts = tmp, project = s,
                                      min.cells = 0, min.features = 0, ...)
    if (use_soupx) {
      # Run SoupX for ambient RNA correction
      soup <- SoupX::load10X(dataDir = paste0(data_dir, sample_meta$path[i],
                                              "/", tenx_dir),
                             keepDroplets = TRUE, ...)
      # Estimate contamination
      if (!is.null(plot_dir)) {
        png(paste0(plot_dir, "soupX_", s, ".png"), width = 7, height = 5,
            res = 150, units = "in")
        soup <- SoupX::autoEstCont(sc = soup, forceAccept = TRUE, ...)
        dev.off()
      } else {
        soup <- SoupX::autoEstCont(sc = soup, forceAccept = TRUE,
                                   doPlot = FALSE, ...)
      }
      # Adjust counts
      gex <- SoupX::adjustCounts(sc = soup, ...)
    } else {
      # No ambient RNA correction, reduction counts matrix
      gex <- Seurat::GetAssayData(object = tmp, slot = "counts")
    }

    ##
    # Create final Seurat object and add metadata information
    seu[[s]] <- Seurat::CreateSeuratObject(
      counts = gex, project = s, min.cells = min.cells,
      min.features = min.features, ...)
    # Rename cells to distinguish across datasets
    seu[[s]] <- Seurat::RenameCells(object = seu[[s]], add.cell.id = s)
    # Compute Mitochondrial percentage
    # TODO: Hack assuming either human OR mouse, and since for one of them we
    # will have all 0s, OR (|) will give us the same information we need.
    seu[[s]]$percent.mito <- Seurat::PercentageFeatureSet(seu[[s]],
                                                          pattern = "^MT-|^mt-")

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
      scr <- scrub$Scrublet(reticulate::r_to_py(
        Seurat::GetAssayData(object = seu[[s]], slot = "counts"))$T$tocsc(),
        expected_doublet_rate = expected_doublet_rate)
      doublet_results <- reticulate::py_to_r(scr$scrub_doublets(...))
      doublet_score <- doublet_results[[1]]
      names(doublet_score) <- colnames(seu[[s]])
      doublet_prediction <- doublet_results[[2]]
      names(doublet_prediction) <- colnames(seu[[s]])
      # Add scrublet output
      seu[[s]] <- Seurat::AddMetaData(seu[[s]], metadata = doublet_score,
                                      col.name = "doublet_score")
      seu[[s]] <- Seurat::AddMetaData(seu[[s]], metadata = doublet_prediction,
                                      col.name = "doublet_prediction")
    }
  }
  # We have a single sample, so removing from list
  if (length(seu) == 1) { seu <- seu[[1]] }
  # Return generated Seurat object or list of Seurat objects
  return(seu)
}


#' @rdname spatial_create_seurat_object
#'
#' @title Create Seurat object based on spatial sample metadata.
#'
#' @description This function creates Seurat object for each spatial
#' (e.g. Visium) sample in the metadata file.
#'
#' @param ... Additional named parameters passed to Seurat.
#' @inheritParams run_qc_pipeline
#'
#' @return List of Seurat objects equal the number of samples in
#' the sample metadata file. If a single sample, returns a Seurat object.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
spatial_create_seurat_object <- function(
    data_dir, sample_meta = NULL, sample_meta_filename = NULL,
    meta_colnames = c("donor", "condition", "pass_qc"), tenx_dir = "outs", ...) {
  pass_qc = sample = Barcode <- NULL
  # Check if metadata object or filename is given as input
  sample_meta <- .get_metadata(sample_meta = sample_meta,
                               sample_meta_filename = sample_meta_filename)

  # Create Seurat object for each sample
  seu <- list()
  for (i in 1:NROW(sample_meta)) {
    s <- sample_meta$sample[i]

    seu[[s]] <- Seurat::Load10X_Spatial(
      data.dir = paste0(data_dir, sample_meta$path[i], "/", tenx_dir),
      slice = s, filter.matrix = TRUE, ...)
    Seurat::Idents(seu[[s]]) <- s
    seu[[s]]@meta.data$orig.ident <- Seurat::Idents(seu[[s]])

    # Rename cells to distinguish across datasets
    seu[[s]] <- Seurat::RenameCells(object = seu[[s]], add.cell.id = s)

    # If manual filtering of spots occured, e.g. through Loupe browser,
    # read file and subset Seurat object accordingly
    fname <- paste0(data_dir, sample_meta$path[i], "/", tenx_dir, "/", s, "_filtered.csv")
    if (file.exists(fname)) {
      df <- read.csv(fname) |>
        dplyr::filter(pass_qc == "pass") |>
        dplyr::mutate(Barcode = paste0(s, "_", Barcode))
      seu[[s]] <- subset(seu[[s]], cells = df$Barcode)
    }

    # Compute Mitochondrial percentage
    # TODO: Hack assuming either human OR mouse, and since for one of them we
    # will have all 0s, OR (|) will give us the same information we need.
    seu[[s]]$percent.mito <- Seurat::PercentageFeatureSet(
      seu[[s]], pattern = "^MT-|^mt-")

    # Add required sample metadata information
    seu[[s]]$sample <- as.factor(s)
    seu[[s]]$donor <- as.factor(sample_meta$donor[i])
    seu[[s]]$condition <- as.factor(sample_meta$condition[i])

    # Add user defined sample metadata information
    for (m in meta_colnames) {
      seu[[s]][[m]] <- as.factor(sample_meta[[m]][i])
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
#' @param ... Additional parameters.
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
qc_filter_seurat_object <- function(seu, nfeat_thresh, mito_thresh, ...){
  UseMethod("qc_filter_seurat_object")
}

# Default function for the generic function 'qc_filter_seurat_object'
qc_filter_seurat_object.default <- function(seu, nfeat_thresh,
                                            mito_thresh, ...){
  stop("Object 'seu' should be either list or Seurat object!")
}


#' @rdname qc_filter_seurat_object
#'
#' @export
qc_filter_seurat_object.list <- function(seu, nfeat_thresh, mito_thresh, ...) {
  # Perform actual QC filtering
  for (s in names(seu)) {
    seu[[s]] <- qc_filter_seurat_object.Seurat(
      seu = seu[[s]], nfeat_thresh=nfeat_thresh, mito_thresh=mito_thresh, ...)
  }
  return(seu)
}

#' @rdname qc_filter_seurat_object
#'
#' @export
qc_filter_seurat_object.Seurat <- function(seu, nfeat_thresh, mito_thresh, ...) {
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
  seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize",
                               scale.factor = 10000)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst",
                                      nfeatures = n_hvgs, ...)
  # PCA reduction
  if (!is.null(npcs)) {
    seu <- Seurat::ScaleData(seu, ...)
    seu <- Seurat::RunPCA(seu, features = Seurat::VariableFeatures(seu),
                          npcs = npcs, ...)
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
  # Number of dimensions to perform UMAP
  dims <- dims[dims <= NCOL(seu@reductions[[reduction]])]
  # Run UMAP
  seu <- Seurat::RunUMAP(seu, dims = dims, seed.use = seed,
                         reduction = reduction, ...)
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
  # Number of dimensions to find neighbors
  dims <- dims[dims <= NCOL(seu@reductions[[reduction]])]
  # Identify neighbors
  seu <- Seurat::FindNeighbors(seu, dims = dims, reduction = reduction, ...)
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
  # Identify neighbors
  seu <- Seurat::FindClusters(seu, resolution = resolution,
                              random.seed = random.seed, ...)
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
  # Identify neighbors
  markers <- Seurat::FindAllMarkers(seu, random.seed = 1, ...)
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
compute_module_score <- function(seu, features, name, ctrl = 100,
                                 seed = 1, ...) {
  # Required for the really low quality samples so they pass without errors
  ctrl <- ifelse(NCOL(seu) <= ctrl, NCOL(seu) - 5, ctrl)
  features <- list(markers = features[features %in% rownames(seu)])[1]
  if (length(features$markers) == 0) {
    message("No features to compute module score.")
    return(seu)
  }
  seu <- Seurat::AddModuleScore(object = seu, features = features,
                                ctrl = ctrl, name = name, seed = seed, ...)
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
  # Order metadata by sample name
  sample_meta <- sample_meta |> dplyr::arrange(sample)

  # Return sample metadata
  return(sample_meta)
}

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
.get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
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
    if (is.character(meta_vals)) {
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
.spatial_plot_dims <- function(feat_len) {
  if (!is.numeric(feat_len)) {
    message("Invalid argument:  Setting 'feat_len' for plotting to 2")
    feat_len <- 2
  } else if (feat_len == 0){
    feat_len <- 2
  }
  if (feat_len == 1) {
    width = 6
    height = 7
    ncols = feat_len
  } else if (feat_len <= 3 & feat_len > 1) {
    width = 6 * feat_len
    height = 7
    ncols = feat_len
  } else {
    width = 24
    height = 7 * ceiling(feat_len / 4)
    ncols = 4
  }
  return(list(width = width, height = height, ncols = ncols))
}



.internal_col_pal <- function() {
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
