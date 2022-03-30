
#' @title Tailored module score calculation
#'
#' @description This function adapts the AddModuleScore Seurat function.
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
compute_module_score <- function(seu, features, name, ctrl = 100, seed = 1, ...) {
  ctrl <- ifelse(NCOL(seu) <= ctrl, NCOL(seu) - 10, ctrl)
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
      embedding = Seurat::Embeddings(object, reduction = "pca"),
      error = function(e) {
        if (verbose) { message("Harmony needs PCA. Trying to run PCA now.") }
        tryCatch(
          object = Seurat::RunPCA(object, assay = assay.use, verbose = verbose),
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
  embedding <- embedding[, dims.use]
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
.as_factor_metadata_seurat <- function(seu) {
  for (meta_ids in colnames(seu@meta.data)) {
    meta_vals <- seu@meta.data[[meta_ids]]
    if (is.character(meta_vals)) {
      meta_levels <- gtools::mixedsort(unique(meta_vals))
      seu@meta.data[[meta_ids]] <- factor(meta_vals, levels = meta_levels)
    }
  }
  return(seu)
}

# General function for obtaining consistent plotting dimensions
.plot_dims <- function(feat_len) {
  if (feat_len <= 3) {
    width = 6 * feat_len
    height = 6
    ncols = feat_len
  } else {
    width = 24
    height = 6 * ceiling(feat_len / 4)
    ncols = 4
  }
  return(list(width = width, height = height, ncols = ncols))
}
