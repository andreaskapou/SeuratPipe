
#' @rdname run_harmony_pipeline
#'
#' @title (Internal) pipeline for Harmony integration
#'
#' @description This is an internal pipeline which wraps all analysis steps in
#' a function. This is meant to be used as a guide and in each analysis, users
#' should write their own function (possibly adapting the code of this
#' function). This is the reason this function is not explicitly 'exported',
#' i.e. even if you load the library, you cannot directly call the function,
#' but you should call seurat_pipeline:::run_harmony_pipeline. Elements in 'io'
#' and 'opts' are not checked, so function might return errors if not all
#' required values are not defined.
#'
#' @param seu_obj Seurat object (required).
#' @param io Input/output settings defined in main function.
#' @param opts Options settings defined in main function.
#' @param ... Additional named parameters passed to Seurat's or Harmony
#' functions.
#'
#' @return None.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
run_harmony_pipeline <- function(seu_obj, io, opts, ...) {
  pcs_remove_name <- ""
  # Iterate over the number of principal components
  for (npcs in opts$npcs) {
    # Define number of dimensions to use
    dims.use <- seq(from = 1, to = npcs, by = 1)
    if (!is.null(opts$pcs_to_remove)) {
      pcs_remove_name <- "_manual"
      dims.use <- dims.use[-opts$pcs_to_remove]
    }

    obj_name <- paste0(opts$int_method, "_npcs", npcs, pcs_remove_name)
    out_dir <- paste0(io$out_dir, "/", obj_name, "/")
    if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}

    # Preprocessing, integration and QC
    qc_dir <- paste0(out_dir, "/qc/")
    dir.create(qc_dir, recursive = TRUE)
    if (!file.exists(paste0(io$out_dir, "/int_", obj_name, ".rds"))) {
      seu <- harmony_analysis(seu = seu_obj, npcs = npcs, dims.use = dims.use,
                              plot_dir = qc_dir, n_hvgs = opts$n_hvgs,
                              max.iter.harmony = opts$max.iter.harmony,
                              assay = opts$assay, seed = opts$seed,
                              fig.res = opts$fig.res, ...)
      opts$npcs <- npcs
      # Total number of dimensions used for Harmony
      opts$harmony_dims <- length(dims.use)
      # Store final object prior to exploring the data
      saveRDS(list(seu = seu, io = io, opts = opts),
              file = paste0(io$out_dir, "int_", obj_name, ".rds"))
    } else{
      obj <- readRDS(paste0(io$out_dir, "int_", obj_name, ".rds"))
      seu <- obj$seu
      rm(obj)
      gc(verbose = FALSE)
    }

    # Explore the data
    tmp <- lapply(X = opts$ndim, FUN = function(ndim) {
      if (ndim > opts$harmony_dims) { return(0) }
      # Create directories
      dim_dir <- paste0(out_dir, "/dim", ndim, "/")
      qc_dir <- paste0(dim_dir, "/qc/")
      dir.create(qc_dir, recursive = TRUE)
      # Number of dimensions to perform UMAP and clustering
      ndim <- min(ndim, opts$harmony_dims)
      # Run UMAP
      seu <- Seurat::RunUMAP(seu, dims = 1:ndim, seed.use = opts$seed,
                             reduction = "harmony", ...)
      # QC plots
      dimred_qc_plots(seu = seu, reductions = c("umap", "pca"),
                      metadata_to_plot = opts$meta_to_plot,
                      qc_to_plot = opts$qc_to_plot, plot_dir = qc_dir,
                      max.cutoff = opts$max.cutoff,
                      legend.position = "top",
                      col_pal = opts$col_pal, dims_plot = opts$dims_plot,
                      pt.size = opts$pt.size, fig.res = opts$fig.res)

      # Compute module scores for modules of interest
      module_dir <- paste0(dim_dir, "/modules/")
      dir.create(module_dir, recursive = TRUE)
      seu <- module_score_analysis(seu = seu, modules_group = opts$modules_group,
                                   plot_dir = module_dir, reduction = "umap",
                                   max.cutoff = opts$max.cutoff, ctrl = opts$ctrl,
                                   seed = opts$seed, fig.res = opts$fig.res,
                                   legend.position = "top", col_pal = opts$col_pal,
                                   dims_plot = opts$dims_plot, ...)
      # Different clustering resolutions and DGE
      cl_dir <- paste0(dim_dir, "/clusters/")
      dir.create(cl_dir, recursive = TRUE)
      seu <- cluster_analysis(seu = seu, ndim = ndim, res = opts$res,
                              logfc.threshol = opts$logfc.threshold,
                              min.pct = opts$min.pct, only.pos = TRUE,
                              topn = opts$topn_genes, plot_dir = cl_dir,
                              plot_cluster_markers = opts$plot_cluster_markers,
                              modules_group = opts$modules_group,
                              cluster_reduction = "harmony",
                              plot_reduction = "umap", max.cutoff = opts$max.cutoff,
                              cl.legend.position = "right", legend.position = "top",
                              col_pal = opts$col_pal, pt.size = opts$pt.size,
                              label = opts$label, label.size = opts$label.size,
                              assay = opts$assay, ctrl = opts$ctrl, seed = opts$seed,
                              fig.res = opts$fig.res, ...)
    })
  }
}
