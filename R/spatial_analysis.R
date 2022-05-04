#' @rdname run_spatial_qc_pipeline
#'
#' @title Spatial QC pipeline
#'
#' @description This function implements all the analysis steps for perfoming
#' QC. These include: 1. reading all sample information from metadata object/file
#' and generating one Seurat object per sample. 2. Performs SoupX (ambient
#' RNA removal) and Scrublet (doublet detection) if user defines the
#' corresponding parameters. 3. Filter Seurat object according to QC
#' criteria 4. Generate correspond QC plots.
#'
#' @param alpha Controls opacity of spots. Provide as a vector specifying the
#' min and max range of values (between 0 and 1).
#' @param pt.size.factor Scale the size of the spots.
#' @param max.cutoff Vector of maximum cutoff values for each feature, may
#' specify quantile in the form of 'q##' where '##' is the quantile
#' (eg, 'q1', 'q10').
#' @param min.cutoff Vector of minimum cutoff values for each feature, may
#' specify quantile in the form of 'q##' where '##' is the quantile
#' (eg, 'q1', 'q10').
#' @param spatial_col_pal Continuous colour palette to use from viridis package to
#' colour spots on tissue, default "inferno".
#' @param crop Crop the plot in to focus on points plotted. Set to FALSE to
#' show entire background image.
#' @inheritParams run_qc_pipeline
#'
#' @return List of Seurat objects as the length of the number of samples in
#' the sample metadata file. If a single sample, return a Seurat object
#' instead of a list.
#'
#' @author C.A.Kapourani \email{C.A.Kapourani@@ed.ac.uk}
#'
#' @export
run_spatial_qc_pipeline <- function(
    data_dir, sample_meta, sample_meta_filename = NULL, nfeat_thresh = 500,
    mito_thresh = 5, meta_colnames = c("donor", "condition", "pass_qc"),
    out_dir = NULL, qc_to_plot = c("nFeature_Spatial", "nCount_Spatial", "percent.mito"),
    alpha = c(0.1, 0.9), pt.size.factor = 1.1, max.cutoff = "q98", min.cutoff = NA,
    spatial_col_pal = "inferno", crop = FALSE,
    tenx_dir = "outs", obj_filename = "seu_qc", ...) {

  # Store all parameters for reproducibility
  opts <- c(as.list(environment()), list(...))

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

  # Create Seurat object
  seu <- spatial_create_seurat_object(
    data_dir = data_dir, sample_meta = opts$sample_meta, sample_meta_filename = NULL,
    meta_colnames = meta_colnames, tenx_dir = tenx_dir, ...)

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

    spat_plot_dim <- .spatial_plot_dims(feat_len = length(qc_to_plot))
    pdf(paste0(plot_dir, "spatial_preqc.pdf"), width = spat_plot_dim$width,
        height = spat_plot_dim$height, useDingbats = FALSE)
    for (s in names(seu)) {
      print(spatial_feature_plot(
        seu[[s]], features = qc_to_plot, alpha = alpha,
        pt.size.factor = pt.size.factor, ncol = spat_plot_dim$ncols,
        max.cutoff = max.cutoff, min.cutoff = min.cutoff, crop = crop,
        col_pal = spatial_col_pal, legend.position = "top", title = s, ...))
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
                              ncol = spat_plot_dim$ncols, pt.size = 0.05, ...))
      }
      dev.off()

      spat_plot_dim <- .spatial_plot_dims(feat_len = length(qc_to_plot))
      pdf(paste0(plot_dir, "spatial_qc.pdf"), width = spat_plot_dim$width,
          height = spat_plot_dim$height, useDingbats = FALSE)
      for (s in names(seu)) {
        print(spatial_feature_plot(seu[[s]], features = qc_to_plot,
                                   alpha = alpha, pt.size.factor = pt.size.factor,
                                   ncol = spat_plot_dim$ncols, max.cutoff = max.cutoff,
                                   min.cutoff = min.cutoff,
                                   crop = crop, col_pal = spatial_col_pal,
                                   legend.position = "top", title = s, ...))
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
