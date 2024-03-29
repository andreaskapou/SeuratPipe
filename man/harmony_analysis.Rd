% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{harmony_analysis}
\alias{harmony_analysis}
\title{Analysis steps for Harmony integration}
\usage{
harmony_analysis(
  seu,
  batch_id = "sample",
  npcs,
  dims.use = NULL,
  plot_dir = NULL,
  n_hvgs = 3000,
  max.iter.harmony = 50,
  seed = 1,
  fig.res = 200,
  ...
)
}
\arguments{
\item{seu}{Seurat object (required).}

\item{batch_id}{Name of batch to try and remove with data integration
(required). Can also be a vector if multiple batch information are present.
Should be a column name in Seurat 'meta.data'. Default is "sample".
This parameter is called 'group.by.vars' in Harmony.}

\item{npcs}{Number of principal components.}

\item{dims.use}{Vector with PCs to use for Harmony integration, e.g. 1:50.}

\item{plot_dir}{Directory to save generated plots. If NULL, plots are
not saved.}

\item{n_hvgs}{Number of highly variable genes (HVGs) to compute, which will
be used as input to PCA.}

\item{max.iter.harmony}{Maximum number of iterations for Harmony integration.}

\item{seed}{Set a random seed, for reproducibility.}

\item{fig.res}{Figure resolution in ppi (see 'png' function).}

\item{...}{Additional named parameters passed to RunHarmony and other Seurat
processing functions, such as RunPCA and ScaleData.}
}
\value{
Updated Seurat object with integrated and processed data.
}
\description{
This function implements all the analysis steps for performing
Harmony integration on a Seurat object. These include, 1. merge all samples
in a single Seurat object (if a list of Seurat objects is provided) 2.
data normalisation, 3. identification of HVGs, 4. Scaling of expression
values and computing PCA, 5. running Harmony integration, 8. Generating
plots.
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
