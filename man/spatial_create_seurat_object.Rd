% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_analysis.R
\name{spatial_create_seurat_object}
\alias{spatial_create_seurat_object}
\title{Create Seurat object based on spatial sample metadata.}
\usage{
spatial_create_seurat_object(
  data_dir,
  sample_meta = NULL,
  sample_meta_filename = NULL,
  meta_colnames = c("donor", "condition", "pass_qc"),
  tenx_dir = "outs",
  ...
)
}
\arguments{
\item{data_dir}{Parent directory where all sample 10x files are stored.
Think of it as project directory.}

\item{sample_meta}{Sample metadata information in a Data.frame like object.
Columns should at least contain 'sample', 'donor', 'condition' and 'pass_qc'.}

\item{sample_meta_filename}{Filename of sample metadata information, same as 'meta'
parameter above. User should provide one of 'meta' or 'meta_filename'.}

\item{meta_colnames}{Sample metadata column names to store in Seurat metadata.}

\item{tenx_dir}{Name of 10x base directory, e.g. with outputs after running
cellranger. Default 'premrna_outs', i.e. assumes single-nuclei RNA-seq.}

\item{...}{Additional named parameters passed to Seurat.}
}
\value{
List of Seurat objects equal the number of samples in
the sample metadata file. If a single sample, returns a Seurat object.
}
\description{
This function creates Seurat object for each spatial
(e.g. Visium) sample in the metadata file.
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
