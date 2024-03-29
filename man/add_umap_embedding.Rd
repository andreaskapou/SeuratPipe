% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{add_umap_embedding}
\alias{add_umap_embedding}
\title{Add UMAP embedding in Seurat object}
\usage{
add_umap_embedding(seu, embedding)
}
\arguments{
\item{seu}{Seurat object}

\item{embedding}{UMAP embedding. If a character, the function will assume a
filename is given and will read the corresponding UMAP embedding. Otherwise,
it assumes a UMAP embedding matrix is provided as input containing only the
embeddings as columns.}
}
\description{
Add UMAP embedding in existing Seurat object. This is the case when running
the pipeline and then you want to attach the stored UMAP embedding, stored
as CSV file, instead of recomputing the UMAP.
}
