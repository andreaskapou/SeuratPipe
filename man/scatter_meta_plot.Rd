% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_plot.R
\name{scatter_meta_plot}
\alias{scatter_meta_plot}
\title{Tailored scatter plot of metadata}
\usage{
scatter_meta_plot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  pt.size = 1.2
)
}
\arguments{
\item{seu}{Seurat object (required).}

\item{features}{Vector of features in metadata to plot.}

\item{pt.size}{Point size.}
}
\value{
A ggplot2 list with all combinations of elements in features
}
\description{
This function creates a scatter plot to assess feature
to feature relationships
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
