% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_plot.R
\name{dim_plot_tailored}
\alias{dim_plot_tailored}
\title{Tailored dim plot}
\usage{
dim_plot_tailored(
  seu,
  group.by,
  reduction = "umap",
  legend.position = "right",
  col_pal = NULL,
  label = FALSE,
  label.size = 7,
  pt.size = 1.4,
  pt.shape = 21,
  pt.stroke = 0.05,
  pt.alpha = 1,
  dims_plot = c(1, 2),
  ...
)
}
\arguments{
\item{seu}{Seurat object}

\item{group.by}{Name of metadata column to group (color) cells by (required).}

\item{reduction}{Dimensionality reduction to use.}

\item{legend.position}{Position of legend, default "right" (set to "none"
for clean plot).}

\item{col_pal}{Continuous colour palette to use, default "RdYlBu".}

\item{label}{Whether to label the groups in 'reduction' space.}

\item{label.size}{Sets size of labels.}

\item{pt.size}{Adjust point size for plotting.}

\item{pt.shape}{Adjust point shape for plotting.}

\item{pt.stroke}{Stroke value for each point.}

\item{pt.alpha}{Adjust alpha value for each point.}

\item{dims_plot}{Dimensions to plot, must be a two-length numeric vector
specifying x- and y-dimensions.}

\item{...}{Additional parameters passed to ggplot2::geom_point.}
}
\value{
A ggplot2 object.
}
\description{
This function generates the same plot as \code{dim_plot},
although it focuses on a single group and generates slightly better
looking plot.
}
\author{
C.A.Kapourani \email{C.A.Kapourani@ed.ac.uk}
}
