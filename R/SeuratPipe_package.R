#' @title \code{SeuratPipe}: Streamlining Seurat analysis
#' @description SeuratPipe contains common analysis, processing and
#' plotting functions with an attempt to streamline single-cell genomics
#' analysis using the Seurat package.
#' @docType package
#' @name SeuratPipe
#'
#' @return seurat_pipeline main package documentation.
#'
#'
#' @author C.A.Kapourani \email{kapouranis.andreas@@gmail.com}
#'
#' @rawNamespace importFrom(magrittr, "%>%")
#' @rawNamespace importFrom(rlang, "%||%")
#' @importFrom grDevices png dev.off dev.list colorRampPalette
#' @importFrom stats filter median
#' @importFrom utils head read.csv write.csv
#' @import ggplot2 BiocStyle
#'
.datatable.aware <- TRUE
NULL
#> NULL


.onLoad <- function(libname = find.package("SeuratPipe"),
                    pkgname = "SeuratPipe"){
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(
      # sample file names from taxstats
      c(# we use the magrittr pipe
        "."
      )
    )
  invisible()
}
