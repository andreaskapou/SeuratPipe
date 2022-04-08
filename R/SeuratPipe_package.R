#' @title \code{SeuratPipe}: Streamlining Seurat analysis
#' @description SeuratPipe contains common analysis, processing and
#' plotting functions with an attempt to streamline single-cell genomics
#' analysis using the Seurat package.
#' @docType package
#' @name SeuratPipe
#'
#' @return SeuratPipe main package documentation.
#'
#'
#' @author C.A.Kapourani \email{kapouranis.andreas@@gmail.com}
#'
#' @rawNamespace importFrom(magrittr, "%>%")
#' @rawNamespace importFrom(rlang, "%||%")
#' @importFrom grDevices png pdf dev.off dev.list colorRampPalette
#' @importFrom stats filter median
#' @importFrom rlang fn_fmls
#' @importFrom utils head read.csv write.csv combn
#' @import ggplot2 BiocStyle
#'
NULL
#> NULL

# global reference to scrub (will be initialized in .onLoad)
scrub <- NULL

.onLoad <- function(libname = find.package("SeuratPipe"),
                    pkgname = "SeuratPipe"){
  # use superassignment to update global reference to scrub
  scrub <<- reticulate::import(module = "scrublet", convert = FALSE,
                               delay_load = TRUE)

  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(
      # sample file names from taxstats
      c(# we use the magrittr pipe
        "."
      )
    )

  col_pal <- c(
    "firebrick", # red
    "skyblue2", # lightblue
    "#CAB2D6", # lt purple
    "#FB9A99", # lt pink
    "#FDBF6F", # lt orange
    "green4", # darkgreen
    "darkgrey", # darkgrey
    "steelblue4",
    "#FF7F00", # orange
    "dodgerblue", # deepblue
    "maroon",
    "#6A3D9A", # purple
    "khaki2", "yellow4", "orchid1", "blue1", "deeppink1", "gold1", "ivory4",
    "darkturquoise", "palegreen2", "green1", "black", "darkorange4", "yellow3"
  )

  assign("discrete_col_pal", col_pal, envir = topenv())
  base::invisible()
}
