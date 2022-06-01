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

  # https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
  col_pal <- c(
    "indianred", # red
    "#6699CB",
    "#FDBF6F", # lt orange
    "#CAB2D6", # lt purple
    "#FB9A99", # lt pink
    "tan3",
    "darkolivegreen4", # darkgreen
    "darkgrey", # darkgrey
    "skyblue2", # lightblue
    "mediumpurple1",
    "darkseagreen3",
    "khaki2", "ivory3",
    "steelblue4",
    "#6A3D9A", # purple
    "seagreen", "orchid1", "blue1", "deeppink1", "gold1",
    "darkturquoise", "darkorange4", "#FF7F00",
    "dodgerblue", "yellow3", "mediumorchid1", "firebrick4",
    "wheat4", "maroon", "grey30", "red2", "burlywood2", "cyan",
    "darkolivegreen2", "yellowgreen"
  )

  assign("discrete_col_pal", col_pal, envir = topenv())
  base::invisible()
}
