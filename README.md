# SeuratPipe

The aim of the `SeuratPipe` package is to remove technical details from the user, enabling wet lab researchers to have an 'initial look' at the data themselves.

## Installation
```R
# install.packages("remotes")
remotes::install_github("andreaskapou/SeuratPipe")
```

### Scrublet Python installation
Scrublet is a method for doublet detection, Github page: [https://github.com/swolock/scrublet](https://github.com/swolock/scrublet). 

The method is implemented in Python, hence when installing `SeuratPipe`, the Scrublet package cannot be installed automatically. We can try to install and run Scrublet within R using the `reticulate` package. I wrote a helper function to ease installation of Scrublet. In R type:
```R
SeuratPipe::install_scrublet()
# type ?install_scrublet for more details
```

If you have any issues with installation, you can still use the SeuratPipe package, although you have to set `use_scrublet = FALSE`, i.e. do not perform any doublet detection in your analysis.

## QC analysis pipeline
`SeuratPipe` assumes the user will have a __sample metadata__ file/data frame which contains all the metadata for each scRNA-seq sample. It also requires the following column names to be present: `sample`, `donor`, `path`, `condition` and `pass_qc` (logical `TRUE` or `FALSE`).

We can then perform all QC pipeline steps using a single wrapper function `run_qc_pipeline` as shown below. This function will also generate all the required QC plots. 

```{R
seu <- run_qc_pipeline(...)
# type ?run_qc_pipeline for argument details
```

## Clustering analysis pipeline
Here we show how to perform the most common analysis pipeline in Seurat for cell type identification. 
The wrapper function performs the following steps: 

1. Data processing, e.g. Log normalisation and scaling.
2. Run PCA. 
3. Perform UMAP and clustering. 
4. Generate plots / outputs in corresponding directories.

```{R
seu <- run_cluster_pipeline(...)
# type ?run_cluster_pipeline for argument details
```

## Harmony integration pipeline
The Harmony integration pipeline is almost identical to `run_cluster_pipeline`, with the main difference of the function accepting a list of Seurat objects (each containing information for a single sample) to perform data integration.
```{R
seu <- run_harmony_pipeline(...)
# type ?run_harmony_pipeline for argument details
```

# Acknowledgments
I would like to thank John Wilson-Kanamori, Jordan Portman and Nick Younger. The `SeuratPipe` package has adapted many of the analyses pipelines that were developed by the computational biology group in Neil Henderson's lab at the University of Edinburgh.

This package was supported by funding from the University of Edinburgh and Medical Research Council (core grant to the MRC Institute of Genetics and Cancer) for the Cross-Disciplinary Fellowship (XDF) programme.
