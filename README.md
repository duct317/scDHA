# scDHA
The scDHA software package can perform cell segregation through unsupervised learning, dimension reduction and visualization, cell classification, and time-trajectory inference on single-cell RNA sequencing data. 

# How to install:
- The package can be installed from this repository.
- Install devtools: `utils::install.packages('devtools')`
- Install the package using: `devtools::install_github('duct317/scDHA')`  
  Or, install with manual and vignette: `devtools::install_github('duct317/scDHA', build_manual = T, build_vignettes = T)`
- If necessary, install miniconda: `reticulate::install_miniconda(force = T)`
- Install tensorflow and keras in python using: `keras::install_keras(tensorflow = "1.10.0")`
- For more information about installation of keras, please visit https://keras.rstudio.com/

# To run the Goolam example:
- Load the package: `library(scDHA)`
- Load SingleCellExperiment package to read dataset: `library(SingleCellExperiment)`
- Load Goolam dataset: `data('Goolam'); data <- t(assay(Goolam)); label <- Goolam$cell_type1`
- Log transform the data: `data <- log2(data + 1)`
- Generating clustering result: `result <- scDHA(data, seed = 1)`
- The clustering result can be found here: `cluster <- result$cluster`
- Calculating adjusted Rand Index using mclust package: `mclust::adjustedRandIndex(cluster,label)`
- A detailed tutorial on how to use scDHA package is available at http://scdha.tinnguyen-lab.com/  
  Or, a vignette in R Notebook format is available [here](https://github.com/duct317/scDHA/blob/master/vignettes/Example.Rmd)

# How to use the package for new data 
To use our package for new data, the package includes these functions:  
- scDHA: main function, doing dimension reuction and clustering. The input is a matrix with rows as samples and columns as genes.
- scDHA.w: plot the normalized weight variances to select suitable cutoff for gene filtering (optional).
- scDHA.vis: visualization. The input is demension reduction output.
- scDHA.pt: generating pseudotime. The input is demension reduction output.
- scDHA.class: classification new data using available one. The inputs consist of train data matrix, train data label and new data matrix. 
- The result is reproducible by setting seed for these functions.
- More detail about parameters for each function could be found in the manual.

# Citation:
Duc Tran, Hung Nguyen, Bang Tran, Carlo La Vecchia, Hung N. Luu, Tin Nguyen (2021). Fast and precise single-cell data analysis using a hierarchical autoencoder. <i>Nature Communications</i>, <b>12</b>, 1029. doi: 10.1038/s41467-021-21312-2 ([link](https://www.nature.com/articles/s41467-021-21312-2)) 
