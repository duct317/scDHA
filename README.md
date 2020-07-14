# scDHA
The scDHA software package conducts cell segregation through unsupervised learning, dimension reduction and visualization, cell classification, and time-trajectory inference. scDHA currently supports Linux.
# How to install
- The package can be installed from this repository.
- Install devtools: `utils::install.packages('devtools')`
- Install the package using: `devtools::install_github('duct317/scDHA')`
- Install keras package version 2.2.4.1 using: `devtools::install_version('keras', version = '2.2.4.1', repos = 'https://cran.rstudio.com/')`
- Install tensorflow and keras in python using: `keras::install_keras(tensorflow = "1.10.0")`
- For more information about installation of keras, please visit https://keras.rstudio.com/

# To run the Goolam example:
- Download the `goolam.rds` file.
- Load the package: `library(scDHA)`
- Load SingleCellExperiment package to read dataset: `library(SingleCellExperiment)`
- Load Goolam dataset: `data <- t(assay(readRDS('goolam.rds'))); label <- readRDS('goolam.rds')$cell_type1`
- Log transform the data: `data <- log2(data + 1)`
- Generating clustering result: `result <- scDHA(data, seed = 1)`
- The clustering result can be found here: `cluster <- result$cluster`
- Calculating adjusted Rand Index using mclust package: `mclust::adjustedRandIndex(cluster,label)`
