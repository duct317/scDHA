# scDHA
The scDHA software package conducts cell segregation through unsupervised learning, dimension reduction and visualization, cell classification, and time-trajectory inference. scDHA currently supports macOS and Linux.
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

# How to use the package for new data 
To use our package for new data, the package includes these functions:
- scDHA: main function, doing dimension reuction and clustering. The input is a matrix with rows as samples and columns as genes.
- scDHA.w: plot the normalized weights to select suitable cutoff for gene filtering.
- scDHA.vis: visualization. The input is demension reduction output.
- scDHA.pt: generating pseudotime. The input is demension reduction output.
- scDHA.clas: classification new data using available one. The inputs consist of train data matrix, train data label and new data matrix. 
- The result is reproducible by setting seed for these functions.
- More detail about parameters for each function could be found in the manual.

# Citation:
Duc Tran, Hung Nguyen, Bang Tran, Carlo La Vecchia, Hung N. Luu, Tin Nguyen (2019). Fast and precise single-cell data analysis using hierarchical autoencoder. bioRxiv, 799817. doi: 10.1101/799817 ([link](https://www.biorxiv.org/content/10.1101/799817v2)) 
