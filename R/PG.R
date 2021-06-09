#' @importFrom igraph graph.data.frame cluster_louvain membership modularity
phenograph <- function(data, k=30){
  if(is.data.frame(data))
    data <- as.matrix(data)
  
  if(!is.matrix(data))
    stop("Wrong input data, should be a data frame of matrix!")
  
  if(k<1){
    stop("k must be a positive integer!")
  }else if (k > nrow(data)-2){
    stop("k must be smaller than the total number of points!")
  }
  
  t1 <- system.time(neighborMatrix <- find_neighbors(data, k=k+1)[,-1])
  t2 <- system.time(links <- jaccard_coeff(neighborMatrix))
  
  links <- links[links[,1]>0, ]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  t3 <- system.time(g <- graph.data.frame(relations, directed=FALSE))

  t4 <- system.time(community <- igraph::cluster_louvain(g))

  return(list(g, community))
}

find_neighbors <- function(data, k){
  nearest <- find_nn(data, k)$idx
  return(nearest)
}

#' @import RcppAnnoy
find_nn <- function(X, k, 
                    metric = "euclidean",
                    n_trees = 50, search_k = 2 * k * n_trees,
                    tmpdir = tempdir(),
                    ret_index = FALSE,
                    verbose = FALSE) {
  
  res <- annoy_nn(X,
                  k = k,
                  metric = metric,
                  n_trees = n_trees, search_k = search_k,
                  tmpdir = tmpdir,
                  ret_index = ret_index,
                  verbose = verbose)
  
  
  res
}


annoy_nn <- function(X, k = 10,
                     metric = "euclidean",
                     n_trees = 50, search_k = 2 * k * n_trees,
                     tmpdir = tempdir(),
                     grain_size = 1,
                     ret_index = FALSE,
                     verbose = FALSE) {
  ann <- annoy_build(X,
                     metric = metric, n_trees = n_trees,
                     verbose = verbose
  )
  
  res <- annoy_search(X,
                      k = k, ann = ann, search_k = search_k,
                      tmpdir = tmpdir,
                      grain_size = grain_size, verbose = verbose
  )
  
  nn_acc <- sum(res$idx == 1:nrow(X)) / nrow(X)
  
  res <- list(idx = res$idx, dist = res$dist, recall = nn_acc)
  if (ret_index) {
    res$index <- ann
  }
  res
}

annoy_build <- function(X, metric = "euclidean", n_trees = 50,
                        verbose = FALSE) {
  nr <- nrow(X)
  nc <- ncol(X)
  
  ann <- create_ann(metric, nc)
  
  # Add items
  for (i in 1:nr) {
    ann$addItem(i - 1, X[i, ])
  }
  
  # Build index
  ann$build(n_trees)
  
  ann
}

#' @importFrom methods new
# create RcppAnnoy class from metric name with ndim dimensions
create_ann <- function(name, ndim) {
  ann <- switch(name,
                cosine =  methods::new(RcppAnnoy::AnnoyAngular, ndim),
                manhattan = methods::new(RcppAnnoy::AnnoyManhattan, ndim),
                euclidean = methods::new(RcppAnnoy::AnnoyEuclidean, ndim),
                hamming = methods::new(RcppAnnoy::AnnoyHamming, ndim),
                stop("BUG: unknown Annoy metric '", name, "'")
  )
  ann
}

# Search a pre-built Annoy index for neighbors of X
annoy_search <- function(X, k, ann,
                         search_k = 100 * k,
                         tmpdir = tempdir(),
                         grain_size = 1,
                         verbose = FALSE) {

  res <- annoy_search_serial(
    X = X, k = k, ann = ann,
    search_k = search_k,
    verbose = verbose
  )
  
  # Convert from Angular to Cosine distance
  if (methods::is(ann, "Rcpp_AnnoyAngular")) {
    res$dist <- 0.5 * (res$dist * res$dist)
  }
  
  res
}

annoy_search_serial <- function(X, k, ann,
                                search_k = 100 * k,
                                verbose = FALSE) {
  nr <- nrow(X)
  
  idx <- matrix(nrow = nr, ncol = k)
  dist <- matrix(nrow = nr, ncol = k)
  for (i in 1:nr) {
    res <- ann$getNNsByVectorList(X[i, ], k, search_k, TRUE)
    if (length(res$item) != k) {
      stop(
        "search_k/n_trees settings were unable to find ", k,
        " neighbors for item ", i
      )
    }
    idx[i, ] <- res$item
    dist[i, ] <- res$distance
  }
  list(idx = idx + 1, dist = dist)
}
