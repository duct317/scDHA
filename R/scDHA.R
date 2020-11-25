#' @import keras tensorflow doParallel
#' @importFrom parallel makeCluster clusterEvalQ stopCluster
#' @importFrom matrixStats colSums2 rowSds rowMeans2 rowMaxs rowMins colSds colMins
#' @importFrom stats predict
#' @importFrom foreach %dopar% foreach
#' @importFrom stats quantile
#' @importFrom igraph membership
#' @title scDHA
#' @description The main function to perform dimension deduction and clustering.
#' @param data Gene expression matrix, with rows represent samples and columns represent genes.
#' @param k Number of clusters, leave as default for auto detection. Has no effect when \code{do.clus = False}.
#' @param sparse Boolen variable indicating whether data is a sparse matrix. The input must be a non negative sparse matrix.
#' @param n Number of genes to keep after feature selection step.
#' @param ncores Number of processor cores to use.
#' @param gen_fil Boolean variable indicating whether to perform scDHA gene filtering before performing dimension deduction and clustering.
#' @param do.clus Boolean variable indicating whether to perform scDHA clustering. If \code{do.clus = False}, only dimension deduction is performed.
#' @param sample.prob Probability used for classification application only. Leave this parameter as default, no user input is required.
#' @param seed Seed for reproducibility.
#' @return List with the following keys:
#' \itemize{
#' \item cluster - A numeric vector containing cluster assignment for each sample. If \code{do.clus = False}, this values is always \code{NULL}.
#' \item latent - A matrix representing compressed data from the input data, with rows represent samples and columns represent latent variables.
#' }
#' @export


scDHA <- function(data = data, k = NULL, sparse = F, n = 5000, ncores = 15L, gen_fil = T, do.clus = T, sample.prob = NULL, seed = NULL) {
  K = 3
  if(sparse) {
    if(min(data) < 0) stop("The input must be a non negative sparse matrix.")
    if(max(data) - min(data) > 100)
    {
      data@x <- log2(data@x + 1)
    }
    data <- normalize_data_sparse(data)
  }  else {
    if(max(data) - min(data) > 100)
    {
      if(min(data) < 0)
      {
        if(nrow(data) == ncol(data))
        {
          data <- t(data)
          data <- data - matrixStats::rowMins(data)
          data <- t(data)
        } else {
          data <- data - matrixStats::colMins(data)
        }
      }
      data <- log2(data + 1)
    }
    data <- normalize_data_dense(data)
  } 
  gc()
  if(nrow(data) >= 50000)
  {
    res <- scDHA.big(data = data, k = k, K = K, n = n, ncores = ncores, gen_fil = gen_fil, do.clus = do.clus, sample.prob = sample.prob, seed = seed)
  } else {
    res <- scDHA.small(data = data, k = k, K = K, n = n, ncores = ncores, gen_fil = gen_fil, do.clus = do.clus, sample.prob = sample.prob, seed = seed)
  }
  res
}

gene.filtering <- function(data.list, original_dim, batch_size, ncores.ind, wdecay, seed)
{
  or <- list()
  cl <- parallel::makeCluster(3, outfile = "/dev/null")
  registerDoParallel(cl, cores = 3)
  parallel::clusterEvalQ(cl,{
    library(scDHA)
  })
  or <- foreach(i = seq(3)) %dopar%
    {
      if (is.null(seed))
      {
        config <- list()
        config$intra_op_parallelism_threads <- ncores.ind
        config$inter_op_parallelism_threads <- ncores.ind
        session_conf <- do.call(tf$ConfigProto, config)
        sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
        k_set_session(session = sess)
      } else {
        set.seed((seed+i))
        use_session_with_seed((seed+i))
      }
      
      data.tmp <- data.list[[i]]
      batch_size <-round(nrow(data.tmp)/50)
      
      x <- layer_input(shape = c(original_dim))
      h <- layer_dense(x, 50, kernel_constraint = constraint_nonneg())
      x_decoded_mean <- layer_dense(h, original_dim)
      vae <- keras_model(x, x_decoded_mean)
      vae %>% compile(optimizer =   tf$contrib$opt$AdamWOptimizer(wdecay,1e-3), loss = 'mse')
      
      
      his <- vae %>% fit(
        data.tmp, data.tmp,
        shuffle = TRUE,
        epochs = 10,
        batch_size = batch_size,
        verbose = 0
      )
      
      W <- get_weights(get_layer(vae, index = 2))[[1]]
      Wsd <- rowSds(W)
      Wsd[is.na(Wsd)] <- 0
      Wsd <- (Wsd-min(Wsd))/(max(Wsd)-min(Wsd))
      Wsd
    }
  parallel::stopCluster(cl)
  or <- rowMeans2(as.matrix(data.frame(or)))
  or
}

latent.generating <- function(da, or.da, batch_size, K, ens, epsilon_std, lr, beta, 
                              original_dim_reduce, latent_dim, intermediate_dim, epochs, wdecay, ncores.ind, sample.prob, seed)
{
  latent <- list()
  cl <- parallel::makeCluster(K, outfile = "/dev/null")
  doParallel::registerDoParallel(cl, cores = K)
  parallel::clusterEvalQ(cl,{
    library(scDHA)
  })
  latent <- foreach(i = 1:K) %dopar% {
    if (is.null(seed))
    {
      config <- list()
      config$intra_op_parallelism_threads <- ncores.ind
      config$inter_op_parallelism_threads <- ncores.ind
      session_conf <- do.call(tf$ConfigProto, config)
      sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
      k_set_session(session = sess)
    } else {
      set.seed((seed+i))
      use_session_with_seed((seed+i))
    }
    
    if(nrow(or.da) >= 50000)
    {
      da <- as.matrix(or.da[sample.int(nrow(or.da), 5000, replace = F, prob = sample.prob),])
      batch_size <- round(nrow(da)/50)
    } else if(nrow(or.da) > 2000)
    {
      da <- as.matrix(or.da[sample.int(nrow(or.da), 2000, replace = F, prob = sample.prob),])
      batch_size <- round(nrow(da)/50)
    }
    
    x <- layer_input(shape = c(original_dim_reduce))
    
    sampling <- function(arg){
      z_mean <- arg[, 1:(latent_dim)]
      z_var <- arg[, (latent_dim + 1):(2 * latent_dim)]
      
      epsilon <- k_random_normal(
        shape = c(k_shape(z_mean)[[1]]),
        mean=0.,
        stddev=epsilon_std
      )
      
      z_mean + k_sqrt(z_var)*epsilon
    }
    
    if(nrow(da) > 500)
    {
      h <- layer_dense(x , intermediate_dim, activation = NULL)
      h <- layer_batch_normalization(h, batch_size = 100)
    } else {
      h <- layer_dense(x , intermediate_dim, activation = "selu")
    }
    
    z_mean <- layer_dense(h, latent_dim)
    
    z_var <- layer_dense(h, latent_dim, activation = "softmax")
    
    z <- list()
    h_decoded <- list()
    x_decoded_mean <- list()
    for (e in 1:2) {
      z[[e]] <- layer_concatenate(list(z_mean, z_var)) %>%
        layer_lambda(sampling)
      h_decoded[[e]] <- layer_dense(z[[e]], units = intermediate_dim, activation = "selu")
      x_decoded_mean[[e]] <- layer_dense(h_decoded[[e]], original_dim_reduce)
    }
    
    da.list <- list()
    for (e in 1:2) {
      da.list[[e]] <- da
    }
    
    # end-to-end autoencoder
    vae <- keras_model(x, x_decoded_mean)
    
    # encoder, from inputs to latent space
    encoder <- keras_model(x, z_mean)
    
    #WU
    vae %>% compile(optimizer =   tf$contrib$opt$AdamWOptimizer(wdecay,lr[1]), loss = 'mae')
    
    his <- vae %>% fit(
      da, da.list,
      shuffle = TRUE,
      epochs = epochs[1],
      batch_size = batch_size,
      verbose = 0
    )
    
    vae_loss <- function(x, x_pred){
      xent_loss <- loss_mean_squared_logarithmic_error(x, x_pred) 
      kl_loss <- -0.5*k_mean(1 + k_log(z_var) - k_square(z_mean) - z_var, axis = -1L)
      loss <- xent_loss*beta + kl_loss
      loss
    }
    
    #Train
    vae %>% compile(optimizer = tf$contrib$opt$AdamWOptimizer(wdecay,lr[2]), loss = vae_loss)
    
    his <- vae %>% fit(
      da, da.list,
      shuffle = TRUE,
      epochs = epochs[2],
      batch_size = batch_size,
      verbose = 0
    )
    
    predict_on_sparse <- function(data, model, latent_dim)
    {
      folds <- round(seq(1, nrow(data), length.out = round(nrow(data)/10000)))
      tmp <- matrix(ncol = latent_dim, nrow = nrow(data))
      
      for (i in 2:length(folds)) {
        tmp[folds[i-1]:folds[i], ] <- model %>% predict(as.matrix(data[folds[i-1]:folds[i], ]))
      }
      
      tmp
    }
    
    latent.tmp <- list()
    for (ite in 1:ens) {
      his <- vae %>% fit(
        da, da.list,
        shuffle = TRUE,
        epochs = epochs[2]+ite,
        batch_size = batch_size,
        verbose = 0,
        initial_epoch = (epochs[2]+ite-1)
      )
      
      if(nrow(or.da) > 20000)
      {
        tmp <- predict_on_sparse(or.da, encoder, latent_dim)
      } else {
        tmp <- encoder %>% predict(as.matrix(or.da))
      }
      latent.tmp[[length(latent.tmp)+1]] <- tmp
    }
    latent.tmp
  }
  parallel::stopCluster(cl)
  latent
}

scDHA.small <- function(data = data, k = NULL, K = 3, n = 5000, ncores = 15L, gen_fil = T, do.clus = T, sample.prob = NULL, seed = NULL) {
  set.seed(seed)
  
  ncores.ind <- as.integer(max(1,floor(ncores/K)))
  original_dim <- ncol(data)
  wdecay <- 1e-6
  batch_size <- max(round(nrow(data)/50),2)
  epochs <- c(10,20)
  
  lr <- c(5e-4, 5e-4, 5e-4)
  gen_fil <- (gen_fil & (ncol(data) > n))
  n <- ifelse(gen_fil, min(n, ncol(data)), ncol(data))
  epsilon_std <- 0.25
  beta <- 250
  ens  <-  3
  
  intermediate_dim = 64
  latent_dim = 15
  
  #Feature selection
  if (gen_fil) {
    data.list <- lapply(seq(3), function(i) {
      if(!is.null(seed)) set.seed((seed+i))
      if(nrow(data) > 2000)
      {
        ind <- sample.int(nrow(data), 2000, replace = F, prob = sample.prob)
      } else {
        ind <- seq(nrow(data))
      }
      data.tmp <- as.matrix(data[ind,])
      data.tmp
    })
    
    or <- gene.filtering(data.list = data.list, original_dim = original_dim, batch_size = batch_size, ncores.ind = ncores.ind, wdecay = wdecay, seed = seed)
    
    keep <- which(or > quantile(or,(1-min(n,original_dim)/original_dim)))
    
    da <- data[,keep]
    original_dim_reduce <- ncol(da)
    or.da <- da
  } else {
    da <- data
    original_dim_reduce <- ncol(da)
    or.da <- da
  }
  
  #Latent generation
  latent <- latent.generating(da, or.da, batch_size, K, ens, epsilon_std, lr, beta, 
                              original_dim_reduce, latent_dim, intermediate_dim, epochs, wdecay, ncores.ind, sample.prob, seed)
  
  gc(reset = TRUE)
  
  latent.save <- latent
  latent <- list()
  for (i in 1:K) {
    latent <- c(latent, latent.save[[i]])
  }
  
  if (do.clus) {
    result <- list()
    cl <- parallel::makeCluster(ncores, outfile = "/dev/null")
    doParallel::registerDoParallel(cl, cores = ncores)
    parallel::clusterEvalQ(cl,{
      library(scDHA)
    })
    result$all <- foreach(x = latent) %dopar% {
      set.seed(seed)
      if (nrow(x) > 5000)
      {
        cluster <- clus.big(x , k = k)
      } else {
        cluster <- clus(x, k = k)
      }
      cluster
    }
    parallel::stopCluster(cl)
    
    set.seed(seed)
    
    final <- clustercom2(result)
    g.en <- latent[[which.max(sapply(result$all, function(x) adjustedRandIndex(x,final)))]]
    final <- as.numeric(factor(final))
    
    list( cluster = final,
          latent = g.en,
          all.latent = latent,
          keep = colnames(or.da),
          all.res = result$all)
  } else {
    list( all.latent = latent,
          keep = colnames(or.da))
  }
  
}

scDHA.big <- function(data = data, k = NULL, K = 3, n = 5000, ncores = 15L, gen_fil = T, do.clus = T, sample.prob = NULL, seed = NULL) {
  set.seed(seed)
  
  ncores.ind <- as.integer(max(1,floor(ncores/K)))
  original_dim <- ncol(data)
  wdecay <- 1e-4
  batch_size <- max(round(nrow(data)/50),2)
  epochs <- c(10,20)
  
  lr <- rep(1e-3, 3) 
  gen_fil <- (gen_fil & (ncol(data) > n))
  n <- ifelse(gen_fil, min(n, ncol(data)), ncol(data))
  epsilon_std <- 0.25
  beta <- 250
  ens  <-  3
  
  intermediate_dim = 64
  latent_dim = 25
  
  #Feature selection
  if (gen_fil) {
    data.list <- lapply(1:3, function(i) {
      if(!is.null(seed)) set.seed((seed+i))
      ind <- sample.int(nrow(data), 5000, replace = F, prob = sample.prob)
      data.tmp <- as.matrix(data[ind,])
      data.tmp
    })
    
    or <- gene.filtering(data.list = data.list, original_dim = original_dim, batch_size = batch_size, ncores.ind = ncores.ind, wdecay = wdecay, seed = seed)
    
    keep <- which(or > quantile(or,(1-min(n,original_dim)/original_dim)))
    
    da <- 0
    or.da <- data[,keep]
    original_dim_reduce <- ncol(or.da)
  } else {
    da <- 0
    or.da <- data
    original_dim_reduce <- ncol(or.da)
  }
  rm(data)
  gc(reset = T)
  
  #Latent generation
  latent <- latent.generating(da, or.da, batch_size, K, ens, epsilon_std, lr, beta, 
                              original_dim_reduce, latent_dim, intermediate_dim, epochs, wdecay, ncores.ind, sample.prob, seed)
  
  gc(reset = TRUE)
  
  latent.save <- latent
  latent <- list()
  for (i in 1:K) {
    latent <- c(latent, latent.save[[i]])
  }
  
  if (do.clus) {
    result <- list()
    cl <- parallel::makeCluster(ncores, outfile = "/dev/null")
    doParallel::registerDoParallel(cl, cores = ncores)
    parallel::clusterEvalQ(cl,{
      library(scDHA)
    })
    
    idx <- sample.int(nrow(latent[[1]] ), size = 5000)
    result$all <- foreach(x = latent) %dopar% {
      set.seed(seed)
      cluster <- clus.big(x, k = k, n = 5000, nmax = 15)
      
      cluster <- as.numeric(factor(cluster))
      
      cluster
    }
    parallel::stopCluster(cl)
    
    set.seed(seed)
    
    final <- clustercom2(result)
    g.en <- latent[[which.max(sapply(result$all, function(x) adjustedRandIndex(x,final)))]]
    final <- as.numeric(factor(final))

    list( cluster = final,
          latent = g.en,
          all.latent = latent,
          keep = colnames(or.da),
          all.res = result$all)
  } else {
    list( all.latent = latent,
          keep = colnames(or.da))
  }
  
}

#' @importFrom graphics plot
#' @title scDHA.w
#' @description This function will plot a graph with normalized weights of all genes so user can select the appropriate number of genes to keep.
#' @param data Gene expression matrix, with rows represent samples and columns represent genes.
#' @param sparse Boolen variable indicating whether data is a sparse matrix. The input must be a non negative sparse matrix.
#' @param ncores Number of processor cores to use.
#' @param seed Seed for reproducibility.
#' @return A plot with normalized weights of all genes.
#' @examples
#' \donttest{
#' #Generate weight variances for each genes
#' weight_variance <- scDHA.w(data, seed = 1)
#'
#' #Plot weight variances for top 5,000 genes
#' plot(weight_variance, xlab = "Genes", ylab = "Normalized Weight Variance", xlim=c(1, 5000))
#'
#' #Plot the change of weight variances for top 5,000 genes
#' weight_variance_change <- weight_variance[-length(weight_variance)] - weight_variance[-1] 
#' plot(weight_variance_change, xlab = "Genes", ylab = "Weight Variance Change", xlim=c(1, 5000))
#' }
#' @export
scDHA.w <- function(data = data, sparse = F, ncores = 15L, seed = NULL) {
  K = 3
  sample.prob = NULL
  do.clus = T
  gen_fil = T
  k = NULL
  if(sparse) {
    if(min(data) < 0) stop("The input must be a non negative sparse matrix.")
    if(max(data) - min(data) > 100)
    {
      data@x <- log2(data@x + 1)
    }
    data <- normalize_data_sparse(data)
  }  else {
    if(max(data) - min(data) > 100)
    {
      if(min(data) < 0)
      {
        if(nrow(data) == ncol(data))
        {
          data <- t(data)
          data <- data - matrixStats::rowMins(data)
          data <- t(data)
        } else {
          data <- data - matrixStats::colMins(data)
        }
      }
      data <- log2(data + 1)
    }
    data <- normalize_data_dense(data)
  } 
  gc()
  if(nrow(data) >= 50000)
  {
    res <- scDHA.big.w(data = data, k = k, K = K, ncores = ncores, gen_fil = gen_fil, do.clus = do.clus, sample.prob = sample.prob, seed = seed)
  } else {
    res <- scDHA.small.w(data = data, k = k, K = K, ncores = ncores, gen_fil = gen_fil, do.clus = do.clus, sample.prob = sample.prob, seed = seed)
  }
  w <- res
  w <- (w - min(w))/(max(w) - min(w))
  w <- sort(w, decreasing = T)
  
  plot(w, xlab = "Genes", ylab = "Normalized Weight Variance")
  
  w
}

scDHA.small.w <- function(data = data, k = NULL, K = 3, ncores = 15L, gen_fil = T, do.clus = T, sample.prob = NULL, seed = NULL) {
  set.seed(seed)
  
  ncores.ind <- as.integer(max(1,floor(ncores/K)))
  original_dim <- ncol(data)
  wdecay <- 1e-6
  batch_size <- max(round(nrow(data)/50),2)
  epochs <- c(10,20)
  
  lr <- c(5e-4, 5e-4, 5e-4)
  gen_fil <- (gen_fil & (ncol(data) > 5000))
  n <- ifelse(gen_fil, min(5000, ncol(data)), ncol(data))
  epsilon_std <- 0.25
  beta <- 250
  ens  <-  3
  
  intermediate_dim = 64
  latent_dim = 15
  
  #Feature selection

  data.list <- lapply(seq(3), function(i) {
    if(!is.null(seed)) set.seed((seed+i))
    if(nrow(data) > 2000)
    {
      ind <- sample.int(nrow(data), 2000, replace = F, prob = sample.prob)
    } else {
      ind <- seq(nrow(data))
    }
    data.tmp <- as.matrix(data[ind,])
    data.tmp
  })
  
  or <- gene.filtering(data.list = data.list, original_dim = original_dim, batch_size = batch_size, ncores.ind = ncores.ind, wdecay = wdecay, seed = seed)

  or
  
}

scDHA.big.w <- function(data = data, k = NULL, K = 3, ncores = 15L, gen_fil = T, do.clus = T, sample.prob = NULL, seed = NULL) {
  set.seed(seed)
  
  ncores.ind <- as.integer(max(1,floor(ncores/K)))
  original_dim <- ncol(data)
  wdecay <- 1e-4
  batch_size <- max(round(nrow(data)/50),2)
  epochs <- c(10,20)
  
  lr <- rep(1e-3, 3) #c(5e-4, 5e-4, 5e-4)
  gen_fil <- (gen_fil & (ncol(data) > 5000))
  n <- ifelse(gen_fil, min(5000, ncol(data)), ncol(data))
  epsilon_std <- 0.25
  beta <- 250
  ens  <-  3
  
  intermediate_dim = 64
  latent_dim = 25
  
  #Feature selection

  data.list <- lapply(1:3, function(i) {
    if(!is.null(seed)) set.seed((seed+i))
    ind <- sample.int(nrow(data), 5000, replace = F, prob = sample.prob)
    data.tmp <- as.matrix(data[ind,])
    data.tmp
  })
  
  or <- gene.filtering(data.list = data.list, original_dim = original_dim, batch_size = batch_size, ncores.ind = ncores.ind, wdecay = wdecay, seed = seed)

  or
  
}


#' @title Goolam
#'
#' @description Goolam dataset in SingleCellExperiment format, include scRNA-seq data and cell type information.
"Goolam"



