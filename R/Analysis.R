#' @importFrom matrixStats colMeans2 rowSums2
#' @title scDHA visulization
#' @description  Generating 2D embeded data for visulation.
#' @param sc Embedding object produced by the \code{scDHA} function.
#' @param ncores number of processor cores to use.
#' @param seed Seed for reproducibility.
#' @return a list with the following keys:
#' \itemize{
#' \item pred - A matrix representing the 2D projection of single-cell data, where rows represent samples and columns represent latent components.
#' }
#' @export


scDHA.vis <- function(sc = sc, ncores = 15L, seed = NULL) {
  set.seed(seed)

  cal.dis <- function(arg) {
    a <- arg
    a2 <- k_square(a)
    a2sum <- k_sum(a2, axis = -1)
    dis <- k_expand_dims(a2sum) + a2sum
    ab <- k_dot(a, k_transpose(a))
    final.dis <- dis - 2 * ab
    final.dis <- k_sqrt(final.dis+1e-6)
    final.dis
  }
  
  tmp.data <- sc$latent
  
  if(nrow(tmp.data > 50000)) {
    tmp <- find_nn(tmp.data, 11L)
    D <- tmp$dist[, 2:11]
    D[D==0] <- min(D[D!=0])
    D <- 1/D
    C <- tmp$idx[, 2:11]
    g <- list(D = D, C = C)
  } else g <- calG(tmp.data, 10L)
  
  rand.res <- matrix(nrow = nrow(tmp.data), ncol = 20)
  
  for (i in 1:nrow(rand.res)) {
    rand.res[i, ] <- as.numeric(find.vec(i, g, 4, 5, 0.05, 0.95))
  }
  
  data.en <- matrix(nrow = nrow(tmp.data), ncol = ncol(tmp.data))
  
  for (i in 1:nrow(data.en)) {
    tmp <- colMeans2(tmp.data[c(rand.res[i, ], rep(i, 20)), ])
    data.en[i, ] <- tmp
  }
  
  
  idx.main <- 1:nrow(data.en)
  
  tmp.idx <- createFolds(idx.main, k = ceiling(nrow(data.en)/5000)  )
  if(length(tmp.idx) > 1) {
    for (i in 1:(length(tmp.idx) - 1) ) {
      tmp.idx[[i]] <- c(tmp.idx[[i]], sample(tmp.idx[[i+1]], ceiling(length(tmp.idx[[i]])/10 ) ) )
    }
  }
  
  data.list <- list()
  counter <- 1
  for (idx in tmp.idx) {
    data <- data.en[idx, ]
    sim.mat <- 1 - cor(t(data))
    
    sim.mat <- log(sim.mat + 1)

    tmp.cluster <- (kmeans(data, 15, nstart = 100, iter.max = 1e6))$cluster
    y1 <- to_categorical(as.numeric(tmp.cluster) - 1)
    
    data.list[[as.character(counter)]]$data <- data
    data.list[[as.character(counter)]]$y1 <- y1
    data.list[[as.character(counter)]]$sim.mat <- sim.mat
    counter <- counter + 1
    
  }
  
  cl <- parallel::makeCluster(1, outfile = "/dev/null")
  registerDoParallel(cl, cores = 1)
  parallel::clusterEvalQ(cl,{
    library(scDHA)
  })
  pred <- foreach(i = 1:1) %dopar%
    {
      if (is.null(seed))
      {
        config <- list()
        config$intra_op_parallelism_threads <- ncores
        config$inter_op_parallelism_threads <- ncores
        session_conf <- do.call(tf$ConfigProto, config)
        sess <- tf$Session(graph = tf$get_default_graph(), config = session_conf)
        k_set_session(session = sess)
      } else {
        set.seed((seed))
        use_session_with_seed((seed))
      }
      
      eps <- 0.1

      x <- layer_input(ncol(data.en))
      x1 <- layer_gaussian_noise(x, eps)
      h1 <- layer_dense(x1, 1024, activation = "elu") %>% layer_dropout(0.5)
      h2 <- layer_dense(h1, 32, activation = "sigmoid")
      z <- layer_dense(h2, 2) %>% layer_batch_normalization(center = F, momentum = 0.9)
      z1 <- layer_lambda(z, cal.dis, name = "z1")
      
      x_ <- layer_dense(z, ncol(data.list[[1]]$y1), activation = "softmax", name = "x_")
      
      model <- keras_model(x, list(z1, x_))
      model_out <- keras_model(x, z)
      
      model_loss <- function(x_or, x_pred) {
        x_pred <- k_log(x_pred + 1)
        x_pred <- (x_pred - k_mean(x_pred, axis = 1))/k_std(x_pred, axis = 1)
        x_pred <- k_transpose(x_pred)
        x_pred <- x_pred*k_abs(x_pred)/2
        xent_loss <- k_exp( -x_pred )
        xent_loss <- xent_loss*k_cast(k_greater(x_or, 0), k_floatx())
        tmp <- k_sum(xent_loss, axis = -1)
        xent_loss <- k_transpose( k_transpose(xent_loss)/tmp )
        loss <- loss_kullback_leibler_divergence(x_or, xent_loss)
        loss
      }
      
      
      model %>% compile(optimizer = optimizer_adam(5e-3, epsilon = 1e-4), loss = list(model_loss, loss_categorical_crossentropy), loss_weight = c(10, 1)) 
      
      max.ite <- ceiling(50/(length(data.list)**(2/3)) )
      
      for (ite in 1:max.ite) {
        for (subdata in data.list) {
          data <- subdata$data
          y1 <- subdata$y1
          sim.mat <- subdata$sim.mat
          
          idxs <- 1:nrow(data)
          idxs <- createFolds(idxs, min(10, ceiling(nrow(data)/ 50)))
          
          for (idx in idxs) {
            
            tmp.sim.mat <- t(scale(t(sim.mat[idx, idx])))
            sim.mat.dis <- tmp.sim.mat/2*abs(tmp.sim.mat)
            sim.mat.dis.exp <- exp(-sim.mat.dis)
            diag(sim.mat.dis.exp) <- 0
            his <- model %>% train_on_batch(data[idx, ], list(sim.mat.dis.exp / rowSums2(sim.mat.dis.exp), y1[idx, ]))
            
          }
          
        }
        
      }
      
      data <- data.en
      
      pred <- model_out %>% predict(data)
      pred

    }
  parallel::stopCluster(cl)
  
  data <- data.en
  pred <- pred[[1]]
  
  

  resol <- 20
  x.seq <- seq(min(pred[, 1]), max(pred[, 1]) + 0.01, length.out = resol)
  y.seq <- seq(min(pred[, 2]), max(pred[, 2]) + 0.01, length.out = resol)
  
  count <- NULL
  idx <- list()
  
  for (i in 1:(resol - 1)) {
    for (j in 1:(resol - 1)) {
      idx1 <- which(pred[, 1] >= x.seq[i] & pred[, 1] < x.seq[i + 1])
      idx2 <- which(pred[, 2] >= y.seq[j] & pred[, 2] < y.seq[j + 1])
      
      count <- c(count, length(intersect(idx1, idx2)))
      idx[[length(count)]] <- intersect(idx1, idx2)
    }
  }
  
  top.pos <- order(count, decreasing = T)[1]
  top.pos.coor <- matrix(colMeans2(data[idx[[top.pos]], ]), nrow = 1)
  
  for (i in order(count, decreasing = T)) {
    if (length(idx[[i]]) == 0) break
    
    if (length(idx[[i]]) > 1) {
      tmp.coor <- colMeans2(data[idx[[i]], ])
    } else {
      tmp.coor <- data[idx[[i]], ]
    }
    
    if (sum(cor(tmp.coor, t(top.pos.coor)) > 0.9) > 3) next
    
    top.pos <- c(top.pos, i)
    top.pos.coor <- rbind(top.pos.coor, tmp.coor)
    
    if (length(top.pos) >= 50) break
  }
  
  
  pred.tmp <- pred
  nontop.pos <- unique(Reduce(c, lapply((1:length(count))[-top.pos], function(i) idx[[i]])))
  top.pos.real <- unique(Reduce(c, lapply(top.pos, function(i) idx[[i]])))
  
  if (length(nontop.pos) > 0) {
    nontop.top.sim <- cor(t(matrix(data[nontop.pos, ], nrow = length(nontop.pos))), t(top.pos.coor))
    
    for (i in 1:nrow(nontop.top.sim)) {
      idx.tmp <- which.max(nontop.top.sim[i, ])
      idx.tmp <- idx[[top.pos[idx.tmp]]]
      if (length(idx.tmp) > 1) {
        sim.mat.tmp <- cor(data[nontop.pos[i], ], t(data[idx.tmp, ]))
        idx.tmp <- idx.tmp[order(sim.mat.tmp, decreasing = T)[1:min(5, length(idx.tmp)) ]]
        pred.tmp[nontop.pos[i], ] <- colMeans2(pred[idx.tmp, ] )
      } else {
        pred.tmp[nontop.pos[i], ] <- pred[idx.tmp, ] 
      }
    }
  }
  
  sc$pred <- pred.tmp
  
  sc
}


#' @importFrom igraph graph_from_adjacency_matrix mst distances
#' @importFrom stats quantile
#' @title scDHA pseudo time inference
#' @description Inferring pseudo-time data.
#' @param sc Embedding object, produced by \code{scDHA} function.
#' @param start.point Starting point of the trajectory.
#' @param ncores Number of processor cores to use.
#' @param seed Seed for reproducibility.
#' @return List with the following keys:
#' \itemize{
#' \item pt - Pseudo-time values for each sample.
#' }
#' @export
scDHA.pt <- function(sc = sc, start.point = 1, ncores = 15L, seed = NULL) {
  lat.idx <- which(sapply(sc$all.res, function(x) adjustedRandIndex(x, sc$cluster)) > 0.75)
  tmp.list <- lapply(lat.idx, function(i) sc$all.latent[[i]])
  
  t.final <- matrix(ncol = length(tmp.list), nrow = nrow(tmp.list[[1]]))
  counter <- 1
  
  set.seed(seed)
  
  for (x in tmp.list) {
    data <- x
    n <- nrow(data)
    adj <- 1 - cor(t(data))
    g <- graph_from_adjacency_matrix(adj, weighted = T, mode = "undirected")
    g <- mst(g)
    
    dis <- distances(g)
    
    dis[is.infinite(dis)] <- -1

    result <- start.point
    
    t <- dis[result,]
    
    for (cl in unique(sc$all.res[[lat.idx[counter]]])) {
      idx <- which(sc$all.res[[lat.idx[counter]]] == cl)
      tmp <- t[idx]
      tmp.max <- max(tmp)
      tmp.min <- min(tmp)
      tmp <- (quantile(tmp, 0.75) - quantile(tmp, 0.25)) * (tmp - min(tmp)) / (max(tmp) - min(tmp)) + quantile(tmp, 0.25)
      t[idx] <- tmp
    }
    
    t.final[, counter] <- t
    counter <- counter + 1
  }
  
  t.final <- rowMeans2(t.final)
  sc$pt <- t.final
  sc
}



#' @title scDHA classification
#' @description  Perform classification of new data based on available data.
#' @param train Expression matrix of available data, with rows represent samples and columns represent genes.
#' @param train.label A vector containing label for each sample in training data.
#' @param test Expression matrix new data for classification, with rows represent samples and columns represent genes.
#' @param ncores Number of processor cores to use.
#' @param seed Seed for reproducibility.
#' @return A vector contain classified labels for new data.
#' @export
scDHA.class <- function(train = train, train.label = train.label, test = test, ncores = 15L, seed = NULL) {
  data <- rbind(train, test)
  prob <- c(rep(1/nrow(train), nrow(train)), rep(1/nrow(test), nrow(test)))
  sc <- scDHA(data, ncores = ncores, do.clus = F, sample.prob = prob, seed = seed)
  latent <- sc$all.latent
  
  cl <- parallel::makeCluster(ncores, outfile = "/dev/null")
  doParallel::registerDoParallel(cl, cores = ncores)

  result <- foreach(x = latent) %dopar% {
    set.seed(seed)
    train.p <- x[1:length(train.label), ]
    test.p <- x[(length(train.label) + 1):nrow(x), ]
    
    tmp <- train.p
    tmp1 <- test.p
    
    clus.tmp <- train.label
    
    dis.tmp <- 1 - cor(t(tmp1), t(tmp))
    
    res <- rep(0, nrow(test.p))
    
    for (i in 1:nrow(dis.tmp)) {
      tmp2 <- order(dis.tmp[i, ])[1:10]
      tmp3 <- clus.tmp[tmp2]
      if (as.numeric(getmode1(tmp3)[2]) > 0) res[i] <- getmode1(tmp3)[1] else res[i] <- -1
    }
    res
  }
  parallel::stopCluster(cl)
  
  result <- data.frame(result)
  
  result <- apply(result, 1, function(x) getmode1(x)[1])
  
  result
}






















