#' @importFrom matrixStats colMeans2 rowSums2
#' @title scDHA visulization
#' @description  Generating 2D embeded data for visulation.
#' @param sc Embedding object produced by the \code{scDHA} function.
#' @param method Visualization method to use. It can be "UMAP" or "scDHA". The default setting is "UMAP".
#' @param ncores Number of processor cores to use.
#' @param seed Seed for reproducibility.
#' @return a list with the following keys:
#' \itemize{
#' \item pred - A matrix representing the 2D projection of single-cell data, where rows represent samples and columns represent latent components.
#' }
#' @examples
#' \donttest{
#' library(scDHA)
#' #Load example data (Goolam dataset)
#' data('Goolam'); data <- t(Goolam$data); label <- as.character(Goolam$label)
#' #Log transform the data 
#' data <- log2(data + 1)
#' if(torch::torch_is_installed()) #scDHA need libtorch installed
#' {
#'   #Generate clustering result, the input matrix has rows as samples and columns as genes
#'   result <- scDHA(data, ncores = 2, seed = 1)
#'   #Generate 2D representation, the input is the output from scDHA function
#'   result <- scDHA.vis(result, ncores = 2, seed = 1)
#'   #Plot the representation of the dataset, different colors represent different cell types
#'   plot(result$pred, col=factor(label), xlab = "scDHA1", ylab = "scDHA2")
#' }
#' }
#' @export

scDHA.vis <- function(sc = sc, method = "UMAP", ncores = 10L, seed = NULL)
{
  if(method == "UMAP")
  {
    set.seed(seed)
    data <- sc$latent
    pred <- uwot::umap(data, n_threads = ncores)
    sc$pred <- pred
    return(sc)
  } else if (method == "scDHA") {
    return(scDHA.vis.old(sc = sc, ncores = ncores, seed = seed))
  } else {
    stop("'method' should be one of 'UMAP', or 'scDHA'")
  }
}

scDHA.vis.old <- function(sc = sc, ncores = 10L, seed = NULL) {
  set.seed(seed)

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
      if (!is.null(seed))
      {
        set.seed((seed+i))
        torch_manual_seed((seed+i))
      } 
      
      torch::torch_set_num_threads(ifelse(nrow(data.en) < 1000, 1, 2))
      RhpcBLASctl::blas_set_num_threads(1)

      model <- scDHA_model_vis(ncol(data.en), ncol(data.list[[1]]$y1))

      model_loss <- function(x_or, x_pred) {
        x_pred <- torch_log(x_pred + 1)
        x_pred <- (x_pred - torch_mean(x_pred, dim = 1, keepdim = TRUE))/torch_std(x_pred, dim = 1, keepdim = TRUE, unbiased = FALSE)
        x_pred <- x_pred$t()
        x_pred <- x_pred*torch_abs(x_pred)/2
        xent_loss <- torch_exp( -x_pred )
        xent_loss <- xent_loss*torch_greater(x_or, 0)
        tmp <- torch_sum(xent_loss, dim = 2)
        xent_loss <- ( xent_loss$t()/tmp )$t()
        
        x_or <- torch_clamp(x_or, 1e-7, 1)
        xent_loss <- torch_clamp(xent_loss, 1e-7, 1)
        loss <- torch_sum(x_or*torch_log(x_or/xent_loss), dim = 2)
        loss
      }
      
      optimizer <- optim_adam(model$parameters, 5e-3, eps = 1e-4)
      max.ite <- ceiling(50/(length(data.list)**(2/3)) )
      
      for (ite in 1:max.ite) {
        for (subdata in data.list) {
          data <- subdata$data
          y1 <- subdata$y1
          sim.mat <- subdata$sim.mat
          
          idxs <- 1:nrow(data)
          idxs <- createFolds(idxs, min(10, ceiling(nrow(data)/ 50)))
          
          optimizer$zero_grad()
          for (idx in idxs) {
            
            tmp.sim.mat <- t(scale(t(sim.mat[idx, idx])))
            sim.mat.dis <- tmp.sim.mat/2*abs(tmp.sim.mat)
            sim.mat.dis.exp <- exp(-sim.mat.dis)
            diag(sim.mat.dis.exp) <- 0
            
            output <- model(torch_tensor(data[idx, ]))
            loss <- torch_mean(model_loss(torch_tensor(sim.mat.dis.exp / rowSums2(sim.mat.dis.exp)), output[[1]])) + nnf_binary_cross_entropy(output[[2]], torch_tensor(y1[idx, ]), reduction = "mean")
            if(as.numeric(torch_isnan(loss)) == 0)
            {
              loss$backward()
              optimizer$step()
            }
            optimizer$zero_grad()
          }
        }
      }
      
      data <- data.en
      
      model$eval()
      with_no_grad({
        pred <- as.matrix(model$encode_latent(torch_tensor(data))) 
      })
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
  
  top.pos <- order(count, decreasing = TRUE)[1]
  top.pos.coor <- matrix(colMeans2(data[idx[[top.pos]], ]), nrow = 1)
  
  for (i in order(count, decreasing = TRUE)) {
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
        idx.tmp <- idx.tmp[order(sim.mat.tmp, decreasing = TRUE)[1:min(5, length(idx.tmp)) ]]
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
#' @examples
#' \donttest{
#' library(scDHA)
#' #Load example data (Goolam dataset)
#' data('Goolam'); data <- t(Goolam$data); label <- as.character(Goolam$label)
#' #Log transform the data 
#' data <- log2(data + 1)
#' if(torch::torch_is_installed()) #scDHA need libtorch installed
#' {
#'   #Generate clustering result, the input matrix has rows as samples and columns as genes
#'   result <- scDHA(data, ncores = 2, seed = 1)
#'   #Cell stage order in Goolam dataset
#'   cell.stages <- c("2cell", "4cell", "8cell", "16cell", "blast")
#'   #Generate pseudo-time for each cell, the input is the output from scDHA function
#'   result <- scDHA.pt(result, start.point = 1, ncores = 2, seed = 1)
#'   #Calculate R-squared value 
#'   r2 <- round(cor(result$pt, as.numeric(factor(label, levels = cell.stages)))^2, digits = 2)
#' }
#' }
#' @export
scDHA.pt <- function(sc = sc, start.point = 1, ncores = 10L, seed = NULL) {
  RhpcBLASctl::blas_set_num_threads(min(ncores, 4))
  lat.idx <- which(sapply(sc$all.res, function(x) adjustedRandIndex(x, sc$cluster)) > 0.75)
  if(length(lat.idx) == 0) lat.idx <- which(sapply(sc$all.res, function(x) adjustedRandIndex(x, sc$cluster)) > 0.5)
  tmp.list <- lapply(lat.idx, function(i) sc$all.latent[[i]])
  
  if(nrow(tmp.list[[1]]) <= 5000)
  {
    set.seed(seed)
    all.res <- sc$all.res
    tmp.list.or <- tmp.list
  } else {
    set.seed(seed)
    idx.all <- sample.int(nrow(tmp.list[[1]]), 4900)
    if(start.point %in% idx.all)
    {
      start.point <- which(idx.all == start.point)
    } else {
      idx.all <- c(start.point, idx.all)
      start.point <- 1
    }
    all.res <- sc$all.res
    tmp.list.or <- tmp.list
    
    tmp.list <- lapply(tmp.list, function(x) x[idx.all,])
    all.res <- lapply(all.res, function(x) x[idx.all])
  }
  
  t.final <- matrix(ncol = length(tmp.list), nrow = nrow(tmp.list[[1]]))
  counter <- 1
  
  for (x in tmp.list) {
    data <- x
    n <- nrow(data)
    adj <- 1 - cor(t(data))
    g <- graph_from_adjacency_matrix(adj, weighted = TRUE, mode = "undirected")
    g <- mst(g)
    
    dis <- distances(g)
    
    dis[is.infinite(dis)] <- -1
    
    result <- start.point
    
    t <- dis[result,]
    
    for (cl in unique(all.res[[lat.idx[counter]]])) {
      idx <- which(all.res[[lat.idx[counter]]] == cl)
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
  
  if(nrow(tmp.list.or[[1]]) > 5000)
  {
    tmp <- rep(0, nrow(tmp.list.or[[1]]))
    tmp[idx.all] <- t.final
    
    latent <- sc$latent
    latent.t.final <- latent[idx.all, ]
    
    for (i in seq(nrow(tmp.list.or[[1]]))[-idx.all]) {
      d <- 1 - cor(latent[i, ], t(latent.t.final))
      tmp[i] <- mean(t.final[order(d)[1:3]])
    }
    t.final <- tmp
  }
  
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
#' @examples
#' \donttest{
#' library(scDHA)
#' #Load example data (Goolam dataset)
#' data('Goolam'); data <- t(Goolam$data); label <- as.character(Goolam$label)
#' #Log transform the data 
#' data <- log2(data + 1)
#' #Split data into training and testing sets
#' set.seed(1)
#' idx <- sample.int(nrow(data), size = round(nrow(data)*0.75))
#' train.x <- data[idx, ]; train.y <- label[idx]
#' test.x <- data[-idx, ]; test.y <- label[-idx]
#' if(torch::torch_is_installed()) #scDHA need libtorch installed
#' {
#'   #Predict the labels of cells in testing set
#'   prediction <- scDHA.class(train = train.x, train.label = train.y, test = test.x, 
#'                             ncores = 2, seed = 1)
#'   #Calculate accuracy of the predictions
#'   acc <- round(sum(test.y == prediction)/length(test.y), 2)
#'   print(paste0("Accuracy = ", acc))
#' }
#' }
#' @export
scDHA.class <- function(train = train, train.label = train.label, test = test, ncores = 10L, seed = NULL) {
  data <- rbind(train, test)
  prob <- c(rep(1/nrow(train), nrow(train)), rep(1/nrow(test), nrow(test)))
  sc <- scDHA(data, ncores = ncores, do.clus = FALSE, sample.prob = prob, seed = seed)
  latent <- sc$all.latent
  
  cl <- parallel::makeCluster(ncores, outfile = "/dev/null")
  doParallel::registerDoParallel(cl, cores = ncores)
  
  x <- NULL
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






















