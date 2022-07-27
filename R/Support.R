#' @import doParallel
#' @importFrom stats kmeans
#' @importFrom foreach %dopar% foreach %do%
#' @importFrom methods is

clus <- function(data, k = NULL, nmax = 10) #, ncore = 2
{
  if (is.null(k)) k <- nclusterPar(data, nmax = nmax) #, ncore = ncore
  if (nrow(data) < 1e3 & (k < 5))
  {
    k <- kmeans(data, k, nstart = 1e3, iter.max = 1e3)
    k$cluster
  } else
  {
    kknn <- try(specClust(data, k, nn = 7))
    while (is(kknn, "try-error")) {
      k <- k +1
      kknn <- try(specClust(data, k, nn = 7))
    }
    kknn$cluster
  }

}

clus.big <- function(data, k = NULL, n = 2e3, nmax = 10) #, ncore = 2
{
  ind <- sample.int(nrow(data), n)
  ind1 <- (1:nrow(data))[-ind]
  tmp <- data[ind,]
  tmp1 <- data[-ind,]

  clus.tmp <- clus(tmp,k, nmax = nmax)
  
  nn.tmp <- matrix(ncol = 10, nrow = nrow(tmp1))
  if(nrow(tmp1) > 10e3)
  {
    folds <- round(seq(1, nrow(tmp1), length.out = ceiling(nrow(tmp1)/10e3)))
  } else {
    folds <- c(1, nrow(tmp1)) 
  }
  for (i in 2:length(folds)) {
    dis.tmp <- 1 - cor(t(tmp1[folds[i-1]:folds[i], ]), t(tmp))
    for (j in 1:nrow(dis.tmp)) {
      nn.tmp[folds[i-1] - 1 + j,   ] <- order(dis.tmp[j,])[1:10]
    }
  }
  
  res <- rep(0, nrow(data))
  res[ind] <- clus.tmp
  
  for (i in 1:nrow(nn.tmp)) {
    tmp2 <- nn.tmp[i, ]
    tmp3 <- clus.tmp[tmp2]
    res[ind1[i]] <- getmode(tmp3)
    
  }
  res
  
}

clus.louvain <- function(data)
{
  if(nrow(data) <= 10e3)
  {
    n <- nrow(data)
    nn <- 30
    clus.tmp <- membership(phenograph(data, nn)[[2]])
    while (length(clus.tmp) != n) {
      nn <- nn + 10
      clus.tmp <- membership(phenograph(data, nn)[[2]])
    }
    res <- clus.tmp
  } else {
    n <- 9900
    ind <- sample.int(nrow(data), n)
    ind1 <- (1:nrow(data))[-ind]
    tmp <- data[ind,]
    tmp1 <- data[-ind,]
    
    nn <- 30
    clus.tmp <- membership(phenograph(tmp, nn)[[2]])
    while (length(clus.tmp) != n) {
      nn <- nn + 10
      clus.tmp <- membership(phenograph(tmp, nn)[[2]])
    }
    
    nn.tmp <- matrix(ncol = 10, nrow = nrow(tmp1))
    if(nrow(tmp1) > 10e3)
    {
      folds <- round(seq(1, nrow(tmp1), length.out = ceiling(nrow(tmp1)/10e3)))
    } else {
      folds <- c(1, nrow(tmp1)) 
    }
    for (i in 2:length(folds)) {
      dis.tmp <- 1 - cor(t(tmp1[folds[i-1]:folds[i], ]), t(tmp))
      for (j in 1:nrow(dis.tmp)) {
        nn.tmp[folds[i-1] - 1 + j,   ] <- order(dis.tmp[j,])[1:10]
      }
    }
    
    res <- rep(0, nrow(data))
    res[ind] <- clus.tmp
    
    for (i in 1:nrow(nn.tmp)) {
      tmp2 <- nn.tmp[i, ]
      tmp3 <- clus.tmp[tmp2]
      res[ind1[i]] <- getmode(tmp3)
    }
  }
  
  res
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

nclusterPar <- function(data, nmax = 10) 
{
  j <- NULL
  result <- foreach (j = 1:10) %do% {
    set.seed(j)
    idx <- sample(1:nrow(data),min(500, nrow(data)))
    to.test <- matrix(0, nrow = nmax, ncol = 3)
    e <- 3
    for (i in 2:nmax) {

      kknn <- try(specClust(data[idx,], i, nn = 7))
      if (!is(kknn, "try-error"))
      {
        to.test[i,1] <- kknn$betweenss/kknn$totss
        to.test[i,2] <- kknn$tot.withinss
      } else {
        e <- i + 2
      }


    }
    for (i in e:nmax) {
      to.test[i,3] <- (to.test[i,2]-to.test[i-1,2])/to.test[i-1,2]
    }
    c(which.max(to.test[,1]), which.max(to.test[,3]) )
  }
  result <- t(data.frame(result))
  floor(mean(result)+0.5)

}

clustercom2 <- function(result)
{
  test <- matrix(0, ncol = length(result$all), nrow = length(result$all))
  for (i in seq(length(result$all))) {
    for (j in seq(length(result$all))) {
      if (i != j)
        test[i, j] <- adjustedRandIndex(result$all[[i]], result$all[[j]])
    }
  }
  for (i in seq(length(result$all))) {
    test[i, i] <- mean(test[-i, i])
  }
  found <- FALSE
  if (sum(test < 0.7) > 0) {
    i <- 2
  } else {
    i <- 1
  }
  while (!found) {
    k <- kmeans(test, i, nstart = 100, iter.max = 5e3)$cluster
    max <- 0
    for (c in unique(k)) {
      score <- mean(test[which(k == c), which(k == c)])
      if (score > max & length(which(k == c)) > 1) {
        max <- score
        idx <- which(k == c)
      }
    }
    if (max > 0.8)
      found <- TRUE
    if (i > 3)
      found <- TRUE
    i <- i + 1
  }
  res <- t(data.frame(result$all))
  res <- res[idx, ]
  cl.max <- floor(mean(apply(res, 1, function(x) length(unique(x)))) + 0.5)
  res <- data.frame(result$all)
  da <- apply(res, 1, paste, collapse = "#")
  indUnique <- which(!duplicated(da))
  indAll <- match(da, da[indUnique])
  da <- res[indUnique, ]
  if (length(indUnique) < 1e3)
  {
    n <- 1e3 - length(indUnique)
    idx <- sample(1:nrow(res), n ,replace = TRUE)
    tmp <- res[idx,]
    da <- rbind(da,tmp)
  }
  test <- wMetaC(da, (cl.max + 1), hmethod = "ward.D")
  test$finalC[indAll]
}

adjustedRandIndex <- function (x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

find.vec <- function(node, graph, w, l, p, q) {
  result <- matrix(ncol = l, nrow = w)
  pre.node <- NULL
  c.node <- rep(node, w)
  for (j in 1:l) {
    neighbor <- graph[[2]][c.node, ]
    weight <- graph[[1]][c.node, ]
    weight <- weight * q / rowSums2(weight)
    if (is.null(pre.node)) {
      pre.node <- c.node
      c.node <- sample(neighbor[1, ], w, prob = weight[1, ], replace = TRUE)
    } else {
      tmp.node <- c.node
      c.node <- sapply(1:w, function(x) {
        sample(c(pre.node[x], neighbor[x, ]), 1, prob = c(p, weight[x, ]))
      })
      pre.node <- tmp.node
    }
    result[, j] <- c.node
  }
  
  
  result <- as.vector(result)
  result
}

#' @useDynLib scDHA
#' @importFrom Rcpp evalCpp 
#' @importFrom RcppParallel RcppParallelLibs
calG <- function(data, k = 10, distance = 2)
{
  tmp <- DC_para(data, k = k)
  
  D <- sqrt(tmp$D)
  
  C <- tmp$C
  
  D[D==0] <- min(D[D!=0])
  
  D <- 1/D
  
  
  list(D = D, C = C)
  
}


getmode1 <- function(v) {
  uniqv <- unique(v)
  c(uniqv[which.max(tabulate(match(v, uniqv)))], max(tabulate(match(v, uniqv))))
}


to_categorical <- function(cluster) {
  if(min(cluster) == 0) cluster <- cluster + 1
  tmp <- matrix(0, ncol = max(cluster), nrow = length(cluster))
  unique_cl <- seq(max(cluster))
  for (i in 1:length(unique_cl)) {
    tmp[which(cluster == unique_cl[i]) , i] <- 1
  }
  tmp
}