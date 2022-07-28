#' @import cluster
#' @importFrom utils combn
#' @import clusterCrit
#' @importFrom stats as.dist cutree hclust median
#' @import Matrix

wMetaC <- function(nC, hmethod, enN.cluster, minN.cluster, maxN.cluster, sil.thre,
                   height.Ntimes) {
  # This is to obtain the weight matrix for each cluster solution for following
  # meta-clustering
  N <- nrow(nC)  #number of points
  C <- ncol(nC)  #number of clustering methods/times; or K
  AA <- Reduce("+", apply(nC, 2, getA))  #sum over obtained matrices; execute along the column and then matrix sum
  AA <- AA/C
  indA <- Matrix::which(AA != 0, arr.ind = TRUE)  #find non-zero indices of AA
  nd <- vapply(AA[indA], function(x) x * (1 - x), numeric(1))
  newAA <- sparseMatrix(i = indA[, 1], j = indA[, 2], x = nd, dims = c(N, N))
  w0 <- 4/N * Matrix::rowSums(newAA)  #the weight for each point
  e <- 0.01
  w1 <- (w0 + e)/(1 + e)  #adjusted point weight
  x <- as.vector(vapply(seq(C), function(i) {
    paste(nC[, i], "_", i, sep = "")
  }, rep("a", nrow(nC))))  #convert the matrix (N*C) to vector (concatenating them)
  newnC <- matrix(x, nrow = N, byrow = FALSE)  #reshape a vector to a matrix; by column
  R <- unique(x)  #all unique labels
  allC <- length(R)  #number of all unique labels
  cb <- combn(allC, 2)  #all possible combinations (n*(n-1)/2)
  # alls = apply(cb, 2, getss, R = R, x = x, w1 = w1)  #calculate the weight s for all combinations, slow
  allk <- lapply(R, function(r) {
    d <- as.numeric(unlist(strsplit(r, "_")))
    which(nC[,d[2]] == d[1])
  })
  alls <- apply(cb, 2, getss1, R = R, x = x, w1 = w1, allk = allk)  #calculate the weight s for all combinations
  S0 <- sparseMatrix(i = cb[1, ], j = cb[2, ], x = alls, dims = c(allC, allC))  #triangle part of the S
  S <- S0 + t(S0) + diag(allC)
  if (missing(sil.thre)) {
    sil.thre <- 0
  }
  hres <- get_opt_hclust(S, hmethod, N.cluster = enN.cluster, minN.cluster, maxN.cluster,
                         sil.thre, height.Ntimes)  #solely using the silhouette index as the criteria
  tf <- hres$f
  v <- hres$v
  
  newnC[] <- vapply(newnC, function(q) tf[match(q, R)], numeric(1))  #apply to every element; reorganizing the clusters for different results
  finalC <- apply(newnC, 1, function(d) names(sort(table(d), decreasing = TRUE)[1]))  #find the most repeated elements for each row
  N.cluster <- length(unique(finalC))  #note that the number of clusters for meta-clustering is not determined by previous selection, but by the unique number in the final round.
  perc <- 0.5
  if (N.cluster == 1) {
    # better not to have only one cluster
    finalC <- apply(newnC, 1, function(d) {
      x <- sort(table(d), decreasing = TRUE)[seq(2)]
      n0 <- length(x[1])
      if (x[2] >= n0 * perc) {
        y <- names(x[2])
      } else {
        y <- names(x[1])
      }
      return(y)
    })
    N.cluster <- length(unique(finalC))
  }
  # For ease of visualization
  uC <- unique(finalC)  #unique clusters
  
  y0 <- apply(newnC, 1, function(q) {
    t <- rep(0, N.cluster)
    for (i in c(seq(N.cluster))) {
      t[i] <- length(which(q %in% uC[i]))
    }
    return(t)
  })#need to reorganize before counting
  y0 <- t(y0)#transpose
  x0 <- matrix(0, nrow = N, ncol = N.cluster)
  tw <- 0.5
  # print(uC)
  for (i in seq(N)) {
    xind <- which(finalC[i] == uC)
    x0[i, xind] <- 1  #the correct clustering result
    allind <- which(y0[i, ] != 0)  #all the counts
    diffind <- setdiff(allind, xind)  #some other counts which are not the correct cluster
    if (length(diffind) != 0) {
      x0[i, diffind] <- tw * y0[i, diffind]/y0[i, xind]  #use a reduced weight
    }
  }
  out <- list()  #declare
  out$finalC <- finalC
  out$x0 <- x0
  return(out)
}

getA <- function(rowColor) {
  # This is to obtain the weighted co-association matrix for clustering solution
  # rowColor
  N <- length(rowColor)  #number of points
  
  L <- levels(factor(rowColor))
  
  # find indices for each cluster, then all combinations of indices
  tmp <- vapply(L, function(k) {
    r <- which(rowColor %in% k)
    expand.grid(r, r)
  }, c(list(1), list(1)))
  # reshape to the indices
  allind <- matrix(unlist(t(tmp)), ncol = 2, byrow = FALSE)  #need transpose
  A <- sparseMatrix(i = allind[, 1], j = allind[, 2], x = 1, dims = c(N, N))  #non-zero entries
  return(A)
}

getss1 <- function(pind, R, x, w1, allk) {
  # This is to get the element of S
  intset <- intersect(allk[[pind[1]]], allk[[pind[2]]])  #set intersection
  ss <- 0
  if (length(intset) != 0) {
    uset <- union(allk[[pind[1]]], allk[[pind[2]]])  #set union
    ss <- sum(w1[intset])/sum(w1[uset])
  }
  return(ss)
}

getss <- function(pind, R, x, w1) {
  # This is to get the element of S
  pairk <- lapply(pind, getnewk, R = R, x = x, N = length(w1))  #run for two indices
  intset <- intersect(unlist(pairk[1]), unlist(pairk[2]))  #set intersection
  ss <- 0
  if (length(intset) != 0) {
    uset <- union(unlist(pairk[1]), unlist(pairk[2]))  #set union
    ss <- sum(w1[intset])/sum(w1[uset])
  }
  return(ss)
}

getnewk <- function(k, R, x, N) {
  # This is to get the original index of the sample
  k1 <- which(x %in% R[k])  #find samples with k-th cluster
  d1 <- unlist(strsplit(R[k], "_"))  #the name contains only two parts; get the numbering part
  d <- as.numeric(Matrix::tail(d1, n = 1))  #the last element of the split arrays
  newk1 <- k1 - (d - 1) * N  #the index
  return(newk1)
}

get_opt_hclust <- function(mat, hmethod, N.cluster, minN.cluster, maxN.cluster, sil.thre,
                           height.Ntimes) {
  # if no agglomeration method for hierarchical clustering is provided
  if (missing(hmethod) || is.null(hmethod)) {
    hmethod <- "ward.D"  #the default hierarchical clustering agglomeration method is 'ward.D'
  }
  # if no minimum number of clusters is provided
  if (missing(minN.cluster) || is.null(minN.cluster)) {
    minN.cluster <- 2  #by default, we try the minimum number of clusters starting from 2
  }
  # if no maximum number of clusters is provided
  if (missing(maxN.cluster) || is.null(maxN.cluster)) {
    maxN.cluster <- 40  #by default, we try the maximum number of clusters as large as 40 or the number of cells minus 1, whichever is smaller.
  }
  # if no threshold for the maximum Silhouette index is provided
  if (missing(sil.thre) || is.null(sil.thre)) {
    sil.thre <- 0.35  #by default, we use 0.35 to determine whether we use Silhouette index as the criteria to determine the optimal number of clusters
  }
  # if no threshold for the height difference is provided
  if (missing(height.Ntimes) || is.null(height.Ntimes)) {
    height.Ntimes <- 2  #by default, we select the first height which is (height.Ntimes) times larger than the immediate consecutive height
  }
  # just use simple criteria to determine whether they are feature vectors or
  # similarity matrix, and then we use different ways to measure the distance
  if (Matrix::isSymmetric(mat)) {
    # symmmetric matrix
    d <- as.dist(1 - mat)
    flag1 <- 1
  } else {
    d <- as.dist(1 - cor(t(mat)))
    flag1 <- 0
  }
  h <- hclust(d, method = hmethod)  #ward to ward.D
  # if N.cluster is given, we simply use the given N.cluster for hierarchical
  # clustering
  #     if (!missing(N.cluster) && is.numeric(N.cluster)) {
  if (is.numeric(N.cluster)) {
    if (!is.numeric(N.cluster)) {
      stop("The given N.cluster is not a numeric!")
    }
    if (N.cluster%%1 != 0) {
      stop("The given N.cluster is not an integer!")
    }
    if (N.cluster < 2) {
      stop("The given N.cluster is less than 2, which is not suitable for clustering!")
    }
    
    v <- cutree(h, k = N.cluster)  #for different numbers of clusters
    f <- v  #the optimal clustering results
    sil <- silhouette(v, d)
    msil <- median(sil[, 3])
    ch0 <- intCriteria(data.matrix(mat), as.integer(v), "Calinski_Harabasz")
    CHind <- unlist(ch0, use.names = FALSE)  #convert a list to a vector/value
    optN.cluster <- N.cluster
  }
  hres <- list()
  hres$f <- f#optimal clustering results
  hres$v <- v#different numbers of clustering results
  hres$maxsil <- max(msil)
  hres$msil <- msil
  hres$CHind <- CHind
  hres$height <- h$height
  hres$optN.cluster <- optN.cluster
  return(hres)
}