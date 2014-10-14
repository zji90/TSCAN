#' exprmclust
#' 
#' Perform model-based clustering on expression values
#'
#' By default, this function first uses principal component analysis (PCA) to reduce dimensionality of original data. 
#' It then performs model-based clustering on the transformed expression values.
#' A minimum-spanning-tree is constructed to link the cluster centers.
#' The clustering results will be used for TSCAN ordering.
#'
#' @param data The raw single_cell data, which is a numeric matrix or data.frame. Rows represent genes/features and columns represent single cells.
#' @param clusternum A numeric vector specifying all possible cluster numbers. The best cluster number will be picked using BIC.
#' @param modelNames model to be used in model-based clustering. By default "ellipsoidal, varying volume, shape, and orientation" is used.
#' @param reduce Whether to perform the PCA on the expression data.
#' @return if more than one cluster detected, a list containing
#' \itemize{
#'    \item pcareduceres Numeric matrix containing the transformed expression values after PCA.
#'    \item MSTtree igraph object which is the result of constructing MST.
#'    \item clusterid A named vector specifying which cluster the cells belong to.
#'    \item clucenter Numeric matrix of the cluster centers.
#' }
#' if only one cluster detected, a list containing
#' \itemize{
#'    \item pcareduceres Numeric matrix containing the transformed expression values after PCA.
#' }
#' @export
#' @importFrom mclust Mclust mclustBIC
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @references Fraley, C., & Raftery, A. E. (2002). Model-based clustering, discriminant analysis, and density estimation. Journal of the American Statistical Association, 97(458), 611-631.
#' @examples
#' data(lpsdata)
#' procdata <- preprocess(lpsdata)
#' exprmclust(procdata)

exprmclust <- function(data, clusternum = 1:9, modelNames="VVV", reduce = T) { 
      set.seed(12345)
      if (reduce) {
            sdev <- prcomp(t(data),scale=T)$sdev[1:20]
            x <- 1:20
            optpoint <- which.min(sapply(2:10, function(i) {
                  x2 <- pmax(0,x-i)
                  sum(lm(sdev~x+x2)$residuals^2)
            }))
            pcadim = optpoint + 1
            
            tmpdata <- t(apply(data,1,scale))
            colnames(tmpdata) <- colnames(data)
            tmppc <- prcomp(t(tmpdata),scale=T)
            pcareduceres <- t(tmpdata) %*% tmppc$rotation[,1:pcadim]      
      } else {
            pcareduceres <- t(data)
      }
      
      res <- suppressWarnings(Mclust(pcareduceres,G=clusternum,modelNames="VVV"))
      if (res$G != 1) {
            clusterid <- apply(res$z,1,which.max)
            clucenter <- matrix(0,ncol=ncol(pcareduceres),nrow=res$G)
            for (cid in 1:res$G) {
                  clucenter[cid,] <- colMeans(pcareduceres[names(clusterid[clusterid==cid]),,drop=F])
            }
            
            dp <- as.matrix(dist(clucenter))                                            
            gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
            dp_mst <- minimum.spanning.tree(gp)
            list(pcareduceres=pcareduceres,MSTtree=dp_mst,clusterid=clusterid,clucenter=clucenter)      
      } else {
            warning("Only one cluster detected!")
            list(pcareduceres=pcareduceres)
      }
      
}