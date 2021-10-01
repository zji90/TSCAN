#' guided_tscan
#' 
#' Guided pseudotime construction using TSCAN
#'
#' This function should be used when users have already performed PCA or other dimension reductions, and cell clustering. Users also need to know which cell clusters should be linked in what order, which depends either on prior biological knowledge or obtained by guided_MST in an unsupervised fashion. This function fits well with any standard analysis pipeline such as Seurat, and can even accompany other pseudotime analysis methods such as diffusion maps. The biggest advantage is its flexibility. Users can supply any dimension reduction result, any cell clustering, and any order cell clusters should be linked. But it's users' responsibility to make sure the constructed pseudotime is biologically meaningful and interpreted in the correct way. 
#' 
#' @param dr A numeric matrix of reduced dimensions. Can be from PCA, t-sne, UMAP, diffusion map, and many others. Each row is a cell and each column is a reduced dimension. Row names should be the cell names.
#' @param clu A named vector of cell clustering. Names are the cell names. Need to correspond to dr.
#' @param cluorder A vector indicating the order of cell clusters to be linked. Can be a subset of all clusters. Values in cluorder must be included in clu.
#' @return A character vector of cell pseudotime ordering.
#' @export
#' @author Zhicheng Ji <zhicheng.ji@@duke.edu>
#' @examples
#' # see vignette

guided_tscan <- function(dr,clu,cluorder) {
      n <- names(clu)
      clu <- as.character(clu)
      names(clu) <- n
      cm <- t(sapply(unique(clu),function(i) colMeans(dr[clu==i,])))
      
      left <- right <- list()
      for (i in 1:(length(cluorder))) {
            if (i==1) {
                  right[[cluorder[i]]] <- names(clu)[clu==cluorder[i]]
            } else if (i==length(cluorder)) {
                  left[[cluorder[i]]] <- names(clu)[clu==cluorder[i]]
            } else {
                  tmpcell <- names(clu)[clu==cluorder[i]]
                  leftdist <- colSums((t(dr[tmpcell,]) - cm[cluorder[i-1],])^2)
                  rightdist <- colSums((t(dr[tmpcell,]) - cm[cluorder[i+1],])^2)
                  left[[cluorder[i]]] <- tmpcell[leftdist <= rightdist]
                  right[[cluorder[i]]] <- tmpcell[leftdist > rightdist]
            }
      }
      
      ord <- NULL
      for (i in 1:(length(cluorder)-1)) {                  
            difvec <- cm[cluorder[i+1],] - cm[cluorder[i],]
            tmppos <- dr[c(right[[cluorder[i]]],left[[cluorder[i+1]]]),,drop=F] %*% difvec
            pos <- as.vector(tmppos)
            names(pos) <- row.names(tmppos)
            ord <- c(ord,names(sort(pos)))  
      }
      ord
}
