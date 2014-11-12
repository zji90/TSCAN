#' TSCANorder
#' 
#' Construct TSCAN order after exprmclust
#'
#' This function takes the exact output of exprmclust function and construct TSCAN order by mapping all cells onto the path that connects cluster centers.
#' Users can also specify their own path.
#' 
#' @param mclustobj The exact output of the \code{\link{exprmclust}} function.
#' @param MSTorder A numeric vector specifying the order of clusters.
#' @param orderonly Only return the ordering. State or pseudotime information will not be returned
#' @param flip whether to flip the ordering
#' @param listbranch whether to list the ordering results of all possible branches
#' @return if orderonly = F, a vector of ordered cell names. if orderonly = T, a data frame of ordered cell names, cell states and pseudotime. 
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(lpsdata)
#' procdata <- preprocess(lpsdata)
#' lpsmclust <- exprmclust(procdata)
#' TSCANorder(lpsmclust)

TSCANorder <- function(mclustobj,MSTorder = NULL,orderonly=T,flip=F,listbranch=F) {
      set.seed(12345)
      if (is.null(MSTorder)) {
            allsp <- shortest.paths(mclustobj$MSTtree)
            longestsp <- which(allsp == max(allsp), arr.ind = T)
            MSTorder <- get.shortest.paths(mclustobj$MSTtree,longestsp[1,1],longestsp[1,2])$vpath[[1]]        
            if (flip) {
                  MSTorder <- rev(MSTorder)
            }
      }      
      internalorderfunc <- function(internalorder) {
            trimclu <- setdiff(as.vector(V(mclustobj$MSTtree)),internalorder)
            clucenter <- mclustobj$clucenter
            row.names(clucenter) <- paste0("clu",1:nrow(clucenter))
            clusterid <- mclustobj$clusterid 
            trimcell <- names(clusterid[clusterid %in% trimclu])
            pcareduceres <- mclustobj$pcareduceres
            pcareduceres <- pcareduceres[setdiff(row.names(pcareduceres),trimcell),]
            
            TSCANorder <- NULL
            
            for (i in 1:length(internalorder)) {
                  if (i == 1) {
                        currentcluid <- internalorder[i]
                        nextcluid <- internalorder[i + 1]
                        currentclucenter <- clucenter[currentcluid,]
                        nextclucenter <- clucenter[nextcluid,]
                        difvec <- nextclucenter - currentclucenter
                        tmppos <- pcareduceres[names(clusterid[clusterid==currentcluid]),] %*% difvec
                        pos <- as.vector(tmppos)
                        names(pos) <- row.names(tmppos)
                        TSCANorder <- c(TSCANorder,names(sort(pos)))                  
                  } else if (i == length(internalorder)) {
                        currentcluid <- internalorder[i]
                        lastcluid <- internalorder[i - 1]
                        currentclucenter <- clucenter[currentcluid,]
                        lastclucenter <- clucenter[lastcluid,]
                        difvec <- currentclucenter - lastclucenter
                        tmppos <- pcareduceres[names(clusterid[clusterid==currentcluid]),] %*% difvec
                        pos <- as.vector(tmppos)
                        names(pos) <- row.names(tmppos)
                        TSCANorder <- c(TSCANorder,names(sort(pos)))   
                  } else {
                        currentcluid <- internalorder[i]
                        nextcluid <- internalorder[i + 1]
                        lastcluid <- internalorder[i - 1]
                        currentclucenter <- clucenter[currentcluid,]
                        nextclucenter <- clucenter[nextcluid,]
                        lastclucenter <- clucenter[lastcluid,]
                        clupoints <- names(clusterid[clusterid==currentcluid])
                        distlast <- rowSums((pcareduceres[clupoints,]-lastclucenter)^2)
                        distnext <- rowSums((pcareduceres[clupoints,]-nextclucenter)^2)
                        lastpoints <- names(which(distlast < distnext))
                        nextpoints <- names(which(distlast >= distnext))
                        
                        difvec <- currentclucenter - lastclucenter
                        tmppos <- pcareduceres[lastpoints,] %*% difvec
                        pos <- as.vector(tmppos)
                        names(pos) <- row.names(tmppos)
                        TSCANorder <- c(TSCANorder,names(sort(pos)))  
                        
                        difvec <- nextclucenter - currentclucenter
                        tmppos <- pcareduceres[nextpoints,] %*% difvec
                        pos <- as.vector(tmppos)
                        names(pos) <- row.names(tmppos)
                        TSCANorder <- c(TSCANorder,names(sort(pos)))  
                  }
            }            
            if (orderonly) {
                  TSCANorder      
            } else {
                  datadist <- dist(mclustobj$pcareduceres)
                  distmat <- as.matrix(datadist)
                  alldist <- sapply(1:(length(TSCANorder)-1), function(x) {
                        distmat[TSCANorder[x],TSCANorder[x+1]]
                  })
                  ptime <- c(0,cumsum(alldist))
                  ptime <- ptime/max(ptime) * 100
                  data.frame(sample_name=TSCANorder,State=clusterid[TSCANorder],Pseudotime=ptime,stringsAsFactors = F)            
            }
      }
      
      if (listbranch) {
            branchnode <- setdiff(as.vector(V(mclustobj$MSTtree)),MSTorder)
            allres <- list()
            allres[["backbone"]] <- internalorderfunc(MSTorder)  
            for (tmpnode in branchnode) {
                  tmporder <- get.shortest.paths(mclustobj$MSTtree,MSTorder[1],tmpnode)$vpath[[1]] 
                  allres[[paste("branch:",paste(tmporder,collapse = ','))]] <- internalorderfunc(tmporder)
            }
            allres
      } else {
            internalorderfunc(MSTorder)            
      }
      
}

