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
      if (!is.null(MSTorder) & length(MSTorder) == 1) {
            stop("MSTorder is not a path!")
      }
      set.seed(12345)
      clucenter <- mclustobj$clucenter
      row.names(clucenter) <- paste0("clu",1:nrow(clucenter))
      clusterid <- mclustobj$clusterid             
      pcareduceres <- mclustobj$pcareduceres            
      adjmat <- get.adjacency(mclustobj$MSTtree)
      if (is.null(MSTorder)) {
            orderinMST <- 1
            allsp <- shortest.paths(mclustobj$MSTtree)
            longestsp <- which(allsp == max(allsp), arr.ind = T)
            MSTorder <- get.shortest.paths(mclustobj$MSTtree,longestsp[1,1],longestsp[1,2])$vpath[[1]]        
            if (flip) {
                  MSTorder <- rev(MSTorder)
            }
      } else {
            edgeinMST <- sapply(1:(length(MSTorder)-1),function(i) {
                  adjmat[MSTorder[i],MSTorder[i+1]]
            })
            if (sum(edgeinMST==0) > 0) {
                  orderinMST <- 0
            } else {
                  orderinMST <- 1
            }
      }          
      internalorderfunc <- function(internalorder,MSTinout) {
            TSCANorder <- NULL
            
            for (i in 1:(length(internalorder)-1)) {                  
                  currentcluid <- internalorder[i]
                  nextcluid <- internalorder[i + 1]
                  currentclucenter <- clucenter[currentcluid,]
                  nextclucenter <- clucenter[nextcluid,]
                  
                  currentreduceres <- pcareduceres[clusterid==currentcluid,]
                  if (MSTinout) {
                        connectcluid <- as.numeric(names(which(adjmat[currentcluid,] == 1)))      
                  } else {
                        if (i == 1) {
                              connectcluid <- nextcluid      
                        } else {
                              connectcluid <- c(nextcluid,internalorder[i - 1])
                        }                        
                  }
                  
                  cludist <- sapply(connectcluid, function(x) {                              
                        rowSums(sweep(currentreduceres,2,clucenter[x,],"-")^2)
                  })
                  mindistid <- apply(cludist,1,which.min)
                  
                  edgecell <- names(which(mindistid == which(connectcluid == nextcluid)))
                  
                  difvec <- nextclucenter - currentclucenter
                  tmppos <- pcareduceres[edgecell,] %*% difvec
                  pos <- as.vector(tmppos)
                  names(pos) <- row.names(tmppos)
                  TSCANorder <- c(TSCANorder,names(sort(pos)))  
                  
                  nextreduceres <- pcareduceres[clusterid==nextcluid,]     
                  if (MSTinout) {
                        connectcluid <- as.numeric(names(which(adjmat[nextcluid,] == 1)))
                  } else {
                        if (i == length(internalorder)-1) {
                              connectcluid <- currentcluid      
                        } else {
                              connectcluid <- c(currentcluid,internalorder[i + 2])
                        }                        
                  }
                  
                  cludist <- sapply(connectcluid, function(x) { 
                        rowSums(sweep(nextreduceres,2,clucenter[x,],"-")^2)
                  })
                  mindistid <- apply(cludist,1,which.min)
                  
                  edgecell <- names(which(mindistid == which(connectcluid == currentcluid)))
                  
                  difvec <- nextclucenter - currentclucenter
                  tmppos <- pcareduceres[edgecell,] %*% difvec
                  pos <- as.vector(tmppos)
                  names(pos) <- row.names(tmppos)
                  TSCANorder <- c(TSCANorder,names(sort(pos)))  
                  
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
      if (!orderinMST) {
            internalorderfunc(MSTorder,0)            
      } else {
            if (listbranch) {
                  branchnode <- setdiff(as.vector(V(mclustobj$MSTtree)),MSTorder)
                  allres <- list()
                  allres[[paste("backbone",paste(MSTorder,collapse = ','))]] <- internalorderfunc(MSTorder,1)  
                  for (tmpnode in branchnode) {
                        tmporder <- get.shortest.paths(mclustobj$MSTtree,MSTorder[1],tmpnode)$vpath[[1]] 
                        allres[[paste("branch:",paste(tmporder,collapse = ','))]] <- internalorderfunc(tmporder,1)
                  }
                  allres
            } else {
                  internalorderfunc(MSTorder,1)            
            }      
      }
}

