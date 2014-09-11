#' TSPpseudotime
#' 
#' Construct pseudotime using Travelling Salesman Problem (TSP) algorithm
#'
#' This function first uses principal component analysis (PCA) to reduce dimensionality of original data.
#' If not specified by the user, the optimal dimension will be automatically selected by fitting 
#' a set of continuous piecewise regression to the standard deviations of first 20 principal components and choose 
#' the one with the smallest residual sum of squares. The distance matrix between cells are calculated based on the reduced data.
#' Then the function uses nearest insertion algorithm to construct a TSP path as a suboptimal solution. Because TSP path does not 
#' have a definite starting/ending point, users can specify a starting point, otherwise a random cell will be chosen as the starting point.
#' Users can also use the expression value of a gene to determine the optimal starting point. The gene expression value must change monotonically over the true biological process. 
#' K-means clustering will be used to determine the different stages of cell during the biological process.
#'
#' @param data The raw single_cell data, which is a numeric matrix or data.frame. Rows represent genes/features and columns represent single cells.
#' @param dim Either "auto" or a numeric value specifying the reduced dimenionality of PCA. If "auto" the optimal dimension will be chosen automatically. 
#' @param statenum Numeric value specifiying number of cell states in K-means clustering.
#' @param scale Whether to scale gene expressions across all cells to have zero mean and unit variance. Normally this argument should be TRUE.
#' @param startpoint Manually specify the starting point of TSP ordering. Should be one of the column names of data. Omitted when geneexpr is not null.
#' @param flip Logical value specifying whether to flip the ordering.
#' @param geneexpr Gene expression values used to determine optimal starting point. Gene expression values should change monotonically in the true biological process. Names should agree exactly with column names of data (can be of different order).
#' @param exprtrend Trend of gene expression values. Either increasing or decreasing.
#' @param maxtime Numeric value to specify the pseudotime of the ending point.
#' @param kmeansiter Number of iterations of K-means clustering. The function will automatically pick an optimal clustering.
#' @return a list containing
#' \itemize{
#'    \item pseudotime Data frame containing the constructed pseudotime information. First column: cell name. Second column: cell states. Third column: Pseudotime.
#'    \item reduceres Matrix storing the gene expression data after dimension reduction using PCA.
#' }
#' @export
#' @import TSP
#' @import HSMMSingleCell
#' @import Biobase
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @references Rosenkrantz, D. J., Stearns, R. E., & Lewis, II, P. M. (1977). An analysis of several heuristics for the traveling salesman problem. SIAM journal on computing, 6(3), 563-581.
#' @examples
#' library(HSMMSingleCell)
#' library(Biobase)
#' data(HSMM)
#' HSMMdata <- exprs(HSMM)
#' HSMMdata <- HSMMdata[,grep("T0|72",colnames(HSMMdata))]
#' procdata <- preprocess(HSMMdata)
#' #Choose MYOG gene expression as marker gene
#' MYOGexpr <- log2(HSMMdata["ENSG00000122180.4",]+1)
#' HSMMpseudotime <- TSPpseudotime(procdata, geneexpr = MYOGexpr, dim = 2)

TSPpseudotime <- function(data, dim = "auto", statenum = 3, scale = TRUE, startpoint = NULL, flip = FALSE, geneexpr = NULL, exprtrend = "increasing", maxtime = 100, kmeansiter = 10) {    
      
      if (!is.null(geneexpr)) {
            if (!identical(sort(names(geneexpr)),sort(colnames(data)))) {
                  stop("Names of geneexpr does not agree with names given in data!")
            } else {
                  if (exprtrend == "increasing") {
                        geneexpr <- geneexpr
                  } else if (exprtrend == "decreasing") {
                        geneexpr <- - geneexpr
                  } else {
                        stop("Invalid exprtrend argument!")
                  }
            }
      }
      
      genescore <- function(scoreorder) {
            sum(sapply(1:(length(scoreorder)-1),function(x) {
                  sum(geneexpr[scoreorder[(x+1):length(scoreorder)]] - geneexpr[scoreorder[x]])
            }))
      }
      
      if (scale) {
            scaledata <- t(apply(data,1,scale))      
            colnames(scaledata) <- colnames(data)
            pcres <- prcomp(t(scaledata), scale = T)
      } else {
            scaledata <- data
            pcres <- prcomp(t(scaledata), scale = F)
      }      
      

      if (dim == "auto") {
            sdev <- pcres$sdev[1:20]
            x <- 1:20
            rddim <- which.min(sapply(2:10, function(i) {
                  x2 <- pmax(0,x-i)
                  sum(lm(sdev~x+x2)$residuals^2)
            })) + 1
            message(paste("The selected optimal dimension is:",rddim))
      } else if (dim > 1) {
            rddim <- as.integer(dim)
      } else {
            stop('dim should either be "auto" or an integer greater than 1.')
      }
      
      reduceres <- t(scaledata) %*% pcres$rotation[,1:rddim] 
      dp <- as.matrix(dist(reduceres)) 
      
      set.seed(12345)
      dataTSP <- TSP(dp)
      datasolve <- solve_TSP(dataTSP)
      TSPorder <- labels(datasolve)
      
      if (!is.null(geneexpr)) {
            if (!is.null(startpoint))
                  message("Omitting startpoint argument. Using gene expression to select optimal starting point.")
            TSPres <- NULL
            for (flip in c(TRUE:FALSE)) {
                  for (start in 1:length(TSPorder)) {
                        if (start == 1) {
                              order <- TSPorder
                        } else {
                              order <- TSPorder[c(start:(length(TSPorder)),1:(start-1))]
                        }
                        if (flip) {
                              order <- rev(order)
                        }
                        TSPres <- rbind(TSPres,c(start,flip,genescore(order)))      
                  }
            }
            optstart <- TSPres[which.max(TSPres[,3]),1]
            optflip <- TSPres[which.max(TSPres[,3]),2]
            if (optstart == 1) {
                  finalorder <- TSPorder
            } else {
                  finalorder <- TSPorder[c(optstart:(length(TSPorder)),1:(optstart-1))]
            }
            if (optflip) {
                  finalorder <- rev(finalorder)
            }            
      } else {
            if (!is.null(startpoint)) {
                  selectstart <- which(TSPorder == startpoint)
                  if (selectstart == 1) {
                        finalorder <- TSPorder
                  } else {
                        finalorder <- TSPorder[c(selectstart:(length(TSPorder)),1:(selectstart-1))]
                  }
                  if (flip) {
                        finalorder <- rev(finalorder)
                  }                     
            } else {
                  message("Neither startpoint nor geneexpr is given. Starting point is chosen randomly.")
                  finalorder <- TSPorder
            }
      }
      alldist <- sapply(1:(length(finalorder)-1), function(x) {
            dp[finalorder[x],finalorder[x+1]]
      })
      ptime <- c(0,cumsum(alldist))
      ptime <- ptime/max(ptime) * maxtime
      optkmeansstate <- NULL
      optkmeanscost <- NULL
      for (i in 1:kmeansiter) {
            tmpstate <- kmeans(reduceres,centers = statenum)$cluster[finalorder]
            tmpcost <- sum(sapply(1:statenum,function(s) {
                  instatecell <- names(tmpstate[tmpstate==s])
                  sum(dp[instatecell,instatecell])                  
            }))
            if (i == 1){
                  optkmeansstate <- tmpstate
                  optkmeanscost <- tmpcost
            } else {
                  if (tmpcost < optkmeanscost) {
                        optkmeansstate <- tmpstate
                        optkmeanscost <- tmpcost                        
                  } 
            }
      }
      state <- optkmeansstate
      pseudotime <- data.frame(sample_name = finalorder, State = state, Pseudotime = ptime, stringsAsFactors = F)
      row.names(pseudotime) <- NULL
      list(pseudotime,reduceres)
      
}
