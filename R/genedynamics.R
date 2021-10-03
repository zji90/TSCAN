#' genedynamics
#' 
#' plot gene expression values of individual genes along pseudotime
#'
#' This function plots the expression values of individual genes against given pseudotime
#'
#' @param geneexpr A named gene expression vector. Names are cell names.
#' @param order A character vector of pseudotime cell ordering. Must be within the column names of \code{data}.
#' @param k degree of gam model. Smaller k makes the curve smoother.
#' @return ggplot2 object.
#' @export
#' @import ggplot2 mgcv
#' @author Zhicheng Ji <zhicheng.ji@@duke.edu>
#' @examples
#' data(lpsdata)
#' procdata <- preprocess(lpsdata)
#' lpsmclust <- exprmclust(procdata)
#' lpsorder <- TSCANorder(lpsmclust,orderonly=TRUE,flip=TRUE)
#' #Choose STAT1 gene expression to plot
#' STAT2expr <- log2(lpsdata["STAT2",]+1)
#' singlegeneplot(STAT2expr, lpsorder)

genedynamics <- function(geneexpr, order, k=3) {
      pt <- 1:length(order)
      geneexpr <- geneexpr[order]
      pred <- fitted.values(mgcv::gam(geneexpr ~ s(pt,k=3)))
      ggplot() + geom_point(data=data.frame(pt=pt,e=geneexpr),aes(x=pt,y=e),color='royalblue',alpha=0.5) + geom_line(data=data.frame(x=pt,y=pred),aes(x=pt,y=y),color='orange',size=2) + theme_classic() + xlab('pseudotime') + ylab('gene expression')
}
