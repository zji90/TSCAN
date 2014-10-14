#' singlegeneplot
#' 
#' plot expression values of individual genes against pseudotime axis
#'
#' This function plots the expression values of individual genes against given pseudotime
#'
#' @param geneexpr The gene expression values. Names should agree with the pseudotime information.
#' @param TSCANorder The output of function \code{\link{TSCANorder}}.
#' @param cell_size Size of cells in the plot.
#' @return ggplot2 object.
#' @export
#' @import ggplot2 mgcv
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(lpsdata)
#' procdata <- preprocess(lpsdata)
#' lpsmclust <- exprmclust(procdata)
#' lpsorder <- TSCANorder(lpsmclust,orderonly=FALSE,flip=TRUE)
#' #Choose STAT1 gene expression to plot
#' STAT2expr <- log2(lpsdata["STAT2",]+1)
#' singlegeneplot(STAT2expr, lpsorder)

singlegeneplot <- function(geneexpr, TSCANorder, cell_size = 2) {
      Pseudotime <- NULL #To overcome No visible binding for global variable Note in R CMD check
      geneexpr <- geneexpr[TSCANorder[,1]]
      exprdata <- cbind(TSCANorder, geneexpr)
      exprdata$State <- factor(exprdata$State)
      exprdata$predict <- fitted.values(mgcv::gam(geneexpr ~ s(Pseudotime,k=3),data=exprdata))
      q <- ggplot(aes(Pseudotime, geneexpr), data = exprdata)
      q <- q + geom_point(aes_string(color = "State"), size = I(cell_size))
      q <- q + geom_line(aes(Pseudotime, predict), data = exprdata)
      q <- q + ylab("Expression") + xlab("Pseudotime")
      q <- q + theme(strip.background = element_rect(colour = "white", fill = "white")) + 
            theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            theme(panel.background = element_rect(fill = "white"))
      q
}
