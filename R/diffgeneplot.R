#' singlegeneplot
#' 
#' plot expression values of individual genes against pseudotime axis
#'
#' This function plots the expression values of individual genes against given pseudotime
#'
#' @param geneexpr The gene expression values. Names should agree with the pseudotime information.
#' @param pseudotime The pseudotime information. It is typically the first element of the return value of function \code{\link{TSPpseudotime}}.
#' @param cell_size Size of cells in the plot.
#' @return ggplot2 object.
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @seealso \code{\link{TSPpseudotime}} for examples

singlegeneplot <- function(geneexpr, pseudotime, cell_size = 2) {
      geneexpr <- geneexpr[pseudotime[,1]]
      exprdata <- cbind(pseudotime, geneexpr)
      exprdata$State <- factor(exprdata$State)
      exprdata$predict <- fitted.values(gam(geneexpr ~ s(Pseudotime),data=exprdata))
      q <- ggplot(aes(Pseudotime, geneexpr), data = exprdata)
      q <- q + geom_point(aes_string(color = "State"), size = I(cell_size))
      q <- q + geom_line(aes(Pseudotime, predict), data = exprdata)
      q <- q + ylab("Expression") + xlab("Pseudotime")
      q <- q +  theme(strip.background = element_rect(colour = "white", fill = "white")) + 
            theme(panel.border = element_blank(), axis.line = element_line()) + 
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            theme(panel.background = element_rect(fill = "white"))
      q
}
