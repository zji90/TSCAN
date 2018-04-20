#' orderscore
#' 
#' Calculate pseudotemporal ordering scores for orders
#'
#' This function calculates pseudotemporal ordering scores (POS) based on the sub-population information and order information given by users.
#' Cells should come from at least two cell sub-populations. These sub-population should be coded as 0,1,2,...
#'
#' @param subpopulation Data frame with two columns. First column: cell names. Second column: sub-population codes.
#' @param orders A list with various length containing pseudotime orderings.
#' @return a numeric vector of calculated POS.
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(lpsdata)
#' procdata <- preprocess(lpsdata)
#' subpopulation <- data.frame(cell = colnames(procdata), sub = ifelse(grepl("Unstimulated",colnames(procdata)),0,1), stringsAsFactors = FALSE)
#' lpsmclust <- exprmclust(procdata)
#' #Comparing default TSCAN ordering and tuned TSCAN ordering
#' order1 <- TSCANorder(lpsmclust,orderonly = T)
#' order2 <- TSCANorder(lpsmclust, c(1,2,3),orderonly = T)
#' orders <- list(order1,order2)
#' orderscore(subpopulation, orders)

orderscore <- function(subpopulation, orders) {
      subinfo <- subpopulation[,2]
      names(subinfo) <- subpopulation[,1]
      scorefunc <- function(order) {
            scoreorder <- subinfo[order]
            optscoreorder <- sort(scoreorder)
            optscore <- sum(sapply(1:(length(optscoreorder)-1),function(x) {
                  sum(optscoreorder[(x+1):length(optscoreorder)] - optscoreorder[x])
            })) 
            sum(sapply(1:(length(scoreorder)-1),function(x) {
                  sum(scoreorder[(x+1):length(scoreorder)] - scoreorder[x])
            })) / optscore
      }      
      sapply(orders,scorefunc)      
}
