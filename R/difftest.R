#' difftest
#' 
#' testing differentially expressed genes
#'
#' This function tests whether a gene is significantly expressed given pseudotime ordering. Likelihood ratio test is performed to compare a generalized additive model (GAM) with a constant fit to get the p-values. The p-values are adjusted for multiple testing by fdr.
#'
#' @param data The raw single_cell data, which is a numeric matrix or data.frame. Rows represent genes/features and columns represent single cells.
#' @param order A character vector of pseudotime cell ordering. Must be within the column names of \code{data}.
#' @param df Numeric value specifying the degree of freedom used in the GAM model.
#' @return Data frame containing pvalues and FDRs of testing differential expression.
#' @export
#' @import mgcv
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(lpsdata)
#' procdata <- preprocess(lpsdata)
#' lpsorder <- TSCANorder(exprmclust(procdata),orderonly=TRUE)
#' diffval <- difftest(procdata[1:10,],lpsorder)
#' #Selected differentially expressed genes under qvlue cutoff of 0.05
#' row.names(diffval)[diffval$qval < 0.05]

difftest <- function(data, order, df = 3) {   
      ptime <- 1:length(order)
      pval <- apply(data[,order,drop=F], 1, function(x) {
            if (sum(x) == 0) {
                  1
            } else {
                  model <- mgcv::gam(x~s(ptime,k=3))
                  pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual,lower.tail = F)                  
            }             
      })      
      fdr <- p.adjust(pval, method = "fdr")
      res <- data.frame(pval = pval,FDR = fdr)
      res[order(res[,2],res[,1]),]
}
