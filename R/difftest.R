#' difftest
#' 
#' testing differentially expressed genes
#'
#' This function tests whether a gene is significantly expressed given pseudotime ordering. Generalized additive model (GAM) with user-specified
#' degrees of freedoms is compared with a constant fit to get the p-values. The p-values are adjusted for multiple testing using fdr to gain qvalues.
#'
#' @param data The raw single_cell data, which is a numeric matrix or data.frame. Rows represent genes/features and columns represent single cells.
#' @param pseudotime The pseudotime information. It is typically the first element of the return value of function \code{\link{TSPpseudotime}}.
#' @param df Numeric value specifying the degree of freedom used in the GAM model.
#' @return Data frame containing pvalues and qvalues of testing differentially expression.
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @seealso \code{\link{TSPpseudotime}} for examples

difftest <- function(data, pseudotime, df = 3) {      
      row.names(pseudotime) <- pseudotime[,1]
      ptime <- pseudotime[colnames(data),3]
      pval <- apply(data, 1, function(x) {
            anova(gam(x~s(ptime,k=df)),lm(x~1))$s.table[4]            
      })      
      qval <- p.adjust(pval, method = "fdr")
      data.frame(pval = pval,qval = qval)
}
