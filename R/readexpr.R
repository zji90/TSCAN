#' readexpr
#' 
#' read in the expression data and convert into a matrix
#'
#' This function reads in the expression data and convert into a matrix.
#'
#' @param file The path to the file. The file contains the gene expression matrix. Note that the first column must be gene names. The first row (header) must be the cell names.
#' @param sep how the file is separated. For csv file, it is "," separated. For txt file, it is usually "\t" or " " separated.
#' @return Matrix of gene expression.
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' \dontrun{
#' exprdata <- readexpr("expr.txt",sep="\t")
#' }

readexpr <- function(file, sep=",") {
  data <- read.table(file,as.is=T,header=T,sep=sep,row.names = 1)
  row.names(data) <- data[,1]
  as.matrix(data[,-1])
}
