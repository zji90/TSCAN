#' Sinlge-cell RNA-seq data for human skeletal muscle myoblasts (HSMM)
#'
#' The dataset contains 11135 rows and 118 columns. Each row represent a gene and each column represent a single cell. This dataset is a subset of single-cell RNA-seq data for HSMM provided by HSMMsinglecell package. Only cells collected at time point 0 and time point 72 are retained for the purpose of demonstration. Genes which have raw expression values of greater than zero in more than 25% cells are retained. For the original dataset please refer to HSMMSingleCell package.
#'
#' @format A data frame with 11135 rows and 118 variables
#' @references Cole Trapnell and Davide Cacchiarelli et al (2014): The dynamics and regulators of cell fate decisions are revealed by pseudo-temporal ordering of single cells. Nature Biotechnology
#' @source \url{http://www.bioconductor.org/packages/devel/data/experiment/html/HSMMSingleCell.html}
#' @name HSMMdata
NULL
