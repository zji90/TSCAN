#' Compute column means of a matrix based on a grouping variable
#'
#' A utility that is equivalent to \code{\link{rowsum}}, but computes the mean instead of the sum of each subset of rows.
#'
#' @param x A numeric matrix or matrix-like object.
#' @param group A vector or factor specifying the group assignment for each row of \code{x}.
#' @param ... Further arguments to pass to \code{\link{rowsum}}, e.g., \code{reorder}.
#'
#' @return A numeric matrix with one row per level of \code{group},
#' where the value for each column contains the average value across the subset of rows corresponding that level.
#'
#' @author Aaron Lun
#' 
#' @examples
#' x <- matrix(runif(100), ncol = 5)
#' group <- sample(1:8, 20, TRUE)
#' (xmean <- rowmean(x, group))
#'
#' @export
#' @importFrom DelayedArray rowsum DelayedArray
rowmean <- function(x, group, ...) {
    x <- rowsum(DelayedArray(x), group, ...) 
    tab <- table(group)
    x / as.integer(tab[rownames(x)])
}
