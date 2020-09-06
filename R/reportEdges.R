#' Report edge coordinates
#'
#' Provides the coordinates of the start and end of every edge.
#' This is mostly useful for plotting purposes in \code{\link{segments}} or the equivalent \pkg{ggplot2} functionality.
#'
#' @inheritParams createClusterMST
#' @param mst A \link{graph} object containing a MST, typically the output of \code{createClusterMST(x, clusters)}.
#' @param combined Logical scalar indicating whether a single data.frame of edge coordinates should be returned.
#'
#' @return 
#' A data.frame containing the start and end coordinates of segments representing all the edges in \code{mst}.
#' If \code{combined=FALSE}, a list of two data.frames is returned where corresponding rows represent the start and end coordinates of the same edge.
#'
#' @details
#' \code{x} may also be \code{NULL}, in which case the center coordinates are obtained 
#' from the \code{coordinates} vertex attribute of \code{mst}.
#'
#' @author Aaron Lun
#' @references
#' Ji Z and Ji H (2016).
#' TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.
#' \emph{Nucleic Acids Res.} 44, e117
#'
#' @examples
#' # Re-using the example from ?createClusterMST.
#' example(createClusterMST, ask=FALSE) 
#'
#' # Plotting the MST on top of existing visualizations:
#' edges <- reportEdges(x=NULL, mst, combined=FALSE)
#' plot(cells[,1], cells[,2], col=clusters)
#' segments(edges$start$dim1, edges$start$dim2, edges$end$dim1, 
#'      edges$end$dim2, lwd=5)
#'
#' # Use with coordinates other than those used to make the MST:
#' shifted.cells <- cells + 10
#'
#' shift.edges <- reportEdges(shifted.cells, mst, 
#'     clusters=clusters, combined=FALSE)
#' plot(shifted.cells[,1], shifted.cells[,2], col=clusters)
#' segments(shift.edges$start$dim1, shift.edges$start$dim2, 
#'     shift.edges$end$dim1, shift.edges$end$dim2, lwd=5)
#'
#' # Also works for ggplot2:
#' df <- data.frame(shifted.cells, cluster=factor(clusters))
#' shift.edges2 <- reportEdges(shifted.cells, mst, clusters=clusters)
#' 
#' library(ggplot2)
#' ggplot(df) +
#'    geom_point(aes(x=X1, y=X2, color=cluster)) + 
#'    geom_line(data=shift.edges2, mapping=aes(x=dim1, y=dim2, group=edge))
#' 
#' @name reportEdges
NULL

#################################################

#' @importFrom Matrix which
#' @importFrom igraph V
.connect_cluster_mst <- function(x, mst, clusters, combined=TRUE) {
    pairs <- which(mst[] > 0, arr.ind=TRUE)
    pairs <- pairs[pairs[,1] > pairs[,2],,drop=FALSE]

    vertices <- V(mst)
    vnames <- names(vertices)
    group <- paste0(vnames[pairs[,1]], "--", vnames[pairs[,2]])

    if (is.null(x)) {
        x <- do.call(rbind, vertices$coordinates)
        rownames(x) <- vnames
    } else if (!is.null(clusters)) {
        x <- rowmean(x, clusters)
    } 
    if (!identical(sort(rownames(x)), sort(vnames))) {
        stop("cluster names in 'x' or 'clusters' should be identical to those in 'mst'")
    }
    if (is.null(colnames(x))) {
        colnames(x) <- sprintf("dim%i", seq_len(ncol(x)))
    }

    L <- vnames[pairs[,1]]
    R <- vnames[pairs[,2]]
    L <- data.frame(edge=group, x[L,,drop=FALSE])
    R <- data.frame(edge=group, x[R,,drop=FALSE])
    rownames(L) <- rownames(R) <- NULL

    if (combined) {
        rbind(L, R)
    } else {
        list(start=L, end=R)
    }
}

#################################################

#' @export
#' @rdname reportEdges
setGeneric("reportEdges", function(x, ...) standardGeneric("reportEdges"))

#' @export
#' @rdname reportEdges
setMethod("reportEdges", "ANY", .connect_cluster_mst)

#' @export
#' @rdname reportEdges
#' @importFrom Matrix t
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("reportEdges", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .connect_cluster_mst(t(assay(x, assay.type)), ...)
})

#' @export
#' @rdname reportEdges
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("reportEdges", "SingleCellExperiment", function(x, ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .connect_cluster_mst(reducedDim(x, use.dimred), ...)
    } else {
        callNextMethod()
    }
})
