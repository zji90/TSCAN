#' Map cells to edges
#'
#' Map each cell to the closest edge on the MST, reporting also the distance to the corresponding vertices.
#'
#' @param x A numeric matrix of coordinates where each row represents a cell/sample and each column represents a dimension
#' (usually a PC or another low-dimensional embedding, but features or genes can also be used).
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object
#' containing such a matrix in its \code{\link{assays}}, as specified by \code{assay.type}.
#' This will be transposed prior to use.
#'
#' Alternatively, for \linkS4class{SingleCellExperiment}s, this matrix may be extracted from its \code{\link{reducedDims}},
#' based on the \code{use.dimred} specification.
#' In this case, no transposition is performed.
#' @inheritParams createClusterMST
#' @param clusters A factor-like object of the same length as \code{nrow(x)},
#' specifying the cluster identity for each cell in \code{x}.
#' This can also be \code{NULL}, see details.
#' @param mst A \link{graph} object containing a MST,
#' constructed from the same coordinate space as the values in \code{x} (e.g., same PC space, same set of features).
#'
#' @details
#' For each cluster, we consider all edges of the MST involving that cluster.
#' Each cell of that cluster is then mapped to the closest of these edges (where proximity is defined by Euclidean distance).
#' The identity of and distance from each ends of the edge is reported;
#' this can be useful for downstream pseudo-time calculations or to subset cells by those lying on a particular edge.
#' 
#' If \code{clusters=NULL}, each cell can be mapped to \emph{any} edge of the MST.
#' This is useful if the \code{mst} was constructed from a different set of cells than those in \code{x},
#' allowing us to effectively project new datasets onto an existing MST.
#' Note, however, that the new \code{x} must lie in the same coordinate space as the \code{x} used to make \code{mst}.
#'
#' Some cells may simply be mapped to the edge endpoints.
#' This manifests as values of zero for the distances from either end of the edge.
#' For analyses focusing on a specific edge, it may be advisable to filter out such cells 
#' as their edge assignments are arbitrarily assigned and they do not contribute to any transitional process along the edge.
#'
#' @return 
#' A \linkS4class{DataFrame} with one row per row of \code{x}, containing the fields:
#' \itemize{
#' \item \code{left.cluster}, the cluster on one end of the edge to which the cell was assigned.
#' \item \code{right.cluster}, the cluster on the other end of the edge to which the cell was assigned.
#' \item \code{left.distance}, the distance to the cluster centroid on one end.
#' \item \code{right.distance}, the distance to the cluster centroid on the other end.
#' }
#' Note that the sum of the distances will always yield the edge length.
#' 
#' @author Aaron Lun
#' @references
#' Ji Z and Ji H (2016).
#' TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis.
#' \emph{Nucleic Acids Res.} 44, e117
#'
#' @examples
#' # Mocking up a Y-shaped trajectory.
#' centers <- rbind(c(0,0), c(0, -1), c(1, 1), c(-1, 1))
#' rownames(centers) <- seq_len(nrow(centers))
#' clusters <- sample(nrow(centers), 1000, replace=TRUE)
#' cells <- centers[clusters,]
#' cells <- cells + rnorm(length(cells), sd=0.5)
#'
#' # Creating the MST first:
#' mst <- createClusterMST(cells, clusters=clusters)
#' plot(mst)
#'
#' # Mapping cells to the MST:
#' mapping <- mapCellsToEdges(cells, mst, clusters=clusters)
#' head(mapping)
#'
#' # Also works with some random extra cells:
#' extras <- matrix(rnorm(1000), ncol=2)
#' emapping <- mapCellsToEdges(extras, mst, clusters=NULL)
#' head(emapping)
#'
#' @seealso
#' \code{\link{quickPseudotime}}, a wrapper to quickly perform these calculations.
#' @name mapCellsToEdges
NULL

#################################################

#' @importFrom S4Vectors DataFrame
#' @importFrom igraph V
#' @importFrom Matrix which t
.map_cells_to_mst <- function(x, mst, clusters, columns=NULL) {
    # Getting the start and end of every edge.
    pairs <- which(mst[] > 0, arr.ind=TRUE)
    pairs <- pairs[pairs[,1] > pairs[,2],,drop=FALSE]
    vertices <- V(mst)
    vnames <- names(vertices)

    distance.to.edge <- matrix(Inf, nrow(x), nrow(pairs))
    left.gap <- right.gap <- matrix(NA_real_, nrow(x), nrow(pairs))
    if (!is.null(columns)) {
        x <- x[,columns,drop=FALSE]
    }

    # Computing distance of each point from each edge.
    for (i in seq_len(nrow(pairs))) {
        L <- pairs[i,1]
        R <- pairs[i,2]
        edge.left <- vertices$coordinates[[L]]
        edge.right <- vertices$coordinates[[R]]

        delta <- edge.right - edge.left
        edge.len <- sqrt(sum(delta^2))
        delta <- delta/edge.len

        if (is.null(clusters)) {
            candidates <- seq_len(nrow(x))
        } else {
            candidates <- which(clusters %in% vnames[c(L, R)])
        }

        centered <- t(t(x[candidates,,drop=FALSE]) - edge.left)
        proj <- as.numeric(centered %*% delta)
        proj <- pmax(0, pmin(proj, edge.len))
        mapped <- outer(proj, delta)
        dist <- sqrt(rowSums((centered - mapped)^2))

        distance.to.edge[candidates,i] <- dist
        left.gap[candidates,i] <- proj
        right.gap[candidates,i] <- edge.len - proj
    }

    # Mapping onto the closest valid edge per cell.
    m <- max.col(-distance.to.edge, ties.method="first")
    keep <- cbind(seq_along(m), m)
    DataFrame(
        left.cluster=vnames[pairs[m,1]],
        right.cluster=vnames[pairs[m,2]],
        left.distance=left.gap[keep],
        right.distance=right.gap[keep],
        row.names=rownames(x)
    )
}

#################################################

#' @export
#' @rdname mapCellsToEdges
setGeneric("mapCellsToEdges", function(x, ...) standardGeneric("mapCellsToEdges"))

#' @export
#' @rdname mapCellsToEdges
setMethod("mapCellsToEdges", "ANY", .map_cells_to_mst)

#' @export
#' @rdname mapCellsToEdges
#' @importFrom Matrix t
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("mapCellsToEdges", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .create_cluster_mst(t(assay(x, assay.type)), ...)
})

#' @export
#' @rdname mapCellsToEdges
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("mapCellsToEdges", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"), ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .create_cluster_mst(reducedDim(x, use.dimred), clusters=clusters, ...)
    } else {
        callNextMethod(x, clusters=clusters, ...)
    }
})
