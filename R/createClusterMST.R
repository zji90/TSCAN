#' Minimum spanning trees on cluster centroids
#'
#' Build a MST where each node is a cluster centroid and 
#' each edge is weighted by the Euclidean distance between centroids.
#' This represents the most parsimonious explanation for a particular trajectory
#' and has the advantage of being directly intepretable with respect to any pre-existing clusters.
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
#'
#' Alternatively, if \code{clusters=NULL}, a numeric matrix of coordinates for cluster centroids,
#' where each row represents a cluster and each column represents a dimension 
#' Each row should be named with the cluster name.
#' This mode can also be used with assays/matrices extracted from SummarizedExperiments and SingleCellExperiments. 
#' @param ... For the generic, further arguments to pass to the specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method
#' (if \code{use.dimred} is specified) or the ANY method (otherwise).
#' @param clusters A factor-like object of the same length as \code{nrow(x)},
#' specifying the cluster identity for each cell in \code{x}.
#' If \code{NULL}, \code{x} is assumed to already contain coordinates for the cluster centroids.
#' @param columns A character, logical or integer vector specifying the columns of \code{x} to use.
#' If \code{NULL}, all provided columns are used by default. 
#' @param outgroup A logical scalar indicating whether an outgroup should be inserted to split unrelated trajectories.
#' Alternatively, a numeric scalar specifying the distance threshold to use for this splitting.
#' @param outscale A numeric scalar specifying the scaling to apply to the median distance between centroids
#' to define the threshold for outgroup splitting.
#' Only used if \code{outgroup=TRUE}.
#' @param assay.type An integer or string specifying the assay to use from a SummarizedExperiment \code{x}.
#' @param use.dimred An integer or string specifying the reduced dimensions to use from a SingleCellExperiment \code{x}.
#'
#' @section Introducing an outgroup:
#' If \code{outgroup=TRUE}, we add an outgroup to avoid constructing a trajectory between \dQuote{unrelated} clusters.
#' This is done by adding an extra row/column to the distance matrix corresponding to an artificial outgroup cluster,
#' where the distance to all of the other real clusters is set to \eqn{\omega/2}.
#' Large jumps in the MST between real clusters that are more distant than \eqn{\omega} will then be rerouted through the outgroup,
#' allowing us to break up the MST into multiple subcomponents by removing the outgroup.
#'
#' The default \eqn{\omega} value is computed by constructing the MST from the original distance matrix,
#' computing the median edge length in that MST, and then scaling it by \code{outscale}.
#' This adapts to the magnitude of the distances and the internal structure of the dataset
#' while also providing some margin for variation across cluster pairs.
#' Alternatively, \code{outgroup} can be set to a numeric scalar in which case it is used directly as \eqn{\omega}.
#'
#' @section Confidence on the edges:
#' For the MST, we obtain a measure of the confidence in each edge by computing the distance gained if that edge were not present.
#' Ambiguous parts of the tree will be less penalized from deletion of an edge, manifesting as a small distance gain.
#' In contrast, parts of the tree with clear structure will receive a large distance gain upon deletion of an obvious edge.
#'
#' For each edge, we divide the distance gain by the length of the edge to normalize for cluster resolution.
#' This avoids overly penalizing edges in parts of the tree involving broad clusters
#' while still retaining sensitivity to detect distance gain in overclustered regions.
#' As an example, a normalized gain of unity for a particular edge means that its removal
#' requires an alternative path that increases the distance travelled by that edge's length.
#'
#' The normalized gain is reported as the \code{"gain"} attribute in the edges of the MST from \code{\link{createClusterMST}}.
#' Note that the \code{"weight"} attribute represents the edge length.
#' 
#' @return A \link{graph} object containing an MST computed on \code{centers}.
#' Each node corresponds to a cluster centroid and has a numeric vector of coordinates in the \code{coordinates} attribute.
#' The edge weight is set to the Euclidean distance and the confidence is stored as the \code{gain} attribute.
#'
#' @author Aaron Lun
#'
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
#' # Creating the MST:
#' mst <- createClusterMST(cells, clusters)
#' plot(mst)
#'
#' # We could also do it on the centers:
#' mst2 <- createClusterMST(centers, clusters=NULL)
#' plot(mst2)
#'
#' # Works if the expression matrix is in a SE:
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(t(cells), colData=DataFrame(group=clusters))
#' mst3 <- createClusterMST(se, se$group, assay.type=1)
#' plot(mst3)
#' 
#' @name createClusterMST
NULL

#################################################

#' @importFrom igraph graph.adjacency minimum.spanning.tree delete_vertices E V V<-
#' @importFrom stats median dist
.create_cluster_mst <- function(x, clusters, outgroup=FALSE, outscale=3, columns=NULL) {
    if (!is.null(columns)) {
        x <- x[,columns,drop=FALSE]                
    }

    if (!is.null(clusters)) {
        x <- rowmean(x, clusters)
    }

    centers <- as.matrix(x)
    dmat <- dist(centers)
    dmat <- as.matrix(dmat)

    if (!isFALSE(outgroup)) {
        if (!is.numeric(outgroup)) {
            g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
            mst <- minimum.spanning.tree(g)
            med <- median(E(mst)$weight)
            outgroup <- med * outscale
        }

        old.d <- rownames(dmat)

        # Divide by 2 so rerouted distance between cluster pairs is 'outgroup'.
        dmat <- rbind(cbind(dmat, outgroup/2), outgroup/2) 
        diag(dmat) <- 0

        special.name <- strrep("x", max(nchar(old.d))+1L)
        rownames(dmat) <- colnames(dmat) <- c(old.d, special.name)
    }

    g <- graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
    mst <- minimum.spanning.tree(g)
    mst <- .estimate_edge_confidence(mst, g)

    if (!isFALSE(outgroup)) {
        mst <- delete_vertices(mst, special.name)
    }

    # Embed vertex coordinates for downstream use.
    coord.list <- vector("list", nrow(centers))
    names(coord.list) <- rownames(centers)
    for (r in rownames(centers)) {
        coord.list[[r]] <- centers[r,]
    }
    V(mst)$coordinates <- coord.list[names(V(mst))]

    mst 
}

#' @importFrom igraph minimum.spanning.tree E E<- ends get.edge.ids delete.edges
.estimate_edge_confidence <- function(mst, g) {
    edges <- E(mst)
    ends <- ends(mst, edges)
    reweight <- numeric(length(edges))

    for (i in seq_along(edges)) {
        id <- get.edge.ids(g, ends[i,])        
        g.copy <- delete.edges(g, id)
        mst.copy <- minimum.spanning.tree(g.copy)
        reweight[i] <- sum(E(mst.copy)$weight)
    }

    W <- edges$weight
    total <- sum(W)
    offset <- min(W)
    E(mst)$gain <- (reweight - total)/(W + offset/1e8)
    mst
}

#################################################

#' @export
#' @rdname createClusterMST
setGeneric("createClusterMST", function(x, ...) standardGeneric("createClusterMST"))

#' @export
#' @rdname createClusterMST
setMethod("createClusterMST", "ANY", .create_cluster_mst)

#' @export
#' @rdname createClusterMST
#' @importFrom Matrix t
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("createClusterMST", "SummarizedExperiment", function(x, ..., assay.type="logcounts") {
    .create_cluster_mst(t(assay(x, assay.type)), ...)
})

#' @export
#' @rdname createClusterMST
#' @importFrom SingleCellExperiment reducedDim
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("createClusterMST", "SingleCellExperiment", function(x, clusters=colLabels(x, onAbsence="error"), ..., use.dimred=NULL) {
    if (!is.null(use.dimred)) {
        .create_cluster_mst(reducedDim(x, use.dimred), clusters=clusters, ...)
    } else {
        callNextMethod(x, clusters=clusters, ...)
    }
})
