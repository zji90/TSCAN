#' Compute pseudotimes from the MST
#'
#' Compute a pseudotime for each cell lying on each path through the MST from a given starting node.
#'
#' @param mapping A \linkS4class{DataFrame} of MST-mapping information for each cell,
#' usually obtained by running \code{\link{mapCellsToEdges}} with the per-cell coordinate matrix and \code{mst}.
#' @param mst A \link{graph} object containing a MST from \code{\link{createClusterMST}}.
#' @param start A character vector specifying the starting node from which to compute pseudotimes in each component of \code{mst}.
#' Defaults to an arbitrarily chosen node of degree 1 or lower in each component.
#'
#' @details
#' The pseudotimes are returned as a matrix where each row corresponds to cell in \code{x} 
#' and each column corresponds to a path through the MST from \code{start} to all nodes of degree 1.
#' (If \code{start} is itself a node of degree 1, then paths are only considered to all other such nodes.)
#' This format is inspired by that from the \pkg{slingshot} package and provides a compact representation of branching events.
#'
#' Each branching event in the MST results in a new path and thus a new column in the pseudotime matrix.
#' For any given row in this matrix, entries are either \code{NA} or they are identical.
#' This reflects the fact that multiple paths will share a section of the MST for which the pseudotimes are the same.
#' 
#' The starting node in \code{start} is \emph{completely arbitrarily chosen} by \code{orderClusterMST},
#' as directionality is impossible to infer from the expression matrix alone.
#' However, it is often possible to use prior biological knowledge to pick an appropriate cluster as the starting node.
#'
#' @return 
#' A numeric matrix containing the pseudotimes of all cells (rows) across all paths (columns) through \code{mst}.
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
#' # Creating the MST and mapping the cells.
#' mst <- createClusterMST(cells, clusters=clusters)
#' mapping <- mapCellsToEdges(cells, mst, clusters=clusters)
#'
#' # Obtaining pseudo-time orderings.
#' ordering <- orderCells(mapping, mst)
#' unified <- rowMeans(ordering, na.rm=TRUE)
#' plot(cells[,1], cells[,2], col=topo.colors(21)[cut(unified, 21)], pch=16)
#'
#' @seealso
#' \code{\link{quickPseudotime}}, a wrapper to quickly perform these calculations.
#'
#' @export
#' @importFrom igraph V degree adjacent_vertices components E get.edge.ids
orderCells <- function(mapping, mst, start=NULL) {
    comp <- components(mst)$membership
    by.comp <- split(names(comp), comp)
    if (is.null(start)) {
        # Choosing one of the terminal nodes within each component.
        candidates <- names(V(mst)[degree(mst) <= 1])
        start <- vapply(by.comp, function(b) intersect(b, candidates)[1], "")
    } else {
        start <- as.character(start)
        for (b in by.comp) {
            if (length(intersect(b, start))!=1) {
                stop("'start' must have one cluster in each component of 'mst'")
            }
        }
    }

    collated <- list()
    latest <- start
    nstarts <- length(latest)
    parents <- rep(NA_character_, nstarts)
    progress <- rep(list(rep(NA_real_, nrow(mapping))), nstarts)
    cumulative <- numeric(nstarts)

    # Rolling through the tree. This is effectively a recursive algorithm
    # being flattened into an iterative one for simplicity.
    while (length(latest)) {
        new.latest <- new.parents <- character(0)
        new.progress <- list()
        new.cumulative <- numeric(0)

        for (i in seq_along(latest)) {
            curnode <- latest[i]
            curparent <- parents[i]
            all.neighbors <- names(adjacent_vertices(mst, curnode, mode="all")[[1]])

            all.children <- setdiff(all.neighbors, curparent)
            if (length(all.children)==0) {
                collated[[curnode]] <- progress[[i]]
                next
            }

            cum.dist <- cumulative[i]
            if (!is.na(curparent)) {
                edge.id <- get.edge.ids(mst, c(curnode, curparent))
                cum.dist <- cum.dist + E(mst)$weight[edge.id]
            }

            collected.progress <- list()
            for (child in all.children) {
                sofar <- progress[[i]] # yes, the 'i' here is deliberate.

                cur.cells.1 <- mapping$left.cluster == curnode & mapping$right.cluster == child
                sofar[cur.cells.1] <- mapping$left.distance[cur.cells.1] + cum.dist

                cur.cells.2 <- mapping$right.cluster == curnode & mapping$left.cluster == child
                sofar[cur.cells.2] <- mapping$right.distance[cur.cells.2] + cum.dist

                collected.progress[[child]] <- sofar
            }

            new.latest <- c(new.latest, all.children)
            new.parents <- c(new.parents, rep(curnode, length(all.children)))
            new.progress <- c(new.progress, collected.progress)
            new.cumulative <- c(new.cumulative, rep(cum.dist, length(all.children)))
        }

        latest <- new.latest
        parents <- new.parents
        progress <- new.progress
        cumulative <- new.cumulative
    }
    
    do.call(cbind, collated)
}
