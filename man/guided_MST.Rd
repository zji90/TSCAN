% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/guided_MST.R
\name{guided_MST}
\alias{guided_MST}
\title{guided_MST}
\usage{
guided_MST(dr, clu)
}
\arguments{
\item{dr}{A numeric matrix of reduced dimensions. Can be from PCA, t-sne, UMAP, diffusion map, and many others. Each row is a cell and each column is a reduced dimension. Row names should be the cell names.}

\item{clu}{A named vector of cell clustering. Names are the cell names. Need to correspond to dr.}
}
\value{
A minimum-spanning-tree plot generated by igraph package
}
\description{
Minimum spanning tree to guide pseudotime construction
}
\details{
This function should be used when users have already performed PCA or other dimension reductions, and cell clustering. This function will plot a minimum-spanning-tree that helps users to determine which cell clusters should be linked in what order. The result can be fed into guided_tscan function to get the actual pseudotime.
}
\examples{
# see vignette
}
\author{
Zhicheng Ji <zhicheng.ji@duke.edu>
}
