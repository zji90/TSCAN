% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exprmclust.R
\name{exprmclust}
\alias{exprmclust}
\title{exprmclust}
\usage{
exprmclust(
  data,
  clustermethod = "mclust",
  clusternum = 2:9,
  modelNames = "VVV",
  reduce = T,
  cluster = NULL
)
}
\arguments{
\item{data}{The raw single_cell data, which is a numeric matrix or data.frame. Rows represent genes/features and columns represent single cells.}

\item{clustermethod}{Either 'mclust' (model-based clustering) or 'kmeans' (k-means clustering). If 'kmeans', clusternum must be specified.}

\item{clusternum}{For mclust, an integer vector specifying all possible cluster numbers. The best cluster number will be picked using BIC. The minimum value should be two. For kmeans, the number of clusters.}

\item{modelNames}{model to be used in model-based clustering. By default "ellipsoidal, varying volume, shape, and orientation" is used.}

\item{reduce}{Whether to perform the PCA on the expression data.}

\item{cluster}{A vector of user specified clustering results. Must be integers starting from 1 with no gap. THe name of the vector should be the same as the column names of the data. If not null then clusternum and modelNames will both be ignored.}
}
\value{
if more than one cluster detected, a list containing
\itemize{
   \item pcareduceres Numeric matrix containing the transformed expression values after PCA.
   \item MSTtree igraph object which is the result of constructing MST.
   \item clusterid A named vector specifying which cluster the cells belong to.
   \item clucenter Numeric matrix of the cluster centers.
}
if only one cluster detected, a list containing
\itemize{
   \item pcareduceres Numeric matrix containing the transformed expression values after PCA.
}
}
\description{
Perform model-based clustering on expression values
}
\details{
By default, this function first uses principal component analysis (PCA) to reduce dimensionality of original data. 
It then performs model-based clustering on the transformed expression values.
A minimum-spanning-tree is constructed to link the cluster centers.
The clustering results will be used for TSCAN ordering.
}
\examples{
data(lpsdata)
procdata <- preprocess(lpsdata)
exprmclust(procdata)

userclust <- sample(1:2,ncol(lpsdata),replace = T)
names(userclust) <- colnames(procdata)
exprmclust(procdata,cluster=userclust)
}
\references{
Fraley, C., & Raftery, A. E. (2002). Model-based clustering, discriminant analysis, and density estimation. Journal of the American Statistical Association, 97(458), 611-631.
}
\author{
Zhicheng Ji, Hongkai Ji <zji4@zji4.edu>
}
