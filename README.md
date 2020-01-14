TSCAN: Tools For Single-Cell ANalysis
====

## Overview
TSCAN is a software tool developed to better support in silico pseudo-Time reconstruction in Single-Cell RNA-seq ANalysis. A major advantage of TSCAN over many other similar methods is that it allows users to supply their own dimension reduction and cell clustering results to TSCAN (see exprmclust function), and tune the pseudotime ordering along the tree. This high flexibility is critical for scRNA-seq data where there is no universal best way of processing data.

TSCAN uses a cluster-based minimum spanning tree (MST) approach to order cells. Cells are first grouped into clusters and an MST is then constructed to connect cluster centers. Pseudo-time is obtained by projecting each cell onto the tree, and the ordered sequence of cells can be used to study dynamic changes of gene expression along the pseudo-time. Clustering cells before MST construction reduces the complexity of the tree space. This often leads to improved cell ordering. It also allows users to conveniently adjust the ordering based on prior knowledge.

## TSCAN Online User Interface
TSCAN user interface can be directly launched online without installing any software package: https://zhiji.shinyapps.io/TSCAN. 

## TSCAN Installation

TSCAN software can be installed via Github.
Users should have R installed on their computer before installing TSCAN. R can be downloaded here: http://www.r-project.org/.
To install the latest version of TSCAN package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("zji90/TSCAN")
```
To launch user interface after installation, run following commands in R:
```{r }
library(TSCAN)
TSCANui()
```
For users with R programming experience, command line tools are also available in TSCAN R package. Please check the manual package included in the package for details.

## TSCAN Datasets
The example datasets for TSCAN paper can be downloaded on Github: https://github.com/zji90/TSCANdata

## Citation
Please cite the following paper:
Zhicheng Ji and Hongkai Ji. TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis. (2016) Nucleic Acids Research, 44(13):e117.


## Contact the Author
Author: Zhicheng Ji, Hongkai Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zji4@jhu.edu)

Or open a new issue on this Github page
