TSCAN: Tools for Single-Cell ANalysis
====
      
## Overview
TSCAN is a novel software tool to facilitate pseudo-time reconstruction in Single-Cell RNA-seq Analysis. TSCAN uses a cluster-based minimum spanning tree (MST) approach to order cells. Cells are first grouped into clusters and an MST is then constructed to connect cluster centers. Pseudo-time is obtained by projecting each cell onto the tree, and the ordered sequence of cells can be used to study dynamic changes of gene expression along the pseudo-time. Clustering cells before MST construction reduces the complexity of the tree space. This often leads to improved pseudo-time reconstruction. It also brings convenience for users to interactively adjust the analysis by incorporating prior knowledge. TSCAN has a graphical user interface (GUI) to support data visualization and user interaction.

## TSCAN Online User Interface
TSCAN user interface can be directly launched online without installing any software package: https://zhiji.shinyapps.io/TSCAN. PLEASE NOTE: Currently the online version only allows one concurrent user. If the online user interface shows "please wait" for a long time, probably another user is using the online interface and please come back at another time. Users are recommended to install TSCAN on their own computers with following procedures.

## TSCAN Installation

TSCAN software can be installed via Github (recommended) and Bioconductor. 
Users should have R installed on their computer before installing TSCAN. R can be downloaded here: http://www.r-project.org/.

### Install  via Github (Recommended)
To install the latest version of TSCAN package via Github, run following commands in R:
```{r }
if (!require("devtools"))
      install.packages("devtools")
devtools::install_github("TSCAN","zji90")
```
To launch user interface after installation, run following commands in R:
```{r }
library(TSCAN)
TSCANui()
```
For users with R programming experience, command line tools are also available in TSCAN R package. Please check the manual package included in the package for details.

### Install TSCAN via Bioconductor
TSCAN can also be installed via Bioconductor. Note that the TSCAN package may not be most up-to-dated on Bioconductor. To install TSCAN via Bioconductor, run the following commands in R:
```{r }
source("http://bioconductor.org/biocLite.R")
biocLite("TSCAN")
```

## TSCAN Demonstration Video
For users who are not familiar with TSCAN, here is a demonstration video on Youtube for a quick walk-through: https://www.youtube.com/watch?v=zdcBAVe1GBE

## TSCAN User Manual
The user manual for TSCAN command line tools can be viewed in R and is also available at http://www.bioconductor.org/packages/release/bioc/vignettes/TSCAN/inst/doc/TSCAN.pdf. The documentation of all R functions in TSCAN is available at http://www.bioconductor.org/packages/release/bioc/manuals/TSCAN/man/TSCAN.pdf.

## Contact the Author
Author: Zhicheng Ji, Hongkai Ji

Report bugs and provide suggestions by sending email to:
      
Maintainer: Zhicheng Ji (zji4@jhu.edu)

Or open a new issue on this Github page
