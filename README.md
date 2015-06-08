TSCAN: Tools For Single-Cell ANalysis
====

## Overview
Given single-cell RNA-seq data and true experiment time of cells or pseudo-time cell ordering, TSCAN provides convenient functions for users to assign genes into different gene expression patterns such as constant, monotone increasing and increasing then decreasing. TSCAN then performs GO enrichment analysis to analysis the functional roles of genes with same or similar patterns.

## TSCAN Online User Interface
TSCAN user interface can be directly launched online without installing any software package: https://zhiji.shinyapps.io/TSCAN. PLEASE NOTE: Currently the online version only allows one concurrent user. If the online user interface shows "please wait" for a long time, probably another user is using the online interface and please come back at another time. Users are recommended to install TSCAN on their own computers with following procedures.

## TSCAN Installation

TSCAN software can be installed via Github.
Users should have R installed on their computer before installing GSCA. R can be downloaded here: http://www.r-project.org/.
To install the latest version of GSCA package via Github, run following commands in R:
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

## Contact the Author
Author: Zhicheng Ji, Hongkai Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zji4@jhu.edu)

Or open a new issue on this Github page
