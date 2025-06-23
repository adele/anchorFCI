## AnchorFCI: Enhancing Causal Discovery with Reliable Anchors

This package provides an implementation of the anchorFCI by Ribeiro et al. (2024), available at <https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2024.1436947/full>, doi: 10.3389/fgene.2024.1436947

AnchorFCI is an extension of the FCI algorithm designed to improve robustness and discovery power in causal discovery by strategically selecting reliable anchors, while leveraging their known non-ancestral relationships.

It operates on two sets of variables: the first set contains the variables of interest, while the second comprises variables that are not caused by any from the first. 

While this structure is beneficial for various applications, it is especially well-suited for datasets involving phenotypic, clinical, and sociodemographic variables (the first set), alongside genetic variables, such as SNPs (the second set), which are recognized as not being caused by the first set.


### Citation

If you use this work, please cite:

Ribeiro, A. H., Crnkovic, M., Pereira, J. L., Fisberg, R. M., Sarti, 
F. M., Rogero, M. M., Heider, D., & Cerqueira, A. (2024). 
AnchorFCI: Harnessing genetic anchors for enhanced causal discovery of 
cardiometabolic disease pathways. Frontiers in Genetics, 15, 1436947. 
https://doi.org/10.3389/fgene.2024.1436947


```
@article{ribeiro2024anchorfci,
title={AnchorFCI: Harnessing genetic anchors for enhanced causal discovery of cardiometabolic disease pathways}, 
author={Ribeiro, Adèle H. and Crnkovic, Milena and Pereira, Jaqueline L. and Fisberg, 
Regina M. and Sarti, Flavia M. and Rogero, Marcelo M. and Heider, Dominik and Cerqueira, Andressa}, 
journal={Frontiers in Genetics}, 
volume={15}, pages={1436947}, year={2024}, doi={10.3389/fgene.2024.1436947} }
```

### Installation

First, install R (>= 3.5.0) and the following packages:
```r
install.packages(c("FCI.Utils", "pcalg", "igraph", "RBGL", "graph", "doFuture", "gtools", "MXM", "pscl", "DOT", "rsvg"), dependencies=TRUE)
```

Then, install the FCI.Utils R package, available at <https://github.com/adele/FCI.Utils>/

You can install the development version directly from GitHub:

``` r
install.packages("devtools", dependencies=TRUE)
devtools::install_github("adele/FCI.Utils", dependencies=TRUE)
```


Now, you can download the latest tar.gz file with the source code of the anchorFCI R package.

The lastest version is available at <https://github.com/adele/anchorFCI/releases/latest>. 
You can install it with the following command, where `path_to_file` represents the full path and file name of the tar.gz file:

``` r
install.packages(path_to_file, repos=NULL, type="source", dependencies=TRUE)
```


