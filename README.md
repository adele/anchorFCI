## AnchorFCI: an extension of the FCI algorithm designed to improve robustness and discovery power by strategically selecting and integrating reliable anchor variables into the graph, while leveraging known non-ancestral relationships.
  
### Overview

This package provides an implementation of the anchorFCI by Ribeiro et al. (Forthcoming).

AnchorFCI operates on two sets of variables: the first set contains the variables of interest, while the second comprises variables that are not caused by any from the first. 

While this structure is beneficial for various applications, it is especially well-suited for datasets involving phenotypic, clinical, and sociodemographic variables (the first set), alongside genetic variables, such as SNPs (the second set), which are recognized as not being caused by the first set.


### Installation

First, install R (>= 3.5.0) and the following packages:
```r
install.packages(c("FCI.Utils", "pcalg", "igraph", "RBGL", "graph", "doFuture", "gtools", "MXM", "pscl", "DOT", "rsvg"), dependencies=TRUE)
```
You can download the latest tar.gz file with the source code of the IOD R package, available at <https://github.com/adele/anchorFCI/releases/latest>, and install it with the following command, where `path_to_file` represents the full path and file name of the tar.gz file:

``` r
install.packages(path_to_file, repos=NULL, type="source", dependencies=TRUE)
```


