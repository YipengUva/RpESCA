## Change to read only
A new package "RpESCA2" will be released soon to replace this version. And this version will be read-only. In RpESCA2, all the core algorithms will be written in C++ and Rcpp.

## package RpESCA
This R package implements the pESCA (penalized exponential familay simultaneous component analysis) model with various concave penalties. The pESCA model with the group concave (concave L2 norm) penalty is applied in my previous paper "Separating common (global and local) and distinct variation in multiple mixed types data sets", https://arxiv.org/abs/1902.06241.

## Installation

You can install the released version of RpESCA:

``` r
devtools::install_gitlab("uvabda/RpESCA")
```

## to do list
1. write the core algorithm in C++ using Rcpp to further increase the speed

2. add more figures to the vignettes

3. put the package on CRAN

