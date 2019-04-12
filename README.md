# package RpESCA

This package is the R implementation of pESCA models with various concave penalties. The details of the 
algorithms developed in the package can be referred to \url{https://arxiv.org/abs/1902.06241}.

## Installation

You can install the released version of RpESCA from [github](https://github.com) with:

``` r
# a small part of the algorithms are in C++
# make sure you have installed Rcpp package

devtools::install_github("YipengUva/RpESCA")
```

## to do list
1. write the core algorithm in C++ using Rcpp to further increase the speed

2. put it on CRAN
