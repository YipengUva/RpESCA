# package RpESCA

This package is the R implementation of pESCA models with various concave penalties. The details of the 
algorithms developed in the package can be referred to \url{https://arxiv.org/abs/1902.06241}.

## Installation

You can install the released version of RpESCA from [github](https://github.com) with:

``` r
# a small part of the algorithms is in C++
# make sure you have installed Rcpp package

devtools::install_git("https://gitlab.com/uvabda/RpESCA.git")
```

## to do list
1. write the core algorithm in C++ using Rcpp to further increase the speed

2. add more figures to the vignettes

3. put the package on CRAN
