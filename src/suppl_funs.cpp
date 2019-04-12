#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// This is an implementaion of fast verion trace function.
// This function will compute the trace of two matrices.
// [[Rcpp::export]]
double fast_traceC(const arma::mat& X, const arma::mat& Y) {
  int n = X.n_rows, p = X.n_cols;
  if(n>p){
    return arma::trace(X.t()*Y);
  }else{
    return arma::trace(Y*X.t());
  }
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Matrix * diag(vec) * matrix
// [[Rcpp::export]]
arma::mat mat_vec_matC(const arma::mat& X, const arma::colvec& d, const arma::mat& Y) {
	return X* arma::diagmat(d) *Y;
}