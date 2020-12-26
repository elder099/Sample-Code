// [[Rcpp::depends(RcppEigen, bigmemory, BH)]]
#include <RcppEigen.h>
#include <bigmemory/MatrixAccessor.hpp>



// [[Rcpp::export]]
Rcpp::NumericVector prod4(Rcpp::XPtr<BigMatrix> bMPtr, const Rcpp::NumericVector& x) {

  MatrixAccessor<double> macc(*bMPtr);

  int n = bMPtr->nrow();
  int m = bMPtr->ncol();

  Rcpp::NumericVector res(n);
  int i, j;

  for (j = 0; j <= m - 4; j += 4) {
    for (i = 0; i < n; i++) { // unrolling optimization
      res[i] += (x[j] * macc[j][i] + x[j+1] * macc[j+1][i]) +
        (x[j+2] * macc[j+2][i] + x[j+3] * macc[j+3][i]);
    } // The parentheses are somehow important. Try without.
  }
  for (; j < m; j++) {
    for (i = 0; i < n; i++) {
      res[i] += x[j] * macc[j][i];
    }
  }

  return Rcpp::wrap(res);
}
