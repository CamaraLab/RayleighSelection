// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

/* Given an open cover x, it computes the adjacency matrix given by pairwise
 * intersections of the open sets.
 */

// [[Rcpp::export]]
arma::sp_mat adjacencyCpp(List x) {
  List xlist(x);
  int n = xlist.size();
  arma::sp_mat res(n,n);

  for(int i=0; i<n-1; i++) {
    NumericVector loc1 = xlist[i];
    for(int j=i+1; j<n; j++) {
      NumericVector loc2 = xlist[j];
      if (intersect(loc1,loc2).length() > 0) {
        res(i,j) = 1.0;
        res(j,i) = 1.0;
      };
    };
  };
  return res;
}
