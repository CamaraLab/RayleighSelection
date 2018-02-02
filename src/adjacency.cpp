// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

/* Given an open cover x, it computes the adjacency matrix given by pairwise
 * intersections of the open sets.
 */

// [[Rcpp::export]]
arma::sp_mat adjacencyCpp(List x, bool weight) {
  List xlist(x);
  int n = xlist.size();
  arma::sp_mat res(n,n);

  for(int i=0; i<n-1; i++) {

    NumericVector loc1 = xlist[i];

    double loc1_size = (double)loc1.size();

    for(int j=i+1; j<n; j++) {

      NumericVector loc2 = xlist[j];

      double intersection = (double)(intersect(loc1,loc2).size());

      if (intersection > 0) {

        if (weight) {
          double loc2_size = (double)loc2.size();

          double jidx = (intersection/(loc1_size + loc2_size - intersection));

          res(i, j) = 1.0 - jidx;
          res(j, i) = 1.0 - jidx;

        } else {
          res(i,j) = 1.0;
          res(j,i) = 1.0;
        }
      };
    };
  };
  return res;
}
