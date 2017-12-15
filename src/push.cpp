// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

/* Given an open cover x and a function v, it computes the pushforward of v to
 * the open sets of x given by the average. It also performs perm random permutations
 * of the values of v and computes the corresponding pushforwards.
 */

// [[Rcpp::export]]
arma::mat pushCpp(arma::vec v, List x, SEXP perm) {
  List xlist(x);
  arma::vec xv(v);
  arma::vec co = xv;
  int n = xlist.size();
  int perms = as<int >(perm);
  arma::mat res = arma::zeros(perms+1, n);

  for (int i=0; i<perms+1; i++) {
    for (int j=0; j<n; j++) {
      arma::ivec o = xlist[j];
      int u = o.size();
      for (int k=0; k<u; k++) {
        res(i,j) += co(o[k]-1);
      };
      res(i,j) /= u;
    };
  co = arma::shuffle(co);
  };
  return res;
}
