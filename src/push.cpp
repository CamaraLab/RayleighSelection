// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

/* Given an open cover x and a function v, it computes the pushforward of v to
 * the open sets of x given by the average. It also performs perm random permutations
 * of the values of v and computes the corresponding pushforwards.
 */

// [[Rcpp::export]]
List pushCpp(arma::vec v, List x, SEXP perm, arma::sp_mat adjacency) {
  List xlist(x);
  arma::vec xv(v);
  arma::vec co = xv;
  int n = xlist.size();
  int perms = as<int >(perm);
  arma::mat res_vertices = arma::zeros(perms+1, n);
  arma::mat res_edges = arma::zeros(perms+1, n);

  for (int i=0; i<perms+1; i++) {
    for (int j=0; j<n; j++) {

      arma::ivec o1 = xlist[j];
      int u = o1.size();

      for (int k=0; k<u; k++) {
        res_vertices(i,j) += co(o1[k]-1);
      }

      res_vertices(i,j) /= u;

      arma::mat row(adjacency.row(j));
      arma::uvec idxs = arma::find(row != 0);

      if (idxs.size() == 0)
      {
        continue;
      }

      for(arma::uvec::iterator it = idxs.begin(); it != idxs.end(); ++it)
      {
        arma::ivec o2 = xlist[*it];
        arma::ivec edge_cover = arma::intersect(o1, o2);
        int v = edge_cover.size();

        for(int k = 0; k < v; k++)
        {
          res_edges(i, j) += co(edge_cover[k]-1);
        }

        res_edges(i, j) /= v;
      }
    }

    co = arma::shuffle(co);
  };

  return List::create(Named("vertices") = res_vertices,
                      Named("edges") = res_edges);
}
