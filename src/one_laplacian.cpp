// [[Rcpp::depends(RcppArmadillo)]]

#include <map>
#include <vector>
#include <math.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

typedef std::map<std::vector<int>, int> BoundaryMap;

BoundaryMap boundary(arma::rowvec simplex)
{
  int n = simplex.size();
  BoundaryMap bound;

  for(int i = 0; i < n; ++i)
  {
    arma::rowvec simplex_copy(simplex);
    simplex_copy.shed_col(i);
    std::vector<int> key(simplex_copy.begin(), simplex_copy.end());

    bound[key] = pow(-1, i);
  }

  return bound;
}

// [[Rcpp::export]]
arma::mat l1down(arma::mat one_simplices, arma::vec zero_weights, arma::vec one_weights)
{
  int n = one_simplices.n_rows;

  arma::mat l1_down(n, n, arma::fill::zeros);

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      if (i > j)
      {
        continue;
      }
      else if (i == j)
      {
        arma::rowvec simplex = one_simplices.row(i);
        l1_down(i, j) += arma::accu(one_weights(i)/zero_weights.elem(arma::conv_to<arma::uvec>::from(simplex)-1));
      }
      else
      {
        arma::rowvec isimplex = one_simplices.row(i);
        arma::rowvec jsimplex = one_simplices.row(j);
        arma::rowvec zero_simplex = intersect(isimplex, jsimplex);

        if (zero_simplex.size() == 0)
        {
          continue;
        }

        // if the i'th and j'th 1-simplices share a 0-simplex (zero_simplex) compute the boundary of the two
        // 1-simplices
        BoundaryMap ff = boundary(isimplex);
        BoundaryMap ffp = boundary(jsimplex);

        // get the sign of the zero_simplex in the two boundaries,
        // l1_down = sgn(E, d(F))*sgn(E, d(F'))
        std::vector<int> key(zero_simplex.begin(), zero_simplex.end());

        int sgn_e_f = ff[key];
        int sgn_e_fp = ffp[key];

        l1_down(i, j) = sgn_e_f*sgn_e_fp*one_weights(j)/zero_weights(zero_simplex(0)-1);
        l1_down(j, i) = sgn_e_f*sgn_e_fp*one_weights(i)/zero_weights(zero_simplex(0)-1);
      }
    }
  }

  return l1_down;
}
