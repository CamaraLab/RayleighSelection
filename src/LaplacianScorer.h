/* Class to obtain simplex laplacian scores of functions defined on points
 * forming a nerve complex. The class can compute laplacian scores in dimensions
 * 0 and 1. The functions is first pushed to the the 0 (or 1) skeleton by
 * averaging its value over the points associated to each vertex (resp. edge).
 */


#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

class LaplacianScorer{

public:

  LaplacianScorer(const List&, const List&, const arma::sp_mat&, const bool);

  double score(const arma::rowvec&, int);//score a function
  arma::vec shuffle_score(const arma::rowvec&, int, int);//shuffle than score function several times

private:

  arma::sp_mat l0;
  arma::mat l1;
  arma::vec zero_weights;
  arma::vec one_weights;

  unsigned n_vertices, n_edges;

  bool one_forms;

  std::vector<arma::uvec> points_in_vertex;
  std::vector<arma::uvec> points_in_edge;

  double length_const_0, length_const_1;
  };
