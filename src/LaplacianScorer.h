#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

class LaplacianScorer{

public:

  LaplacianScorer(const List&, const List&, const arma::sp_mat&, const bool);

  //score functions
  arma::vec score(const arma::mat&, int);

  //sample scores by shuffling
  arma::mat sample_scores(const arma::mat&, arma::uword, int, int);

  //sample scores together with covariates
  List sample_with_covariate(const arma::mat&, const arma::mat&, arma::uword, int, int);

private:
  // All attributes are vectors with 1 or 2 elements, depending on whether one_forms is true.
  // The first entry refers to 0-d data end the second (if present) refers to 1-d data.

  //laplacians
  std::vector<arma::sp_mat> l;

  //weights
  std::vector<arma::vec> weights;

  //number of simplices
  std::vector<arma::uword> n_simplex;

  //0-based index of points represented by each simplex
  std::vector<std::vector<arma::uvec>> points_in_simplex;

  //length of the constant function (= 1 in all points) in the inner product defined by weights
  std::vector<double> length_const;
  };
