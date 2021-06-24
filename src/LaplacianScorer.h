#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

class LaplacianScorer{

public:

  LaplacianScorer(const List&, const List&, const arma::sp_mat&, const bool);

  //score functions
  arma::vec score(const arma::mat&, int);

  //sample scores by shuffling
  arma::mat sample_scores(const arma::mat&, int, int, int);

  //Sample scores together with covariates
  List sample_with_covariate(const arma::mat&, const arma::mat&, int, int, int);

private:
  //0-d laplacian
  arma::sp_mat l0;

  //1-d laplacian (if present)
  arma::mat l1;

  //weights [0-dim, 1-dim]
  std::vector<arma::vec> weights;

  //number of simplices
  std::vector<int> n_simplex;

  //0-based index of points represented by each simplex
  std::vector<std::vector<arma::uvec>> points_in_simplex;

  //length of the constant function (= 1 in all points) in the inner product defined by weights
  std::vector<double> length_const;
  };
