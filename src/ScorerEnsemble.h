#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "LaplacianScorer.h"
using namespace Rcpp;

class ScorerEnsemble{

public:

  ScorerEnsemble(const List&, const List&, const List&, const bool);

  //score function on all complexes
  arma::mat score(const arma::mat&, int);

  //sample scores of all complexes in tandem
  arma::cube sample_scores(const arma::mat&, arma::uword, int, int);

private:
  //scorers
  std::vector<LaplacianScorer> scorers;

  //number of elements
  int n_scorers;
};
