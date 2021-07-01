//' Computes and samples 0 and 1-dimensional laplacian scorers in ensembles of simplicial complexes.
//'
//' Builds an ensemble of \link{LaplacianScorer} objects that can be used to score functions
//' and to sample scores by shuffling point labels in tandem.
//'
//' @name ScorerEnsemble
//' @field new Constructs a scorer ensemble using lists containing data on the complex
//' and their laplacians.\cr\cr
//' \strong{Use} \code{scorer.ensemple <- new(ScorerEnsemble, comb_laplacian, pts_in_vertex, adjacency, one_forms)}\cr\cr
//' \strong{Parameters}\itemize{
//'  \item \code{comb_laplacian}: List containing outputs of \code{\link{combinatorial_laplacian}}
//'  \item \code{pts_in_vertex}: List of lists of vectors where the i-th vector of the j-th list
//'         contains the index of the points associated to the i-th vertex of the j-th complex
//'  \item \code{adjacency}: List of adjacency matrices for 1-skeleton as a sparse matrix
//'  \item \code{one_forms}: boolean indicating if laplacian of one-forms will be computed}
//'
//' @field score Pushes functions defined by rows of funcs to the dim-skeleton of all complexes by
//' averaging and computes its laplacian score.\cr\cr
//' \strong{Use} \code{scorer$score(funcs, dim)}\cr\cr
//' \strong{Parameters}\itemize{
//'  \item \code{funcs}: functions to be scored as a rows of a dense matrix
//'  \item \code{dim}: dimension of the laplacian (0 or 1)}
//' \strong{Value} Scores of functions as a metrix. The value in position (i,j) is the score of
//'                the i-th function on the j-th complex.
//'
//' @field sample_scores Takes samples of scores by permuting point labels\cr\cr
//' \strong{Use} \code{scorer$sample_scores(funcs, n_perm, dim, n_cores)}\cr\cr
//' \strong{Parameters}\itemize{
//'  \item \code{funcs}: base functions as rows of a dense matrix
//'  \item \code{n_perm}: number of permutations
//'  \item \code{dim}: dimension of the laplacian
//'  \item \code{n_cores}: number of cores to be used, parallelization requires code to be compiled with \code{openmp} }
//' \strong{Value} Dense 3-dimensional array with sampled scores where position (i,j,k) has the score of the j-th
//'                permutation of the i-th function in the k-th complex.
//' @examples
//' library(RayleighSelection)
//'
//' # Create a simplicial complex and compute the associated 0- and 1-laplacians
//' gy.list <- list(
//'   nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5))),
//'   nerve_complex(list(c(1,6,10), c(1,4,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
//' )
//' lout.list <- lapply(gy.list, combinatorial_laplacian, one_forms = TRUE)
//'
//' # Create an associated instance of ScorerEnsemble
//' scorer.ensemble <- new(
//'   ScorerEnsemble,
//'   lout.list,
//'   lapply(gy.list, function(gy) gy$points_in_vertex),
//'   lapply(gy.list, function(gy) gy$adjacency),
//'   TRUE
//' )
//'
//' # Compute the 0-laplacian scores of a a function
//' scorer.ensemble$score(t(as.matrix(c(0,1,1,0,0,0,0,0,0,1))), 0)
//'
//' # Sample scores by shuffling the function
//' scorer$sample_scores(t(as.matrix(c(0,1,1,0,0,0,0,0,0,1))), 10, 0, 1)



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include "ScorerEnsemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif


ScorerEnsemble::ScorerEnsemble(const List& comb_laplacian,
                               const List& pts_in_vertex,
                               const List& adjacency,
                               const bool one_forms){
  n_scorers = comb_laplacian.length();
  for(int i = 0; i < n_scorers; i++){
    scorers.push_back(
      LaplacianScorer(comb_laplacian[i], pts_in_vertex[i], adjacency[i], one_forms)
    );
  }
}

arma::mat ScorerEnsemble::score(const arma::mat& funcs, int dim){
  arma::mat scores(funcs.n_rows, n_scorers);
  for(int i = 0; i < n_scorers; i++) scores.col(i) = scorers[i].score(funcs, dim);
  return scores;
  }

arma::cube ScorerEnsemble::sample_scores(const arma::mat& funcs, arma::uword n_perm,
                                    int dim, int n_cores){
  arma::uword n_funcs = funcs.n_rows;
  arma::cube scores(n_funcs, n_perm, n_scorers);

  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(n_cores)
  #endif
  for(arma::uword i = 0; i < n_funcs; i++){
    arma::rowvec f = funcs.row(i);
    arma::mat shuffled (n_perm, funcs.n_cols);
    for(arma::uword j = 0; j < n_perm; j++) shuffled.row(j) = arma::shuffle(f);
    scores.row(i) = score(shuffled, dim);
  }
  return scores;
}

RCPP_MODULE(mod_scorer_ensemble){
  class_<ScorerEnsemble>("ScorerEnsemble")

  .constructor<List, List, List, bool>()

  .method("score", &ScorerEnsemble::score)
  .method("sample_scores", &ScorerEnsemble::sample_scores)
  ;
}
