// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "LaplacianScorer.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
List rayleight_selectionCpp(const List& comb_lablacian,
                            const List& points_in_vertex,
                            const arma::sp_mat& adjacency,
                            const arma::mat& func, int n_perm,
                            const int num_cores,
                            const bool one_forms = true){
  /*Computes laplacian scores and p-values.
   *
   * Arguments:
   *
   * comb_laplacian: List output by combinatorial_laplacian
   * points_in_vertex: List of points in each vertex
   * adjacency: 1-simplex adjacency matrix (as a sparce matrix)
   * func: Matrix with functions (defined on points). Each row defines a different function
   * n_perm: Number of permutations to be performed when computing p-values
   * num_cores: number of cores to be used. OpenMP is used for multithreading.
   * one_forms: True if, and only if, scores and p-values relative to 1-laplacian
   *            should also be computed.
   *
   *Return: List with scores (R0 and R1) and p-values (p0 and p1)
   *
   */

  LaplacianScorer Scorer(comb_lablacian, points_in_vertex, adjacency, one_forms);

  arma::uword n_funcs = func.n_rows;

  std::vector<double> scores_0(n_funcs);
  std::vector<double> p_values_0(n_funcs);

  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(num_cores)
  #endif
  for(int i = 0; i < n_funcs; i++){
    arma::rowvec f = func.row(i);
    scores_0[i] = Scorer.score(f, 0);
    arma::vec sampled_scores = Scorer.shuffle_score(f, n_perm, 0);
    p_values_0[i] = (double)(arma::accu(sampled_scores <= scores_0[i]))/ n_perm;
  }

  if(not one_forms) return List::create(Named("R0") = scores_0, Named("p0") = p_values_0);

  std::vector<double> scores_1(n_funcs);
  std::vector<double> p_values_1(n_funcs);

  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(num_cores)
  #endif
  for(int i = 0; i < n_funcs; i++){
    arma::rowvec f = func.row(i);
    scores_1[i] = Scorer.score(f, 1);
    arma::vec sampled_scores = Scorer.shuffle_score(f, n_perm, 1);
    p_values_1[i] = (double)(arma::accu(sampled_scores <= scores_1[i]))/ n_perm;
  }

  return List::create(Named("R0") = wrap(scores_0), Named("p0") = wrap(p_values_0),
                      Named("R1") = wrap(scores_1), Named("p1") = wrap(p_values_1));
}
