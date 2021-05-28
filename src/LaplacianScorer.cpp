//' Computes and samples 0 and 1-dimensional laplacian scores.
//'
//' Builds an object to compute and sample values for the combinatorial laplacian score of functions
//' defined on the set of points underlying a given simplex. Takes in data about the simplex and
//' its 0- (and possibly 1-) dimensional laplacian to build a scorer. The scorer pushes functions defined on
//' the points underlying the simplex to the 0- (or 1-) dimensional skeleton by averaging over the points associated
//' to each simplex and evaluates the Rayleigh quotient of the normalized graph laplacean on the pushed function.
//' The scorer can also sample values from the null distribution obtained by shuffling the labels of the underlying points.
//'
//' @name LaplacianScorer
//' @field new Constructs a scorer using  data on the simplex and its associated laplacian.\cr\cr
//' \strong{Use} \code{scorer <- new(LaplacianScorer, comb_laplacian, pts_in_vertex, adjacency, one_forms)}\cr\cr
//' \strong{Parameters}\itemize{
//'  \item \code{comb_laplacian}: output of \code{\link{combinatorial_laplacian}}
//'  \item \code{pts_in_vertex}: list of vectors where the i-th vector contains the index of the points associated to the i-th vertex
//'  \item \code{adjacency}: adjacency matrix for 1-skeleton as a sparse matrix
//'  \item \code{one_forms}: boolean indicating if laplacian of one-forms will be computed}
//'
//' @field score Pushes functions defined by rows of funcs to the dim-skeleton by
//' averaging and computes its laplacian score.\cr\cr
//' \strong{Use} \code{scorer$(funcs, dim)}\cr\cr
//' \strong{Parameters}\itemize{
//'  \item \code{funcs}: functions to be scored as a rows of a dense matrix
//'  \item \code{dim}: dimension of the laplacian (0 or 1)}
//' \strong{Value} Scores of functions as a vector
//'
//' @field sample_scores Takes samples of scores by permuting point labels\cr\cr
//' \strong{Use} \code{scorer$sample_scores(funcs, n_perm, dim, n_cores)}\cr\cr
//' \strong{Parameters}\itemize{
//'  \item \code{funcs}: base functions as rows of a dense matrix
//'  \item \code{n_perm}: number of permutations
//'  \item \code{dim}: dimension of the laplacian
//'  \item \code{n_cores}: number of cores to be used, parallelization requires code to be compiled with \code{openmp} }
//' \strong{Value} Dense matrix with sampled scores where the i-th row has samples for the i-th function
//' @examples
//' library(RayleighSelection)
//'
//' # Create a simplicial complex and compute the associated 0- and 1-laplacians
//' gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
//' lout <- combinatorial_laplacian(gy, one_forms = TRUE)
//'
//' # Create an associated instance of LaplacianScorer
//' scorer <- new(LaplacianScorer, lout, gy$points_in_vertex, gy$adjacency, TRUE)
//'
//' # Compute the 0-laplacian score of a a function
//' scorer$score(t(as.matrix(c(0,1,1,0,0,0,0,0,0,1))), 0)
//'
//' # Sample scores by suffling the function
//' scorer$sample_scores(t(as.matrix(c(0,1,1,0,0,0,0,0,0,1))), 10, 0, 1)
//'

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include "LaplacianScorer.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


LaplacianScorer::LaplacianScorer(const List& comb_lablacian,
                                 const List& pts_in_vertex,
                                 const arma::sp_mat& adjacency,
                                 const bool one_forms){
  //building points in vertex
  points_in_simplex.push_back(std::vector<arma::uvec> ());
  for(const IntegerVector& v: pts_in_vertex){
    points_in_simplex[0].push_back(as<arma::uvec>(v) - 1);
  }
  //building other 0-d attributes
  weights.push_back(as<arma::vec>(comb_lablacian["zero_weights"]));
  length_const.push_back(arma::accu(weights[0]));
  l0 = as<arma::sp_mat>(comb_lablacian["l0"]);
  n_simplex.push_back(pts_in_vertex.length());

  if(one_forms){
    //getting 1-dim laplacian data
    l1 = as<arma::mat>(comb_lablacian["l1up"]) + as<arma::mat>(comb_lablacian["l1down"]);
    weights.push_back(as<arma::vec>(comb_lablacian["one_weights"]));
    length_const.push_back(arma::accu(weights[1]));
    arma::uvec edge_order = as<arma::uvec>(comb_lablacian["adjacency_ordered"]);
    n_simplex.push_back(edge_order.n_elem);

    //building points in each edge
    int i = 0;
    points_in_simplex.push_back( std::vector<arma::uvec>(n_simplex[1]) );
    for (int j=0; j < n_simplex[0]; j++) {
      arma::uvec o1 = points_in_simplex[0][j];
      arma::uvec idxs = arma::find(adjacency.row(j));

      if (idxs.size() == 0) continue;

      for(arma::uvec::iterator it = idxs.begin(); it != idxs.end(); ++it) {
        arma::uvec o2 = points_in_simplex[0][*it];
        arma::uvec edge_cover = arma::intersect(o1, o2);
        //if they do not intersect the edge represents the union
        if (edge_cover.size() == 0) edge_cover = arma::join_cols(o1, o2);
        points_in_simplex[1][edge_order(i)-1] = edge_cover;
        i++;
      }
    }
  }
}


arma::vec LaplacianScorer::score(const arma::mat& funcs, int dim){

  //pushing function
  arma::uword n_funcs = funcs.n_rows;
  arma::mat pushed_func(n_funcs, n_simplex[dim], arma::fill::zeros);
  for(arma::uword i = 0; i < n_funcs; ++i){
    for(arma::uword j = 0; j < n_simplex[dim]; ++j){
      for(const auto p: points_in_simplex[dim][j]) pushed_func(i,j) += funcs(i, p);
      pushed_func(i,j) /= points_in_simplex[dim][j].size();
    }
  }
  //computing component orthogonal to constant (\tilde f)
  pushed_func.each_col() -= pushed_func*weights[dim]/length_const[dim];
  //computing scores
  arma::vec score;
  if(dim == 0){
    score = arma::diagvec(pushed_func*l0*pushed_func.t()) /
      (arma::square(pushed_func)*weights[dim]);
  }else{
    score = arma::diagvec(pushed_func*l1*pushed_func.t()) /
      (arma::square(pushed_func)*weights[dim]);
  }
  score.replace(arma::datum::nan, arma::datum::inf);

  return score;
}

arma::mat LaplacianScorer::sample_scores(const arma::mat& funcs, int n_perm,
                                          int dim, int n_cores){
  arma::uword n_funcs = funcs.n_rows;
  arma::mat scores(n_funcs, n_perm);

  //setting a maximum number of samples to be taken at once
  int max_samp = std::pow(10, 9)/n_simplex[dim];

  #if defined(_OPENMP)
  #pragma omp parallel for num_threads(n_cores)
  #endif
  for(arma::uword i = 0; i < n_funcs; i++){
    arma::rowvec f = funcs.row(i);
    int sampled = 0;

    while(sampled < n_perm){
      int n_samps = std::min(max_samp, n_perm - sampled);

      arma::mat shuffled (n_samps, funcs.n_cols);
      for(int k = 0; k < n_samps; k++) shuffled.row(k) = arma::shuffle(f);

      scores(i, arma::span(sampled, sampled + n_samps - 1)) = score(shuffled, dim).t();
      sampled += n_samps;
      }
  }
  return scores;
}

RCPP_MODULE(mod_laplacian){
  class_<LaplacianScorer>("LaplacianScorer")

  .constructor<List, List, arma::sp_mat, bool>()

  .method("score", &LaplacianScorer::score)
  .method("sample_scores", &LaplacianScorer::sample_scores)
  ;
}
