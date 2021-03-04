#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include "LaplacianScorer.h"

using namespace Rcpp;


LaplacianScorer::LaplacianScorer(const List& comb_lablacian,
                                 const List& pts_in_vertex,
                                 const arma::sp_mat& adjacency,
                                 const bool one_forms):

  one_forms(one_forms),//true if 1-dim laplacian is also to be considered
  n_vertices (pts_in_vertex.length()),
  l0(as<arma::sp_mat>(comb_lablacian["l0"])),
  zero_weights (as<arma::vec>(comb_lablacian["zero_weights"])),
  length_const_0 (arma::accu(zero_weights)){
  //building points_in_vertex
  for(const IntegerVector& v: pts_in_vertex){
    points_in_vertex.push_back(as<arma::uvec>(v) - 1);
  }

  if(one_forms){
    //getting 1-dim laplacian data
    l1 = as<arma::mat>(comb_lablacian["l1up"]) + as<arma::mat>(comb_lablacian["l1down"]);
    one_weights = as<arma::vec>(comb_lablacian["one_weights"]);
    length_const_1 = arma::accu(one_weights);
    arma::uvec edge_order = as<arma::uvec>(comb_lablacian["adjacency_ordered"]);
    n_edges = edge_order.n_elem;

    //building points_in_edge
    int i = 0;
    points_in_edge = std::vector<arma::uvec>(n_edges);
    for (int j=0; j<n_vertices; j++) {
      arma::uvec o1 = points_in_vertex[j];
      arma::uvec idxs = arma::find(adjacency.row(j));

      if (idxs.size() == 0) continue;

      for(arma::uvec::iterator it = idxs.begin(); it != idxs.end(); ++it) {
        arma::uvec o2 = points_in_vertex[*it];
        arma::uvec edge_cover = arma::intersect(o1, o2);
        //if they do not intersect the edge represents the union
        if (edge_cover.size() == 0) edge_cover = arma::join_cols(o1, o2);
        points_in_edge[edge_order(i)-1] = edge_cover;
        i++;
      }
    }
  }
}


double LaplacianScorer::score(const arma::rowvec& func, int dim = 0){
  //push func to skeleton of dimension dim then scores
  assert(dim == 0 or dim == 1);
  if(dim == 0){
    //pushing function
    arma::rowvec pushed_func(n_vertices);
    for(arma::uword j = 0; j < n_vertices; ++j){
      pushed_func.at(j) = arma::accu(func(points_in_vertex[j])) / points_in_vertex[j].size();
    }
    //computing component orthogonal to constant (\tilde f)
    pushed_func -= arma::dot(pushed_func,zero_weights)/length_const_0;
    //computing scores
    double score = arma::dot(pushed_func*l0, pushed_func) /
                      arma::dot(arma::square(pushed_func), zero_weights);
    if(std::isnan(score)) return arma::datum::inf;
    return score;
  }
  assert(one_forms);
  //case dim = 1
  //pushing function
  arma::rowvec pushed_func(n_edges);
  for(arma::uword j = 0; j < n_edges; ++j){
    pushed_func.at(j) = arma::accu(func(points_in_edge[j])) / points_in_edge[j].size();
  }
  //computing component orthogonal to constant (\tilde f)
  pushed_func -= arma::dot(pushed_func,one_weights)/length_const_1;
  //computing scores
  double score = arma::dot(pushed_func*l1, pushed_func) /
    arma::dot(arma::square(pushed_func), one_weights);
  if(std::isnan(score)) return arma::datum::inf;
  return score;
}

arma::vec  LaplacianScorer::shuffle_score(const arma::rowvec& func, int n_perm, int dim = 0){
  //shuffle the function func and pushes to skeleton of dimension dim n_perm times
  //then scores pushed functions

  assert(dim == 0 or dim == 1);

  arma::rowvec f = arma::shuffle(func);

  if(dim == 0){

    arma::mat pushed_func(n_perm, n_vertices);
    //pushing functions
    for(arma::uword i = 0; i < n_perm; ++i){
      for(arma::uword j = 0; j < n_vertices; ++j){
        pushed_func.at(i,j) = arma::accu(f(points_in_vertex[j])) / points_in_vertex[j].size();
      }
      f = arma::shuffle(f);
    }
    assert(one_forms);
    //computing component orthogonal to constant (\tilde f)
    pushed_func.each_col() -= pushed_func*zero_weights/length_const_0;
    //computing scores
    arma::vec score = arma::diagvec(pushed_func*l0*pushed_func.t()) /
                          (arma::square(pushed_func)*zero_weights);

    score.replace(arma::datum::nan, arma::datum::inf);

    return score;
  }
  assert(one_forms);
  //case dim = 1
  //pushing functions
  arma::mat pushed_func(n_perm, n_edges);
  for(arma::uword i = 0; i < n_perm; ++i){
    for(arma::uword j = 0; j < n_edges; ++j){
      pushed_func.at(i,j) = arma::accu(f(points_in_edge[j])) / points_in_edge[j].size();
    }
    f = arma::shuffle(f);
  }
  //computing component orthogonal to constant (\tilde f)
  pushed_func.each_col() -= pushed_func*one_weights/length_const_1;
  //computing scores
  arma::vec score = arma::diagvec(pushed_func*l1*pushed_func.t()) /
                        (arma::square(pushed_func)*one_weights);
  score.replace(arma::datum::nan, arma::datum::inf);
  return score;
}
