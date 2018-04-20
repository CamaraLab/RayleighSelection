// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;


// taken from: https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
std::vector<std::vector<int>> combo(int N, int K)
{
  std::string bitmask(K, 1); // K leading 1's
  bitmask.resize(N, 0); // N-K trailing 0's
  std::vector<std::vector<int>> combos;

  // print integers and permute bitmask
  do {
    std::vector<int> combo;
    for (int i = 0; i < N; ++i) // [0..N-1] integers
    {
      if (bitmask[i]) combo.push_back(i);
    }
    combos.push_back(combo);
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return combos;
}


template<typename T, typename TT>
bool sortbysec(const std::pair<T, TT> &a, const std::pair<T, TT> &b)
{
  return (a.second < b.second);
}

double compute_vertex_index_mean(NumericVector& open_set, IntegerVector& sample_order)
{
  double mean = 0.0;

  for(NumericVector::iterator it = open_set.begin(); it != open_set.end(); ++it)
  {
    IntegerVector::iterator idx_it = std::find(sample_order.begin(), sample_order.end(), *it);
    int idx = std::distance(sample_order.begin(), idx_it);
    mean += (idx+1);
  }

  return mean/((double) open_set.size());
}

int find_index(int vertex,
               std::vector<int>& rank)
{

  std::vector<int>::iterator it = std::find(rank.begin(), rank.end(), vertex+1);
  return std::distance(rank.begin(), it);
}

/* Given an open cover x, it computes the adjacency matrix given by pairwise
 * intersections of the open sets.
 */
// [[Rcpp::export]]
List adjacencyCpp(List x, DataFrame& features, bool weight) {
  List xlist(x);
  int n = xlist.size();
  arma::sp_mat one_simplices(n,n);

  List two_simplices(n);
  two_simplices.fill(arma::sp_mat(n, n));

  IntegerVector sample_order = features["sample"];
  std::vector<std::pair<int, double>> sample_index_mean;

  for(int i = 0; i < n; i++)
  {
    NumericVector cover = xlist[i];
    double mean = compute_vertex_index_mean(cover, sample_order);
    sample_index_mean.push_back(std::pair<int, double>(i+1, mean));
  }

  sort(sample_index_mean.begin(), sample_index_mean.end(), sortbysec<int, double>);
  std::vector<int> rank;
  std::transform(sample_index_mean.begin(), sample_index_mean.end(), back_inserter(rank),
                 (const int& (*)(const std::pair<int, double>&))std::get<0>);

  for(int i=0; i<n-1; i++) {

    NumericVector loc1 = xlist[i];

    double loc1_size = (double)loc1.size();

    for(int j=i+1; j<n; j++) {

      NumericVector loc2 = xlist[j];
      double loc2_size = (double)loc2.size();

      double intersection = (double)(intersect(loc1,loc2).size());

      if (intersection > 0) {
        int idxi = find_index(i, rank);
        int idxj = find_index(j, rank);

        if (idxi > idxj)
        {
          std::swap(idxi, idxj);
        }

        if (weight) {
          double jidx = (intersection/(loc1_size + loc2_size - intersection));
          one_simplices(idxi, idxj) = 1.0 - jidx;
        }
        else {
          one_simplices(i, j) = 1.0;
        }
      };
    };
  };

  for(int i = 0; i < n-1; i++)
  {
    NumericVector cover0 = xlist[rank[i]-1];

    arma::mat row(one_simplices.row(i));
    arma::uvec idxs = arma::find(row != 0);

    if(idxs.size() < 2)
    {
      continue;
    }

    std::vector<std::vector<int>> edge_idxs = combo(idxs.size(), 2);

    for(auto&& e : edge_idxs)
    {
      int j = idxs[e[0]];
      int k = idxs[e[1]];

      if(one_simplices(j, k) == 0 && one_simplices(k, j) == 0)
      {
        continue;
      }

      NumericVector cover1 = xlist[rank[j]-1];
      NumericVector cover2 = xlist[rank[k]-1];

      double intersection = (double)(intersect(intersect(cover0, cover1), cover2).size());

      if (intersection == 0)
      {
        continue;
      }

      Rcout << "Inserting two simplex: "
            << rank[i] << " "
            << rank[j] << " "
            << rank[k] << std::endl;

      arma::sp_mat adjacency = two_simplices(rank[i]-1);
      adjacency(rank[j]-1, rank[k]-1) = 1;
      two_simplices(rank[i]-1) = adjacency;
    }
  }

  return List::create(Named("one_simplices") = one_simplices,
                      Named("two_simplices") = two_simplices,
                      Named("order") = rank);
}
