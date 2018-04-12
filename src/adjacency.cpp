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
  double mean;

  for(NumericVector::iterator it = open_set.begin(); it != open_set.end(); ++it)
  {
    IntegerVector::iterator idx_it = std::find(sample_order.begin(), sample_order.end(), *it);
    int idx = std::distance(sample_order.begin(), idx_it);
    mean += (idx+1);
  }

  return mean/((double) open_set.size());
}

double find_index_mean(NumericVector& open_set,
                       int vertex,
                       IntegerVector& sample_order,
                       std::map<int, double>& sample_index_mean)
{
  double mean;
  std::map<int, double>::iterator it = sample_index_mean.find(vertex);

  if (it == sample_index_mean.end())
  {
    mean = compute_vertex_index_mean(open_set, sample_order);
    sample_index_mean.insert(std::pair<int, double>(vertex, mean));
  }
  else
  {
    mean = it->second;
  }

  return mean;
}

/* Given an open cover x, it computes the adjacency matrix given by pairwise
 * intersections of the open sets.
 */
// [[Rcpp::export]]
List adjacencyCpp(List x, DataFrame& features, bool weight) {
  List xlist(x);
  int n = xlist.size();
  arma::sp_mat one_simplices(n,n);
  arma::sp_mat upper_tri_one_simplices(n, n);
  List two_simplices(n);
  two_simplices.fill(arma::sp_mat(n, n));

  IntegerVector sample_order = features["sample"];
  std::map<int, double> sample_index_mean;

  for(int i=0; i<n-1; i++) {

    NumericVector loc1 = xlist[i];

    double loc1_size = (double)loc1.size();

    for(int j=i+1; j<n; j++) {

      NumericVector loc2 = xlist[j];
      double loc2_size = (double)loc2.size();

      double intersection = (double)(intersect(loc1,loc2).size());

      if (intersection > 0) {

        double cover1 = find_index_mean(loc1, i, sample_order, sample_index_mean);
        double cover2 = find_index_mean(loc2, j, sample_order, sample_index_mean);

        if (weight) {
          double jidx = (intersection/(loc1_size + loc2_size - intersection));

          if(cover1 < cover2)
          {
            one_simplices(j, i) = 1.0 - jidx;
          }
          else
          {
            one_simplices(i, j) = 1.0 - jidx;
          }

        } else {

          if(cover1 < cover2)
          {
            one_simplices(j, i) = 1.0;
          }
          else
          {
            one_simplices(i, j) = 1.0;
          }
        }

        upper_tri_one_simplices(i, j) = 1.0;
      };
    };
  };

  for(int i = 0; i < n-1; i++)
  {
    NumericVector cover0 = xlist[i];

    arma::mat row(upper_tri_one_simplices.row(i));

    arma::uvec idxs = arma::find(row != 0);
    arma::uvec subset = arma::find(idxs > i);
    idxs = idxs(subset);

    if(idxs.size() < 2)
    {
      continue;
    }

    std::vector<std::vector<int>> edge_idxs = combo(idxs.size(), 2);

    for(auto&& k : edge_idxs)
    {
      if(one_simplices(idxs[k[0]], idxs[k[1]]) == 0 && one_simplices(idxs[k[1]], idxs[k[0]]) == 0)
      {
        continue;
      }

      NumericVector cover1 = xlist[idxs[k[0]]];
      NumericVector cover2 = xlist[idxs[k[1]]];

      double intersection = (double)(intersect(intersect(cover0, cover1), cover2).size());

      if (intersection == 0)
      {
        continue;
      }

      double cover0_mean = find_index_mean(cover0, i, sample_order, sample_index_mean);
      double cover1_mean = find_index_mean(cover1, idxs[k[0]], sample_order, sample_index_mean);
      double cover2_mean = find_index_mean(cover2, idxs[k[1]], sample_order, sample_index_mean);

      std::vector<std::pair<int, double>> vertices;
      vertices.push_back(std::make_pair(i, cover0_mean));
      vertices.push_back(std::make_pair(idxs[k[0]], cover1_mean));
      vertices.push_back(std::make_pair(idxs[k[1]], cover2_mean));

      sort(vertices.begin(), vertices.end(), sortbysec<int, double>);

      Rcout << "Inserting two simplex: "
            << vertices[0].first+1 << " "
            << vertices[1].first+1 << " "
            << vertices[2].first+1 << std::endl;

      arma::sp_mat adjacency = two_simplices(vertices[0].first);
      adjacency(vertices[1].first, vertices[2].first) = 1;
      two_simplices(vertices[0].first) = adjacency;
    }
  }

  return List::create(Named("one_simplices") = one_simplices,
                      Named("two_simplices") = two_simplices);
}
