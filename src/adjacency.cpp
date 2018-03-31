// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

template<typename T, typename TT>
bool sortbysec(const std::pair<T, TT> &a, const std::pair<T, TT> &b)
{
  return (a.second < b.second);
}

// taken from: https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
template<typename T>
std::vector<T> combo(const T& c, int k)
{
  std::vector<T> combos;

  int n = c.size();
  int combo = (1 << k) - 1;       // k bit sets
  while (combo < 1<<n) {

    T comb;
    int n = c.size();
    for (int i = 0; i < n; ++i) {
      if ((combo >> i) & 1)
        comb.push_back(c[i]);
    }
    combos.push_back(comb);

    int x = combo & -combo;
    int y = combo + x;
    int z = (combo & ~y);
    combo = z / x;
    combo >>= 1;
    combo |= y;
  }

  return combos;
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
  arma::cube two_simplices(n, n, n, arma::fill::zeros);

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
      };
    };
  };

  // generate all possible two simplices and for each possible two simplex generate a list of its faces
  // check that its faces are contained in the adjacency matrix of one simplices, if so, compute the
  // triple intersection of the covers
  std::vector<int> vertices(n);
  std::iota(vertices.begin(), vertices.end(), 0);
  std::vector<std::vector<int>> two_simplices_candidates = combo(vertices, 3);

  for(auto&& it : two_simplices_candidates)
  {
    std::vector<std::vector<int>> faces = combo(it, 2);
    bool contains_faces = true;

    for(auto&& itt : faces)
    {
      if(one_simplices(itt[0], itt[1]) == 0 && one_simplices(itt[1], itt[0]) == 0)
      {
        contains_faces = false;
        break;
      }
    }

    if (!contains_faces)
    {
      continue;
    }

    int v0 = it[0];
    int v1 = it[1];
    int v2 = it[2];
    NumericVector cover0 = xlist[v0];
    NumericVector cover1 = xlist[v1];
    NumericVector cover2 = xlist[v2];

    double intersection = (double)(intersect(intersect(cover0, cover1), cover2).size());

    if (intersection > 0)
    {
      double cover0_mean = find_index_mean(cover0, v0, sample_order, sample_index_mean);
      double cover1_mean = find_index_mean(cover1, v1, sample_order, sample_index_mean);
      double cover2_mean = find_index_mean(cover2, v2, sample_order, sample_index_mean);

      std::vector<std::pair<int, double>> vertices;
      vertices.push_back(std::make_pair(v0, cover0_mean));
      vertices.push_back(std::make_pair(v1, cover1_mean));
      vertices.push_back(std::make_pair(v2, cover2_mean));

      sort(vertices.begin(), vertices.end(), sortbysec<int, double>);
      two_simplices(vertices[0].first, vertices[1].first, vertices[2].first) = 1;
    }
  }

  return List::create(Named("one_simplices") = one_simplices,
                      Named("two_simplices") = two_simplices);
}
