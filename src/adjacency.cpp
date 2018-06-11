// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;


// taken from:
// https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
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


/* Given an open set and a vector (feature_order) which assigns an order to the features contained in the open
 * set this function returns the mean of the orders of the features in the open set, this is used to
 * assign an order to the open sets (ie the vertices of the nerve complex)
 *                                  vertices = < 1, 2, 3, 4, 5 >
 * ex: open_set = < 2, 4, 5 >, feature_order = < 3, 1, 5, 2, 3 >
       computex_vertex_index_mean(<2, 4, 5 >, < 3, 1, 5, 2, 3 >) returns (1 + 2 + 3)/3
 */
double compute_vertex_index_mean(NumericVector& open_set, IntegerVector& feature_order)
{
  double mean = 0.0;

  for(NumericVector::iterator it = open_set.begin(); it != open_set.end(); ++it)
  {
    IntegerVector::iterator idx_it = std::find(feature_order.begin(), feature_order.end(), *it);
    int idx = std::distance(feature_order.begin(), idx_it);
    mean += (idx+1);
  }

  return mean/((double) open_set.size());
}

/* Given a vertex and a vector which assigns an order to all vertices return the order
 * ex: vertices = < 1, 2, 3, 4 >
                    |  |  |  |   => vertex order is 3 < 4 < 1 < 2
           rank = < 3, 4, 1, 2 >

   find_rank(4, rank) returns 2
 */
int find_rank(int vertex,
               std::vector<int>& rank)
{
  std::vector<int>::iterator it = std::find(rank.begin(), rank.end(), vertex+1);
  return std::distance(rank.begin(), it);
}

/* Given an open cover x, it computes the adjacency matrix given by pairwise
 * intersections of the open sets.
 */
// [[Rcpp::export]]
List adjacencyCpp(List x, IntegerVector& feature_order, bool weight) {
  List xlist(x);
  int n = xlist.size();
  arma::sp_mat one_simplices(n,n);

  List two_simplices(n);
  two_simplices.fill(arma::sp_mat(n, n));

  // for a given vertex i in the nerve complex get the associated open set and compute the mean of the
  // indices (ie the orders) of the features in the open set, feature_index_mean will then be a vector of
  // pairs < i, mean >
  std::vector<std::pair<int, double>> feature_index_mean;

  for(int i = 0; i < n; i++)
  {
    NumericVector cover = xlist[i];
    double mean = compute_vertex_index_mean(cover, feature_order);
    feature_index_mean.push_back(std::pair<int, double>(i+1, mean));
  }

  // sort all pairs < i, mean > by the mean and return the i's to get the order of the vertices in the complex
  sort(feature_index_mean.begin(), feature_index_mean.end(), sortbysec<int, double>);
  std::vector<int> rank;
  std::transform(feature_index_mean.begin(), feature_index_mean.end(), back_inserter(rank),
                 (const int& (*)(const std::pair<int, double>&))std::get<0>);

  // build the adjacency matrix for the one simplices by computing pairwise intersections
  for(int i=0; i<n-1; i++) {

    NumericVector loc1 = xlist[i];

    double loc1_size = (double)loc1.size();

    for(int j=i+1; j<n; j++) {

      NumericVector loc2 = xlist[j];
      double loc2_size = (double)loc2.size();

      double intersection = (double)(intersect(loc1,loc2).size());

      if (intersection > 0) {
        // the adjancency matrix is directed and will go from the vertex with smaller rank to larger rank
        // note the rows and columns of the adjancency matrix will be labeled according to rank, ie if the
        // rank is 3 < 4 < 1 < 2 then one_simplices will be:
        //  |3 4 1 2
        // -|-------
        // 3|0 . . .
        // 4|0 0 . .  where only . will be non-zero
        // 1|0 0 0 .
        // 2|0 0 0 0

        int idxi = find_rank(i, rank);
        int idxj = find_rank(j, rank);

        if (idxi > idxj)
        {
          std::swap(idxi, idxj);
          std::swap(i, j);
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
    // the outer for loop loops over the rows in the one simplices adjancency matrix, but
    // since the one simplices adjacency matrix is labeled by the rank we need to do some
    // indirection for indexing into the open cover to get the associated open set with each vertex
    // ex: if the rank is 3 < 4 < 1 < 2 then when i = 0, we want to get the open set associated with vertex 3
    NumericVector cover0 = xlist[rank[i]-1];

    // next we want to see how many 0-simplices the given vertex is adjancent to. if it is not adjancent to
    // at least 2 other 0-simplices then it cannot belong to a 2-simplex
    arma::mat row(one_simplices.row(i));
    arma::uvec idxs = arma::find(row != 0);

    if(idxs.size() < 2)
    {
      continue;
    }

    // for a 2-simplex that contains vertex rank[i] there will be one 1-simplex that does not contain
    // rank[i] we want to generate all such 1-simplices and check that they are contained in the complex
    // ex: for the 2-simplex <1, 2, 3> we need to check that the edge <2, 3> exists since from the
    // above checks we already know <1, 2> and <1, 3> exist
    std::vector<std::vector<int>> edge_idxs = combo(idxs.size(), 2);

    for(auto&& e : edge_idxs)
    {
      int j = idxs[e[0]];
      int k = idxs[e[1]];

      if(one_simplices(j, k) == 0 && one_simplices(k, j) == 0)
      {
        continue;
      }

      // if the faces of the 2-simplex are contained, compute the triple intersection
      // if it is non-zero insert the two-simplex
      NumericVector cover1 = xlist[rank[j]-1];
      NumericVector cover2 = xlist[rank[k]-1];

      double intersection = (double)(intersect(intersect(cover0, cover1), cover2).size());

      if (intersection == 0)
      {
        continue;
      }

      // Rcout << "Inserting two simplex: "
      //      << rank[i] << " "
      //      << rank[j] << " "
      //      << rank[k] << std::endl;

      // two_simplices is a list of sparse matrices. for a given 2-simplex <i, j, k> the first index i is
      // the i'th sparse matrix in the list, the second indices (j, k) are the row and columns for that matrix
      arma::sp_mat adjacency = two_simplices(rank[i]-1);
      adjacency(rank[j]-1, rank[k]-1) = 1;
      two_simplices(rank[i]-1) = adjacency;
    }
  }

  return List::create(Named("one_simplices") = one_simplices,
                      Named("two_simplices") = two_simplices,
                      Named("order") = rank);
}
