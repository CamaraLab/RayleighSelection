#' Generates a k-nearest neighbor complex from a distance matrix.
#'
#' Takes a distance matrix, computes a k-nearest neighbor graph, and returns
#' a simplicial complex that can be used by other functions in the
#' \code{RayleighSelection} package.
#'
#' @param dist a distance matrix.
#' @param k nunber of nearest neighbors
#' @return An object of the class \code{simplicial}. The class \code{simplicial} inherits from
#' the class \code{igraph}.
#' @examples
#' library(RayleighSelection)
#' # Load MNIST dataset
#' data("mnist")
#'
#' # Compute a distance matrix for a subset of the MNIST dataset
#' mnist_test <- mnist[,1:800]
#' mnist_test_distances <- (1.0 - cor(mnist_test))
#'
#' # Compute nearest neighbor complex with k = 5
#' gg <- knn_complex(mnist_test, 5)
#'
#' # Plot the skeleton of the nearest neighbor graph colored by the intensity of the 500th pixel
#' plot_skeleton(gg, k=as.numeric(mnist_test[500,]), spring.constant = 1.8)
#'
#' @export
#'
knn_complex <- function(dist, k)
{
  kneighbors <- kNN(dist, k)
  id <- kneighbors$id
  nodes <- dim(id)[1]
  rows <- 1:nodes
  adjacency <- matrix(0, nrow=nodes, ncol=nodes)

  for (col in 1:k)
  {
    adjacency[cbind(rows, id[,col])] <- 1
    adjacency[cbind(id[,col], rows)] <- 1
  }

  return (graph_to_complex(adjacency))
}
