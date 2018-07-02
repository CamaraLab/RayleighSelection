#' Generates a Vietoris-Rips complex from a distance matrix.
#'
#' Takes a distance matrix and returns its Vietoris-Rips complex at a given level.
#' The output is a clique complex that can be used by other functions in the
#' \code{RayleighSelection} package.
#'
#' @param dist a distance matrix.
#' @param epsilon filtration distance at which the complex is generated
#' @return An object of the class \code{simplicial}. The class \code{simplicial} inherits from
#' the class \code{igraph}.
#' @examples
#' library(RayleighSelection)
#' # Load MNIST example dataset
#' data("mnist")
#'
#' # Compute a distance matrix for a subset of the MNIST dataset using only pixels with high variance
#' mnist_test <- mnist[,1:800]
#' mnist_test_top <- mnist_test[apply(mnist_test, 1, var) > 0.9,]
#' mnist_test_distances <- (1.0 - cor(mnist_test_top))
#'
#' # Compute Vietoris-Rips complex at filtration distance 0.5
#' gg <- vr_complex(mnist_test_distances, 0.5)
#'
#' # Plot the skeleton of the Vietoris-Rips complex colored by the intensity of the 500th pixel
#' plot_skeleton(gg, k=as.numeric(mnist_test[500,]))
#'
#' @export
#'
vr_complex <- function(dist, epsilon)
{
  # Builds adjacency matrix and igraph object
  idx <- (dist <= epsilon)
  return(graph_to_complex(idx))
}
