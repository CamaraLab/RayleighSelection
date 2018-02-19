#' Generates a Vietoris-Rips complex from a distance matrix.
#'
#' Takes a distance matrix and returns its Vietoris-Rips complex at a given level.
#' The output is a simplicial complex that can be used by other functions in the
#' \code{RayleighSelection} package.
#'
#' @param dist a distance matrix.
#' @param epsilon filtration distance at which the complex is generated
#' @param t parameter of the exponential weights \code{exp(-d_ij**2/t)}. By default it is set to infinity.
#' @return An object of the class \code{simplicial}. The class \code{simplicial} inherits from
#' the class \code{igraph}.
#' @examples
#' library(RayleighSelection)
#' # Load pre-processed LFW dataset (aligned, cropped, and normalized)
#' data("lfw")
#'
#' # Compute a distance matrix for a subset of the LFW dataset using only pixels with high variance
#' library(dimRed)
#' lfw_test <- lfw[1:800,1:700]
#' lfw_test_top <- lfw_test[apply(lfw_test, 1, var) > 0.9,]
#' lfw_test_distances <- (1.0 - cor(lfw_test_top))
#'
#' # Compute Vietoris-Rips complex at filtration distance 0.5 and weight parameter t = 5.0
#' gg <- vr_complex(lfw_test_distances, 0.5, t=5.0)
#'
#' # Plot the skeleton of the Vietoris-Rips complex colored by the intensity of the 500th pixel
#' plot_skeleton(gg, k=as.numeric(lfw_test[500,]), spring.constant = 0.02, seed = 3)
#'
#' @export
#'
vr_complex <- function(dist, epsilon, t = Inf)
{
  # Builds adjacency matrix and igraph object
  idx <- (dist <= epsilon)
  adjacency <- idx*dist
  adjacency <- exp(-adjacency**2/t)*idx
  diag(adjacency) <- 0
  g2 <- graph.adjacency(adjacency, mode='undirected', weighted=TRUE)

  # Enriches the class with information about the open cover and adjacency matrix
  g2$points_in_vertex <- seq(1, sqrt(length(lfw_test_distances)), 1)
  g2$adjacency <- get.adjacency(g2, sparse = TRUE)

  # Decorates the graph to include node sizes, color, etc.
  V(g2)$size <- 4.0
  V(g2)$label <- ""
  V(g2)$frame.color <- "black"

  class(g2) <-  c('simplicial', 'igraph')
  return(g2)
}
