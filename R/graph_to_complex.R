#' Generates a simplicial complex from a graph.
#'
#' Takes the adjacency matrix of a graph and returns a simplicial complex that can be used
#' by other functions in the \code{RayleighSelection} package.
#'
#' @param adjacency a weighted adjacency matrix.
#' @return An object of the class \code{simplicial}. The class \code{simplicial} inherits from
#' the class \code{igraph}.
#' @examples
#' library(RayleighSelection)
#' # Load MNIST dataset
#' data("mnist")
#'
#' # Compute a correlation matrix for a subset of the LFW dataset using only pixels with high variance
#' mnist_test <- mnist[,1:800]
#' mnist_test_top <- mnist_test[apply(mnist_test, 1, var) > 0.9,]
#' mnist_test_distances <- cor(mnist_test_top)
#'
#' # Compute a simplifical complex
#' gg <- graph_to_complex(mnist_test_distances)
#'
#' # Plot the skeleton of the simplicial complex colored by the intensity of the 500th pixel
#' plot_skeleton(gg, k=as.numeric(mnist_test[500,]), spring.constant = 0.02, seed = 3)
#'
#' @export
#'
graph_to_complex <- function(adjacency)
{
  # Builds igraph object
  diag(adjacency) <- 0
  g2 <- graph.adjacency(adjacency, mode='undirected', weighted=TRUE)

  # Enriches the class with information about the open cover and adjacency matrix
  g2$points_in_vertex <- seq(1, sqrt(length(adjacency)), 1)
  g2$order <- seq(1, sqrt(length(adjacency)), 1)
  g2$adjacency <- get.adjacency(g2, sparse = TRUE)
  g2$adjacency[lower.tri(g2$adjacency)] <- 0
  #g2$one_simplices <- g2$adjacency
  g2$two_simplices <- list()
  siz <- sqrt(length(g2$adjacency))
  for (i in 1:siz){
    g2$two_simplices[[i]] <- Matrix(0, siz, siz, sparse=TRUE)
  }
  for (m in suppressWarnings(cliques(g2,3,3))) {
    g2$two_simplices[[as.numeric(m[1])]][as.numeric(m[2]), as.numeric(m[3])] <- 1
  }

  # Decorates the graph to include node sizes, color, etc.
  V(g2)$size <- 4.0
  V(g2)$label <- ""
  V(g2)$frame.color <- "black"

  class(g2) <-  c('simplicial', 'igraph')
  return(g2)
}
