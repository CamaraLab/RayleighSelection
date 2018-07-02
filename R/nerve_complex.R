#' Generates the nerve complex of an open cover.
#'
#' Takes the open cover of a finite set and returns its nerve
#' complex. The output is a nerve complex that can be used by other functions in the
#' \code{RayleighSelection} package. The orientation of simplices is based on the order of the
#' points in the point cloud.
#'
#' @param open_cover a list of open sets with \code{open_cover[[i]]} a vector of indices of points
#' corresponding to the \code{i}-th open set in the open cover.
#' @return An object of the class \code{simplicial}. The class \code{simplicial} inherits from
#' the class \code{igraph}.
#' @examples
#' # Example 1
#' library(RayleighSelection)
#' gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
#' plot_skeleton(gy,k=c(0,1,1,0,0,0,0,0,0,1))
#'
#'
#' # Example 2: MNIST dataset
#' data("mnist")
#'
#' # Compute reduced representation using Laplacian eigenmap of pixels with high variance
#' library(dimRed)
#' leim <- LaplacianEigenmaps()
#' mnist_top <- mnist[apply(mnist, 1, var) > 10000,]
#' emb <- leim@fun(as(t(mnist_top), "dimRedData"), leim@stdpars)
#'
#' # Compute Mapper representation using the Laplacian eigenmap as an auxiliary function and correlation
#' # distance as metric
#' library(TDAmapper)
#' mnist_distances <- (1.0 - cor(mnist_top))
#' m2 <- mapper2D(distance_matrix = mnist_distances,
#'                filter_values = list(emb@data@data[,1], emb@data@data[,2]),
#'                num_intervals = c(30,30),
#'                percent_overlap = 35,
#'                num_bins_when_clustering = 10);
#'
#' # Compute the nerve complex
#' gg <- nerve_complex(m2$points_in_vertex)
#'
#' # Plots simplicial complex colored by the value of the 301th pixel
#' plot_skeleton(gg, k=mnist[301,])
#'
#' @export
#'
nerve_complex <- function(open_cover) {
  # Builds adjacency matrix and igraph object
  kal <- which(lapply(open_cover, length) %in% 1)
  open_cover[kal] <- NULL
  open_cover <- unique(open_cover)

  feature_order <- 1:max(unlist(open_cover))
  complex <- adjacencyCpp(open_cover, feature_order)
  adjacency <- Matrix(complex$one_simplices, sparse = TRUE)

  sorted_order <- order(as.numeric(complex$order))
  adjacency <- adjacency[sorted_order, sorted_order]

  g2 <- graph.adjacency(adjacency, mode="directed", weighted='NULL')

  # Decorates the graph to include node sizes, color, etc.
  V(g2)$size <- log2(as.numeric(lapply(open_cover, length))+1)
  V(g2)$label <- ""
  V(g2)$frame.color <- "black"

  # Enriches the class with information about the open cover and adjacency matrix
  g2$points_in_vertex <- open_cover
  g2$adjacency <- get.adjacency(g2, sparse = TRUE)
  g2$one_simplices <- complex$one_simplices
  g2$two_simplices <- complex$two_simplices
  g2$order <- complex$order

  class(g2) <- c('simplicial', 'igraph')
  return(g2)
}
