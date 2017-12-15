#' Generates the nerve complex of an open cover.
#'
#' Takes the open cover of a finite set and returns its nerve
#' complex. The output is a simplicial complex that can be used by other functions in the
#' \code{RayleighSelection} package.
#'
#' @param open_cover a list of open sets with \code{open_cover[[i]]} a vector of indices of points
#' corresponding to the \code{i}-th open set in the open cover.
#' @return An object of the class \code{simplicial}. The class \code{simplicial} inherits from
#' the class \code{igraph}.
#'
#' @export
#'
nerve_complex <- function(open_cover) {
  # Builds adjacency matrix and igraph object
  adjacency <- Matrix(adjacencyCpp(open_cover), sparse = TRUE)
  diag(adjacency)<-1
  g2 <- graph.adjacency(adjacency, mode="undirected")

  # Prunes vertices that only one element and simplifies the graph
  kal <- which(lapply(open_cover, length) %in% 1)
  open_cover[kal] <- NULL
  g2 <- delete_vertices(g2, kal)
  g2 <- contract(g2, match(open_cover, unique(open_cover)))
  g2 <- simplify(g2)
  open_cover <- unique(open_cover)

  # Decorates the graph to include node sizes, color, etc.
  V(g2)$size <- log2(as.numeric(lapply(open_cover, length))+1)*2.0
  V(g2)$label <- ""
  V(g2)$frame.color <- "black"

  # Enriches the class with information about the open cover and adjacency matrix
  g2$points_in_vertex <- open_cover
  g2$adjacency <- get.adjacency(g2, sparse = TRUE)
  class(g2) <- c('simplicial', 'igraph')
  return(g2)
}
