#' @useDynLib RayleighSelection
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import data.table
#' @import Matrix
#' @import igraph
#' @import parallel
#' @import ggplot2
#' @import ForceAtlas2
#' @import dbscan

# This is a wrapper of the function pushCpp() only used for internal purposes.
push <- function(lo, g2, pushforward) {
  if (!missing(lo)) {
    lo[is.na(lo)] <- 0
    if (sum(abs(lo))!=0) {
      r1 <- pushCpp(as.numeric(lo), g2$points_in_vertex, 0, g2$adjacency, FALSE)
      return(log2(1.0 + ((r1$vertices-min(r1$vertices))/(max(r1$vertices)-min(r1$vertices)))))
    } else {
      return(0)
    }
  } else {
    return(0)
  }
}

