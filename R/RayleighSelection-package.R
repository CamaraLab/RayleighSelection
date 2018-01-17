#' @useDynLib RayleighSelection
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import data.table
#' @import Matrix
#' @import igraph
#' @import parallel
#' @import TDAmapper
#' @import grid
#' @import ggplot2

# This is a wrapper of the function pushCpp() only used for internal purposes.
push <- function(lo, g2, pushforward) {
  if (!missing(lo)) {
    lo[is.na(lo)] <- 0
    if (sum(abs(lo))!=0) {
      r1 <- log2(1.0 + pushCpp(as.numeric(lo), g2$points_in_vertex, 0))
      return((r1-min(r1))/(max(r1)-min(r1)))
    } else {
      return(0)
    }
  } else {
    return(0)
  }
}


# Multiple plot function used only for internal purposes (taken from http://www.cookbook-r.com).
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
