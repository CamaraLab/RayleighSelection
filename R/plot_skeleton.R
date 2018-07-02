#' Plots a simplicial complex.
#'
#' Plots the 0-, 1- and 2-simplices of a simplicial complex and colors the vertices according
#' to the value of one or more functions with support on the underlying set of points. The complex
#' is visualized using a force-directed graph layout acting on the 1-skeleton. Node sizes are
#' proportional to the number of points in the corresponding open set.
#'
#' @param g2 an object of the class \code{simplicial} containing the nerve complex.
#' @param r a numeric vector or matrix specifying one or more functions with
#' support on the set of points that will be mapped to the red channel. Each column
#' corresponds to a point and each row specifies a different function. If more than one
#' function is specified, the average of all rows is considered.
#' @param g similar to \code{r} but mapped to the green channel.
#' @param b similar to \code{r} but mapped to the blue channel.
#' @param k similar to \code{r} but mapped to a rainbow palette. Incompatible with
#' the parameters \code{r}, \code{g}, and \code{b}.
#' @param pushforward pushforward function that maps \code{r}, \code{g}, \code{b},
#' and \code{k} to functions with support on the simplices. By default is set to
#' \code{mean}. Other choices may slow down the computation substantially.
#' @param seed integer specifying the seed used to initialize the force-directed layout. By
#' default is set to 10.
#' @param iterations number of iterations used in the Force Atlas 2 layout. By default is set
#' to 1500. A larger value may be required for optimal visualization of large graphs.
#' @param file if specified, exports the 1-skeleton to graphviz DOT \code{file}
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
plot_skeleton <- function(g2, r, g, b, k, pushforward = mean, seed = 10, iterations = 1500, file="") {
  set.seed(seed)
  if (!missing(k)) {
    k1 <- push(k, g2, pushforward)
    V(g2)$color = rainbow(256, end=0.7)[256-round(k1*255)]
  } else {
    k1 <- push(r, g2, pushforward)
    k2 <- push(g, g2, pushforward)
    k3 <- push(b, g2, pushforward)
    V(g2)$color = rgb(1-(k2+k3)/2,1-(k1+k3)/2,1-(k1+k2)/2)
  }

  p = layout.forceatlas2(as.undirected(g2), directed=TRUE, iterations=iterations, plotstep=0)
  plot(norm_coords(p),cex=0.1,xaxt='n', yaxt='n', ann=FALSE, bty='n')
  idxs <- which(lapply(g2$two_simplices, any) == TRUE, arr.ind=T)
  for (m in idxs) {
    qwy <- which(g2$two_simplices[[m]]!=0,arr.ind = T)
    for (k in 1:nrow(qwy)) {
      polygon(norm_coords(p)[c(m, as.numeric(qwy[k,1]), as.numeric(qwy[k,2])),], col = rgb(0.75,0.85,0.95),
            border = rgb(0.75,0.85,0.95))
    }
  }
  plot(g2, layout = p, edge.arrow.size=0.4, add = TRUE)
  if (file != "") {
    V(g2)$fillcolor <- V(g2)$color
    V(g2)$fixedsize <- TRUE
    V(g2)$height <- V(g2)$size
    V(g2)$width <- V(g2)$size
    g2 <- delete_vertex_attr(g2, 'label')
    g2 <- delete_vertex_attr(g2, 'color')
    g2 <- delete_vertex_attr(g2, 'size')
    g2 <- delete_vertex_attr(g2, 'frame.color')
    write_graph(g2, file=file, format="dot")
  }
}
