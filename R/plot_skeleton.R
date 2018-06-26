#' Plots the 1-skeleton of a simplicial complex.
#'
#' Plots the 1-skeleton graph of a nerve complex and colors the vertices according
#' to the value of one or more functions with support on the underlying set of points. The 1-skeleton
#' is visualized using a force-directed graph layout with node sizes proportional to the
#' number of points in the corresponding open set.
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
#' library(RayleighSelection)
#' # Load pre-processed MNIST test dataset
#' data("mnist")
#'
#' # Compute reduced representation using Laplacian eigenmap of pixels with high variance
#' library(dimRed)
#' leim <- LaplacianEigenmaps()
#' mnist_top <- lfw[apply(mnist, 1, var) > 10000,]
#' emb <- leim@fun(as(t(mnist_top), "dimRedData"), leim@stdpars)
#'
#' # Compute Mapper representation using the Laplacian eigenmap as an auxiliary function and correlation
#' # distance as metric
#' library(TDAmapper)
#' mnist_distances <- (1.0 - cor(mnist_top))
#' m2 <- mapper2D(distance_matrix = mnist_distances,
#'                filter_values = list(emb@data@data[,1], emb@data@data[,2]),
#'                num_intervals = c(50,50),
#'                percent_overlap = 35,
#'                num_bins_when_clustering = 10);
#'
#' # Compute the nerve complex
#' gg <- nerve_complex(m2$points_in_vertex)
#'
#' # Plot the skeleton of the nerve complex colored by the intensity of the 201st pixel
#' plot_skeleton(gg, k=as.numeric(mnist[201,]))
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

  plot(g2, layout = layout.forceatlas2(as.undirected(g2), directed=TRUE, iterations=iterations, plotstep=0),
       edge.arrow.size=0.4)
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
