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
#' @param spring.constant number specifying the value of the spring constant used in the
#' force-directed layout. By default is set to 0.3.
#' @param file if specified, exports the 1-skeleton to graphviz DOT \code{file}
#' @examples
#' library(RayleighSelection)
#' # Load pre-processed LFW dataset (aligned, cropped, and normalized)
#' data("lfw")
#'
#' # Compute reduced representation using Laplacian eigenmap of pixels with high variance
#' library(dimRed)
#' leim <- LaplacianEigenmaps()
#' lfw_top <- lfw[apply(lfw, 1, var) > 0.9,]
#' emb <- leim@fun(as(t(lfw_top), "dimRedData"), leim@stdpars)
#'
#' # Compute Mapper representation using the Laplacian eigenmap as an auxiliary function and correlation
#' # distance as metric
#' library(TDAmapper)
#' lfw_distances <- (1.0 - cor(lfw_top))
#' m2 <- mapper2D(distance_matrix = lfw_distances,
#'                filter_values = list(emb@data@data[,1], emb@data@data[,2]),
#'                num_intervals = c(40,40),
#'                percent_overlap = 30,
#'                num_bins_when_clustering = 10);
#'
#' # Compute the nerve complex
#' gg <- nerve_complex(m2$points_in_vertex)
#'
#' # Plot the skeleton of the nerve complex colored by the intensity of the 1063rd pixel
#' plot_skeleton(gg, k=as.numeric(lfw[1063,]), spring.constant = 3.0, seed = 3)
#'
#' @export
#'
plot_skeleton <- function(g2, r, g, b, k, pushforward = mean, seed = 10, spring.constant = 0.3, file="") {
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
  plot(g2, layout = layout_with_graphopt(g2, spring.constant = spring.constant))
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
