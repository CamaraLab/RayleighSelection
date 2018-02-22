#' Evaluates the significance of a feature across the parameter space of Mapper.
#'
#' A natural way to generate a nerve complex from a point cloud is to use the Mapper algorithm,
#' implemented in the package \code{TDAmapper}. The function \code{scan_mapper2D} asseses the
#' significance of a feature in the nerve complexes that result across a region of the
#' parameter space of Mapper.
#'
#' @param f a numeric vector specifying the function with support on the set of points
#' whose significance will assesed.
#' @param distance_matrix a matrix of pairwise dissimilarities associated to the high
#' dimensional point cloud.
#' @param filter_values a list of two numeric vectors containing the values of the two
#' filter functions in the point cloud.
#' @param num_intervals_min a list of two integers specifying the minimum value for the \code{num_intervals}
#' parameter of \code{mapper2D} to be used in the scan.
#' @param num_intervals_max a list of two integers specifying the maximum value for the \code{num_intervals}
#' parameter of \code{mapper2D} to be used in the scan.
#' @param num_intervals_step an integer specifying the number of steps for the \code{num_intervals}
#' parameter of \code{mapper2D} to be used in the scan. By default it is set to 10.
#' @param percent_overlap_min a number between 0 and 100 specifying the minimum value for the
#' \code{percent_overlap} parameter of \code{mapper2D} to be used in the scan.
#' @param percent_overlap_max a number between 0 and 100 specifying the maximum value for the
#' \code{percent_overlap} parameter of \code{mapper2D} to be used in the scan.
#' @param percent_overlap_steps an integer specifying the number of steps for the \code{percent_overlap}
#' parameter of \code{mapper2D} to be used in the scan. By default it is set to 10.
#' @param num_bins_when_clustering a positive integer specifying the value of the \code{num_bins_when_clustering}
#' parameter of \code{mapper2D} to be used in the scan. By default it is set to 10.
#' @param weight specificies whether to weight 1-simplicies by Jaccard distance
#' @param shift real number specifying a shift that is added to \code{f}. Shifts are useful to
#' reduce the statistical power and rank very significant features without the need of using a
#' too large number of permutations. By default is set to 0.
#' @param num_perms a positive integer specifying the number of permutations used to build the null distribution for each
#' feature. By default it is set to 1000.
#' @param seed integer specifying the seed used to initialize the generator of permutations.
#' By default is set to 10.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#' @param significance_threshold level of significance to be indicated in the output.
#'
#' @return A data frame listing the p-value, q-value, and Rayleigh quotient score of the
#' feature \code{f} for each configuration of parameters of \code{mapper2D}. The distribution
#' of p-values and Rayleigh quotient scores across the space of parameters of \code{mapper2D}
#' are visualized as heatmaps.
#' @examples
#' library(RayleighSelection)
#' # Load pre-processed MNIST test dataset
#' data("mnist")
#'
#' # Compute reduced representation using Laplacian eigenmap of pixels with high variance
#' library(dimRed)
#' leim <- LaplacianEigenmaps()
#' mnist_top <- mnist[apply(mnist, 1, var) > 10000,]
#' emb <- leim@fun(as(t(lfw_top), "dimRedData"), leim@stdpars)
#'
#' # Compute correlation distance using pixels with high variance
#' mnist_distances <- (1.0 - cor(mnist_top))
#'
#' # Evaluate the significance of the R score of the 201th pixel accross a section of the
#' # Mapper parameter space
#' scan_mapper2D(mnist[201,],
#'               mnist_distances,
#'               list(emb@data@data[,1], emb@data@data[,2]),
#'               c(25,25), c(55,55), 10, 20, 50, 10)
#'
#' @export
#'
scan_mapper2D <- function(f, distance_matrix, filter_values, num_intervals_min,
                          num_intervals_max, num_intervals_steps = 10, percent_overlap_min,
                          percent_overlap_max, percent_overlap_steps = 10,
                          num_bins_when_clustering = 10, weight = TRUE, shift = 0.0, num_perms = 1000,
                          seed = 10, num_cores = 1, significance_threshold = 0.05) {
  yy <- seq(percent_overlap_min, percent_overlap_max, length.out = percent_overlap_steps)
  xx1 <- seq(num_intervals_min[1], num_intervals_max[1], length.out = num_intervals_steps)
  xx2 <- seq(num_intervals_min[2], num_intervals_max[2], length.out = num_intervals_steps)
  r <- mapply(c, xx1, xx2)

  # Each worker evaluates a slice of the parameter space with fix value of the num_intervals parameter
  worker <- function(m) {
    th <- NULL
    th$p <- NULL
    th$R <- NULL
    th$interval_index <- NULL
    th$interval1 <- NULL
    th$interval2 <- NULL
    th$overlap <- NULL
    for (percen in yy) {
      # We use sink to avoid the verbosity of mapper2D. There must be better ways to do this.
      {sink("/dev/null"); m2 <- mapper2D(distance_matrix = distance_matrix, filter_values = filter_values, num_intervals = round(r[,m]), percent_overlap = percen, num_bins_when_clustering = num_bins_when_clustering); sink();}
      gg <- nerve_complex(m2$points_in_vertex, weight = weight)
      pr <- rayleigh_selection(gg, f, shift = shift, num_perms = num_perms)
      th$p <- c(th$p, pr$p)
      th$R <- c(th$R, pr$R)
      th$interval1 <- c(th$interval1, round(r[,m][1]))
      th$interval2 <- c(th$interval2, round(r[,m][2]))
      th$interval_index <- c(th$interval_index, m)
      th$overlap <- c(th$overlap, percen)
    }
    return(th)
  }
  resul <- mclapply(as.list(1:ncol(r)), worker, mc.cores = num_cores)

  # Puts the results in a data frame.
  tth <- resul[[1]]
  for (m in 2:ncol(r)) {
    tth$p <- c(tth$p, resul[[m]]$p)
    tth$R <- c(tth$R, resul[[m]]$R)
    tth$interval1 <- c(tth$interval1, resul[[m]]$interval1)
    tth$interval2 <- c(tth$interval2, resul[[m]]$interval2)
    tth$interval_index <- c(tth$interval_index, resul[[m]]$interval_index)
    tth$overlap <- c(tth$overlap, resul[[m]]$overlap)
  }
  tth <- as.data.frame(tth)

  # Visulizes the results
  multiplot(ggplot(data = tth, aes(x = tth$interval_index, y = tth$overlap, fill = tth$R))
            + geom_raster(interpolate = TRUE)
            + labs(title = "", x = "Intervals", y = "Overlap (%)", fill = "R")
            + theme(legend.position="bottom"),
            ggplot(data = tth, aes(x = tth$interval_index, y = tth$overlap, fill = log10(tth$p+(1.0/num_perms))))
            + geom_raster(interpolate = TRUE)
            + scale_fill_gradient2(low = rgb(1,0,0), mid = rgb(1,0.5,0), high = rgb(1,1,1), midpoint = log10((1.0/num_perms))/2)
            + labs(title = "", x = "Intervals", y = "Overlap (%)", fill = "p-value (log10)")
            + theme(legend.position="bottom")
            + stat_contour(aes(z = tth$p), breaks = significance_threshold, colour = "yellow"), layout = matrix(c(1,2), ncol=2, byrow=TRUE))
  return(tth)
}

