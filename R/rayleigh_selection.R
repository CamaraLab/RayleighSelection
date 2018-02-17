#' Ranks features using Rayleigh quotients.
#'
#' Given a nerve complex and a set of features consisting on functions with support on
#' the set of points underlying the complex, it asseses the significance of each feature
#' in the simplicial complex by computing the Rayleigh quotient of the function with
#' respect to the combinatorial Laplace operator of the complex and comparing it with
#' the null distribution which results from reshufling the values of the function across
#' the point cloud.
#'
#' @param g2 an object of the class \code{simplicial} containing the nerve complex.
#' @param f a numeric vector or matrix specifying one or more functions with support on
#' the set of points whose significance will be assesed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param shift real number specifying a shift that is added to \code{f}. Shifts are useful to
#' reduce the statistical power and rank very significant features without the need of using a
#' too large number of permutations. By default is set to 0.
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param seed integer specifying the seed used to initialize the generator of permutations.
#' By default is set to 10.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#'
#' @return Returns a data frame with the value of the Rayleigh quotient score, its p-values, and
#' its value adjusted for multiple hypotheis testing using Benjamini-Hochberg procedure for each
#' feature.
#' @examples
#' library(RayleighSelection)
#' # Load pre-processed MNIST test dataset
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
#'                num_intervals = c(50,50),
#'                percent_overlap = 35,
#'                num_bins_when_clustering = 10);
#'
#' # Compute the nerve complex
#' gg <- nerve_complex(m2$points_in_vertex)
#'
#' # Compute R score, p-value, and q-value for the pixels 201st to 205th
#' rayleigh_selection(gg, mnist[201:205,])
#'
#' @export
#'
rayleigh_selection <- function(g2, f, shift = 0.0, num_perms = 1000, seed = 10, num_cores = 1) {
  # Check class of f
  if (class(f) != 'data.frame') {
    f <- as.data.frame(f)
  }

  # Degree matrix
  dd <- rowSums(g2$adjacency)

  # Scalar Laplace operator
  col <- diag(dd)-g2$adjacency

  # Evaluates R and p for a feature fo
  cornel <- function(fo) {
    kk<-pushCpp(as.numeric(fo), g2$points_in_vertex, num_perms)
    qlom <- rowSums(dd*kk^2)
    if (sum(abs(qlom))==0.0) {
      qt <- qlom
    } else {
      qt <- rowSums((kk%*%col)*kk)/qlom
    }
    ph <- NULL
    ph$R <- qt[1]
    ph$p <- (sum(qt<=qt[1])-1.0)/num_perms
    return(ph)
  }

  # Each worker evaluates R anp p for a set fu of features
  worker <- function(fu) {
    qh <- NULL
    qh$R <- NULL
    qh$p <- NULL
    for (m in row.names(fu)) {
      d <- cornel(fu[m,])
      qh$R <- rbind(qh$R, d$R)
      qh$p <- rbind(qh$p, d$p)
    }
    return(data.frame(qh, row.names = row.names(fu)))
  }
  if (num_cores == 1 || nrow(f) == 1) {
    qqh <- worker(f+shift)
  } else {
    # If more than one core then split the features in num_cores parts accordingly
    wv <- floor(nrow(f)/num_cores)
    wr <- nrow(f) - wv*num_cores
    work <- list()
    if (wr>0) {
      for (m in 1:wr) {
        work[[m]] <- (f[(1+(m-1)*(wv+1)):(m*(wv+1)),] + shift)
      }
      for (m in (wr+1):num_cores) {
        work[[m]] <- (f[(1+wr+(m-1)*wv):(wr+m*wv),] + shift)
      }
    } else {
      for (m in 1:num_cores) {
        work[[m]] <- (f[(1+(m-1)*wv):(m*wv),] + shift)
      }
    }
    reul <- mclapply(work, worker, mc.cores = num_cores)
    qqh <- reul[[1]]
    for (m in 2:num_cores) {
      qqh <- rbind(qqh, reul[[m]])
    }
  }

  # Adjust for multiple hypothesis testing
  qqh$q <- p.adjust(qqh$p, method = 'BH')
  return(qqh)
}
