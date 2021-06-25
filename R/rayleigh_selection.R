#' Ranks features using the Combinatorial Laplacian Score for 0- and 1-forms.
#'
#' Given a nerve or a clique complex and a set of features consisting of functions with support on
#' the set of points underlying the complex, it asseses the significance of each feature
#' in the simplicial complex by computing its scalar and vectorial Combinatorial Laplacian
#' Score and comparing it with the null distribution that results from reshufling many times the values of
#' the function across the point cloud. For nerve complexes, feature functions induce 0- and
#' 1-forms in the complex by averaging the function across the points associated to 0- and 1-simplices
#' respectively. For clique complexes, feature functions are directly 0-forms in the complex and 1-forms
#' are obtained by averaging the function across the two vertices connected by each edge.
#'
#' The calculation of p-values can be optimized by iteratively doubling the number of samples of the
#' null distribution until convergence is reached. Two version of this iteraction scheme are implemented.
#' In the fist one, a p-value is considered convergent if there are at least 10 samples of the null
#' distribution that do not exceed the associated Combinatorial Laplacian Score. In the second one, a p-value is considered
#' convergent of the condition above holds, and, in case there are less than 10 small samples a generalized
#' Pareto distribution (GPD) is used to approximate a p-value. A p-value obtained by a GPD is considered
#' convergent if the relative variation is small in the last 3 iteractions and the quartiles
#' of the approximation are relatively close.
#'
#'
#' @param g2 an object of the class \code{simplicial} containing the nerve or clique complex.
#' @param f a numeric vector or matrix specifying one or more functions with support on
#' the set of points whose significance will be assessed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param seed integer specifying the seed used to initialize the generator of permutations.
#' By default is set to 10.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#' @param mc.preschedule when set to TRUE parallel compulations are prescheduled, see
#' \link[parallel]{mclapply}. Only has effect if \code{num_cores} > 1 and the code is not being run on Windows.
#' By default is set to TRUE.
#' @param one_forms when set TRUE the Combinatorial Laplacian Score for 1-forms is
#' also computed. By default is set to FALSE.
#' @param weights when set to TRUE it takes 2-simplices into account when computing weights.
#' By default is set to FALSE.
#' @param optimize.p string indicating the type of optimization used for computing p-values.
#' Must have value \code{NULL} for no optimization, \code{"perm"} for optimizing the calculation of
#' p-values using only permutations, or \code{"gpd"} for using a permutations and GPD in optimizing p-value calculation.
#' By default is set to \code{NULL}.
#' @param max_perm maximum number of permutations to be used when computing p-values, only
#' relevant when \code{optimize.p} is set to \code{"perm"} or \code{"gpd"}. By default is set to 16000.
#' @param pow positive number indicating the power to which the samples of the null distribution and the associated
#' score are to be transformed before computing a GPD approximation (only used when
#' \code{optimize.p} is set to \code{"gdp"}).
#' @param nextremes vector of integers with the candidate number of extremes for fitting a GDP.
#' Only used when \code{optimize.p} is set to \code{"gdp"}. By default is set to
#' \code{c(seq(50, 250, 25), seq(300, 500, 50), seq(600, 1000, 100))}.
#' @param alpha level of FDR control for choosing the number of extremes. Only used when
#' \code{optimize.p} is set to \code{"gdp"}. By default is set to 0.15.
#'
#'
#' @details When computing a p-value using a GPD, only null distribution samples in the first quartile are considered.
#' The Combinatorial Laplacian Score and associated null distribution samples are transformed by the function
#' \deqn{f(x) = (1 - (x - loc)/scale)^pow}
#' where \eqn{loc} is the first quartile of the null distribution and \eqn{scale} is the first quartile minus
#' the 5%-quantile. A number of extremes for fitting a GPD is chosen using the ForwardStop p-value adjustment, and
#' quartiles for the p-value estimates are obtained by sampling GDP parameters form a multivariate normal distribution.
#'
#' @return Returns a data frame with the value of the Combinatorial Laplacian Score for 0- and 1-forms,
#' the p-values, and the q-values computed using Benjamini-Hochberg procedure. If \code{optimize.p} is set
#' to \code{"perm"} or \code{"gpd"} then then number of samples at which convergence of p-values was obtained is
#' also returned.
#'
#' @examples
#' # Example 1
#' library(RayleighSelection)
#' gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
#' rayleigh_selection(gy,t(as.data.frame(c(0,1,1,0,0,0,0,0,0,1))), one_forms = TRUE)
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
#' # Compute R score, p-value, and q-value for the pixels 301st to 305th
#' rayleigh_selection(gg, mnist[301:305,], one_forms = TRUE)
#'
#' @export
#'

rayleigh_selection <- function(g2, f, num_perms = 1000, seed = 10, num_cores = 1,
                               mc.preschedule = TRUE, one_forms = FALSE, weights = FALSE,
                               optimize.p = NULL, max_perms = 16000, pow = 1,
                               nextremes = c(seq(50, 250, 25), seq(300, 500, 50), seq(600, 1000, 100)),
                               alpha = 0.15){
  # Check class of f
  if (!is(f,'matrix') && !is(f,'Matrix')) {
    if (is(f,'numeric')) {
      f <- t(as.matrix(f))
    } else {
      f <- as.matrix(f)
    }
  }

  if(!is.null(optimize.p) && optimize.p != "perm" && optimize.p != "gpd"){
    optimize.p <- NULL
    warning(
      "optimize.p must be either NULL, 'perm' or 'gpd'. Proceding with no p-value optimization."
      )
  }


  lout <- combinatorial_laplacian(g2, one_forms, weights)
  scorer <- new(LaplacianScorer,lout, g2$points_in_vertex, g2$adjacency, one_forms)


  dims <- if(one_forms) c(0,1) else 0 # dimension to be considered

  out <- list()
  is.windows <- Sys.info()['sysname'] == "Windows"
  use.gpd <- !is.na(optimize.p) && optimize.p == "gpd"
  use.mclapply  <- !is.windows && nrow(f) > 1 && num_cores > 1

  # Splitting functions for use with mclapply
  if(use.mclapply) f.list <- split(f, row(f))

  set.seed(seed)

  for(d in dims){
    # Names of columns
    Rd <- sprintf("R%d", d)
    nd <- sprintf("n%d.conv", d)
    pd <- sprintf("p%d", d)
    qd<-sprintf("q%d", d)

    R <- as.numeric(scorer$score(f, d)) # computing scores
    out[[Rd]] <- R

    if(is.null(optimize.p)){
      # No optimization for p values is done
      if(use.mclapply){
        # Using mclapply to parallelize sampling
        worker <- function(func) scorer$sample_scores(t(as.matrix(func)), num_perms, d, 1)

        samp.list <- parallel::mclapply(f.list, worker,
                                        mc.cores =  num_cores,
                                        mc.preschedule = mc.preschedule)
        samp <- samp.list[[1]]
        for(k in 2:nrow(f)) samp <- rbind(samp, samp.list[[k]])
      }else{
        # Paralellization (if any) is done in C++ code with OMP
        samp <- scorer$sample_scores(f, num_perms, d, num_cores)
      }
      out[[pd]] <- apply(samp <= R, 1, sum) / num_perms # computing p-values
    }else{
      # Optimization for p-values is done
      if(use.mclapply){
        # Using mcmapply to parallelize p-value approximation
        work <- function(score, func){
          return(
            optim.p(score, t(as.matrix(func)), scorer, d, use.gpd, num_perms,
                    max_perms, n.cores = 1, nextremes, alpha)
          )
        }

        p.vals <- parallel::mcmapply(work, R, f.list,
                                     mc.cores = num_cores,
                                     mc.preschedule = mc.preschedule)

        out[[pd]] <- as.numeric(p.vals["p",])
        out[[nd]] <- as.numeric(p.vals["n.conv",])
      }else{
        # Paralellization (if any) is done in C++ with OMP
        p.vals <- optim.p(R, f, scorer, d, use.gpd, num_perms,
                          max_perms, n.cores = num_cores, nextremes, alpha)
        out[[pd]] <- p.vals$p
        out[[nd]] <- p.vals$n.conv
      }
    }
    # Adjusting p-values with the Benjamini-Hochberg procedure
    out[[qd]] <- p.adjust(out[[pd]], method = 'BH')
  }
  return(as.data.frame(out))
}
