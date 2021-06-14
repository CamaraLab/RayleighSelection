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
#' @param one_forms when set FALSE only the Combinatorial Laplacian Score for 0-forms is
#' computed. By default is set to TRUE.
#' @param weights when set to TRUE it takes 2-simplices into account when computing weights.
#' By default is set to FALSE.
#' @param optimize.p string indicating the type of optimization used for computing p-values.
#' Must have value \code{NULL} for no optimization, \code{"perm"} for optimizing the calculation of
#' p-values using only permutations, or \code{"gpd"} for using a permutations and GPD in optimizing p-value calculation.
#' @param max_perm maximum number of permutations to be used when computing p-values, only
#' relevant when \code{optimize.p} is set to \code{"perm"} or \code{"gpd"}.
#' @param pow positive number indicating the power to which the samples of the null distribution and the associated
#' score are to be transformed before computing a GPD approximation (only used when
#' \code{optimize.p} is set to \code{"gdp"}).
#' @param nextremes vector of integers with the candidate number of extremes for fitting a GDP
#' (only used when \code{optimize.p} is set to \code{"gdp"}).
#' @param alpha level of FDR control for fitting a GPD (only used when \code{optimize.p} is set to \code{"gdp"}).
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
#' rayleigh_selection(gy,t(as.data.frame(c(0,1,1,0,0,0,0,0,0,1))))
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
#' rayleigh_selection(gg, mnist[301:305,])
#'
#' @export
#'
rayleigh_selection <- function(g2, f, num_perms = 1000, seed = 10, num_cores = 1,
                               one_forms = TRUE, weights = FALSE,
                               optimize.p = NULL, max_perms = 10000,
                               pow = 1,
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

  stopifnot( is.null(optimize.p) || optimize.p == "perm" || optimize.p == "gpd" )

  gpd.approx <- function(samples, R, n.samp.quar = 250){
    # Computes gpd approximations of p-values and quartiles of these approximations
    #
    # Arguments
    #  samples: samples of scores obtained by permutation
    #  R: scores to have p-values computed
    #  n.samp.quar: number of samples taken when computing quartiles
    #
    # Return
    #  data frame with the natural log of the p-values and the first ant third quartiles of the approximation

    out <- data.frame(log.p = rep(NA, length(R)),
                      log.q1 = rep(NA, length(R)),
                      log.q3 = rep(NA, length(R)))

    # Only consider values in the lowers quantile for fitting a gpd
    num.ext <- sort(nextremes[nextremes*4 < ncol(samples)])
    if(length(num.ext) == 0) return(out)

    for(i in 1:length(R)){
      # Transform samples and scores
      samp.quant <- quantile(samples[i,], probs = c(0.05, 0.25), names = FALSE)
      samp.loc <- samp.quant[2]
      samp.scale <- samp.quant[2] - samp.quant[1]
      samp <- ( 1 - (samples[i, samples[i,] <= samp.loc] - samp.loc)/samp.scale )**pow
      score <- ( 1 - (R[i] - samp.loc)/samp.scale)^pow

      # Using ForwardStop to select the number of extreme samples
      gpd.test <- tryCatch(
                    eva::gpdSeqTests(samp, nextremes = num.ext, method = "ad"),
                    error = function(e) NULL)
      if(is.null(gpd.test)) next
      if(all(gpd.test$ForwardStop > alpha)){
        n.selected <- num.ext[length(num.ext)]
      }else if(any(gpd.test$ForwardStop > alpha)){
        if(gpd.test$ForwardStop[1] <= alpha) next
        n.selected <- num.ext[which(gpd.test$ForwardStop <= alpha)[1] - 1]
      }else next

      # Fitting gpd and computing (log of) p-value
      fit <- eva::gpdFit(samp, nextremes = n.selected)
      loc <- unname(min(fit$exceedances))
      scale <- unname(fit$par.ests[1])
      shape <- unname(fit$par.ests[2])
      if(scale <= 0 | ( (score - loc) / scale )*shape <= -1 ) next
      if(shape == 0){
        log.p <- (score - loc) / scale
      }else{
        log.p <- (-1/shape) * log1p( ( (score - loc) / scale)*shape )
      }
      if(is.infinite(log.p) | is.na(log.p)) next

      # Sampling parameters to compute quartiles of p-value
      coef <- MASS::mvrnorm(n = n.samp.quar, mu = fit$par.ests, Sigma = fit$varcov)
      degenerate <- (coef[,1] <= 0) | ( ( (score - loc)/coef[,1] )*coef[,2] <= -1)
      coef <- coef[!degenerate,]
      coef.1 <- coef[coef[,2] != 0,]
      coef.2 <- coef[coef[,2] == 0,]
      log.p.samples <- (-1/coef.1[,2]) * log1p( ( (score - loc) / coef.1[,1] )*coef.1[,2] )
      log.p.samples <- append(log.p.samples, (score - loc) / coef.2[,1])
      log.p.samples <- append(log.p.samples, rep(-Inf, sum(degenerate)))
      q <- quantile(log.p.samples, probs = c(0.25, 0.75), names = FALSE)
      out$log.p[i] <- log.p
      out$log.q1[i] <- q[1]
      out$log.q3[i] <- q[2]
    }
    return(out)
  }

  optim.p <- function(R, dim, use.gpd){
    # Optimizes the calculation of p-values by permutation
    #
    # Arguments
    #  R: score
    #  dim: dimension of forms (0 or 1)
    #  use.gpd: Boolean indicating if a generalized pareto is to be used
    #
    # Return
    #  data frame with p-values and number of samples where convergence was achieved

    n_funcs <- nrow(f)
    idx <- 1:n_funcs # index of non-convergent functions
    p <- rep(NA, n_funcs)
    n.conv <- rep(NA, n_funcs)

    # if the score is positive infinity the p-value is one
    p[is.infinite(R) & R > 0] <- 1
    n.conv[is.infinite(R) & R > 0] <- 0
    idx<- idx[!(is.infinite(R) & R > 0)]

    n_samps <- num_perms
    samp <- scorer$sample_scores(f[idx, ,drop = F], n_samps, dim, num_cores)
    # setting up matrix to keep track of gpd p-values
    if(use.gpd) log.p.gpd <- matrix(rep(NA, 3*length(idx)), nrow = length(idx))
    while(TRUE){
      # if there are enough (at least 10) small the p-value converges
      num.le <- apply(samp <= R[idx], 1, sum)
      conv <- num.le >= 10
      p[idx[conv]] <- num.le[conv]/n_samps
      n.conv[idx[conv]] <- n_samps
      idx <- idx[!conv]
      samp <- samp[!conv, , drop = F]

      if(use.gpd & any(!conv)){
        # try gpd approximation
        log.p.gpd <- log.p.gpd[!conv, , drop = F]
        gpd.out <- gpd.approx(samp, R[idx])

        #condition 0: all values are finite
        cond0 <- apply(log.p.gpd, 1, function(x) all(is.finite(x))) &
          is.finite(gpd.out$log.q1) & is.finite(gpd.out$log.q3) &  is.finite(gpd.out$log.p)
        #condition 1: relative variation of log p values is small
        cond1 <- apply(abs(log.p.gpd/ gpd.out$log.p - 1) < 0.1, 1, all)
        #condition 3: quartiles are close to predicted values
        cond2 <- gpd.out$log.q1/gpd.out$log.p < 1.1 & gpd.out$log.q3/gpd.out$log.p > 0.9

        conv <- cond0 & cond1 & cond2
        conv[is.na(conv)] <- FALSE
        p[idx[conv]] <- exp(gpd.out$log.p[conv])
        n.conv[idx[conv]] <- n_samps
        idx <- idx[!conv]
        samp <- samp[!conv, , drop = F]
        log.p.gpd <- cbind(log.p.gpd[!conv, , drop = F], gpd.out$log.p[!conv])
        if(2*n_samps <= max_perms) log.p.gpd <- log.p.gpd[, 2:ncol(log.p.gpd)]
      }
      if(n_samps < max_perms){
        # get new samples
        num.new.samps = min(n_samps, max_perms - n_samps)
        samp <- cbind(samp, scorer$sample_scores(f[idx, , drop = F], num.new.samps, dim, num_cores))
        n_samps <- n_samps + num.new.samps
      }else{
        break
      }
    }
    p[idx] <- num.le[!conv]/n_samps
    return(data.frame(p = p, n.conv = n.conv))
  }

  lout <- combinatorial_laplacian(g2, one_forms, weights)

  scorer <- new(LaplacianScorer,lout, g2$points_in_vertex, g2$adjacency, one_forms)

  R0 <- scorer$score(f, 0)
  out <- data.frame(R0 = R0, row.names = row.names(f))

  if(is.null(optimize.p)){
    samp <- scorer$sample_scores(f, num_perms, 0, num_cores)
    out$p0 <- apply(samp <= out$R0, 1, sum) / num_perms
  }else{
  p.vals <- optim.p(R0, 0, optimize.p == "gpd")
  out$p0 <- p.vals$p
  out$n0.conv <- p.vals$n.conv
  }
  out$q0 <- p.adjust(out$p0, method = 'BH')

  if(!one_forms) return(out)

  out$R1 <- scorer$score(f, 1)[,1]
  if(is.null(optimize.p)){
    samp <- scorer$sample_scores(f, num_perms, 1, num_cores)
    out$p1 <- apply(samp <= out$R1, 1, sum) / num_perms
  }else{
  p.vals <- optim.p(out$R1, 1, optimize.p == "gpd")
  out$p1 <- p.vals$p
  out$n1.conv <- p.vals$n.conv
  }
  out$q1 <- p.adjust(out$p1, method = 'BH')
  return(out)

}

#needed to load module
Rcpp::loadModule("mod_laplacian", TRUE)
