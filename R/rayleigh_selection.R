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
#' @param one_forms when set TRUE the Combinatorial Laplacian Score for 1-forms is
#' also computed. By default is set to FALSE.
#' @param weights when set to TRUE it takes 2-simplices into account when computing weights.
#' By default is set to FALSE.
#'
#' @return Returns a data frame with the value of the Combinatorial Laplacian Score for 0- and 1-forms,
#' the p-values, and the q-values computed using Benjamini-Hochberg procedure.
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
rayleigh_selection <- function(g2, f, num_perms = 1000, seed = 10, num_cores = 1, one_forms = FALSE,
                               weights = FALSE) {
  # Check class of f
  if (!is(f,'matrix') && !is(f,'Matrix')) {
    if (is(f,'numeric')) {
      f <- t(as.matrix(f))
    } else {
      f <- as.matrix(f)
    }
  }

  # If system is Windows, can't use mclapply
  if (Sys.info()['sysname'] == "Windows" && num_cores > 1) {
    stop("Windows does not support parallel::mclapply, please set num_cores to 1")
  }

  lout <- combinatorial_laplacian(g2, one_forms, weights)
  col <- lout$l0
  zero_weights <- lout$zero_weights
  one_weights <- lout$one_weights
  diji <- lout$adjacency_ordered
  if (one_forms) {
    l1_up <- lout$l1up
    l1_down <- lout$l1down
  }

  # Evaluates R and p for a feature fo
  cornel <- function(fo) {
    kmn<-pushCpp(as.numeric(fo), g2$points_in_vertex, num_perms, g2$adjacency, one_forms)
    kk <- kmn$vertices
    kk <- kk-matrix(rep(kk%*%zero_weights/sum(zero_weights),dim(kk)[2]),dim(kk))
    qlom <- rowSums(t(zero_weights*t(kk^2)))
    if (sum(abs(qlom))==0.0) {
      qt <- rep(Inf,length(qlom))
    } else {
      qt <- rowSums((kk%*%col)*kk)/qlom
      qt[is.nan(qt)] <- Inf
    }
    ph <- NULL
    ph$R0 <- qt[1]
    ph$p0 <- (sum(qt<=qt[1])-1.0)/num_perms
    if (one_forms) {
      kkv <- kmn$edges[,order(diji)]
      kkv <- kkv-matrix(rep(kkv%*%one_weights/sum(one_weights),dim(kkv)[2]),dim(kkv))
      qlomv <- rowSums(t(one_weights*t(kkv^2)))
      if (sum(abs(qlomv))==0.0) {
        qtv <- rep(Inf,length(qlomv))
      } else {
        qtv <- rowSums((kkv%*%(l1_up+l1_down))*kkv)/qlomv
        qtv[is.nan(qtv)] <- Inf
      }
      ph$R1 <- qtv[1]
      ph$p1 <- (sum(qtv<=qtv[1])-1.0)/num_perms
    }
    return(ph)
  }

  # Each worker evaluates R anp p for a set fu of features
  worker <- function(fu) {
    qh <- NULL
    qh$R0 <- NULL
    qh$p0 <- NULL
    if (one_forms) {
      qh$R1 <- NULL
      qh$p1 <- NULL
    }
    for (m in 1:nrow(fu)) {
      d <- cornel(fu[m,])
      qh$R0 <- rbind(qh$R0, d$R0)
      qh$p0 <- rbind(qh$p0, d$p0)
      if (one_forms) {
        qh$R1 <- rbind(qh$R1, d$R1)
        qh$p1 <- rbind(qh$p1, d$p1)
      }
    }
    return(data.frame(qh, row.names = row.names(fu)))
  }

  if (num_cores > nrow(f)) {
    num_cores <- nrow(f)
  }

  if (num_cores == 1 || nrow(f) == 1) {
    qqh <- worker(f)
  } else {
    # If more than one core then split the features in num_cores parts accordingly
    wv <- floor(nrow(f)/num_cores)
    wr <- nrow(f) - wv*num_cores
    work <- list()
    if (wr>0) {
      for (m in 1:wr) {
        work[[m]] <- (f[(1+(m-1)*(wv+1)):(m*(wv+1)),])
      }
      for (m in (wr+1):num_cores) {
        work[[m]] <- (f[(1+wr+(m-1)*wv):(wr+m*wv),])
      }
    } else {
      for (m in 1:num_cores) {
        work[[m]] <- (f[(1+(m-1)*wv):(m*wv),])
      }
    }
    reul <- mclapply(work, worker, mc.cores = num_cores)
    qqh <- reul[[1]]
    for (m in 2:num_cores) {
      qqh <- rbind(qqh, reul[[m]])
    }
  }

  # Adjust for multiple hypothesis testing
  qqh$q0 <- p.adjust(qqh$p0, method = 'BH')
  if (one_forms) {
    qqh$q1 <- p.adjust(qqh$p1, method = 'BH')
  }

  if (one_forms) {
    return(qqh[,c(1,2,5,3,4,6)])
  } else {
    return(qqh)
  }
}
