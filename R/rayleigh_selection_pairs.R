#' Ranks pairs of features using the Combinatorial Laplacian Score for 0- and 1-forms.
#'
#' Given a nerve or a clique complex, a set of features consisting of functions with support on
#' the set of points underlying the complex, and a list of pairs of features,
#' it asseses the significance of each pair of features
#' in the simplicial complex by computing its scalar and vectorial Combinatorial Laplacian
#' Score and comparing it with the null distribution that results from reshufling many times the values of
#' the function across the point cloud. For nerve complexes, feature functions induce 0- and
#' 1-forms in the complex by averaging the function across the points associated to 0- and 1-simplices
#' respectively. For clique complexes, feature functions are directly 0-forms in the complex and 1-forms
#' are obtained by averaging the function across the two vertices connected by each edge.
#'
#' @param g2 an object of the class \code{simplicial} containing the nerve or clique complex.
#' @param f a numeric vector or matrix specifying one or more functions with support on
#' the set of points whose significance will be assesed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param f_pairs a 2 column matrix where each row specifes the indices or names
#' of a pair of points on which the Comb. Lap. score will be computed
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
#'
#' # Example 1
#' gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
#' features <- rbind(c(0,1,1,0,0,0,0,0,0,1), c(1,1,0,1,1,0,0,1,1,0))
#' row.names(features) <- c("f1","f2")
#' pairs <- pairs <- matrix(c("f1","f1","f2","f2","f1","f2"), ncol=2, byrow=T)
#' rayleigh_selection_pairs(gy,features,pairs)


rayleigh_selection_pairs <- function(g2, f, f_pairs, num_perms = 1000, seed = 10, num_cores = 1, one_forms = TRUE,
                               weights = FALSE) {
  # Check class of f
  if (class(f) != 'data.frame') {
    f <- as.data.frame(f)
  }

  # Check class of f_pairs
  if (class(f_pairs) == 'list') {
    f_pairs <- matrix(unlist(f_pairs), ncol=2, byrow=T)
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

  # Evaluates R and p for a pair of features fo
  cornel <- function(fo) {
    kmn1 <- pushCpp(as.numeric(f[fo[1],]), g2$points_in_vertex, num_perms, g2$adjacency, one_forms)
    kmn2 <- pushCpp(as.numeric(f[fo[2],]), g2$points_in_vertex, num_perms, g2$adjacency, one_forms)
    kk1 <- kmn1$vertices
    kk2 <- kmn2$vertices
    kk1 <- kk1-matrix(rep(kk1%*%zero_weights/sum(zero_weights),dim(kk1)[2]),dim(kk1))
    kk2 <- kk2-matrix(rep(kk2%*%zero_weights/sum(zero_weights),dim(kk2)[2]),dim(kk2))
    qlom <- rowSums(t(zero_weights*t(kk1*kk2)))
    if (sum(abs(qlom))==0.0) {
      qt <- rep(Inf,length(qlom))
    } else {
      qt <- rowSums((kk1%*%col)*kk2)/qlom
      qt[is.nan(qt)] <- Inf
    }
    ph <- NULL
    ph$R0 <- qt[1]
    ph$p0 <- (sum(qt<=qt[1])-1.0)/num_perms
    if (one_forms) {
      kkv1 <- kmn1$edges[,order(diji)]
      kkv2 <- kmn2$edges[,order(diji)]
      kkv1 <- kkv1-matrix(rep(kkv1%*%one_weights/sum(one_weights),dim(kkv1)[2]),dim(kkv1))
      kkv2 <- kkv2-matrix(rep(kkv2%*%one_weights/sum(one_weights),dim(kkv2)[2]),dim(kkv2))
      qlomv <- rowSums(t(one_weights*t(kkv1*kkv2)))
      if (sum(abs(qlomv))==0.0) {
        qtv <- rep(Inf,length(qlomv))
      } else {
        qtv <- rowSums((kkv1%*%(l1_up+l1_down))*kkv2)/qlomv
        qtv[is.nan(qtv)] <- Inf
      }
      ph$R1 <- qtv[1]
      ph$p1 <- (sum(qtv<=qtv[1])-1.0)/num_perms
    }
    return(ph)
  }

  # Each worker evaluates R anp p for a set fu of pairs of features
  worker <- function(fu) {
    qh <- NULL
    qh$R0 <- NULL
    qh$p0 <- NULL
    if (one_forms) {
      qh$R1 <- NULL
      qh$p1 <- NULL
    }
    for (i in 1:nrow(fu)) {
      d <- cornel(fu[i,])
      qh$R0 <- rbind(qh$R0, d$R0)
      qh$p0 <- rbind(qh$p0, d$p0)
      if (one_forms) {
        qh$R1 <- rbind(qh$R1, d$R1)
        qh$p1 <- rbind(qh$p1, d$p1)
      }
    }
    return(data.frame(qh))
  }

  if (num_cores > nrow(f_pairs)) {
    num_cores <- nrow(f_pairs)
  }

  if (num_cores == 1 || nrow(f_pairs) == 1) {
    qqh <- worker(f_pairs)
  } else {
    # If more than one core then split the features in num_cores parts accordingly
    wv <- floor(nrow(f_pairs)/num_cores)
    wr <- nrow(f_pairs) - wv*num_cores
    work <- list()
    if (wr>0) {
      for (m in 1:wr) {
        work[[m]] <- (f_pairs[(1+(m-1)*(wv+1)):(m*(wv+1)),])
      }
      for (m in (wr+1):num_cores) {
        work[[m]] <- (f_pairs[(1+wr+(m-1)*wv):(wr+m*wv),])
      }
    } else {
      for (m in 1:num_cores) {
        work[[m]] <- (f_pairs[(1+(m-1)*wv):(m*wv),])
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

  # Add feature columns
  qqh <- cbind(f = f_pairs[,1], g = f_pairs[,2], qqh)

  if (one_forms) {
    return(qqh[,c(1,2,3,4,7,5,6,8)])
  } else {
    return(qqh)
  }
}
