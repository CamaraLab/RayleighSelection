#' Ranks pairs of features according to the Adjacency Score
#'
#' Uses the Laplacian of a nerve or clique complex to assess the significant colocalization of
#' given pairs of features defined over the points in that complex. Significance is computed by
#' comparing to a null distribution created by reshuffling the features across the points. Pairs
#' of features are permuted in the same way.
#'
#' @param adj_matrix a binary adjacency matrix for the points over which the
#'  features of f are defined
#' @param f a numeric vector or matrix specifying one or more functions with support on
#' the set of points whose significance will be assesed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param f_pairs a 2 column matrix where each row specifies the row names
#'  of two functions or features to be assessed
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#'
#' @return Returns a data frame with the value of the Adjacency Score,
#' the p-values, and the q-values computed using Benjamini-Hochberg procedure.
#'
#' @export
#'
adjacency_score <- function(adj_matrix, f, f_pairs, num_perms = 1000, num_cores = 1) {

  # Check class of f
  if (class(f) != 'matrix') {
    f <- as.matrix(f)
  }

  if(!isSymmetric(adj_matrix)) {
    adj_sym <- 1*((adj_matrix+t(adj_matrix)) > 0)
  } else {
    adj_sym <- adj_matrix
  }

  zero_weights <- rowSums(adj_sym)
  col <- diag(zero_weights) - adj_sym

  permutations <- NULL
  if (num_perms > 0) {
    permutations <- t(mcmapply(function(x) sample(1:ncol(f)), 1:num_perms, mc.cores=num_cores))
  }
  permutations <- rbind(1:ncol(f), permutations)
  # Permute each feature, result is a list of matrices where each matrix corresponds to all the permutations for each feature
  perm_f <- mclapply(1:nrow(f), function(i) as(t(sapply(1:nrow(permutations), function(j) f[i,][permutations[j,]])), class(f)), mc.cores=num_cores)
  names(perm_f) <- row.names(f)

  # Evaluates R and p for a pair of features feature fo
  cornel <- function(fo) {
    kk <- perm_f[[fo[1]]]
    kk <- kk-matrix(rep(kk%*%zero_weights/sum(zero_weights),dim(kk)[2]),dim(kk))
    kk2 <- perm_f[[fo[2]]]
    kk2 <- kk2-matrix(rep(kk2%*%zero_weights/sum(zero_weights),dim(kk2)[2]),dim(kk2))
    qlom <- abs(rowSums((kk%*%diag(zero_weights))*kk2))
    if (sum(qlom)==0.0) {
      qt <- rep(Inf,length(qlom))
    } else {
      qt <- rowSums((kk%*%col)*kk2)/qlom
      qt[is.nan(qt)] <- Inf
    }
    ph <- NULL
    ph$R0 <- qt[1]
    ph$p0 <- (sum(qt<=qt[1])-1.0)/num_perms
    return(ph)
  }

  # Each worker evaluates R anp p for a set fu of pairs of features
  worker <- function(fu) {
    qh <- NULL
    qh$R0 <- NULL
    qh$p0 <- NULL
    for (m in 1:nrow(fu)) {
      d <- cornel(fu[m,])
      qh$R0 <- rbind(qh$R0, d$R0)
      qh$p0 <- rbind(qh$p0, d$p0)
    }
    return(data.frame(qh, row.names = row.names(fu)))
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

  return(qqh)
}
