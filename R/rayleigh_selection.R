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
#' @param shift real number specifying a shift that is added to \code{f}. By default is set to 0.
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param seed integer specifying the seed used to initialize the generator of permutations.
#' By default is set to 10.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#' @param adjacency when set to \code{TRUE} the adjacency matrix is used instead of the Laplacian,
#' as in Rizvi, Camara, et al. Nat. Biotechnol. 35 (2017). By default is set to \code{FALSE}.
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
rayleigh_selection <- function(g2, f, shift = 0.0, num_perms = 1000, seed = 10, num_cores = 1,
                               adjacency = FALSE, L1 = FALSE) {
  # Check class of f
  if (class(f) != 'data.frame') {
    f <- as.data.frame(f)
  }

  siz <- sqrt(length(g2$adjacency))

  adj_sym <- g2$adjacency+t(g2$adjacency)
  if (adjacency) {
    dd <- rep(1, siz)
    col <- -(adj_sym)
  } else {
    # Scalar Laplace operator
    dd <- rowSums(adj_sym)
    col <- diag(dd)-adj_sym
  }


  # Evaluates R and p for a feature fo
  cornel <- function(fo) {
    kk<-pushCpp(as.numeric(fo), g2$points_in_vertex, num_perms)
    kk <- kk-matrix(rep(kk%*%dd/sum(dd),dim(kk)[2]),dim(kk))
    qlom <- rowSums(dd*kk^2)
    if (sum(abs(qlom))==0.0) {
      qt <- rep(Inf,length(qlom))
    } else {
      qt <- rowSums((kk%*%col)*kk)/qlom
      qt[is.nan(qt)] <- Inf
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

  # If L1 compute also L_1 Laplacian
  if (L1) {
    zero_simplices = as.data.frame(matrix(1:siz, siz, 1))

    ## get the non-zero indices of the entries in the adjacency matrix, sort them by the row index, and
    ## convert the row index into the vertex using the order of the vertices, since the adjacency matrix returned
    ## by adjacencyCpp will have rows and columns labeled by the order
    one_simplices_idx <- data.frame(which(as.matrix(g2$one_simplices != 0), arr.ind=T))
    one_simplices_idx <- one_simplices_idx[with(one_simplices_idx, order(row)), ]
    one_simplices <- data.frame(t(apply(one_simplices_idx, 1, function(x) g2$order[unlist(x)])))
    rownames(one_simplices) <- NULL

    ## the boundary function takes a n-simplex and a dataframe containing all (n-1)-simplices in the complex
    ## it then looks for the row in the faces dataframe which equals the n-simplex without the i'th vertex and
    ## stores the sign in the $sign column of the faces
    boundary <- function(simplex, faces) {
      n <- length(simplex)
      faces["sign"] <- 0
      for(i in 1:n)
      {
        row_idx <- which(apply(faces[,1:n-1,drop=F], 1, function(r) all(r == simplex[-i]) == TRUE), arr.ind=T)
        faces[row_idx, "sign"] <- (-1)**(i+1)
      }
      faces <- faces[faces$sign != 0, ]
      return(faces)
    }

    n <- dim(one_simplices)[1]
    l1_up <- matrix(0, n, n)
    l1_down <- matrix(0, n, n)

    ## since complex$two_simplices is a list of sparse matrices, we first get the i'th indices of all the
    ## 2-simplices <i, j, k> by checking which matrices in the list are non-zero
    idxs <- which(lapply(g2$two_simplices, any) == TRUE, arr.ind=T)

    for(idx in idxs)
    {
      ## given the i'th vertex (idx) of the two-simplices, find all the j and k vertices by looking for non-zero
      ## entries of the i'th matrix in the list of matrices, create the vecter of two_simplices <i, j, k>,
      ## and compute their boundaries
      edges <- which(as.matrix(g2$two_simplices[[idx]] != 0), arr.ind=T)
      two_simplices <- t(apply(edges, 1, function(edge) c(idx, edge)))
      bound <- apply(two_simplices, 1, function(two_simplex) {boundary(two_simplex, one_simplices)})

      lapply(bound, function(two_simplex_boundary) {
        ## for each two simplex boundary computed above, compute l1_up
        for(i in rownames(two_simplex_boundary))
        {
          ## row_sign = sgn(F, d(Fbar))
          row_sign <- two_simplex_boundary[i, "sign"]
          idxi <- as.numeric(i)
          for(j in rownames(two_simplex_boundary))
          {
            if(i > j)
            {
              next
            }
            idxj <- as.numeric(j)
            if(i == j)
            {
              ## need to use <<- since this is a closure
              l1_up[idxi, idxi] <<- l1_up[idxi, idxi] + 1
            }
            else
            {
              ## col_sign = sgn(F', d(Fbar))
              col_sign <- two_simplex_boundary[j, "sign"]
              l1_up[idxi, idxj] <<- row_sign*col_sign
              l1_up[idxj, idxi] <<- row_sign*col_sign
            }
          }
        }
      })
    }

    for(i in 1:nrow(one_simplices))
    {
      for(j in 1:nrow(one_simplices))
      {
        if(i > j)
        {
          next
        }
        else if(i == j)
        {
          l1_down[i, j] <- l1_down[i, j] + length(one_simplices[i,])
        }
        else
        {
          zero_simplex <- intersect(one_simplices[i,], one_simplices[j, ])
          if(length(zero_simplex) == 0)
          {
            next
          }

          ## if the i'th and j'th 1-simplices share a 0-simplex (zero_simplex) compute the boundary of the two
          ## 1-simplices
          f <- boundary(one_simplices[i,], zero_simplices)
          fp <- boundary(one_simplices[j, ], zero_simplices)

          ## get the sign of the zero_simplex in the two boundaries,
          ## l1_down = sgn(E, d(F))*sgn(E, d(F'))
          row <- which(f[,1,drop=F] == zero_simplex[1,], arr.ind=T)[1]
          sgn_e_f <- f[row,]$sign
          row <- which(fp[,1,drop=F] == zero_simplex[1,], arr.ind=T)[1]
          sgn_e_fp <- fp[row, ]$sign

          l1_down[i, j] <- sgn_e_f*sgn_e_fp
          l1_down[j, i] <- sgn_e_f*sgn_e_fp
        }
      }
    }

    print(l1_up + l1_down)

  }

  return(qqh)
}
