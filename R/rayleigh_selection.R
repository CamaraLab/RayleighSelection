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
#' the set of points whose significance will be assesed in the simplicial complex. Each
#' column corresponds to a point and each row specifies a different function.
#' @param num_perms number of permutations used to build the null distribution for each
#' feature. By default is set to 1000.
#' @param seed integer specifying the seed used to initialize the generator of permutations.
#' By default is set to 10.
#' @param num_cores integer specifying the number of cores to be used in the computation. By
#' default only one core is used.
#' @param one_forms when set FALSE only the Combinatorial Laplacian Score for 0-forms is
#' computed. By default is set to TRUE.
#'
#' @return Returns a data frame with the value of the Combinatorial Laplacian Score for 0- and 1-forms,
#' the p-values, and the q-values computed using Benjamini-Hochberg procedure.
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
rayleigh_selection <- function(g2, f, num_perms = 1000, seed = 10, num_cores = 1, one_forms = TRUE,
                               weights = TRUE) {
  # Check class of f
  if (class(f) != 'data.frame') {
    f <- as.data.frame(f)
  }

  siz <- sqrt(length(g2$adjacency))

  zero_simplices <- as.data.frame(matrix(1:siz, siz, 1))

  ## get the non-zero indices of the entries in the adjacency matrix, sort them by the row index, and
  ## convert the row index into the vertex using the order of the vertices, since the adjacency matrix returned
  ## by adjacencyCpp will have rows and columns labeled by the order
  one_simplices_idx <- data.frame(which(as.matrix(g2$one_simplices != 0), arr.ind=T), row.names = NULL)
  one_simplices_idx <- one_simplices_idx[with(one_simplices_idx, order(row)), ]
  one_simplices <- data.frame(t(apply(one_simplices_idx, 1, function(x) g2$order[unlist(x)])))
  rownames(one_simplices) <- NULL
  adjo <- matrix(rep(0,siz**2),c(siz,siz))
  for (ik in 1:nrow(one_simplices)) {
    adjo[one_simplices[ik,1],one_simplices[ik,2]] <- ik
  }
  diji <- as.numeric(t(adjo))
  diji <- diji[diji != 0]

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
      row_idx <- which(apply(faces[,1:n-1,drop=F], 1, function(r) all(r == rev(simplex[-i])) == TRUE), arr.ind=T)
      faces[row_idx, "sign"] <- (-1)**i
    }
    faces <- faces[faces$sign != 0, ]
    return(faces)
  }

  n <- dim(one_simplices)[1]

  ## since complex$two_simplices is a list of sparse matrices, we first get the i'th indices of all the
  ## 2-simplices <i, j, k> by checking which matrices in the list are non-zero
  idxs <- which(lapply(g2$two_simplices, any) == TRUE, arr.ind=T)

  # compute weights for 1-simplices
  if (weights) {
    one_weights <- rep(0,nrow(one_simplices))
    for(idx in idxs)
    {
      edges <- which(as.matrix(g2$two_simplices[[idx]] != 0), arr.ind=T)
      two_simplices <- t(apply(edges, 1, function(edge) c(idx, edge)))
      bound <- apply(two_simplices, 1, function(two_simplex) {boundary(two_simplex, one_simplices)})
      for (mkj in bound) {
        one_weights[as.numeric(rownames(mkj))] <- (one_weights[as.numeric(rownames(mkj))] + 1)
      }
    }
    one_weights[one_weights==0] <- 1

    # compute weights of 0-simplices
    zero_weights <- rep(0,siz)
    for(i in 1:nrow(one_simplices))
    {
      zero_weights[one_simplices[i,1]] <- zero_weights[one_simplices[i,1]] + one_weights[i]
      zero_weights[one_simplices[i,2]] <- zero_weights[one_simplices[i,2]] + one_weights[i]
    }
  } else {
    one_weights <- rep(1,nrow(one_simplices))
    zero_weights <- rep(1,siz)
  }

  # Compute L_0 Laplacian
  adj_sym <- g2$adjacency+t(g2$adjacency)
  col <- diag(zero_weights)-adj_sym

  if (one_forms) {
    l1_up <- matrix(0, n, n)
    l1_down <- matrix(0, n, n)

    for(idx in idxs)
    {
      ## given the i'th vertex (idx) of the two-simplices, find all the j and k vertices by looking for non-zero
      ## entries of the i'th matrix in the list of matrices, create the vector of two_simplices <i, j, k>,
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
              l1_up[idxi, idxi] <<- l1_up[idxi, idxi] + 1/one_weights[idxi]
            }
            else
            {
              ## col_sign = sgn(F', d(Fbar))
              col_sign <- two_simplex_boundary[j, "sign"]
              l1_up[idxi, idxj] <<- row_sign*col_sign/one_weights[idxi]
              l1_up[idxj, idxi] <<- row_sign*col_sign/one_weights[idxj]
            }
          }
        }
      })
    }

    # Very slow double loop. It should be recoded in C++
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
          l1_down[i, j] <- l1_down[i, j] + sum(one_weights[i]/zero_weights[as.numeric(one_simplices[i,])])
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
          ff <- boundary(one_simplices[i,], zero_simplices)
          ffp <- boundary(one_simplices[j, ], zero_simplices)

          ## get the sign of the zero_simplex in the two boundaries,
          ## l1_down = sgn(E, d(F))*sgn(E, d(F'))
          row <- which(ff[,1,drop=F] == zero_simplex[1,], arr.ind=T)[1]
          sgn_e_f <- ff[row,]$sign
          row <- which(ffp[,1,drop=F] == zero_simplex[1,], arr.ind=T)[1]
          sgn_e_fp <- ffp[row, ]$sign

          l1_down[i, j] <- sgn_e_f*sgn_e_fp*one_weights[j]/zero_weights[as.numeric(zero_simplex)]
          l1_down[j, i] <- sgn_e_f*sgn_e_fp*one_weights[i]/zero_weights[as.numeric(zero_simplex)]
        }
      }
    }
  }

  # Evaluates R and p for a feature fo
  cornel <- function(fo) {
    kmn<-pushCpp(as.numeric(fo), g2$points_in_vertex, num_perms, g2$adjacency)
    kk <- kmn$vertices
    kk <- kk-matrix(rep(kk%*%zero_weights/sum(zero_weights),dim(kk)[2]),dim(kk))
    qlom <- rowSums(zero_weights*kk^2)
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
      qlomv <- rowSums(one_weights*kkv^2)
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
    for (m in row.names(fu)) {
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
