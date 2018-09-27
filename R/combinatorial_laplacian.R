#' Computes the Combinatorial Laplacian for 0- and 1-forms.
#'
#' Given a nerve or a clique complex, it computes the Combinatorial Laplacians for 0- and 1-forms and
#' the weights of 0- and 1-simplices.
#'
#' @param g2 an object of the class \code{simplicial} containing the nerve or clique complex.
#' @param one_forms when set FALSE only the Combinatorial Laplacian for 0-forms is
#' computed. By default is set to TRUE.
#' @param weights when set to TRUE it takes 2-simplices into account when computing weights.
#' By default is set to FALSE.
#'
#' @return Returns a list with the Combinatorial Laplacians for 0- and 1-forms
#' and the weights for 0- and 1-simplices.
#' @examples
#' library(RayleighSelection)
#' gy <- nerve_complex(list(c(1,4,6,10), c(1,2,7), c(2,3,8), c(3,4,9,10), c(4,5)))
#' combinatorial_laplacian(gy)
#'
#' @export
#'
combinatorial_laplacian <- function(g2, one_forms = TRUE, weights = FALSE) {
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
  if (weights || one_forms) {
    idxs <- which(lapply(g2$two_simplices, any) == TRUE, arr.ind=T)
  }

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
  } else {
    one_weights <- rep(1,nrow(one_simplices))
  }

  # compute weights of 0-simplices
  zero_weights <- rep(0,siz)
  for(i in 1:nrow(one_simplices))
  {
    zero_weights[one_simplices[i,1]] <- zero_weights[one_simplices[i,1]] + one_weights[i]
    zero_weights[one_simplices[i,2]] <- zero_weights[one_simplices[i,2]] + one_weights[i]
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

    l1_down <- l1down(as.matrix(one_simplices), zero_weights, one_weights)
  }

  out <- list(l0 = col, zero_weights = zero_weights)
  if (one_forms) {
    out[['l1up']] <- l1_up
    out[['l1down']] <- l1_down
  }
  return(out)
}
