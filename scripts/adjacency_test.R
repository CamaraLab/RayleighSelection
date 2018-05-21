library('RayleighSelection')

## open_cover = list(c(1, 2, 3), c(2, 3), c(3, 4), c(5, 6))
## open_cover = list(c(5, 6, 3), c(6, 3), c(3, 4), c(1, 2))
open_cover = list(c(1, 2), c(2, 3), c(3, 4), c(2, 4, 5))
n = 5
samples = 1:n
feat = rep(1, n)
features = data.frame(sample = samples, feature = feat)
complex <- adjacencyCpp(open_cover, features, TRUE)
## complex

zero_simplices = as.data.frame(matrix(1:length(open_cover), length(open_cover), 1))
one_simplices_idx <- data.frame(which(as.matrix(complex$one_simplices != 0), arr.ind=T))
one_simplices_idx <- one_simplices_idx[with(one_simplices_idx, order(row)), ]
one_simplices <- data.frame(t(apply(one_simplices_idx, 1, function(x) complex$order[unlist(x)])))
rownames(one_simplices) <- NULL

idxs <- which(lapply(complex$two_simplices, any) == TRUE, arr.ind=T)

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

for(idx in idxs)
{
  edges <- which(as.matrix(complex$two_simplices[[idx]] != 0), arr.ind=T)
  two_simplices <- t(apply(edges, 1, function(edge) c(idx, edge)))
  bound <- apply(two_simplices, 1, function(two_simplex) {boundary(two_simplex, one_simplices)})

  lapply(bound, function(two_simplex_boundary) {
    two_simplex_boundary <- two_simplex_boundary[two_simplex_boundary$sign != 0, ]

    for(i in rownames(two_simplex_boundary))
    {
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
          l1_up[idxi, idxi] <<- l1_up[idxi, idxi] + 1
        }
        else
        {
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

      f <- boundary(one_simplices[i,], zero_simplices)
      fp <- boundary(one_simplices[j, ], zero_simplices)

      row <- which(f[,1,drop=F] == zero_simplex[1,], arr.ind=T)[1]
      sgn_e_f <- f[row,]$sign
      row <- which(fp[,1,drop=F] == zero_simplex[1,], arr.ind=T)[1]
      sgn_e_fp <- fp[row, ]$sign

      l1_down[i, j] <- sgn_e_f*sgn_e_fp
      l1_down[j, i] <- sgn_e_f*sgn_e_fp
    }
  }
}

l1_up
l1_down
l1_up + l1_down
