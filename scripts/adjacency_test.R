library(RayleighSelection)

## open_cover = list(c(1, 2, 3), c(2, 3), c(3, 4), c(5, 6))
## open_cover = list(c(5, 6, 3), c(6, 3), c(3, 4), c(1, 2))
open_cover = list(c(1, 2), c(2, 3), c(3, 4), c(2, 4, 5))
n = 5
samples = 1:n
feat = rep(1, n)
features = data.frame(sample = samples, feature = feat)
complex <- adjacencyCpp(open_cover, features, TRUE)
## complex

one_simplices_idx<- data.frame(which(as.matrix(complex$one_simplices != 0), arr.ind=T))
one_simplices_idx[with(one_simplices_idx, order(row)), ]
one_simplices <- data.frame(t(apply(one_simplices_idx, 1, function(x) complex$order[unlist(x)])))

idxs <- which(lapply(complex$two_simplices, any) == TRUE, arr.ind=T)

boundary <- function(simplex, faces) {
  n <- length(simplex)

  faces["sign"] <- 0

  for(i in 1:n)
  {
    simp <- rep(two_simplex)
    simp[i] <- NA
    simp <- simp[!is.na(simp)]
    row_idx <- which(apply(faces[,1:n-1], 1, function(r) all(r == simp) == TRUE), arr.ind=T)
    faces[row_idx, "sign"] <- (-1)**(i+1)
  }

  return(faces)
}

n <- dim(one_simplices)[1]
l1_up <- matrix(0, n, n)

for(idx in idxs)
{
  edges <- which(as.matrix(complex$two_simplices[[idx]] != 0), arr.ind=T)
  two_simplex <- c(idx, edges)
  bound <- boundary(two_simplex, one_simplices)
  print(bound)
  for(i in 1:n)
  {
    row_sign <- bound$sign[i]

    if(row_sign == 0)
      next

    for(j in i:n)
    {
      if(i == j)
      {
        l1_up[i, i] = l1_up[i, i] + 1
      }
      else
      {
        col_sign <- bound$sign[j]
        l1_up[i, j] = row_sign*col_sign
      }
    }
  }
}

l1_up

