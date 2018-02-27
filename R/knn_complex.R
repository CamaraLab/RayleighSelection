library(dbscan)

knn_complex <- function(dist, k, t = Inf)
{
  kneighbors <- kNN(dist, k)
  id <- kneighbors$id
  nodes <- dim(id)[1]
  rows <- 1:nodes
  adjacency <- matrix(0, nrow=nodes, ncol=nodes)

  for (col in 1:k)
  {
    adjacency[cbind(rows, id[,col])] <- 1
    adjacency[cbind(id[,col], rows)] <- 1
  }

  adjacency <- exp(-dist**2/t)*adjacency

  return (graph_to_complex(adjacency))
}
