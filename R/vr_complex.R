vr_complex <- function(dist, epsilon, weight = TRUE, t = Inf)
{
  idx <- (dist <= epsilon)

  if (!weight)
  {
    weight <- 'NULL'
  }

  adjacency <- idx*dist
  adjacency <- exp(-adjacency**2/t)*idx

  diag(adjacency) <- 0

  g2 <- graph.adjacency(adjacency, mode='undirected', weighted=weight)
  g2$adjacency <- get.adjacency(g2, sparse = TRUE)
  class(g2) <-  c('simplicial', 'igraph')
  return(g2)
}
