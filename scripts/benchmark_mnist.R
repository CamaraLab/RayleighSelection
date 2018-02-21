library(ggplot2)
library(parallel)
library(clue)
library(infotheo)

data("mnist_results_weighted")

sorted_results <- mnist_results_weighted[order(mnist_results_weighted$p),]

benchmark_mnist <- function(sorted_results, q, r, tSNE = FALSE, perplexity = 10, num_cores = 4) {
  data("mnist")

  number <- strsplit(colnames(mnist), '\\.')
  number <- unlist(number)
  number <- unique(number[seq(1, length(number), 2)])
  clusters <- seq(1, length(number))
  names(clusters) <- number

  mnist_labels <- vector(mode = 'integer', length = length(mnist))
  names(mnist_labels) <- colnames(mnist)
  mnist_names <- colnames(mnist)
  for (idx in seq(mnist_names))
  {
    name <- strsplit(mnist_names[[idx]], '\\.')
    mnist_labels[idx] <- clusters[name[[1]][1]]
  }

  bench <- function(q) {
    top_q_pixels <- sorted_results[1:q, ]
    pixel_indices <- as.numeric(substring(row.names(top_q_pixels),2))
    pixels <- data.frame(mnist[pixel_indices, ])

    set.seed(1)
    if (tSNE) {
      tSNE <- Rtsne::Rtsne(t(pixels), perplexity=perplexity)
      plot(tSNE$Y, col=as.numeric(mnist_labels), pch=16, cex = 0.5)
    }

    results <-  data.frame(ac=double(), mi=double())

    for (rep in 1:r)
    {
      km <- kmeans(t(pixels), length(number), iter.max = 100)

      ac <- accuracy(km$cluster, mnist_labels)
      mi <- natstobits(mutinformation(km$cluster, mnist_labels)/max(entropy(km$cluster), entropy(mnist_labels)))

      results[rep, ] <- c(ac, mi)
    }
    return(results)
  }

  results <- mclapply(q, bench, mc.cores = num_cores)
  ac <- lapply(results, "[[", 'ac')
  ac <- matrix(unlist(ac), nrow=r, ncol=length(q), byrow=FALSE)

  mi <- lapply(results, "[[", 'mi')
  mi <- matrix(unlist(mi), nrow=r, ncol=length(q), byrow=FALSE)

  ac_mean <- colMeans(ac)
  mi_mean <- colMeans(mi)
  result <- data.frame(q, ac_mean, mi_mean)

  ac_plot <- ggplot(data=result, aes(x=q, y=ac_mean, group=1)) + geom_line() + geom_point() +
    ggtitle("Accuracy")
  mi_plot <- ggplot(data=result, aes(x=q, y=mi_mean, group=2)) + geom_line() + geom_point() +
    ggtitle("Mutual Information")

  multiplot(ac_plot, mi_plot, cols=2)

  return(list(ac=ac, mi=mi))
}

# eq 6 from He, Cai, Niyogi paper Laplacian Score for Feature Selection
accuracy <- function(clusterA, clusterB) {
  map <- minWeightBipartiteMatching(clusterA, clusterB)
  nom <- sum(clusterA == map[clusterB])
  return (nom/length(clusterB))
}

# taken from: https://www.r-bloggers.com/matching-clustering-solutions-using-the-hungarian-method/
# labels from cluster A will be matched on the labels from cluster B
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
    idsA <- unique(clusteringA)  # distinct cluster ids in a
    idsB <- unique(clusteringB)  # distinct cluster ids in b
    nA <- length(clusteringA)  # number of instances in a
    nB <- length(clusteringB)  # number of instances in b
    if (length(idsA) != length(idsB) || nA != nB) {
        stop("number of cluster or number of instances do not match")
    }

    nC <- length(idsA)
    tupel <- c(1:nA)

    # computing the distance matrix
    assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
    for (i in 1:nC) {
        tupelClusterI <- tupel[clusteringA == i]
        solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
            nA_I <- length(tupelA_I)  # number of elements in cluster I
            tupelB_I <- tupel[clusterIDsB == i]
            nB_I <- length(tupelB_I)
            nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
            return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
        }, clusteringB, tupelClusterI)
        assignmentMatrix[i, ] <- solRowI
    }

    # optimization
    result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
    attr(result, "assignmentMatrix") <- assignmentMatrix
    return(result)
}

