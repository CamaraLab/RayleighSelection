library(clue)
library(infotheo)

data("lfw")
data("lfw_results")

sorted_results <- lfw_results_weighted[order(lfw_results_weighted$p),]

# assign images in lfw to clusters each image from the same person is in one cluster:
# ie:         Alejandro_Toledo.1           Alejandro_Toledo.2
#                          1                            1       ...
#                 Alvaro_Uribe.1               Alvaro_Uribe.2
#                          2                            2
people <- strsplit(colnames(lfw), '\\.')
people <- unlist(people)
people <- unique(people[seq(1, length(people), 2)])
clusters <- seq(1, length(people))
names(clusters) <- people

lfw_labels <- vector(mode = 'integer', length = length(lfw))
names(lfw_labels) <- colnames(lfw)
lfw_names <- colnames(lfw)
for (idx in seq(lfw_names))
{
  name <- strsplit(lfw_names[[idx]], '\\.')
  lfw_labels[idx] <- clusters[name[[1]][1]]
}


benchmark_lfw <- function(q, r, perplexity = 10) {
  top_q_pixels <- sorted_results[1:q, ]
  pixel_indices <- as.numeric(row.names(top_q_pixels))
  pixels <- data.frame(lfw[pixel_indices, ])

  set.seed(1)
  tSNE <- Rtsne::Rtsne(t(pixels), perplexity=perplexity)
  plot(tSNE$Y, col=as.numeric(lfw_labels), pch=16, cex = 0.5)

  results <-  data.frame(ac=double(), mi=double())

  for (rep in 1:r)
  {
    km <- kmeans(t(pixels), length(people), iter.max = 100)

    ac <- accuracy(km$cluster, lfw_labels)
    mi <- natstobits(mutinformation(km$cluster, lfw_labels)/max(entropy(km$cluster), entropy(lfw_labels)))

    results[rep, ] <- c(ac, mi)
  }

  return (results)
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
    require(clue)
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

