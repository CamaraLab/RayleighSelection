#' Benchmarks feature selection algorithms using various datasets.
#'
#' Benchmarks a ranked list of features using a k-means classifier based on the top \code{q} features
#' in the list. The normalized mutual information is used as a metric.
#'
#' @param sorted_results a dataframe with each row corresponding to a different feature. Rows are
#' ranked according to decreasing significance.
#' @param q an integer or a vector of integers specifying the number of features taken into account
#' in the evaluation. If a list of integers is specified, benchmarks are sequentially performed for
#' the different values in the list.
#' @param r an integer specifying the number of times the benchmark is repeated for each value of
#' \code{q}. By default is set to 10.
#' @param num_cores number of cores used for the computation. By default is set to 1.
#' @param dataset a string specifying the dataset used for the benchmark. Supported datasets are
#' \code{gisette} (default), \code{dexter}, and \code{mnist}.
#'
#' @return Returns a matrix with \code{r} rows and \code{q} columns containing the normalized mutual
#' information estimated from the empirical probability distribution. In addition, the average mutual
#' information is ploted as a function of \code{q}.
#' @examples
#' library(RayleighSelection)
#' # Load GISETTE dataset and an example of ranked list produced using \code{rayleigh_selection}
#' data("gisette")
#' data("gisette_results")
#'
#' # Rank pixels according to their p-value
#' sorted_results <- gisette_results[order(gisette_results$p),]
#'
#' # Benchmarks the ranked list of pixels
#' benchmark(sorted_results, seq(500, length(row.names(sorted_results)), 1000))
#'
#' @export
#'
benchmark <- function(sorted_results, q, r = 10, num_cores = 1, dataset = "gisette") {
  if (dataset == "mnist") {
    m <- mnist
  }
  else if (dataset == "gisette") {
    m <- gisette
  }
  else if (dataset == "dexter") {
    m <- dexter
  }
  number <- strsplit(colnames(m), '\\.')
  number <- unlist(number)
  number <- unique(number[seq(1, length(number), 2)])
  clusters <- seq(1, length(number))
  names(clusters) <- number

  m_labels <- vector(mode = 'integer', length = length(m))
  names(m_labels) <- colnames(m)
  m_names <- colnames(m)
  for (idx in seq(m_names))
  {
    name <- strsplit(m_names[[idx]], '\\.')
    m_labels[idx] <- clusters[name[[1]][1]]
  }

  bench <- function(q) {
    top_q_pixels <- sorted_results[1:q, ]
    if (dataset == "dexter") {
      pixel_indices <- as.numeric(row.names(top_q_pixels))
    }
    else {
      pixel_indices <- as.numeric(substring(row.names(top_q_pixels),2))
    }
    pixels <- data.frame(m[pixel_indices, ])
    results <-  data.frame(mi=double())
    for (rep in 1:r)
    {
      km <- kmeans(t(pixels), length(number), iter.max = 100, nstart = 10)
      mi <- mutinformation(km$cluster, m_labels)/max(entropy(km$cluster), entropy(m_labels))
      results[rep, ] <- mi
    }
    return(results)
  }

  results <- mclapply(q, bench, mc.cores = num_cores)

  mi <- lapply(results, "[[", 'mi')
  mi <- matrix(unlist(mi), nrow=r, ncol=length(q), byrow=FALSE)

  mi_mean <- colMeans(mi)
  result <- data.frame(q, mi_mean)

  mi_plot <- ggplot(data=result, aes(x=q, y=mi_mean, group=2)) + geom_line() + geom_point() +
    ggtitle("Mutual Information")
  plot(mi_plot)

  colnames(mi) <- q

  return(mi)
}


