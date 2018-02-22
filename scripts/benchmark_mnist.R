library(ggplot2)
library(parallel)
library(clue)
library(infotheo)

benchmark_mnist <- function(sorted_results, q, r, num_cores = 4) {
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
    results <-  data.frame(mi=double())
    for (rep in 1:r)
    {
      km <- kmeans(t(pixels), length(number), iter.max = 100, nstart = 10)
      mi <- mutinformation(km$cluster, mnist_labels)/max(entropy(km$cluster), entropy(mnist_labels))
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


data("mnist_results_weighted")
sorted_results <- mnist_results_weighted[order(-mnist_results_weighted$var),]
bench_var_euclid <- benchmark_mnist(sorted_results, seq(100, length(row.names(sorted_results)), 20), 20, 10)
save(bench_var_euclid, file = 'data/bench_var_euclid.RData', compress = TRUE)

sorted_results <- mnist_results_weighted[order(mnist_results_weighted$p),]
bench_rayleigh_euclid <- benchmark_mnist(sorted_results, seq(100, length(row.names(sorted_results)), 20), 20, 10)
save(bench_rayleigh_euclid, file = 'data/bench_rayleigh_euclid.RData', compress = TRUE)

data("mnist_results")
sorted_results <- mnist_results[order(mnist_results$p),]
bench_rayleigh_euclid_unweighted <- benchmark_mnist(sorted_results, seq(100, length(row.names(sorted_results)), 20), 20, 10)
save(bench_rayleigh_euclid_unweighted, file = 'data/bench_rayleigh_euclid_unweighted.RData', compress = TRUE)

data("mnist_results_shifted")
sorted_results <- mnist_results_shifted[order(mnist_results_shifted$p),]
bench_rayleigh_euclid_shifted <- benchmark_mnist(sorted_results, seq(100, length(row.names(sorted_results)), 20), 20, 10)
save(bench_rayleigh_euclid_shifted, file = 'data/bench_rayleigh_euclid_shifted.RData', compress = TRUE)
