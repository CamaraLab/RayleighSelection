library(devtools)

load_all("../")

data("lfw")
data("lfw_results")

benchmark_lfw <- function(q, r) {
  sorted_results <- lfw_results[order(-lfw_results$p),]
  top_q_pixels <- sorted_results[1:q, ]
  pixel_indices <- as.numeric(row.names(top_q_pixels))
  pixels <- data.frame(lfw[pixel_indices[1], ])

  for (index in tail(pixel_indices, -1))
  {
    pixels[nrow(pixels) + 1,] <- lfw[index,]
  }

  for (rep in 1:r)
  {
    km <- kmeans(t(pixels), 62)
  }

  return(km)
}
