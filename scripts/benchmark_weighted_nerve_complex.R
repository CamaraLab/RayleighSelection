library(TDAmapper)
library(dimRed)
library(microbenchmark)

data("lfw")

leim <- LaplacianEigenmaps()
lfw_top <- lfw[apply(lfw, 1, var) > 0.9,]
emb <- leim@fun(as(t(lfw_top), "dimRedData"), leim@stdpars)

lfw_distances <- (1.0 - cor(lfw_top))

m2 <- mapper2D(distance_matrix = lfw_distances,
               filter_values = list(emb@data@data[,1], emb@data@data[,2]),
               num_intervals = c(40,40),
               percent_overlap = 30,
               num_bins_when_clustering = 10);

# Compute the nerve complex

microbenchmark(nerve_complex(m2$points_in_vertex))
microbenchmark(nerve_complex(m2$points_in_vertex, weight = FALSE))
