library(devtools)
load_all("/home/seshu/dev/RayleighSelection", reset = TRUE, recompile = TRUE)

open_cover = list(c(1, 2), c(2, 3), c(3, 4), c(2, 4, 5))
n = 5
order = 1:n
complex <- nerve_complex(open_cover)

ones <- matrix(1, nrow=n, ncol=1)
pushCpp(ones, complex$points_in_vertex, 0, complex$adjacency)

identity <- diag(5)
pushCpp(identity, complex$points_in_vertex, 0, complex$adjacency)
