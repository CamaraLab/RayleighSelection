library(devtools)
load_all("/home/seshu/dev/RayleighSelection", reset = TRUE, recompile = TRUE)

open_cover = list(c(1, 2), c(2, 3), c(3, 4), c(2, 4, 5))
complex <- nerve_complex(open_cover)
plot_skeleton(complex)
ones <- matrix(1, nrow=5, ncol=1)
pushCpp(ones, complex$points_in_vertex, 0, complex$adjacency)

pushCpp(c(0,1,0,0,0), complex$points_in_vertex, 10, complex$adjacency)
