# RayleighSelection

```RayleighSelection``` is an R package for feature selection in topological spaces. Features are defined as differential forms on a simplicial complex and their significance is assesed through Rayleigh quotients

## Installation
```
library(devtools)
install_github("CamaraLab/RayleighSelection")
```

## Tutorials
[Nerve complex on toy data](https://github.com/CamaraLab/RayleighSelection/blob/master/examples/plot_nerve_example.md)

Given an open cover and a feature on points, compute the Combinatorial Laplacian scores of that feature on the nerve complex of the cover.

[Vietoris-Rips on cyclic scRNA-seq data](https://github.com/CamaraLab/RayleighSelection/blob/master/examples/vr_cycle_example.md)

Given the PCA results of mouse embryonic cells in two differentiation protocols and an ordering on the cells, create a Vietoris-Rips complex. Compute the Combinatorial Laplacian score of gene expression on either just the 0-forms (fast) or both 0-forms and 1-forms (slow).

[Nerve complex on Mapper representation of MNIST](https://github.com/CamaraLab/RayleighSelection/blob/master/examples/mnist_example.md)

Run Mapper on the MNIST dataset to compute an open cover on the handwriting samples, then compute the Combinatorial Laplacian score of the pixel intensity.