# RayleighSelection

```RayleighSelection``` is an R package for feature selection in topological spaces. Features are defined as differential forms on a simplicial complex and their significance is assesed through Rayleigh quotients

## Installation
```
library(devtools)
install_github("CamaraLab/RayleighSelection")
```

## Tutorials
[Nerve complex on toy data](https://github.com/CamaraLab/RayleighSelection/tree/master/examples/plot_nerve_example.Rmd)
Given an open cover and a feature on that cover, compute the Laplacian scores of that feature on the nerve complex of the cover.

[Vietoris-Rips on cyclic scRNA-seq data](https://github.com/CamaraLab/RayleighSelection/tree/master/examples/vr_cycle_example.Rmd)
Given the PCA results of mouse embryonic cells in two differentiation protocols and an ordering on the cells, create a Vietoris-Rips complex. Compute either just the 0th (fast) or both 0th and 1st (slow) Laplacian scores on the gene expression.

[Nerve complex on Mapper representation of MNIST](https://github.com/CamaraLab/RayleighSelection/tree/master/examples/mnist_example.Rmd)
Run Mapper on the MNIST dataset to compute an open cover on the handwriting samples, then compute the Laplacian score of the pixel intensity.