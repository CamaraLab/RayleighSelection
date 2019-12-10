# RayleighSelection

```RayleighSelection``` is an R package for feature selection in topological spaces. Features are defined as differential forms on a simplicial complex and their significance is assesed through Rayleigh quotients

## Installation
```
library(devtools)
install_github("CamaraLab/RayleighSelection")
```

*Note for Windows users*: There is currently an error with installing this package on Windows (see [Issue #15](https://github.com/CamaraLab/RayleighSelection/issues/15))

Please use our Docker image [camaralab/rayleigh-selection](https://hub.docker.com/r/camaralab/rayleigh-selection) to run an RStudio server (v3.4 or v3.6) with RayleighSelection already installed:

```docker run -d --rm -p 8787:8787 -v "<dir_path>:/home/rstudio/<dir_name>" -e USER=rstudio -e PASSWORD=<password> camaralab/rayleigh-selection```

After running the above command, RStudio should be available at localhost:8787 in your browser with the local directory at \<dir_path\> mounted in Home.

## Tutorials
[Nerve complex on toy data](https://github.com/CamaraLab/RayleighSelection/blob/master/examples/plot_nerve_example.md)

Given an open cover and a feature on points, compute the Combinatorial Laplacian scores of that feature on the nerve complex of the cover.

[Vietoris-Rips on cyclic scRNA-seq data](https://github.com/CamaraLab/RayleighSelection/blob/master/examples/vr_cycle_example.md)

Given the PCA results of mouse embryonic cells in two differentiation protocols and an ordering on the cells, create a Vietoris-Rips complex. Compute the Combinatorial Laplacian score of gene expression on either just the 0-forms (fast) or both 0-forms and 1-forms (slow).

[Nerve complex on Mapper representation of MNIST](https://github.com/CamaraLab/RayleighSelection/blob/master/examples/mnist_example.md)

Run Mapper on the MNIST dataset to compute an open cover on the handwriting samples, then compute the Combinatorial Laplacian score of the pixel intensity.