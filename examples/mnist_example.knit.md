---
title: "MNIST Example"
output: rmarkdown::github_document
---

#### Calculate the Combinatorial Laplacian score of pixels on a nerve complex created by TDA mapper run on the MNIST dataset


```r
library(RayleighSelection)
```

## Load the mnist dataset

```r
data("mnist")
```

## Compute reduced representation
#### Using Laplacian eigenmap of pixels with high variance

```r
library(dimRed)
```

```
## Loading required package: DRR
```

```
## Loading required package: kernlab
```

```
## Loading required package: CVST
```

```
## Loading required package: Matrix
```

```
## 
## Attaching package: 'dimRed'
```

```
## The following object is masked from 'package:stats':
## 
##     embed
```

```
## The following object is masked from 'package:base':
## 
##     as.data.frame
```

```r
leim <- LaplacianEigenmaps()
mnist_top <- mnist[apply(mnist, 1, var) > 10000,]
emb <- leim@fun(as(t(mnist_top), "dimRedData"), leim@stdpars)
```

```
## 2018-10-23 16:16:58: Creating weight matrix
```

```
## 2018-10-23 16:17:23: Eigenvalue decomposition
```

```
## Eigenvalues: 6.178736e-02 4.905515e-02 1.200188e-15
```

```
## 2018-10-23 16:18:04: DONE
```

## Compute Mapper representation
#### Using the Laplacian eigenmap as an auxiliary function and correlation distance as metric

```r
library(TDAmapper)
mnist_distances <- (1.0 - cor(mnist_top))
m2 <- mapper2D(distance_matrix = mnist_distances,
                filter_values = list(emb@data@data[,1], emb@data@data[,2]),
                num_intervals = c(30,30),
                percent_overlap = 35,
                num_bins_when_clustering = 10);
```

## Compute the nerve complex

```r
gg <- nerve_complex(m2$points_in_vertex)
```

## Compute 0-form and 1-form Comb. Lap. scores, p-value, and q-value
#### For the 301st through 305th pixels

```r
rayleigh_selection(gg, mnist[301:305,])
```

```
##             R0    p0      q0       R1    p1     q1
## X301 0.1774026 0.000 0.00000 1.526852 0.279 0.6975
## X302 0.1855035 0.000 0.00000 1.498672 0.183 0.6975
## X303 0.2114068 0.000 0.00000 1.560080 0.513 0.8550
## X304 0.2536190 0.003 0.00375 1.631599 0.720 0.9000
## X305 0.2999599 0.078 0.07800 1.903140 0.955 0.9550
```
