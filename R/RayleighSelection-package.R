# Needed to load C++ module
Rcpp::loadModule("mod_laplacian", TRUE)
Rcpp::loadModule("mod_scorer_ensemble", TRUE)


#' @useDynLib RayleighSelection
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @import data.table
#' @import Matrix
#' @import igraph
#' @import parallel
#' @import ggplot2
#' @import ForceAtlas2
#' @import dbscan
#' @import dimRed
#' @import loe
#' @import RSpectra
#' @import RANN
#' @import TDAmapper
#' @import abind

# dummy function to force imports
.import.dummy <- function() {return(NULL)}

