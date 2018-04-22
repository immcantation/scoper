# Project documentation and imports

#' SCOPe
#'
#' Provides a computational framework for unsupervised identification B cell
#' clones from adaptive immune receptor repertoire sequencing (AIRR-Seq) datasets. 
#' This method is based on spectral clustering of the junction sequences of B cell 
#' receptors (BCRs, Immunoglobulins) that share the same V gene, J gene and 
#' junction length.
#'
#' @name        scope
#' @docType     package
#' @import      alakazam
#' @import      shazam
#' @import      doParallel
#' @import      foreach
#' @import      dplyr
#' @import      ggplot2
#' @import      stringi
#' @importFrom  methods     new
#' @importFrom  stats       density kmeans sd uniroot
#' @importFrom  iterators   icount
#' @importFrom  lazyeval    interp
NULL
