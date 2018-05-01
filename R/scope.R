# Project documentation and imports

#' The SCOPe package
#'
#' Provides a computational framework for unsupervised identification B cell
#' clones from adaptive immune receptor repertoire sequencing (AIRR-Seq) datasets.
#' This method is based on spectral clustering of the junction sequences of B cell
#' receptors (BCRs, Immunoglobulins) that share the same V gene, J gene and
#' junction length.
#'
#' @section Spectral Clustering for clOne Partitioning (SCOPe):
#'
#' \itemize{
#'   \item  \link{defineClonesScope}:      Clustering sequences into clonal groups.
#'   \item  \link{clonesAnalysis}:         Summary statistics and visualization of the
#'                                         clonal clustering results.
#' }
#'
#' @name        scope
#' @docType     package
#' @references
#' \itemize{
#'   \item  Nouri N and Kleinstein SH. A spectral clustering-based method for identifying
#'   clones from high-throughput B cell repertoire sequencing data. Bioinformatics, (in press).
#'  }
#'
#' @import      alakazam
#' @import      shazam
#' @import      doParallel
#' @import      foreach
#' @import      dplyr
#' @import      ggplot2
#' @import      stringi
#' @import      methods     
#' @importFrom  stats       density kmeans sd uniroot
#' @importFrom  iterators   icount
#' @importFrom  lazyeval    interp
NULL
