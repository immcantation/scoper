SCOPer
-------------------------------------------------------------------------------

SCOPer provides a computational framework for B cell clones identification
from adaptive immune receptor repertoire sequencing (AIRR-Seq) datasets. 
Three models are included (identical, hierarchical, and spectral) 
that perform clustering among sequences of BCRs/IGs (B cell receptors/immunoglobulins) 
which share the same V gene, J gene and junction length. SCOPer is part of the 
[Immcantation](http://immcantation.readthedocs.io) analysis framework.

Contact
-------------------------------------------------------------------------------

For help and questions please contact the [Immcantation Group](mailto:immcantation@googlegroups.com)


# Dependencies

**Depends:** ggplot2  
**Imports:** alakazam, shazam, doParallel, foreach, dplyr, Rcpp, seqinr, data.table, stringi, stringr, methods, stats, rlang  
**Suggests:** knitr, rmarkdown, testthat


# Authors

[Nima Nouri](mailto:nima.nouri@yale.edu) (aut, cre)  
[Jason Vander Heiden](mailto:jason.vanderheiden@yale.edu) (ctb)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


# Citing


To cite the scoper package or spectral clustering-based model in publications, please use:

Nouri N, Kleinstein S (2018). “A spectral clustering-based method for identifying clones from high-throughput B
cell repertoire sequencing data.” _Bioinformatics_, i341-i349. doi: 10.1093/bioinformatics/bty235 (URL:
https://doi.org/10.1093/bioinformatics/bty235).

Nouri N, Kleinstein S (2019). “Somatic hypermutation analysis for improved identification of B cell clonal
families from next-generation sequencing data.” _bioRxiv_. doi: 10.1101/788620 (URL:
https://doi.org/10.1101/788620).

To cite the hierarchical clustering-based model in publications, please use:

Gupta N, Adams K, Briggs A, Timberlake S, Vigneault F, Kleinstein S (2017). “Hierarchical clustering can
identify B cell clones with high confidence in Ig repertoire sequencing data.” _The Journal of Immunology_,
2489-2499. doi: 10.4049/jimmunol.1601850 (URL: https://doi.org/10.4049/jimmunol.1601850).


