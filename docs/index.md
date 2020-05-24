SCOPer
-------------------------------------------------------------------------------

SCOPer provides a computational framework for the identification of B cell 
clonal relationships from Adaptive Immune Receptor Repertoire sequencing 
(AIRR-Seq) data. It includes methods for assigning clonal identifiers using
sequence identity, hierarchical clustering, and spectral clustering.
SCOPer is part of the [Immcantation](http://immcantation.readthedocs.io) 
analysis framework.

Contact
-------------------------------------------------------------------------------

For help and questions please contact the 
[Immcantation Group](mailto:immcantation@googlegroups.com)


# Dependencies

**Depends:** ggplot2  
**Imports:** alakazam, shazam, data.table, doParallel, dplyr, foreach, magrittr, methods, Rcpp, rlang, scales, stats, stringi  
**Suggests:** knitr, rmarkdown, testthat


# Authors

[Nima Nouri](mailto:nima.nouri@yale.edu) (aut)  
[Edel Aron](mailto:edel.aron@yale.edu) (ctb)  
[Jason Vander Heiden](mailto:jason.vanderheiden@gmail.com) (aut, cre)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


# Citing


To cite the scoper package or spectral clustering-based model in publications,
please use:

Nouri N, Kleinstein S (2018). “A spectral clustering-based method for identifying
clones from high-throughput B cell repertoire sequencing data.” _Bioinformatics_,
i341-i349. doi: 10.1093/bioinformatics/bty235 (URL:
https://doi.org/10.1093/bioinformatics/bty235).

Nouri N, Kleinstein S (2019). “Somatic hypermutation analysis for improved
identification of B cell clonal families from next-generation sequencing data.”
_bioRxiv_. doi: 10.1101/788620 (URL: https://doi.org/10.1101/788620).

To cite the hierarchical clustering-based model in publications, please use:

Gupta N, Adams K, Briggs A, Timberlake S, Vigneault F, Kleinstein S (2017).
“Hierarchical clustering can identify B cell clones with high confidence in Ig
repertoire sequencing data.” _The Journal of Immunology_, 2489-2499. doi:
10.4049/jimmunol.1601850 (URL: https://doi.org/10.4049/jimmunol.1601850).




# License

AGPL-3
