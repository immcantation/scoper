[![](http://cranlogs.r-pkg.org/badges/grand-total/scoper)](https://www.r-pkg.org/pkg/scoper)
[![](https://cranlogs.r-pkg.org/badges/scoper)](https://www.r-pkg.org/pkg/scoper)
[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

SCOPer
-------------------------------------------------------------------------------

[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)

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
**Imports:** alakazam, shazam, data.table, doParallel, dplyr, foreach, methods, Rcpp, rlang, scales, stats, stringi, tidyr  
**Suggests:** knitr, rmarkdown, testthat


# Authors

[Nima Nouri](mailto:nima.nouri@yale.edu) (aut)  
[Edel Aron](mailto:edel.aron@yale.edu) (ctb)  
[Jason Vander Heiden](mailto:jason.vanderheiden@gmail.com) (aut, cre)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)


# Citing


To cite the scoper package or spectral clustering-based model in
publications, please use:

Nouri N, Kleinstein S (2018). “A spectral clustering-based method for
identifying clones from high-throughput B cell repertoire sequencing
data.” _Bioinformatics_, i341-i349. doi: 10.1093/bioinformatics/bty235
(URL: https://doi.org/10.1093/bioinformatics/bty235).

Nouri N, Kleinstein S (2020). “Somatic hypermutation analysis for
improved identification of B cell clonal families from next-generation
sequencing data.” _PLOS Computational Biology_, *16*(6), e1007977. doi:
10.1371/journal.pcbi.1007977 (URL:
https://doi.org/10.1371/journal.pcbi.1007977).

To cite the hierarchical clustering-based model in publications, please
use:

Gupta N, Adams K, Briggs A, Timberlake S, Vigneault F, Kleinstein S
(2017). “Hierarchical clustering can identify B cell clones with high
confidence in Ig repertoire sequencing data.” _The Journal of
Immunology_, 2489-2499. doi: 10.4049/jimmunol.1601850 (URL:
https://doi.org/10.4049/jimmunol.1601850).




# License

AGPL-3
