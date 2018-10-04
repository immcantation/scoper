SCOPer
-------------------------------------------------------------------------------

Provides a computational framework for unsupervised identification B cell
clones from adaptive immune receptor repertoire sequencing (AIRR-Seq) datasets. 
This method is based on spectral clustering of the junction sequences of B cell 
receptors (BCRs, Immunoglobulins) that share the same V gene, J gene and 
junction length.

Build Instructions
-------------------------------------------------------------------------------

To build from the [source code](http://bitbucket.org/kleinstein/scoper),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2"))
```

To install the latest development code via devtools:

```R
library(devtools)
install_bitbucket("kleinstein/scoper@default")
```

Note, using `install_bitbucket` will not build the documentation. To generate the 
documentation, clone the repository and build as normal. Then run the following 
R commands from the package root:

```R
library(devtools)
install_deps(dependencies=T)
document()
install()
```

Contact
-------------------------------------------------------------------------------

For help and questions please contact the [Immcantation Group](mailto:immcantation@googlegroups.com)


# Dependencies

**Depends:** ggplot2  
**Imports:** alakazam, shazam, doParallel, foreach, dplyr, stringi, methods, stats, iterators, lazyeval  
**Suggests:** knitr, rmarkdown, testthat


# Authors

[Nima Nouri](mailto:nima.nouri@yale.edu) (aut, cre)  
[Jason Vander Heiden](mailto:jason.vanderheiden@yale.edu) (ctb)  
[Steven Kleinstein](mailto:steven.kleinstein@yale.edu) (aut, cph)
