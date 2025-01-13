# Download and Installation

Download
-------------------------------------------------------------------------------

The latest stable release of SCOPer can be downloaded from 
[CRAN](http://cran.rstudio.com/web/packages/scoper) or 
[GitHub](https://github.com/immcantation/scoper/tags).

Installing Released Versions
-------------------------------------------------------------------------------

The simplest way to install SCOPer is via CRAN:

```R
install.packages("scoper")
```

Downloaded source builds from GitHub may be installed in the usual way:

```R
install.packages("scoper_x.y.z.tar.gz", repos=NULL, type="source")
```

Building Development Versions
-------------------------------------------------------------------------------

To build from the [source code](http://github.com/immcantation/scoper),
first install the build dependencies:

```R
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "Rcpp"))
```

To install the latest development code via devtools:

```R
library(devtools)
install_github("immcantation/scoper@master")
```

Note, using `install_github` will not build the documentation. To generate the 
documentation, clone the repository and build as normal using devtools, 
roxygen and knitr:

```R
library(devtools)
install_deps()
document()
build()
install()
```
