WRScpp
======

This package provides C++ sub-routines for several iterative procedures in the [R package `WRS`](https://r-forge.r-project.org/projects/wrs/) for robust statistics by [Rand Wilcox](http://dornsife.usc.edu/cf/labs/wilcox/wilcox-faculty-display.cfm). These sub-routines provide substantial performance boosts.

##Installation

To install this package, three R packages are required:

* `WRS`:
	install.packages("")

* `RcppArmadillo` and `Rcpp`: These two packages allow convenient interface between `R` and `C++`. `RcppArmadillo` requires `Rcpp`, so installing the first one will automatically prompt `R` to install the latter.
	install.packages("RcppArmadillo")     

