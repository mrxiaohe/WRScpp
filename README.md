WRScpp
======

This package provides C++ sub-routines for several iterative procedures in the [R package `WRS`](https://r-forge.r-project.org/projects/wrs/) for robust statistics by [Rand Wilcox](http://dornsife.usc.edu/cf/labs/wilcox/wilcox-faculty-display.cfm). These sub-routines provide substantial performance boosts.

##Installation

####To install this package, three R packages are required:

* `RcppArmadillo` and `Rcpp`: These two packages allow convenient interface between `R` and `C++`. `RcppArmadillo` requires `Rcpp`, so installing the first one will automatically prompt `R` to install the latter.
* `devtools`: This package allows R users to install packages hosted on Github. 


      install.packages( c("RcppArmadillo", "devtools") )



####Next, load `devtools` and install the `WRScpp` binary:

	library("devtools")
	install_github(repo="WRScpp", username="mrxiaohe")  


####Lastly, load `WRScpp`:

	library("WRScpp")




