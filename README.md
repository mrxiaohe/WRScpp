WRScpp
======

This package provides C++ sub-routines for several iterative procedures in the [R package `WRS`](https://r-forge.r-project.org/projects/wrs/) for robust statistics by [Rand Wilcox](http://dornsife.usc.edu/cf/labs/wilcox/wilcox-faculty-display.cfm). These sub-routines provide substantial performance boosts.

##Installation

####To install this package, three `R` packages are required:

* `RcppArmadillo` and `Rcpp`: These two packages allow convenient interface between `R` and `C++`. `RcppArmadillo` requires `Rcpp`, so installing the first one will automatically prompt `R` to install the latter.
* `devtools`: This package allows R users to install packages hosted on Github. 


        install.packages( c("RcppArmadillo", "devtools") )



####Next, load `devtools` and install the `WRScpp` binary:

	library("devtools")
	install_github(repo="WRScpp", username="mrxiaohe")  


####Lastly, load `WRScpp`:

	library("WRScpp")



##Package content:

    ls("package:WRScpp")
 [1] "ancGLOB_pv_C"    "ancGLOB_sub2"    "ancGLOB_sub3"    "aov2depth"       "ddepGMC_C"       "ddepGMC_sub"     "ddepGMC_sub_cen" "fdepthv2_C"      "mgvar_C"         "outmgv_C"       
[11] "outpro_C"        "skip_boot"       "skip_C"          "stsreg_C"        "tshdreg_C"       "tsreg_C"         "tstsreg_C"      
