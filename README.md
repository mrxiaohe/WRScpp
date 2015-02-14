WRScpp
======


This package provides `C++` sub-routines for several iterative procedures in the `R` package `WRS` for robust statistics by [Dr. Rand Wilcox](http://dornsife.usc.edu/cf/labs/wilcox/wilcox-faculty-display.cfm). These `C++` sub-routines can provide substantial performance boosts. The raw code is in [here](https://github.com/mrxiaohe/robustmethods_cplusplus). 

###[Installation]

####To install this package, four R packages are required:

* `WRS`
* `RcppArmadillo` and `Rcpp`: These two packages allow convenient interface between R and C++. `RcppArmadillo` requires `Rcpp`, so installing the first one will automatically prompt R to install the latter.
* `devtools`: This package allows R users to install packages hosted on Github.

        install.packages("WRS", repos="http://R-Forge.R-project.org", type="source")
        install.packages( c("RcppArmadillo", "devtools") )

####Next, load devtools and install the WRScpp binary:

    library("devtools")
    install_github(repo="WRScpp", username="mrxiaohe")  


###[Examples]

Load `WRScpp` and `WRS`:

    library("WRScpp")
    library("WRS")

Also load `rbenchmark` for performance benchmarks:

    #install.packages("rbenchmark")     #If it hasn't been installed yet.
    library("rbenchmark")

Create a mock dataset called dataset1 using the code below:

    dataset1<-structure(list(y=c(-1.49719416897806,  -1.04139301128124,   0.0945236304367192, 2.34517600028403,  -2.01910934921224, 
                                  0.738818483718415, -1.29894622500239,  -2.28749064172037,  -1.83460914213215,   1.7264401742697, 
                                 -0.282651383157718, -0.86805194949488,  -0.864427449770363,  1.3890252146112,   -1.06755713953516, 
                                 -2.32897163610783,   1.66038960720445,   1.61487788113575,   1.0899088587434,    1.59281468275878, 
                                 -3.92950882381207,  -2.86897048491236,   0.756543092818702,  0.75605628995123,  -1.94551222893336, 
                                 -0.521647504402827,  0.750976521143092,  2.05661942206825,  -0.0553654356461042, 1.83487282390805, 
                                 -0.681331876121451, -1.09742432884331,   0.809885352209668, -1.75444757180733,  -1.32646138930983, 
                                  0.93932182852857,   2.69573032525144,   1.19980666244453,   2.66656782068436,  -2.04469844361269, 
                                  2.25443895842938,  -0.287436324399811, -0.646793604946014, -1.12134602060874,  -0.123057506010716, 
                                  3.40147762961413,  -0.105013757259212,  2.57888932230486,  -0.248784082384957, -0.80495910311552, 
                                 -1.62736886754023,  -1.80623698171649,   0.0826863988021116,-1.02443966620053,  -0.624640883344864), 
                           x = c(-0.414242088066754,  0.0475817644136152,-0.815134097927552,  0.578444566236884, -0.459034615425856, 
                                  0.809378901441766, -0.977367705972685, -1.67738162375379,  -0.430291674794372,  1.42019731595199, 
                                 -0.749035091787927, -1.35225060418292,  -0.492196888552811,  0.539778435067201,  0.0459914831995694, 
                                 -1.2521152391345,    1.00360926349214,   1.09644063218488,   0.508112310623724,  1.32898360956295, 
                                 -1.5116989382724,   -1.80654236075974,   1.06240878829949,  -0.091275817909321, -0.570211504099203, 
                                 -0.268748606177556,  0.764654371144826,  1.9589087924292,   -0.906694839411102,  0.450015221346758, 
                                 -0.254560403793998, -1.30082788198172,   1.15863187888594,  -1.61536768356543,  -0.858653684897159, 
                                  1.13914138747698,   0.775779787679006,  0.267099868753253,  1.32191674807275,  -0.23978424525859,
                                  1.93490877730752,   1.08101027023282,  -1.51349689876255,  -0.758066482747475, -0.592268993080326, 
                                  0.676081370171539, -0.342776942447623,  1.0517207133769,   -0.39294982859476,   0.622178484006425, 
                                  6.22656137599703,   6.45244734765708,   8.50469010153436,   6.24821960890362,   7.56808156590792
                                  )), .Names = c("y", "x"), row.names = c(NA, -55L), class = "data.frame")


####1. `tsreg_C()`: Theil-Sen regression estimator

Use Gauss-Seidel algorithm when there is more than one predictor.

(1). Create a scatter plot for this dataset:

    plot( y ~ x, dataset1 )

![plot](http://imageshack.us/a/img542/717/6ts.png)


(2). Run the least squares linear regression using `lm()`, and plot a regression line based on the result.

    model.lm <- lm( y ~ x, data = dataset1 )
    summary(model.lm)

    Coefficients:
                Estimate Std. Error t value Pr(>|t|)
    (Intercept) -0.17700    0.22667  -0.781    0.438
    x            0.13528    0.09723   1.391    0.170

    Residual standard error: 1.617 on 53 degrees of freedom
    Multiple R-squared:  0.03524,   Adjusted R-squared:  0.01703 
    F-statistic: 1.936 on 1 and 53 DF,  p-value: 0.1699

    abline( model.lm, col="red" )

![plot](http://imageshack.us/a/img812/6629/qg63.png)


(3). Run Theil-Sen regression estimator `tsreg_C()`, and plot a regression line based on the result.

    model.tsreg <- with(dataset1, tsreg(x, y))

    model.tsreg[-2]                 #Display results but suppress residuals
    $coef
     Intercept            
    -0.3281116  0.9554134 

    $Strength.Assoc
    [1] 0.715131

    $Explanatory.Power
    [1] 0.5114123


    abline( model.tsreg$coef, col="blue" )

![plot](http://imageshack.us/a/img542/9886/l5z.png)



(4). Compare runtime between `tsreg()` and `tsreg_C()` using `benchmark()`. `tsreg()` is entirely coded in `R`, whereas portions of `tsreg_C()` are coded in `C++`.

Create another, larger dataset

    set.seed(1); dataset2 <- matrix(rnorm(1000), ncol=5)
    head(dataset2)
               [,1]       [,2]       [,3]       [,4]        [,5]
    [1,] -0.6264538  0.4094018  1.0744410 -0.3410670 -1.08690882
    [2,]  0.1836433  1.6888733  1.8956548  1.5024245 -1.82608301
    [3,] -0.8356286  1.5865884 -0.6029973  0.5283077  0.99528181
    [4,]  1.5952808 -0.3309078 -0.3908678  0.5421914 -0.01186178
    [5,]  0.3295078 -2.2852355 -0.4162220 -0.1366734 -0.59962839
    [6,] -0.8204684  2.4976616 -0.3756574 -1.1367339 -0.17794799


    benchmark( replications = 100, 
               tsreg( x=dataset2[, 1:4], y=dataset2[, 5] ),
               tsreg_C( x=dataset2[, 1:4], y=dataset2[, 5] )
             )
                                         test replications elapsed relative user.self sys.self  user.child sys.child
    2 tsreg_C(dataset2[, 1:4], dataset2[, 5])          100   8.834    1.000     8.473    0.388           0         0
    1   tsreg(dataset2[, 1:4], dataset2[, 5])          100  43.232    4.894    38.848    4.520           0         0
    


####2. `tshdreg_C()`: A variation of Theil-Sen regression estimator

This function uses Harrell-Davis estimator rather than the usual sample median. Also, the intercept is taken to be the median of the residuals.

(1). Use this function on dataset1, and plot a regression line based on the result.

    model.tshdreg <- with(dataset1, tshdreg(x, y))
    model.tshdreg[-2]

    $coef
    [1] -0.09980145  0.94965936
    
    $Strength.Assoc
    [1] 0.7108241

    $Explanatory.Power
    [1] 0.5052709


    abline( model.tshdreg$coef, col="green")

![plot](http://imageshack.us/a/img842/784/sa15.png)


(2). Compare runtime between `tshdreg()` and `tshdreg_C()` using `benchmark()`.

    benchmark( replications = 100, 
               tshdreg( x=dataset2[, 1:4], y=dataset2[, 5] ),
               tshdreg_C( x=dataset2[, 1:4], y=dataset2[, 5] )
             )
                                           test replications elapsed relative user.self sys.self user.child sys.child
    2 tshdreg_C(dataset2[, 1:4], dataset2[, 5])          100  27.774    1.000    26.784    0.758          0         0
    1   tshdreg(dataset2[, 1:4], dataset2[, 5])          100 145.014    5.221   140.692    3.050          0         0


####3. `stsreg_C()`: A variation of Theil-Sen regression estimator

Slopes are selected such that some robust measure of variance of residuals is minimized. By default, percentage bend midvariance (see the function `pbvar()`) is minimized.

(1). Apply `stsreg_C()` to dataset1 and plot a regression line based on the result.

    model.stsreg <- with(dataset1, stsreg_C(x, y))
    model.stsreg[-2]

    $coef
    [1] -0.3424785  1.2573545

    $Strength.Assoc
    [1] 0.9411352

    $Explanatory.Power
    [1] 0.8857355

    abline( model.stsreg$coef, col="purple")

![plot](http://imageshack.us/a/img819/1843/gj5.png)




(2). Compare performance between stsreg() and stsreg_C() using system.time().

    system.time( stsreg_C(dataset2[,1:4], dataset2[,5]) )
      user  system elapsed 
    24.679   0.234  25.002 

    system.time( stsreg(dataset2[,1:4], dataset2[,5]) )
       user  system elapsed 
    839.375  25.478 867.326 


####4. `tstsreg_C()`: A modified version of Theil-Sen regression estimator

The function first uses `stsreg_C()` to compute the initial estimate, next uses the computed residuals to determine outliers, lastly does regular Theil-Sen on data with regression outliers removed.

(1). Apply `tstreg_C()` on dataset1; plot a regression line.

    model.tstsreg <- with(dataset1, tstsreg_C(x, y))
    model.tstsreg[-2]

    $coef
     Intercept            
    0.04935162 1.24334983 

    abline(model.tstsreg$coef, col="black")

![plot](http://imageshack.us/a/img825/6753/uqkh.png)


(2). Compare performance between `tstreg()` and `tstreg_C()` using `system.time()`.

    system.time( tstsreg_C(dataset2[,1:4], dataset2[,5]) )
      user  system elapsed 
    25.517   0.281  26.299 
    system.time( tstsreg(dataset2[,1:4], dataset2[,5]) )
       user  system elapsed 
    740.578  23.169 762.316 


[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/mrxiaohe/wrscpp/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

