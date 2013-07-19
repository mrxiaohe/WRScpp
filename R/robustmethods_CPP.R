require("RcppArmadillo")
require("Rcpp")

stsreg_C<-function(x,y,xout=FALSE,outfun=out,iter=10,sc=pbvar,varfun=pbvar,
corfun=pbcor,plotit=FALSE,...){
#
#  Compute Theil-Sen regression estimator
#
#  Use Gauss-Seidel algorithm
#  when there is more than one predictor
#
#
	x<-as.matrix(x)
	xx<-cbind(x,y)
	xx<-elimna(xx)
	x<-xx[,1:ncol(x)]
	x<-as.matrix(x)
	y<-xx[,ncol(x)+1]
	temp<-NA
	x<-as.matrix(x)
	if(xout){
		x<-as.matrix(x)
		flag<-outfun(x,plotit=plotit,...)$keep
		x<-x[flag,]
		y<-y[flag]
		x<-as.matrix(x)
	}
	if(ncol(x)==1){
		temp1<-.Call("stsregp1_C", X=x,Y=y)
		coef<-temp1$coef
		res<-temp1$res
	}
	if(ncol(x)>1){
		temp1<-.Call("stsreg_for", X=x,Y=y, IT=as.integer(iter))
		coef<-c(temp1$alpha, temp1$beta)
		res<-temp1$res
	}
	yhat<-y-res
	stre=NULL
	e.pow<-varfun(yhat)/varfun(y)
	if(!is.na(e.pow)){
		if(e.pow>=1)
			e.pow<-corfun(yhat,y)$cor^2
		e.pow=as.numeric(e.pow)
		stre=sqrt(e.pow)
	}
	list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

tstsreg_C<-function(x,y,sc=pbvar,xout=FALSE,outfun=out,plotit=FALSE,...){
	#
	# Compute a modified Theil-Sen regression estimator.
	# Use s-type initial estimate, eliminate points with
	# outlying residuals, then do regular Theil-Sen
	#
	require("RcppArmadillo")
	require("Rcpp")

	x<-as.matrix(x)
	xx<-cbind(x,y)
	xx<-elimna(xx)
	x<-xx[,1:ncol(x)]
	x<-as.matrix(x)
	y<-xx[,ncol(x)+1]
	x<-as.matrix(x)
	if(xout){
		x<-as.matrix(x)
		flag<-outfun(x,plotit=plotit,...)$keep
		x<-x[flag,]
		y<-y[flag]
		x<-as.matrix(x)
	}
	res=stsreg_C(x,y)$res
	chk<-abs(res-median(res))/mad(res)
	xx<-x[chk<=2,]
	yy<-y[chk<=2]
	temp<-tsreg(xx,yy)
	list(coef=temp$coef,residuals=temp$res)
}

tshdreg_C<- function(x,y,HD=TRUE,xout=FALSE,outfun=out,iter=10,varfun=pbvar,
corfun=pbcor,plotit=FALSE,tol=.0001,...){
	#
	#  Compute Theil-Sen regression estimator
	#
	#  Use back-fitting
	#  when there is more than one predictor
	#  and estimate intercept using Harrel-Davis estimator
	#
	require("RcppArmadillo")
	require("Rcpp")

	x<-as.matrix(x)
	xx<-cbind(x,y)
	xx<-elimna(xx)
	x<-xx[,1:ncol(x)]
	x<-as.matrix(x)
	y<-xx[,ncol(x)+1]
	temp<-NA
	x<-as.matrix(x)
	if(xout){
		x<-as.matrix(x)
		flag<-outfun(x,plotit=plotit,...)$keep
		x<-x[flag,]
		y<-y[flag]
		x<-as.matrix(x)
	}
	if(ncol(x)==1){
		coef<-.Call("tshd_C", X=x, Y=y, hd=as.integer(HD))
		res<-y-coef[2]*x-coef[1]
	}
	if(ncol(x)>1){
		temp1<-.Call("tshdreg_for", X=x,Y=y, IT=as.integer(iter), TOL=tol)
		coef<-c(temp1$alpha, temp1$beta)
		res<-temp1$res
	}
	yhat<-y-res
	stre=NULL
	temp=varfun(y)
	if(temp==0)print('Warning: When computing strength of association, measure of variation=0')
	e.pow=NULL
	if(temp>0){
		e.pow<-varfun(yhat)/varfun(y)
		if(!is.na(e.pow)){
			if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
			e.pow=as.numeric(e.pow)
			stre=sqrt(e.pow)
		}
	}
	res=NULL
	list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}


tsreg_C<-function(x,y,xout=FALSE,outfun=out,iter=10,varfun=pbvar,
corfun=pbcor,plotit=FALSE,WARN=TRUE,HD=FALSE,...){
	#
	#  Compute Theil-Sen regression estimator
	#
	#  Use Gauss-Seidel algorithm
	#  when there is more than one predictor
	#
	#
	require("RcppArmadillo")
	require("Rcpp")

	x<-as.matrix(x)
	xx<-cbind(x,y)
	xx<-elimna(xx)
	x<-xx[,1:ncol(x)]
	x<-as.matrix(x)
	y<-xx[,ncol(x)+1]
	temp<-NA
	x<-as.matrix(x)
	if(xout){
		x<-as.matrix(x)
		flag<-outfun(x,plotit=plotit,...)$keep
		x<-x[flag,]
		y<-y[flag]
		x<-as.matrix(x)
	}
	if(ncol(x)==1){
		temp1<-.Call("tsp1reg_C", X=x, Y=y, HD=as.integer(HD))
		coef<-temp1$coef
		res<-temp1$res
	} 
	if(ncol(x)>1){
		temp1<-.Call("tsreg_for", X=x, Y=y, IT=as.integer(iter), HD=as.integer(HD))
		coef<-c(temp1$alpha,temp1$beta)
		res<-temp1$res
	}
	yhat<-y-res
	stre=NULL
	temp=varfun(y)
	if(temp==0){
		if(WARN)print("Warning: When computing strength of association, measure of variation=0")
	}
	e.pow=NULL
	if(temp>0){
		e.pow<-varfun(yhat)/varfun(y)
		if(!is.na(e.pow)){
			if(e.pow>=1)e.pow<-corfun(yhat,y)$cor^2
			e.pow=as.numeric(e.pow)
			stre=sqrt(e.pow)
		}
	}
	list(coef=coef,residuals=res,Strength.Assoc=stre,Explanatory.Power=e.pow)
}

outmgv_C<-function(x,y=NULL,plotit=TRUE,outfun=outbox,se=TRUE,op=1,
cov.fun=rmba,xlab="X",ylab="Y",SEED=TRUE,STAND=FALSE,...){
	#
	# Check for outliers using mgv method
	#
	# NOTE: if columns of the input matrix are reordered, this can
	# have an effect on the results due to rounding error when calling
	# the R function eigen.
	#
	#  (Argument STAND is included simply to avoid programming issues when outmgv is called by other 
	#  functions.)
	#
	require("RcppArmadillo")
	require("Rcpp")

	if(is.null(y[1]))m<-x
	if(!is.null(y[1]))m<-cbind(x,y)
	m=elimna(m)
	m=as.matrix(m)
	nv=nrow(m)
	temp<-mgvar_C(m,se=se,op=op,cov.fun=cov.fun,SEED=SEED)
	#if(fast)temp<-mgvdep.for(m,se=se)$distance
	temp[is.na(temp)]<-0
	if(ncol(m)==1){
		temp2=outpro(m)
		nout=temp2$n.out
		keep=temp2$keep
		temp2=temp2$out.id
	}
	if(ncol(m)>1){
		if(ncol(m)==2)temp2<-outfun(temp,...)$out.id
		if(ncol(m)>2)temp2<-outbox(temp,mbox=TRUE,gval=sqrt(qchisq(.975,ncol(m))))$out.id
		vec<-c(1:nrow(m))
		flag<-rep(T,nrow(m))
		flag[temp2]<-F
		vec<-vec[flag]
		vals<-c(1:nrow(m))
		keep<-vals[flag]
		if(plotit && ncol(m)==2){
			x<-m[,1]
			y<-m[,2]
			plot(x,y,type="n",xlab=xlab,ylab=ylab)
			flag<-rep(T,length(y))
			flag[temp2]<-F
			points(x[flag],y[flag],pch="*")
			points(x[temp2],y[temp2],pch="o")
		} else if(plotit && ncol(m)==3){
			require("scatterplot3d")
			x1<-m[,1]
			x2<-m[,2]
			y<-m[,3]
			flag<-rep(T,length(y))
			flag[temp2]<-F
			scatter3d<-scatterplot3d(x1, x2, y, type="n")
 			scatter3d$points3d(x1[flag], x2[flag], y[flag], pch="*")
 			scatter3d$points3d(x1[temp2], x2[temp2], y[temp2], pch="+", col="red")
		}
		nout=0
		if(!is.na(temp2[1]))nout=length(temp2)
	}
	list(n=nv,n.out=nout,out.id=temp2,keep=keep)
}

mgvar_C<-function(m,se=FALSE,op=0,cov.fun=covmve,SEED=TRUE){
#
# Find the center of a scatterplot, add point that
# increases the generalized variance by smallest amount
# continue for all points
# return the generalized variance
#  values corresponding to each point.
# The central values and point(s) closest to it get NA
#
# op=0 find central points using pairwise differences
# op!=0 find central points using measure of location
# used by cov.fun
#
# choices for cov.fun include
# covmve
# covmcd
# tbs (Rocke's measures of location
# rmba (Olive's median ball algorithm)
#
require("RcppArmadillo")
require("Rcpp")

	if(op==0)temp<-apgdis(m,se=se)$distance
	if(op!=0)temp<-out(m,cov.fun=cov.fun,plotit=FALSE,SEED=SEED)$dis
	flag<-(temp!=min(temp))
	temp2<-temp
	temp2[!flag]<-max(temp)
	flag2<-(temp2!=min(temp2))
	flag[!flag2]<-F
	if(sum(flag)>0)
		varvec<-.Call("mgvar_while", X=as.numeric(flag), M=m)
	else varvec<-NA
}


fdepthv2_C<-function(m,pts=NA,plotit=TRUE){
#
# Determine depth of points in pts relative to
# points in m
#
# Draw a line between each pair of distinct points
# and determine depth of the projected points.
# The final depth of a point is its minimum depth
# among all projections.
#
# This function is slower than fdepth and requires
# space for a nc by nc matrix, nc=(n^2-n)/2.
# But it allows
# data to have a singular covariance matrix
# and it provides a more accurate approximation of
# halfspace depth.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data and pts=NA
#
#  When plotting,
#  center is marked with a cross, +.
#
require("RcppArmadillo")
require("Rcpp")
	m<-elimna(m) # Remove missing values
	if(!is.na(pts[1]))remm<-m
	if(!is.matrix(m))dep<-unidepth(m)
	if(is.matrix(m)){
		nm<-nrow(m)
		nt<-nm
		nm1<-nm+1
	    if(!is.na(pts[1])){
    		if(ncol(m)!=ncol(pts))
    			stop("Number of columns of m is not equal to number of columns for pts")
			nt<-nm+nrow(pts)
			}
		}
	    if(ncol(m)==1)depth<-unidepth(m)
    	if(ncol(m)>1){
			m<-elimna(m) # Remove missing values
			nc<-(nrow(m)^2-nrow(m))/2
		  #  if(is.na(pts[1]))mdep <- matrix(0,nrow=nc,ncol=nrow(m))
		   # if(!is.na(pts[1])){
			#	mdep <- matrix(0,nrow=nc,ncol=nrow(pts))
			#}
		#ic<-0
		if(is.na(pts[1])) pts=matrix(, 2,2)
		mdep<-t(.Call("fdepthv2_for", M=m, PTS=pts))
		dep<-apply(mdep,2,min)
	}
    	if(ncol(m)==2 &&is.na(pts[1])){
			flag<-chull(m)
			dep[flag]<-min(dep)
		}
	    if(ncol(m)==2){
    		if(is.na(pts[1]) && plotit){
				plot(m, pch="+", cex=0.7)
				x<-m
				temp<-dep
				flag<-(temp>=median(temp))
				xx<-x[flag,]
				xord<-order(xx[,1])
				xx<-xx[xord,]
				temp<-chull(xx)
				xord<-order(xx[,1])
				xx<-xx[xord,]
				temp<-chull(xx)
				lines(xx[temp,], col="red")
				lines(xx[c(temp[1],temp[length(temp)]),], , col="red")
			}
		}
		dep
}


outpro_C<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
xlab="VAR 1",ylab="VAR 2",STAND=FALSE,tr=.2,q=.5,pr=TRUE,...){
#
# Detect outliers using a modification of the
# Stahel-Donoho  projection method.
#
# Determine center of data cloud, for each point,
# connect it with center, project points onto this line
# and use distances between projected points to detect
# outliers. A boxplot method is used on the
# projected distances.
#
# plotit=TRUE creates a scatterplot when working with
# bivariate data.
#
# op=T
# means the .5 depth contour is plotted
# based on data with outliers removed.
#
# op=F
# means .5 depth contour is plotted without removing outliers.
#
#  MM=F  Use interquatile range when checking for outliers
#  MM=T  uses MAD.
#
#  If value for center is not specified,
#  there are four options for computing the center of the
#  cloud of points when computing projections:
#
#  cop=2 uses MCD center
#  cop=3 uses median of the marginal distributions.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)#  cop=7 uses the spatial (L1) median
#
#  args q and tr having are not used by this function. They are included to deal
#  with situations where smoothers have optional arguments for q and tr
#
#  When using cop=2, 3 or 4, default critical value for outliers
#  is square root of the .975 quantile of a
#  chi-squared distribution with p degrees
#  of freedom.
#
#  STAND=T means that marginal distributions are standardized before
#  checking for outliers.
#
#  Donoho-Gasko (Tukey) median is marked with a cross, +.
#
	m<-as.matrix(m)
	if(pr){
		if(!STAND){
			#if(ncol(m)>1)cat("STAND=FALSE. If measures are on different scales,", 
			#					"might want to use STAND=TRUE\n")
		}
	}
	library(MASS)
	m=elimna(m)
	m<-as.matrix(m)
	nv=nrow(m)
	if(ncol(m)==1){
		dis<-(m-median(m,na.rm=TRUE))^2/mad(m,na.rm=TRUE)^2
		dis<-sqrt(dis)
		dis[is.na(dis)]=0
		crit<-sqrt(qchisq(.975,1))
		chk<-ifelse(dis>crit,1,0)
		vec<-c(1:nrow(m))
		outid<-vec[chk==1]
		keep<-vec[chk==0]
	}
	if(ncol(m)>1){
		if(STAND)m=standm(m,est=median,scat=mad)
		if(is.na(gval) && cop==1)gval<-sqrt(qchisq(.95,ncol(m)))
		if(is.na(gval) && cop!=1)gval<-sqrt(qchisq(.975,ncol(m)))
		if(cop==1 && is.na(center[1])){
		if(ncol(m)>2)center<-dmean(m,tr=.5,cop=1)
		if(ncol(m)==2){
			tempd<-NA
			for(i in 1:nrow(m))
			tempd[i]<-depth(m[i,1],m[i,2],m)
			mdep<-max(tempd)
			flag<-(tempd==mdep)
			if(sum(flag)==1)center<-m[flag,]
			if(sum(flag)>1)center<-apply(m[flag,],2,mean)
		}
	}
	if(cop==2 && is.na(center[1])){
		center<-cov.mcd(m)$center
	}
	if(cop==4 && is.na(center[1])){
		center<-cov.mve(m)$center
	}
	if(cop==3 && is.na(center[1])){
		center<-apply(m,2,median)
	}
	if(cop==5 && is.na(center[1])){
		center<-tbs(m)$center
	}
	if(cop==6 && is.na(center[1])){
		center<-rmba(m)$center
	}
	if(cop==7 && is.na(center[1])){
		center<-spat(m)
	}
	outid.flag<-.Call("outpro_for", 
				 M=m, 
				 GVAL=gval, 
				 CENTER=center, 
				 MM=MM
				 )
	idv<-1:nrow(m)
	outid<-idv[outid.flag]
	keep<-idv[!outid.flag]
	if(ncol(m)==2){
		if(plotit){
			plot(m[,1],m[,2],type="n",xlab=xlab,ylab=ylab)
			points(m[keep,1],m[keep,2],pch="*")
			if(length(outid)>0)points(m[outid,1],m[outid,2],pch="o")
				if(op){
					tempd<-NA
					keep<-keep[!is.na(keep)]
					mm<-m[keep,]
					for(i in 1:nrow(mm))tempd[i]<-depth(mm[i,1],mm[i,2],mm)
					mdep<-max(tempd)
					flag<-(tempd==mdep)
					if(sum(flag)==1)center<-mm[flag,]
					if(sum(flag)>1)center<-apply(mm[flag,],2,mean)
					m<-mm
				}
				points(center[1],center[2],pch="+")
				x<-m
				temp<-fdepth(m,plotit=FALSE)
				flag<-(temp>=median(temp))
				xx<-x[flag,]
				xord<-order(xx[,1])
				xx<-xx[xord,]
				temp<-chull(xx)
				xord<-order(xx[,1])
				xx<-xx[xord,]
				temp<-chull(xx)
				lines(xx[temp,])
				lines(xx[c(temp[1],temp[length(temp)]),])
			}
		}
	}
	list(n=nv,n.out=length(outid),out.id=outid,keep=keep)
}



skip_C<-function(m,cop=6,MM=FALSE,op=3,mgv.op=0,outpro.cop=3,...){
#
# m is an n by p matrix
#
# Compute skipped location and covariance matrix
#
# op=1:
# Eliminate outliers using a projection method
# That is, first determine center of data using:
#
# cop=1 Donoho-Gasko median,
# cop=2 MCD,
# cop=3 marginal medians.
#  cop=4 uses MVE center
#  cop=5 uses TBS
#  cop=6 uses rmba (Olive's median ball algorithm)
#
# For each point
# consider the line between it and the center,
# project all points onto this line, and
# check for outliers using
#
# MM=F, a boxplot rule.
# MM=T, rule based on MAD and median
#
# Repeat this for all points. A point is declared
# an outlier if for any projection it is an outlier
#
# op=2 use mgv (function outmgv) method to eliminate outliers
#
# Eliminate any outliers and compute means
#  using remaining data.
# mgv.op=0, mgv uses all pairwise distances to determine center of the data
# mgv.op=1 uses MVE
# mgv.op=2 uses MCD
#
temp<-NA
m<-elimna(m)
if(op==1)temp<-outpro(m,plotit=FALSE,MM=MM,cop=outpro.cop,...)$keep
if(op==2)temp<-outmgv(m,plotit=FALSE,op=mgv.op,...)$keep
if(op==3)temp<-outpro_C(m,plotit=FALSE,MM=MM,cop=outpro.cop,...)$keep
if(op==4)temp<-outmgv_C(m,plotit=FALSE,op=mgv.op,...)$keep
val<-var(m[temp,])
loc<-apply(m[temp,],2,mean)
list(center=loc,cov=val)
}


skip_boot <- function(x, data, cop, MM, nboot){
	if( cop == 1 && ncol(x) > 2 )
		stop("the funciton `dmean()` hasn't been implemented in C++ yet.")
	if( cop != 3 & cop != 6 )
		stop(" `cov.mcd()`, `cov.mve()`, `tbs()`, and `spat()` have not been implemented in C++ yet.")
	bvec <- .Call("skip_boot", x, MM, cop, data) 
	matrix(bvec, nboot, byrow=TRUE)
}



ddepGMC_version2<-function(x,est=skipSPR,na.rm=TRUE,alpha=.05,grp=NA,nboot=500,plotit=TRUE,SEED=TRUE,STAND=FALSE ,CENTER=TRUE,...){
	require("RcppArmadillo")
	extras <- list(...)      #Capture dot arguments for the skip function
	library(parallel)
	if(is.list(x)){
		nv<-NA
		for(j in 1:length(x))
			nv[j]<-length(x[[j]])
			if(var(nv) !=0){
				stop('The groups are stored in list mode and appear to have different sample sizes')
			}
			temp<-matrix(NA,ncol=length(x),nrow=nv[1])
			for(j in 1:length(x))
				temp[,j]<-x[[j]]
			x<-temp
	}
	J<-ncol(x)
	if(!is.na(grp[1])){       #Select the groups of interest
		J<-length(grp)
		for(j in 1:J)
			temp[,j]<-x[,grp[j]]
		x<-temp
	}
	if(na.rm)
		x<-elimna(x) # Remove any rows with missing values.
	bvec<-matrix(0,ncol=J,nrow=nboot)
	hval=NA
	if(SEED)set.seed(2) # set seed of random number generator so that
						# results can be duplicated.
	n<-nrow(x)
	flag=identical(est,skip)
	#if(flag)

	##ADDED THIS IF CONDITIONAL WHEN skip_boot() is chosen.
	if(identical(est, skip_boot)){                                    
		m.arg 		<- match(c("MM", "cop"), names(extras), 0L)
		m.dot 		<- extras[m.arg]
		args.list 	<- modifyList(list(x, data=list(1:nrow(x)-1),  MM= FALSE, nboot=1, cop = 3), m.dot)
		totv 		<- as.vector(do.call("skip_boot", args.list))
	} else
		totv<-skip(x,...)$center 
	gv   <-rep(mean(totv),J)  #Grand mean
	data <-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
	if(!identical(est, skip_boot)){
		data <- listm(t(data))
		#bvec <- mclapply(data,ddepGMC_sub_cen,x,est=est,mc.preschedule=TRUE)
		bvec <- lapply(data,ddepGMC_sub_cen,x,est=est)
		bvec <- t(matl(bvec))
	}else if(identical(est, skip_boot)){
		data 	  <- listm(t(data - 1))
		args.list$data  <- data
		args.list$nboot <- nboot
		bvec 	  <- do.call("skip_boot", args.list)	
		bvec 	  <- matrix(bvec, nboot, byrow=FALSE)
	}
	
	bplus<-nboot+1
	m1<-rbind(bvec,gv)
	center<-totv
	cmat<-var(bvec)
	discen<-mahalanobis( m1, totv, cmat)
	if(plotit && ncol(x)==2){
		plot(bvec,xlab='Group 1',ylab='Group 2')
		temp.dis<-order(discen[1:nboot])
		ic<-round((1-alpha)*nboot)
		xx<-bvec[temp.dis[1:ic],]
		xord<-order(xx[,1])
		xx<-xx[xord,]
		temp<-chull(xx)
		lines(xx[temp,])
		lines(xx[c(temp[1],temp[length(temp)]),])
		abline(0,1)
	}
	sig.level<-sum(discen[bplus]<=discen)/bplus
	list(p.value=sig.level,center=totv,weighted.grand.mean=gv[1])
}

ddepGMC_sub<-function(data,est,x,...){
	v=est(x[data,],...)$center
	v
}

ddepGMC_sub_cen<-function(data,est,x,...){
	v=est(x[data,],...)$center
	v
}

ancGLOB_sub2<-function(bvec,fr1=fr1,fr2=fr2,est=est,SEED=SEED,...){
	p=ancGLOB_sub3(bvec[,1],bvec[,2],bvec[,3],bvec[,4],est=est,SEED=SEED,fr1=fr1,fr2=fr2,plotit=FALSE,...)$p.value
	p
}

aov2depth<-function(x1,x2,est=tmean,nboot=500,SEED=TRUE,...){
#
# 2 by K ANOVA independent groups
#
#   Main effect Factor A only
#
# Strategy: Use depth of zero based on estimated
# differences for each column  of the K levels of Factor B
# That is, testing no main effects for Factor A in 
# a manner that takes into account the pattern of the
# measures of location rather then simply averaging 
# across columns.
#
	if(is.matrix(x1)||is.data.frame(x1))x1=listm(x1)
	if(is.matrix(x2)||is.data.frame(x2))x2=listm(x2)
	J=length(x1)
	if(J!=length(x2))
		stop('x1 and x2 should have same number of groups')
	if(SEED)set.seed(2)
	for(j in 1:J){
		x1[[j]]=na.omit(x1[[j]])
		x2[[j]]=na.omit(x2[[j]])
	}
	n1=mapply(x1,FUN=length)
	n2=mapply(x2,FUN=length)
	bplus=nboot+1
	bvec1=matrix(NA,nrow=nboot,ncol=J)
	bvec2=matrix(NA,nrow=nboot,ncol=J)
	for(j in 1:J){
		data1=matrix(sample(x1[[j]],size=n1[j]*nboot,replace=TRUE),nrow=nboot)
		data2=matrix(sample(x2[[j]],size=n2[j]*nboot,replace=TRUE),nrow=nboot)
		bvec1[,j]=apply(data1,1,est,...)
		bvec2[,j]=apply(data2,1,est,...)
	}

	difb=bvec1-bvec2
	est1=mapply(x1,FUN=est,...)
	est2=mapply(x2,FUN=est,...)
	dif=est1-est2
	m1=var(difb)
	nullvec=rep(0,J)
	difz=rbind(difb,nullvec)
	dis=mahalanobis(difz,dif,m1)
	sig=sum(dis[bplus]<=dis)/bplus
	list(p.value=sig,est1=est1,est2=est2,dif=dif,n1=n1,n2=n2)
}


ancGLOB_sub3<-function(x1,y1,x2,y2,xout=FALSE,outfun=outpro,est=tmean,pcrit=NULL,p.crit=NULL,iter=100,
nboot=500,SEED=TRUE,MC=FALSE,nmin=12,pts=NULL,fr1=1,fr2=1,plotit=TRUE,xlab='X',ylab='Y',LP=TRUE,...){
#
#
	if(SEED)set.seed(2)
	x1<-as.matrix(x1)
	p1<-ncol(x1)+1
	p<-ncol(x1)
	if(p>1)
		stop('Current version is for one independent variable only')
	xy<-cbind(x1,y1)
	xy<-elimna(xy)
	x1<-xy[,1:p]
	y1<-xy[,p1]
	xy<-cbind(x2,y2)
	xy<-elimna(xy)
	x2<-xy[,1:p]
	y2<-xy[,p1]
	if(xout){
		m<-cbind(x1,y1)
		flag<-outfun(x1,plotit=FALSE,...)$keep
		m<-m[flag,]
		x1<-m[,1:p]
		y1<-m[,p1]
		m<-cbind(x2,y2)
		flag<-outfun(x2,plotit=FALSE,...)$keep                            
		m<-m[flag,]                                                        
		x2<-m[,1:p]                                                              
		y2<-m[,p1] 
	}
	N1=length(y1)
	N2=length(y2)
	if(is.null(pts[1])){
		isub<-c(1:5)  # Initialize isub
		test<-c(1:5)
		xorder<-order(x1)
		y1<-y1[xorder]
		x1<-x1[xorder]
		xorder<-order(x2)
		y2<-y2[xorder]
		x2<-x2[xorder]
		n1<-1
		n2<-1
		vecn<-1
		for(i in 1:length(x1))
			n1[i]<-length(y1[near(x1,x1[i],fr1)])
		for(i in 1:length(x1))
			n2[i]<-length(y2[near(x2,x1[i],fr2)])
		for(i in 1:length(x1)){
			vecn[i]<-min(n1[i],n2[i])
		}
		sub<-c(1:length(x1))
		isub[1]<-min(sub[vecn>=nmin])
		isub[5]<-max(sub[vecn>=nmin])
		isub[3]<-floor((isub[1]+isub[5])/2)
		isub[2]<-floor((isub[1]+isub[3])/2)
		isub[4]<-floor((isub[3]+isub[5])/2)
		pts=x1[isub]
		g1=list()
		g2=list()
		for (i in 1:5){
			g1[[i]]<-y1[near(x1,x1[isub[i]],fr1)]
			g2[[i]]<-y2[near(x2,x1[isub[i]],fr2)]
		}
	}
	if(!is.null(pts[1])){
		if(length(pts)<2)stop('Should have at least two points (use the R function ancova)')
		g1=list()
		g2=list()
		for (i in 1:length(pts)){
			g1[[i]]<-y1[near(x1,pts[i],fr1)]
			g2[[i]]<-y2[near(x2,pts[i],fr2)]
		}
	}
	n1=lapply(g1,length)
	res=aov2depth(g1,g2,est=est,SEED=SEED,nboot=nboot, ...)
	if(plotit)
		runmean2g(x1,y1,x2,y2,nboot=nboot,fr=fr1,est=est,xout=xout,...)
	list(p.value=res$p.value,est1=res$est1,est2=res$est2,dif=res$dif,pts=pts,n1=res$n1,n2=res$n2)
}



ancGLOB_pv_version2<-function(n1,n2,est=tmean,fr1=.8,fr2=.8,nboot=500,SEED=TRUE,iter=100, cpp=FALSE, ...){
	require("RcppArmadillo")   #ADDED 
	extras <- list(...)        #ADDED
	if(SEED)
		set.seed(45)
	bvec=list()
	np1=min(c(n1,n2))+1
	nmax=max(c(n1,n2))
	library(parallel)

	for(i in 1:iter){
		bvec[[i]]=rmul(nmax,p=4)
		if(n1!=n2)
			bvec[[i]][np1:nmax,1:2]=NA
	}
	##### I ADDED AN IF STATEMENT HERE SO THAT WHEN cpp IS SET TO `TRUE` the function 
	##### CALLS THE C++ CODE
	if(!cpp){
		prm=mclapply(bvec,ancGLOB_sub2,fr1=fr1,fr2=fr2,est=est,SEED=SEED,...)
	} else {
		if(identical(est, tmean))
			estimator <- 1
		else if(identical(est, median))
			estimator <- 2
		else if(identical(est, mean))
			estimator <- 3
		else 
			stop("Only mean, trimmed mean, and median are implemented for the C++ code so far.")
		dot.names <- names(extras)
		if( sum("tr" != dot.names) >0 )
			stop("The following arguments are not available for the C++ code: ", 
                  paste(dot.names, collapse=", "))
        if(!("tr" %in% dot.names))
        	tr <- 0.2
        else 
        	tr <- extras[["tr"]]
		bvec <- lapply(lapply(bvec, as.vector), na.omit)
		dot.names <- names(extras) 

		###############################  CALLING C++ ROUTINE ###############################
		prm  <- .Call("ancGLOB_sub2_C", bvec, n1, n2, fr1, fr2, NULL, estimator, tr)
		####################################################################################

	}
	prm=as.vector(matl(prm))
	prm=sort(elimna(prm))
	prm
}
