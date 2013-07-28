stsreg_C<-function(x,y,xout=FALSE,outfun=out,iter=10,sc=pbvar,varfun=pbvar,
corfun=pbcor,plotit=FALSE,...){
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
	stsreg_out <- .Call("stsreg_C", X=as.matrix(x), Y=as.matrix(y), IT=as.integer(iter))
	coef <- stsreg_out$coef
	res  <- stsreg_out$res
	yhat <- y-res
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
	tshdreg_output <- .Call("tshdreg_C", as.matrix(x), as.matrix(y), as.integer(iter), as.double(tol), TRUE)
	coef <- tshdreg_output$coef
	res  <- tshdreg_output$res
	yhat <- y - res
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
	if(op==0)
		temp <- apgdis(m,se=se)$distance
	if(op!=0)
		temp <- out(m,cov.fun=cov.fun,plotit=FALSE,SEED=SEED)$dis
	flag  <- (temp!=min(temp))
	temp2 <- temp
	temp2[!flag] <- max(temp)
	flag2 <- (temp2 != min(temp2))
	flag[!flag2] <- F
	if(sum(flag) > 0)
		.Call("mgvar_while", X = flag, M=m)
	else 
		NA
}


fdepthv2_C <- function(m, pts=NA, plotit=TRUE){
	m <- as.matrix( elimna(m) )
	if(!is.na(pts[1])){
    	if(ncol(m) != ncol(pts))
    		stop("Number of columns of m is not equal to number of columns for pts")
	}
	if( ncol(m) == 1 )
	    	depth <- .Call("unidepth_C", as.vector( m ), as.vector( m ))
    if( ncol(m) > 1 ){
    	nm <- nrow( m )
    	if( is.matrix(pts) )
    		m <- rbind( m, pts )
    	mdep <- t(.Call("fdepthv2_C", m, nm ))
		dep<-apply(mdep,2,min)
	}
    if(ncol(m)==2 && is.na(pts[1])){
		flag<-chull(m)
		dep[flag]<-min(dep)
	}
	if(ncol(m)==2 && is.na(pts[1]) && plotit){
		plot(m, pch="+", cex=0.7)
		x <- m
		temp <- dep
		flag <- (temp>=median(temp))
		xx   <-x[flag,]
		xord <-order(xx[,1])
		xx   <-xx[xord,]
		temp <-chull(xx)
		xord <-order(xx[,1])
		xx   <-xx[xord,]
		temp <-chull(xx)
		lines(xx[temp,], col="red")
		lines(xx[c(temp[1],temp[length(temp)]),], , col="red")
	}
	dep
}


outpro_C<-function(m,gval=NA,center=NA,plotit=TRUE,op=TRUE,MM=FALSE,cop=3,
xlab="VAR 1",ylab="VAR 2",STAND=FALSE,tr=.2,q=.5,pr=TRUE,...){
	m<-as.matrix(m)
	if(pr){
		if(!STAND){
			#if(ncol(m)>1)cat("STAND=FALSE. If measures are on different scales,", 
			#					"might want to use STAND=TRUE\n")
		}
	}
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



ddepGMC_C<-function(x, est=skip_boot, na.rm=TRUE, alpha=.05, grp=NA, nboot=500, plotit=TRUE, SEED=TRUE,...){
	extras <- list(...)     
	require(parallel)
	x <- as.matrix(x)
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
	if(!is.na(grp[1])){     
		J<-length(grp)
		for(j in 1:J)
			temp[,j]<-x[,grp[j]]
		x<-temp
	}
	if(na.rm)
		x<-as.matrix(elimna(x))
	bvec<-matrix(0,ncol=J,nrow=nboot)
	hval=NA
	if(SEED)set.seed(2)
						
	n<-nrow(x)
	flag=identical(est,skip)
	#if(flag)

	if(identical(est, skip_boot)){                                    
		m.arg 		<- match(c("MM", "cop"), names(extras), 0L)
		m.dot 		<- extras[m.arg]
		args.list 	<- modifyList(list(x, data=list(1:nrow(x)-1),  MM= FALSE, nboot=1, cop = 3), m.dot)
		totv 		<- as.vector(do.call("skip_boot", args.list))
	} else
		totv <- est(x,...)$center 
	gv   <-rep(mean(totv),J)  #Grand mean
	data <-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
	if(!identical(est, skip_boot)){
		data <- listm(t(data))
		bvec <- mclapply(data,ddepGMC_sub_cen,x,est=est,mc.preschedule=TRUE)
		#bvec <- lapply(data,ddepGMC_sub_cen,x,est=est)
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
		plot(bvec,xlab='Group 1',ylab='Group 2', type="n")
		temp.dis<-order(discen[1:nboot])
		ic<-round((1-alpha)*nboot)
		xx<-bvec[temp.dis[1:ic],]
		points(xx,xlab='Group 1',ylab='Group 2', cex=0.8)
		pts.outside <- bvec[temp.dis[(ic+1):length(temp.dis)], ]
		points(pts.outside, col="red", pch="+")
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
		isub<-c(1:5)  
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
		sub.id <- sub[vecn>=nmin]
		if(length(sub.id) > 0){
			isub[1]<-min(sub[sub.id ])
			isub[5]<-max(sub[sub.id ])
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
	if(length(sub.id) > 0){
		n1=lapply(g1,length)
		res=aov2depth(g1,g2,est=est,SEED=SEED,nboot=nboot, ...)
		if(plotit)
			runmean2g(x1,y1,x2,y2,nboot=nboot,fr=fr1,est=est,xout=xout,...)
		list(p.value=res$p.value,est1=res$est1,est2=res$est2,dif=res$dif,pts=pts,n1=res$n1,n2=res$n2)
	} else 
		list(p.value = NA)
}



ancGLOB_pv_C<-function(n1,n2,est=tmean,fr1=.8,fr2=.8,nboot=500,SEED=TRUE,iter=100, cpp=TRUE, ...){
	options(warn=1)
	extras <- list(...)        
	if(SEED)
		set.seed(45)
	bvec=list()
	np1=min(c(n1,n2))+1
	nmax=max(c(n1,n2))
	require(parallel)

	for(i in 1:iter){
		bvec[[i]]=rmul(nmax,p=4)
		if(n1!=n2)
			bvec[[i]][np1:nmax,1:2]=NA
	}
	if(!cpp){
		prm=mclapply(bvec,ancGLOB_sub2,fr1=fr1,fr2=fr2,est=est,SEED=SEED,...)
	} else {
		cat("Calling C++ sub-routine...\n")
		estimatorlist <- c("tmean", "median", "mean", "hd")
		estimator <- mapply(function(tar, cur){ ae <- all.equal(tar, get(cur)); if( is.logical(ae) ) TRUE else FALSE},
			   				tar=list(est), cur=estimatorlist)
		estimator <- seq_along(estimatorlist)[estimator]
		if( sum(estimator) == 0 )
			stop("Only the following estimators are implemented for the C++ sub-routine right now:",
					paste(estimatorlist, collapse=", "))
		
		dot.names <- names(extras)
		if( estimator == 1 && sum("tr" != dot.names) >0 ){
			dot.names <- dot.names[dot.names != "tr"]
			stop("The following arguments are not available for the C++ code: ", 
                  paste(dot.names, collapse=", "))
        }else if( estimator == 1 && !("tr" %in% dot.names) ){
        	tr <- 0.2
        }else if( estimator == 1){
        	tr <- extras[["tr"]]
        }else if( estimator == 4 && sum("q" != dot.names) >0 ){
			dot.names <- dot.names[dot.names != "q"]
			stop("The following arguments are not available for the C++ code: ", 
                  paste(dot.names, collapse=", "))
		}else if( estimator == 4 && !("q" %in% dot.names) ){
        	q <- 0.5
        }else if( estimator == 4){
        	q <- extras[["q"]]
        }else if( estimator !=1 && estimator !=4 && length(dot.names)>0 ){ 
        	warning("The following args are not compatible with the chosen estimator and are ignored: ", 
        			paste(dot.names, collapse=", "))
        	tr <- 0.2
        }else 
        	tr <- 0.2
		bvec <- lapply(lapply(bvec, as.vector), na.omit)
		prm  <- .Call("ancGLOB_sub2_C", bvec, n1, n2, fr1, fr2, 
					  NULL, estimator, if(estimator==4) q else tr)
	}
	options(warn=0)
	prm=as.vector(matl(prm))
	prm=sort(elimna(prm))
	prm
}

