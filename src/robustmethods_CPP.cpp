#include <R.h>
#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

IntegerVector order(NumericVector x);

template <class RandomAccessIterator, class StrictWeakOrdering>
void sort(RandomAccessIterator first, RandomAccessIterator last, StrictWeakOrdering comp);

struct val_order{
		int order;
		double value;
};

bool compare(const val_order & a, const val_order & b){return (a.value<b.value);}


// arma in; arma out
arma::mat sumbatC1(arma::mat X, arma::colvec T) {
    mat y = X.rows(find(T == 1));
    return y;
}

// rcpp in; arma out
arma::mat submatrix(NumericMatrix X, NumericVector T) {
    mat Xmat(X.begin(), X.nrow(), X.ncol(), false);
    colvec tIdx(T.begin(), T.size(), false); 
    mat y = Xmat.rows(find(tIdx == 1));
    return y;
}

// rcpp in; rcpp out
NumericMatrix sumbatC3(NumericMatrix X, LogicalVector condition) { 
    int n=X.nrow(), k=X.ncol();
    NumericMatrix out(sum(condition),k);
    for (int i = 0, j = 0; i < n; i++) {
        if(condition[i]) {
            out(j,_) = X(i,_);
            j = j+1;
        }
    }
    return(out);
}


NumericVector eigenval(arma::mat M, Function eigen) {
    List output=eigen(M);
    
    return output(0);
}

double gvar_C(arma::mat M){
	arma::mat temp=arma::cov(M);
	arma::vec values=arma::eig_sym(temp);
	double output=arma::as_scalar(prod(values));
	return output;
}

double sum(NumericVector v){
	int n=v.size();
	double Sum=0.0;
	for(int i=0;i<n;i++)
		Sum+=v(i);
		
	return Sum;
}

// [[Rcpp::export]]
//NumericVector mgvar_while(NumericVector flag, NumericMatrix m){
RcppExport SEXP mgvar_while(SEXP X, SEXP M){
	NumericVector flag(X);
	NumericMatrix m(M);
 	int nr=m.nrow();
 	NumericVector f=flag;
 	arma::mat m_sub;
 	NumericVector varvec(nr, 0.0);
	while(sum(f)>0.0){
		int ic=0;
	 	NumericVector chk;
	 	IntegerVector remi;
		for(int i=0; i<nr; i++){
 			if(f(i)==1.0) {
 				ic+=1;
	 			NumericVector f2=ifelse(f==1.0, 0.0, 1.0);
 				m_sub=submatrix(m, f2);
 				arma::rowvec m_i=m(i,_);
 				m_sub.insert_rows(m_sub.n_rows, m_i);
	 			chk.push_back(gvar_C(m_sub));
 				remi.push_back(i);
 			}
		}
		IntegerVector sor=order(chk);
 		int k=remi(sor(0));
		varvec(k)=chk(sor(0));
		f(k)=0.0;
		R_CheckUserInterrupt();
	}
	return varvec;
}

IntegerVector order(NumericVector x){
	int n=x.size();
	std::vector<int> output(n);
	std::vector<val_order> index(n);
	for(int i=0;i<x.size();i++){
		index[i].value=x(i);
		index[i].order=i;
	}
	std::sort(index.begin(), index.end(), compare); 
	for(int i=0;i<x.size();i++){
		output[i]=index[i].order;
	}
	return wrap(output);
}

double median(NumericVector x) {
    double output;
  	int n=x.size();
  	std::vector<double> xcopy(n);
    std::copy(x.begin(), x.end(), xcopy.begin());
    std::sort(xcopy.begin(), xcopy.end());
    if(n%2==0) {
	   output=((xcopy[n/2] + xcopy[n/2-1])/2);
    } else {
      output=xcopy[n/2];
   }
    return output;
}

double mad(NumericVector x, double constant=1.4826){
	//NumericVector xcopy=Rcpp::clone(x);
	int n=x.size();
	double center=median(x);
	NumericVector diff=x-center;
	NumericVector diff_abs=ifelse(diff<0, -diff, diff);
	return median(diff_abs)*constant;
}


double mad(arma::vec x, double constant=1.4826){
	NumericVector X = Rcpp::as<Rcpp::NumericVector>(wrap(x));
	double output = mad(X);
	return output;
}

NumericVector outer_pos(NumericVector x, NumericVector y){
	std::vector<double> output;
	int n=x.size();
	double temp;
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			temp=x[j]-x[i];
			if(temp>0){
				output.push_back((y[j]-y[i])/temp);
			}
		}
	} 
	return Rcpp::wrap(output);
}

double hd_C(NumericVector x, NumericVector q){
	int n=x.size();
	double m1=(n+1.0)*q(0), m2=(n+1.0)*(1.0-q(0));
	double output=0.0;
	IntegerVector i=seq_len(n);
	NumericVector i2(n);
	i2=i*1.0;
	NumericVector w=pbeta(i2*1.0/n, m1, m2) - pbeta((i2-1.0)/n, m1, m2);
	NumericVector xcopy=Rcpp::clone<Rcpp::NumericVector>(x);
	std::sort(xcopy.begin(), xcopy.end());
	for(int j=0; j<n; j++){
		output+=(xcopy(j)*w(j));
	}
	return(output);
}

// [[Rcpp::export]]
RcppExport SEXP tshd_C(SEXP X, SEXP Y, SEXP hd){
	NumericVector x(X), y(Y);
	IntegerVector HD(hd);
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector coef(2);
	if(sum(HD)==0){
		coef(1)=median(v1v2);
	}else{
		coef(1)=hd_C(v1v2,wrap(0.5));
	}
	NumericVector res=y-coef(1)*x;
	coef(0)=hd_C(res,wrap(0.5));
	return coef;
}

NumericVector tshd_C(NumericVector x, NumericVector y, IntegerVector HD){
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector coef(2);
	if(sum(HD)==0){
		coef(1)=median(v1v2);
	}else{
		coef(1)=hd_C(v1v2,wrap(0.5));
	}
	NumericVector res=y-coef(1)*x;
	coef(0)=hd_C(res,wrap(0.5));
	return coef;
}


// [[Rcpp::export]]
RcppExport SEXP pbvar_C(SEXP X, SEXP BETA){
	NumericVector x(X);
	NumericVector beta(BETA);
	int n=x.size();
	double Median=median(x);
	NumericVector w=ifelse((x-Median)<=0, -(x-Median), (x-Median));
	double pbvar;
	std::sort(w.begin(), w.end());
	double omega=w[floor((1.0-beta(0))*(n+0.5))-1];
	if(omega>0.0){
		double len=0;
		NumericVector z(n);
		for(int i=0;i<n;i++){
			if((x(i)-Median)/omega<= -1.0)
				z(i)=-1.0;
			if ((x(i)-Median)/omega>=1.0)
				z(i)=1.0;
			if(fabs((x(i)-Median)/omega) < 1.0){
				z(i)=(x(i)-Median)/omega;
				len+=1.0;
			}
		}
		double z_sq_sum=0;
		for(int i=0;i<n;i++){
			z_sq_sum+=pow(z(i), 2.0);
		}
		pbvar=(double)n*pow(omega, 2)*z_sq_sum/pow(len, 2);
		return wrap(pbvar);
	} else
		return wrap(0.0);
}

double pbvar_C(NumericVector x, double beta){
	int n=x.size();
	double Median=median(x);
	NumericVector w=ifelse((x-Median)<=0, -(x-Median), (x-Median));
	double pbvar;
	std::sort(w.begin(), w.end());
	double omega=w[floor((1.0-beta)*(n+0.5))-1];
	if(omega>0.0){
		double len=0;
		NumericVector z(n);
		for(int i=0;i<n;i++){
			if((x(i)-Median)/omega<= -1.0)
				z(i)=-1.0;
			if ((x(i)-Median)/omega>=1.0)
				z(i)=1.0;
			if(fabs((x(i)-Median)/omega) < 1.0){
				z(i)=(x(i)-Median)/omega;
				len+=1.0;
			}
		}
		double z_sq_sum=0;
		for(int i=0;i<n;i++){
			z_sq_sum+=pow(z(i), 2.0);
		}
		pbvar=(double)n*pow(omega, 2)*z_sq_sum/pow(len, 2);
		return pbvar;
	} else
		return 0.0;
}

// [[Rcpp::export]]
RcppExport SEXP tsp1reg_C(SEXP X, SEXP Y, SEXP hd){
	NumericVector x(X), y(Y);
	IntegerVector HD(hd);
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector coef(2);
	coef(1)=median(v1v2);
	if(sum(HD)==0){
		coef(0)=median(y)-coef(1)*median(x);
	}else{
		NumericVector temp1(1), temp2(1);
		temp1(0)=hd_C(y,wrap(0.5));
		temp2(0)=hd_C(x,wrap(0.5));
		coef(0)=temp1(0)-coef(1)*temp2(0);
	}
	NumericVector res=y-coef(1)*x-coef(0);
	return List::create(_["coef"]=coef, _["res"]=res);
}


// [[Rcpp::export]]

RcppExport SEXP stsregp1_C(SEXP X, SEXP Y){
	NumericVector x(X), y(Y);
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector allvar(n_s);
	double temp;
	for(int i=0; i<n_s;i++){
		temp = pbvar_C(y-v1v2(i)*x, 0.2);
		allvar(i) = temp;
		R_CheckUserInterrupt();
	}
	IntegerVector b1_id=order(allvar);
	NumericVector coef(2);
	coef(1)=v1v2[b1_id(0)];
	coef(0)=median(y)-coef[1]*median(x);
	NumericVector res=y-coef(1)*x-coef(0);
	return List::create(_["coef"]=coef, _["res"]=res);
}


NumericVector stsregp1_coef(SEXP X, SEXP Y){
	NumericVector x(X), y(Y);
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector allvar(n_s);
	for(int i=0; i<n_s;i++){
		NumericVector temp(pbvar_C(wrap(y-v1v2(i)*x), wrap(0.2)));
		allvar(i)=temp(0);
		R_CheckUserInterrupt();
	}
	IntegerVector b1_id=order(allvar);
	NumericVector coef(2);
	coef(1)=v1v2[b1_id(0)];
	coef(0)=median(y)-coef[1]*median(x);
	return coef;
}


NumericVector tsp1reg_C(NumericVector x, NumericVector y, IntegerVector HD){
	int n=x.size();
	IntegerVector ord=order(x);
	NumericVector x_ord(n), y_ord(n);
	for(int i=0; i<n; i++){
		x_ord(i)=x[ord(i)];
		y_ord(i)=y[ord(i)];	
	}
	NumericVector v1v2=outer_pos(x_ord, y_ord);
	int n_s=v1v2.size();	
	NumericVector coef(2);
	coef(1)=median(v1v2);
	if(sum(HD)==0){
		coef(0)=median(y)-coef(1)*median(x);
	}else{
		NumericVector temp1(1), temp2(1);
		temp1(0)=hd_C(y,wrap(0.5));
		temp2(0)=hd_C(x,wrap(0.5));
		coef(0)=temp1(0)-coef(1)*temp2(0);
	}
	NumericVector res=y-coef(1)*x-coef(0);
	return coef;
}

RcppExport SEXP stsreg_for(SEXP X, SEXP Y, SEXP IT){
	NumericMatrix x(X);
	NumericVector y(Y);
	IntegerVector it(IT);
	int ncols=x.ncol();
	int nrows=x.nrow();
	NumericVector temp(ncols);
	for(int i=0; i<ncols;i++){
		NumericVector tempcoef=tsp1reg_C(x(_,i), y, wrap(1));
		temp(i)=tempcoef(1);
	}
	arma::colvec res=(Rcpp::as<arma::colvec>(y))
	              -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp));
	double alpha=median((Rcpp::as<NumericVector>(wrap(res))));
//	NumericVector tempold(temp);
	arma::colvec r(nrows);
	NumericVector tempcol(nrows);
	for(int i=0; i<it(0); i++){
		for(int j=0; j<ncols; j++){
			tempcol=x(_,j);
			r=(Rcpp::as<arma::colvec>(y))
			  -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
		      -alpha + temp(j)*(Rcpp::as<arma::colvec>(tempcol));
		    temp(j)=stsregp1_coef(wrap(tempcol), wrap(r))(1);
		}
		alpha=median(Rcpp::as<NumericVector>(wrap(((Rcpp::as<arma::colvec>(y))
					 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))))));
		R_CheckUserInterrupt();
//		NumericVector tempold(temp);
	}
	res=(Rcpp::as<arma::colvec>(y))
	    -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
	    -alpha;
	return List::create(_["alpha"]=alpha, _["beta"]=temp, _["res"]=res);
}


RcppExport SEXP tshdreg_for(SEXP X, SEXP Y, SEXP IT, SEXP TOL){
	NumericMatrix x(X);
	NumericVector y(Y);
	NumericVector tol(TOL);
	IntegerVector it(IT);
	int ncols=x.ncol();
	int nrows=x.nrow();
	IntegerVector HD(1);
	HD(0)=1;
	NumericVector temp(ncols), tempold(ncols);
	for(int i=0; i<ncols;i++){
		temp(i)=tshd_C(x(_,i), y, HD)(1);
		tempold(i)=temp(i);
	}
	arma::colvec res=(Rcpp::as<arma::colvec>(y))
	                 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp));
	double alpha=hd_C((Rcpp::as<NumericVector>(wrap(res))), wrap(0.5));
	std::vector<double> diff(ncols);
	arma::colvec r(nrows);
	NumericVector tempcol(nrows);
	for(int i=0; i<it(0); i++){
		for(int j=0; j<ncols; j++){
			tempcol=x(_,j);
			r=(Rcpp::as<arma::colvec>(y))
			  -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
		      -alpha + temp(j)*(Rcpp::as<arma::colvec>(tempcol));
		    temp(j)=tshd_C(tempcol, Rcpp::as<Rcpp::NumericVector>(wrap(r)), HD)(1);
			diff[j]=temp(j)-tempold(j);
			if(diff[j]<0.0) diff[j]=-1.0*diff[j];
			tempold(j)=temp(j);
		}
		R_CheckUserInterrupt();
		std::sort(diff.begin(), diff.end());
		if(diff[diff.size()-1]<tol(0)) break;
		alpha=hd_C(Rcpp::as<NumericVector>(wrap(((Rcpp::as<arma::colvec>(y))
					 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))))), wrap(0.5));
		NumericVector tempold(clone(temp));
	}
	res=(Rcpp::as<arma::colvec>(y))
	    -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
	    -alpha;
	return List::create(_["alpha"]=alpha, _["beta"]=temp, _["res"]=res);
}




RcppExport SEXP tsreg_for(SEXP X, SEXP Y, SEXP IT, SEXP hd){
	NumericMatrix x(X);
	NumericVector y(Y);
	IntegerVector it(IT);
	int ncols=x.ncol();
	int nrows=x.nrow();
	IntegerVector HD(hd);
	NumericVector temp(ncols), tempold(ncols);
	for(int i=0; i<ncols;i++){
		temp(i)=tsp1reg_C(x(_,i), y, HD)(1);
		tempold(i)=temp(i);
	}
	arma::colvec res=(Rcpp::as<arma::colvec>(y))
	                 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp));
	double alpha;
	if(sum(HD)==0){
		alpha=median(res);
	} else {
		alpha=hd_C((Rcpp::as<NumericVector>(wrap(res))), wrap(0.5));
	}
	std::vector<double> diff(ncols);
	arma::colvec r(nrows);
	NumericVector tempcol(nrows);
	for(int i=0; i<it(0); i++){
		for(int j=0; j<ncols; j++){
			tempcol=x(_, j);
			r=(Rcpp::as<arma::colvec>(y))
			  -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
		      -alpha + temp(j)*(Rcpp::as<arma::colvec>(tempcol));
		    temp(j)=tsp1reg_C(x(_,j), Rcpp::as<Rcpp::NumericVector>(wrap(r)), HD)(1);
			tempold(j)=temp(j);
		}
		R_CheckUserInterrupt();
		if(sum(HD)==0){
			alpha=median(Rcpp::as<NumericVector>(wrap(((Rcpp::as<arma::colvec>(y))
					 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))))));
		} else {
			alpha=hd_C(Rcpp::as<NumericVector>(wrap(((Rcpp::as<arma::colvec>(y))
					 -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))))), wrap(0.5));
		}
	}
	res=(Rcpp::as<arma::colvec>(y))
	    -(Rcpp::as<arma::mat>(x))*(Rcpp::as<arma::colvec>(temp))
	    -alpha;
	return List::create(_["alpha"]=alpha, _["beta"]=temp, _["res"]=res);
}




double sign(double x)
{
	 double output;
	 if (std::isnan(x)) {
	    	output=NA_REAL;
	} else {
    		output=(((x > 0) ? 1 : ((x == 0)? 0 : -1)));
    	}
    return output;
}


NumericVector unidepth1(NumericVector x, NumericVector y){
	int nx=x.size();
	int ny=y.size();
	NumericVector output(ny);
	double pup, pdown;
	for(int i=0; i<ny; i++){
		double temp1=0, temp2=0;
		for(int j=0; j<nx; j++){
			if(y(i)<=x(j)) temp1+=1;
			if(y(i)<x(j)) temp2+=1;
		}
		pup=temp1/((double)nx);
		pdown=1-temp2/((double)nx);
		if(pup<pdown)
		{	output(i)=pup;  }
		else output(i)=pdown;
	} 
	return output;
}

NumericVector subsetVec(NumericVector A, int start, int end) {
  NumericVector B(end-start+1) ;
  std::copy(A.begin() + start, A.begin() + end+1, B.begin() ) ;
  return B;
}


RcppExport SEXP fdepthv2_for(SEXP M, SEXP PTS){
	NumericMatrix pts(PTS);
	NumericMatrix m(M);
	int nrows=m.nrow();
	//IntegerVector mdep_dim(MDEP_DIM);
	//int mdep_nr=mdep_dim(1), mdep_nc=mdep_dim(0); 
	int mdep_nc=(nrows*nrows-nrows)/2, mdep_nr=nrows;
	if(!R_IsNA(pts(0,0))) mdep_nr=pts.nrow();
	NumericMatrix mdep(mdep_nr, mdep_nc);
	NumericVector dis_temp;
	int ic=0;
	NumericVector B(m.ncol()), BB(m.ncol()), A(m.ncol()), temp(m.ncol());
	double bot;
	for(int iall=0; iall<nrows; iall++){
		for(int i=0; i<nrows; i++){
			R_CheckUserInterrupt();
			if(iall<i){
				ic+=1;
				//NumericVector B=m(i,_)-m(iall,_);
				NumericVector dis;
				B=m(i,_)-m(iall,_);
				BB=B*B;
				//NumericVector BB=B*B;
				//double bot=sum(BB);
				bot=sum(BB);
				if(bot!=0.0){
					if(R_IsNA(pts(0,0))){
						for(int j=0; j<nrows; j++){
							//NumericVector A=m(j,_)-m(iall,_);
							A=m(j,_)-m(iall,_);
							//NumericVector temp=(sum(A*B)*B)/bot;
							temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
							//dis_temp.push_back(sum(temp));
						}
					}
					if(!R_IsNA(pts(0,0))){
						for(int j=0; j<nrows; j++){
							//NumericVector A=m(j,_)-m(iall,_);
							A=m(j,_)-m(iall,_);
							//NumericVector temp=(sum(A*B)*B)/bot;
							temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
						}
						for(int j=0; j<pts.nrow(); j++){
							//NumericVector A=pts(j,_)-m(iall,_);
							A=pts(j,_)-m(iall,_);
							//NumericVector temp=(sum(A*B)*B)/bot;
							temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
						}
					}
					//
					// For ic_th projection, store depths of
					// points in mdep[ic,]
					//
					if(R_IsNA(pts(0,0))){
						NumericVector unidepth_out=unidepth1(dis, dis);
						for(int k=0; k<mdep_nr;k++)
							mdep(k,ic-1)=unidepth_out(k);					
					} else if(!R_IsNA(pts(0,0))){
						NumericVector dis_1=subsetVec(dis, 0, nrows-1);
						NumericVector dis_2=subsetVec(dis, nrows, nrows+pts.nrow()-1);
						NumericVector unidepth_out=unidepth1(dis_1, dis_2);
						for(int k=0; k<mdep_nr;k++)
							mdep(k,ic-1)=unidepth_out(k);					
					}
				}
				if(bot==0.0) {
					for(int k=0; k<mdep_nr;k++)
						mdep(k,ic-1)=0.0;					
				}
			}
		}
	}
	return(mdep);
}
	
RcppExport SEXP fdepthv2_for2(SEXP M, SEXP PTS){
	NumericMatrix pts(PTS);
	NumericMatrix m(M);
	int nrows=m.nrow();
	//IntegerVector mdep_dim(MDEP_DIM);
	//int mdep_nr=mdep_dim(1), mdep_nc=mdep_dim(0); 
	int mdep_nc=(nrows*nrows-nrows)/2, mdep_nr=nrows;
	if(!R_IsNA(pts(0,0))) mdep_nr=pts.nrow();
	NumericMatrix mdep(mdep_nr, mdep_nc);
	NumericVector dis_temp;
	int ic=0;
	for(int iall=0; iall<nrows; iall++){
		for(int i=0; i<nrows; i++){
			R_CheckUserInterrupt();
			if(iall<i){
				ic+=1;
				NumericVector B=m(i,_)-m(iall,_);
				NumericVector dis;
				NumericVector BB=B*B;
				double bot=sum(BB);
				if(bot!=0.0){
					if(R_IsNA(pts(0,0))){
						for(int j=0; j<nrows; j++){
							NumericVector A=m(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
							//dis_temp.push_back(sum(temp));
						}
					}
					if(!R_IsNA(pts(0,0))){
						for(int j=0; j<nrows; j++){
							NumericVector A=m(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
						}
						for(int j=0; j<pts.nrow(); j++){
							NumericVector A=pts(j,_)-m(iall,_);
							NumericVector temp=(sum(A*B)*B)/bot;
							dis.push_back(sign(sum(A*B))*pow(sum(pow(temp,2)), 0.5));
						}
					}
					//
					// For ic_th projection, store depths of
					// points in mdep[ic,]
					//
					if(R_IsNA(pts(0,0))){
						NumericVector unidepth_out=unidepth1(dis, dis);
						for(int k=0; k<mdep_nr;k++)
							mdep(k,ic-1)=unidepth_out(k);					
					} else if(!R_IsNA(pts(0,0))){
						NumericVector dis_1=subsetVec(dis, 0, nrows-1);
						NumericVector dis_2=subsetVec(dis, nrows, nrows+pts.nrow()-1);
						NumericVector unidepth_out=unidepth1(dis_1, dis_2);
						for(int k=0; k<mdep_nr;k++)
							mdep(k,ic-1)=unidepth_out(k);					
					}
				}
				if(bot==0.0) {
					for(int k=0; k<mdep_nr;k++)
						mdep(k,ic-1)=0.0;					
				}
			}
		}
	}
	return(mdep);
}

NumericVector idealf(NumericVector x){
	NumericVector xcopy=Rcpp::clone(x), output(2);
	int n=xcopy.size();
	int j=floor((double)n/4.0 + 5.0/12.0)-1;
	std::sort(xcopy.begin(), xcopy.end());
	double g=n/4.0-j-1+5.0/12.0;
	int k=n-j-1;
	output(0)=(1.0-g)*xcopy(j)+g*xcopy(j+1);
	output(1)=(1.0-g)*xcopy(k)+g*xcopy(k-1);
	return output;
}


RcppExport SEXP mahalanobis(SEXP X, SEXP CENTER, SEXP COV){
	BEGIN_RCPP
	arma::mat x = Rcpp::as<arma::mat>(X);
	arma::mat center = Rcpp::as<arma::mat>(CENTER);
	arma::mat covs = Rcpp::as<arma::mat>(COV);
	int n = x.n_rows, ngrp = x.n_cols;
	arma::vec mdis(n);
	arma::mat x_i=mat(ngrp, 1);
	arma::mat xt = x.t();
	arma::mat Y;
	for (int i=0; i<n; i++){
		x_i = xt.col(i) - center;
		Y = arma::solve( covs, x_i );
		mdis(i) = arma::as_scalar(x_i.t() * Y);
	}
	return(wrap(mdis));
	END_RCPP
}

NumericVector mahalanobis(NumericMatrix X, NumericMatrix CENTER, NumericMatrix COV){
	BEGIN_RCPP
	arma::mat x = Rcpp::as<arma::mat>(X);
	arma::mat center = Rcpp::as<arma::mat>(CENTER);
	arma::mat covs = Rcpp::as<arma::mat>(COV);
	int n = x.n_rows, ngrp = x.n_cols;
	arma::vec mdis(n);
	arma::mat x_i=mat(ngrp, 1);
	arma::mat xt = x.t();
	arma::mat Y;
	for (int i=0; i<n; i++){
		x_i = xt.col(i) - center;
		Y = arma::solve( covs, x_i );
		mdis(i) = arma::as_scalar(x_i.t() * Y);
	}
	return(wrap(mdis));
	END_RCPP
}

 arma::mat mahalanobis(arma::mat x, arma::mat center, arma::mat covs){
	int n = x.n_rows, ngrp = x.n_cols;
	arma::vec mdis(n);
	arma::mat x_i=mat(ngrp, 1);
	arma::mat xt = x.t();
	arma::mat Y;
	for (int i=0; i<n; i++){
		x_i = xt.col(i) - center;
		Y = arma::solve( covs, x_i );
		mdis(i) = arma::as_scalar(x_i.t() * Y);
	}
	return(mdis);
}



RcppExport SEXP rmba(SEXP X, SEXP Csteps){
	//Computes the reweighted MBA estimator
	//Original R code provided by David Olive
	//int p = X.ncol(), n = X.nrow();
	arma::mat x = Rcpp::as<arma::mat>(X);
	IntegerVector csteps(Csteps);
	int n = x.n_rows, p = x.n_cols;
	arma::mat covs = arma::cov(x);
	arma::mat mns = mean(x, 0);
	arma::vec md2 = vec(n);
	double medd2 = 0.0;
	for(int i=0; i<csteps(0); i++){
		md2 = mahalanobis(x, mns.t(), covs);
		medd2 = arma::as_scalar(arma::median(md2));
		mns = arma::mean( x.rows(find(md2 <= medd2 ) ));
		covs = arma::cov( x.rows(find(md2 <= medd2 ) ));	
	}
	arma::mat covb = covs(span::all, span::all);
	arma::mat mnb = mns(span::all, span::all);
	double critb = arma::as_scalar(prod(diagvec(chol(covb))));
	arma::mat covv = eye(p, p);
	arma::mat med = median(x, 0);
	md2 = mahalanobis(x, med.t(), covv);
	medd2 = arma::as_scalar(arma::median(md2));
	mns = arma::mean( x.rows(find(md2 <= medd2 ) ));
	covs = arma::cov( x.rows(find(md2 <= medd2 ) ));	
	for(int i=0; i<csteps(0); i++){
		md2 = mahalanobis(x, mns.t(), covs);
		medd2 = arma::as_scalar(arma::median(md2));
		mns = arma::mean( x.rows(find(md2 <= medd2 ) ));
		covs = arma::cov( x.rows(find(md2 <= medd2 ) ));	
	}
	double crit = arma::as_scalar(prod(diagvec(chol(covs))));
	if( crit < critb ){
		critb = crit;
		covb = covs;
		mnb = mns;
	}
	arma::vec rd2 = mahalanobis(x, mnb.t(), covb);
	double Const = (arma::as_scalar(arma::median(rd2)))/(R::qchisq(0.5, (double)p, 1, 0));
	covb = Const * covb;
	double up = R::qchisq(0.975, (double)p, 1, 0);
	arma::mat rmnb = arma::mean( x.rows(find(rd2 <= up)) );
	arma::mat rcovb = arma::cov( x.rows(find(rd2 <= up)) );
	rd2 = mahalanobis(x, rmnb.t(), rcovb);
	Const = (arma::as_scalar(arma::median(rd2)))/(R::qchisq(0.5, (double)p, 1, 0));
	rcovb = Const * rcovb;
	rd2 = mahalanobis(x, rmnb.t(), rcovb);

	rmnb = arma::mean( x.rows(find(rd2 <= up)) );
	rcovb = arma::cov( x.rows(find(rd2 <= up)) );
	rd2 = mahalanobis(x, rmnb.t(), rcovb);	
	Const = (arma::as_scalar(arma::median(rd2)))/(R::qchisq(0.5, (double)p, 1, 0));
	rcovb = Const * rcovb;
	return List::create( _["center"] = rmnb, _["cov"] = rcovb );
}



List rmba(NumericMatrix X, IntegerVector csteps){
	//Computes the reweighted MBA estimator
	//Original R code provided by David Olive
	//int p = X.ncol(), n = X.nrow();
	arma::mat x 	= Rcpp::as<arma::mat>(X);
	int n 			= x.n_rows, p = x.n_cols;
	arma::mat covs 	= arma::cov(x);
	arma::mat mns 	= mean(x, 0);
	arma::vec md2 	= vec(n);
	double medd2 	= 0.0;
	for(int i=0; i<csteps(0); i++){
		md2 		= mahalanobis(x, mns.t(), covs);
		medd2 		= arma::as_scalar(arma::median(md2));
		mns 		= arma::mean( x.rows(find(md2 <= medd2 ) ));
		covs 		= arma::cov( x.rows(find(md2 <= medd2 ) ));	
	}
	arma::mat covb 	= covs(span::all, span::all);
	arma::mat mnb 	= mns(span::all, span::all);
	double critb 	= arma::as_scalar(prod(diagvec(chol(covb))));
	arma::mat covv 	= eye(p, p);
	arma::mat med 	= median(x, 0);
	md2 			= mahalanobis(x, med.t(), covv);
	medd2 			= arma::as_scalar(arma::median(md2));
	mns 			= arma::mean( x.rows(find(md2 <= medd2 ) ));
	covs 			= arma::cov( x.rows(find(md2 <= medd2 ) ));	
	for(int i=0; i<csteps(0); i++){
		md2 		= mahalanobis(x, mns.t(), covs);
		medd2 		= arma::as_scalar(arma::median(md2));
		mns 		= arma::mean( x.rows(find(md2 <= medd2 ) ));
		covs 		= arma::cov( x.rows(find(md2 <= medd2 ) ));	
	}
	double crit = arma::as_scalar(prod(diagvec(chol(covs))));
	if( crit < critb ){
		critb 		= crit;
		covb 		= covs;
		mnb 		= mns;
	}
	arma::vec rd2 	= mahalanobis(x, mnb.t(), covb);
	double Const 	= (arma::as_scalar(arma::median(rd2)))/(R::qchisq(0.5, (double)p, 1, 0));
	covb 			= Const * covb;
	double up 		= R::qchisq(0.975, (double)p, 1, 0);
	arma::mat rmnb 	= arma::mean( x.rows(find(rd2 <= up)) );
	arma::mat rcovb = arma::cov( x.rows(find(rd2 <= up)) );
	rd2 			= mahalanobis(x, rmnb.t(), rcovb);
	Const 			= (arma::as_scalar(arma::median(rd2)))/(R::qchisq(0.5, (double)p, 1, 0));
	rcovb 			= Const * rcovb;
	rd2 			= mahalanobis(x, rmnb.t(), rcovb);
	rmnb 			= arma::mean( x.rows(find(rd2 <= up)) );
	rcovb 			= arma::cov( x.rows(find(rd2 <= up)) );
	rd2 			= mahalanobis(x, rmnb.t(), rcovb);	
	Const 			= (arma::as_scalar(arma::median(rd2)))/(R::qchisq(0.5, (double)p, 1, 0));
	rcovb 			= Const * rcovb;
	return List::create( _["center"] = rmnb, _["cov"] = rcovb );
}


RcppExport SEXP outpro_for(SEXP M, SEXP GVAL, SEXP CENTER, SEXP MM){
	NumericMatrix m(M);
	int nrows=m.nrow();
	NumericVector gval(GVAL), center(CENTER);
	LogicalVector mm(MM);
	NumericVector outid; 
	LogicalVector flag(nrows, false);
	NumericVector B(m.ncol()), A(m.ncol()), temp(m.ncol()), BB(m.ncol()), temp_idealf(2);

	for(int i=0; i<nrows; i++){
		R_CheckUserInterrupt();
		B  = m(i,_)-center;
		BB = B*B;
		NumericVector dis(nrows);
		double bot = sum(BB); 
		R_CheckUserInterrupt();
		if(bot!=0){
			for(int j=0;j<nrows;j++){
				A = m(j,_)-center;
				temp = sum(A*B)*B/bot;
				dis(j) = pow(sum(temp*temp), 0.5);
			}
			temp_idealf=idealf(dis);
			double cu;
			if(!mm(0))
				cu=median(dis)+gval(0)*(temp_idealf(1)-temp_idealf(0));
			else if(mm(0))
				cu=median(dis)+gval(0)*mad(dis);
			for(int k=0; k<nrows; k++){
				if(dis(k)>cu)
					flag(k)=true;
			}
		}
	}
	return(wrap(flag));
}



double depths1( double m, double j ){
	double dep;
	if( m < j ){
		dep = 0.0;
	} else if( j == 1.0 ){
		dep = m;
	} else if( j == 2.0 ){
		dep = m*(m - 1)/2.0;
	} else if( j == 3.0 ){
		dep = m*(m - 1)*(m - 2)/6.0;
	}
	return dep;
}

RcppExport SEXP depth(SEXP U, SEXP V, SEXP M){
	NumericMatrix m(M); 
	NumericVector UU(U), VV(V), alpha;
	IntegerVector fv;
	double u = UU(0), v = VV(0);
	double nums = 0.0, numh = 0.0, sdep = 0.0, p = acos(-1.0);
	double p2 = p*2.0, eps = 0.000001;
	int n = m.nrow(), nt = 0;
	for( int i=0; i < n; i++ ){
		double dv = ( m(i, 0) - u )*( m(i, 0) - u ) + ( m(i, 1) - v )*( m(i, 1) - v );
		dv = pow(dv, 0.5);
		if( dv <= eps ){
			nt += 1;
		}else{
			double xu = (m(i, 0) - u)/dv, yu = (m(i, 1) - v)/dv;
			if( fabs(xu) > fabs(yu) ){
				if( m(i,0) >= u ){
					alpha.push_back( asin(yu) );
					if( alpha( alpha.size()-1 ) < 0.0 ){
						alpha( alpha.size()-1 ) += p2;
					} 
				}else {
					alpha.push_back( p - asin( yu ) );
				}
			} else {
				if( m(i, 1) >= v ){
					alpha.push_back( acos(xu) );
				} else {
					alpha.push_back( p2 - acos(xu) );
				}
			}
			if( alpha( alpha.size()-1 ) >= (p2 - eps) ){
				alpha( alpha.size()-1 ) = 0.0;
			}
		}
	}
	int nn = n - nt;
	if( nn <= 1 ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		numh += nt;
		return(wrap( numh/n ));
	}
	std::sort( alpha.begin(), alpha.end() );
	double angle = alpha(0) - alpha( alpha.size()-1 ) + p2;
	for( int i=1; i<nn; i++ ){
		double dif = (alpha(i) - alpha(i-1));  
		if( angle < dif )
			angle = dif;
	}
	if( angle > (p + eps) ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		numh += nt;
		return(wrap( numh/n ));
	}
	angle = alpha(0);
	int nu = 0;
	for( int i=0; i<nn; i++ ){
		alpha(i) = alpha(i) - angle;
		if( alpha(i) < (p - eps) ){
			nu += 1;
		}
	}
	if( nu >= nn ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		//if( n >=3 ){
		//	sdep = nums/depths1( double(n), 3.0 );
		//}
		numh += nt;
		return(wrap( numh/n ));
	}

//  Mergesort the alpha with their antipodal angles beta and at the same time update I, F(I), and NBAD.
	int ja = 1, jb = 1, nn2 = nn*2, nbad = 0, I = nu, nf = nn;
	double alphk = alpha(0), betak = alpha(nu) - p;
	IntegerVector fv_id;
	for( int i=0; i<nn2; i++ ){
		double add = alphk + eps;
		if ( add < betak ){
			nf += 1;
			if( ja < nn ){
				ja += 1;
				alphk = alpha(ja - 1);
			} else {
				alphk = p2 + 1.0;
			}
		} else {
			I += 1;
			int nn1 = nn + 1;
			if( I == nn1 ){
				I = 1;
				nf = nf - nn;
			}
			fv.push_back(nf);
			fv_id.push_back(I-1);
			int nfi = nf - I;
			nbad += (int)depths1( (double)nfi, 2.0 );
			if( jb < nn ){
				jb += 1;
				if( (jb + nu) <= nn ){
					betak = alpha( jb + nu - 1 ) - p;
				} else {
					betak = alpha( jb + nu - nn - 1 ) + p;
				}
			} else {
				betak = p2 + 1.0;
			}
		}
	}
	// Computation of numh for halfspace depth.
	int gi = 0;
	ja = 1;
	arma::vec fv2( fv.size() );
	for( int i=0; i < fv.size(); i++ ){
		fv2( fv_id(i) ) = fv(i);
	}
	arma::uvec fv_id2 = sort_index( Rcpp::as<arma::vec>(fv_id) );
	angle = alpha(0);
	int dif = nn - fv2(0);
	NumericVector numh_temp(2);
	numh_temp(0) = fv2(0)*1.0;
	numh_temp(1) = dif*1.0;
	numh = min(numh_temp);
	for( int i=1; i < nn; i++ ){
		double aeps = angle + eps;
		if( alpha(i) <= aeps ){
			ja += 1;
		} else {
			gi += ja;
			ja = 1;
			angle = alpha(i);
		}
		int ki = fv2(i) - gi;
		int nnki = nn - ki;
		if( ki >= nnki ){
			ki = nnki;
		}
		if( numh >= ki )
			numh = ki;
	}
	nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
			depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
			depths1( (double)nt, 3 );
	numh = numh + nt;
	return( wrap(numh/n) );
}


double depth(double u, double v, NumericMatrix m){
	NumericVector alpha;
	IntegerVector fv;
	double nums = 0.0, numh = 0.0, sdep = 0.0, p = acos(-1.0);
	double p2 = p*2.0, eps = 0.000001;
	int n = m.nrow(), nt = 0;
	double xu, yu;
	for( int i=0; i < n; i++ ){
		double dv = ( m(i, 0) - u )*( m(i, 0) - u ) + ( m(i, 1) - v )*( m(i, 1) - v );
		dv = pow(dv, 0.5);
		if( dv <= eps ){
			nt += 1;
		}else{
			xu = (m(i, 0) - u)/dv;
			yu = (m(i, 1) - v)/dv;
			if( fabs(xu) > fabs(yu) ){
				if( m(i,0) >= u ){
					alpha.push_back( asin(yu) );
					if( alpha( alpha.size()-1 ) < 0.0 ){
						alpha( alpha.size()-1 ) += p2;
					} 
				}else {
					alpha.push_back( p - asin( yu ) );
				}
			} else {
				if( m(i, 1) >= v ){
					alpha.push_back( acos(xu) );
				} else {
					alpha.push_back( p2 - acos(xu) );
				}
			}
			if( alpha( alpha.size()-1 ) >= (p2 - eps) ){
				alpha( alpha.size()-1 ) = 0.0;
			}
		}
	}
	int nn = n - nt;
	if( nn <= 1 ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		numh += nt;
		return( numh/n );
	}
	std::sort( alpha.begin(), alpha.end() );
	double angle = alpha(0) - alpha( alpha.size()-1 ) + p2;
	for( int i=1; i<nn; i++ ){
		double dif = (alpha(i) - alpha(i-1));  
		if( angle < dif )
			angle = dif;
	}
	if( angle > (p + eps) ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		numh += nt;
		return( numh/n );
	}
	angle = alpha(0);
	int nu = 0;
	for( int i=0; i<nn; i++ ){
		alpha(i) = alpha(i) - angle;
		if( alpha(i) < (p - eps) ){
			nu += 1;
		}
	}
	if( nu >= nn ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		numh += nt;
		return( numh/n );
	}
//  Mergesort the alpha with their antipodal angles beta and at the same time update I, F(I), and NBAD.
	int ja = 1, jb = 1, nn2 = nn*2, nbad = 0, I = nu, nf = nn;
	double alphk = alpha(0), betak = alpha(nu) - p, add;
	IntegerVector fv_id;
	for( int i=0; i<nn2; i++ ){
		add = alphk + eps;
		if ( add < betak ){
			nf += 1;
			if( ja < nn ){
				ja += 1;
				alphk = alpha(ja - 1);
			} else {
				alphk = p2 + 1.0;
			}
		} else {
			I += 1;
			int nn1 = nn + 1;
			if( I == nn1 ){
				I = 1;
				nf = nf - nn;
			}
			fv.push_back(nf);
			fv_id.push_back(I-1);
			int nfi = nf - I;
			nbad += (int)depths1( (double)nfi, 2.0 );
			if( jb < nn ){
				jb += 1;
				if( (jb + nu) <= nn ){
					betak = alpha( jb + nu - 1 ) - p;
				} else {
					betak = alpha( jb + nu - nn - 1 ) + p;
				}
			} else {
				betak = p2 + 1.0;
			}
		}
	}
	// Computation of numh for halfspace depth.
	int gi = 0;
	ja = 1;
	arma::vec fv2( fv.size() );
	for( int i=0; i < fv.size(); i++ ){
		fv2( fv_id(i) ) = fv(i);
	}
	arma::uvec fv_id2 = sort_index( Rcpp::as<arma::vec>(fv_id) );
	angle = alpha(0);
	int dif = nn - fv2(0);
	NumericVector numh_temp(2);
	numh_temp(0) = fv2(0)*1.0;
	numh_temp(1) = dif*1.0;
	numh = min(numh_temp);
	for( int i=1; i < nn; i++ ){
		double aeps = angle + eps;
		if( alpha(i) <= aeps ){
			ja += 1;
		} else {
			gi += ja;
			ja = 1;
			angle = alpha(i);
		}
		int ki = fv2(i) - gi;
		int nnki = nn - ki;
		if( ki >= nnki ){
			ki = nnki;
		}
		if( numh >= ki )
			numh = ki;
	}
	nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
			depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
			depths1( (double)nt, 3 );
	numh = numh + nt;
	return( numh/n );
}


RcppExport SEXP outpro_C(SEXP M, SEXP GVAL, SEXP CENTER, SEXP MM, SEXP COP, SEXP TR, SEXP Q){
	BEGIN_RCPP
	arma::mat m = Rcpp::as<arma::mat>(M);
	NumericVector gval(GVAL), center(CENTER);
	double tr = Rcpp::as<double>(wrap(TR)), q = Rcpp::as<double>(wrap(Q));
	int cop = Rcpp::as<int>(wrap(COP));
	int nv = m.n_rows, nc = m.n_cols;
	LogicalVector mm(MM), flag(nv, false);  

	IntegerVector outid, keepid; 
	if(nc == 1){
		double Med = median( Rcpp::as<Rcpp::NumericVector>(wrap(m)) );
		double Mad = mad ( Rcpp::as<Rcpp::NumericVector>(wrap(m)) );
		NumericVector dis( nv );
		double crit = pow( R::qchisq( 0.975, 1, 1, 0 ), 0.5 );
		for( int i = 0; i < nv; i++ ){
			dis(i) = ( m(i, 0) - Med )*( m(i, 0) - Med )/( Mad*Mad );
			dis(i) = pow( dis(i), 0.5);
			if( dis(i) > crit ){
				outid.push_back(i);
			} else {
				keepid.push_back(i);
			}
		}
	}
	if( nc > 1 ){
		if( R_IsNA(gval(0)) && cop == 1 )
			gval(0) = pow( R::qchisq(0.95, (double)nc, 1, 0), 0.5 );
		if( R_IsNA(gval(0)) && cop != 1 )
			gval(0) = pow( R::qchisq(0.975, (double)nc, 1, 0), 0.5 );
		if( cop == 1 && R_IsNA(center(0)) ){
			if( nc == 2 ){
				NumericVector tempd( nv );
				for( int i = 0; i < nv; i++ ){
					tempd(i) = depth( m(i, 0), m(i, 1), Rcpp::as<Rcpp::NumericMatrix>(wrap(m)) );
				}
				double mdep = *std::max_element( tempd.begin(), tempd.end() );
				arma::uvec flag = find( Rcpp::as<arma::vec>(tempd) == mdep );
				if( flag.n_elem == 1 ){
					center = Rcpp::as<Rcpp::NumericMatrix>( wrap(m) )(flag(0), _);
				} else if( flag.n_elem > 1 ){
					center = Rcpp::as<Rcpp::NumericVector>( wrap( mean( m.rows(flag) ) ));
				}
			}
		}
		if( cop == 3 && R_IsNA(center(0)) )
			center = Rcpp::as<Rcpp::NumericMatrix>(wrap(arma::median( m )))(0,_);
		if( cop == 6 && R_IsNA(center(0)) ){
			center = rmba( Rcpp::as<Rcpp::NumericMatrix>(wrap( m )), 5)(0);
		}
		NumericMatrix m2 = Rcpp::as<Rcpp::NumericMatrix>( wrap(m) );
		NumericVector B(nc), A(nc), BB(nc), temp(nc), temp_idealf(2);
		double bot;
		for( int i=0; i < nv; i++ ){
			B = m2(i,_) - center;
			BB = B*B;
			NumericVector dis(nv);
			bot = sum(BB); 
			if(bot!=0){
				for( int j=0; j<nv; j++ ){
					A = m2(j,_) - center;
					temp = sum( A*B )*B/bot;
					dis(j)=pow(sum( temp*temp ), 0.5);
				}
				temp_idealf=idealf(dis);
				double cu;
				if( !mm(0) )
					cu = median( dis ) + gval( 0 )*( temp_idealf( 1 ) - temp_idealf( 0 ) );
				else if( mm(0) )
					cu = median( dis ) + gval( 0 )*mad( dis );
				for( int k=0; k < nv; k++ ){
					if( dis(k) > cu )
						flag(k) = true;
				}
			}
		}		
		for( int i=0; i < nv; i++ ){
			if(flag(i))
				outid.push_back(i);
			else
				keepid.push_back(i);
		}
	}
	return List::create(_["outid"] = outid, _["keepid"] = keepid);
	END_RCPP
}




IntegerVector outpro_C(NumericMatrix M, NumericVector gval, NumericVector center, 
	          			LogicalVector mm, int cop, double tr, double q){
	arma::mat m = Rcpp::as<arma::mat>(M);
	int nv = m.n_rows, nc = m.n_cols;
	LogicalVector flag(nv, false);  
	IntegerVector outid, keepid; 
	if(nc == 1){
		double Med = median( Rcpp::as<Rcpp::NumericVector>(wrap(m)) );
		double Mad = mad ( Rcpp::as<Rcpp::NumericVector>(wrap(m)) );
		NumericVector dis( nv );
		double crit = pow( R::qchisq( 0.975, 1, 1, 0 ), 0.5 );
		for( int i = 0; i < nv; i++ ){
			dis(i) = ( m(i, 0) - Med )*( m(i, 0) - Med )/( Mad*Mad );
			dis(i) = pow( dis(i), 0.5);
			if( dis(i) > crit ){
				outid.push_back(i);
			} else {
				keepid.push_back(i);
			}
		}
	}
	if( nc > 1 ){
		if( R_IsNA(gval(0)) && cop == 1 )
			gval(0) = pow( R::qchisq(0.95, (double)nc, 1, 0), 0.5 );
		if( R_IsNA(gval(0)) && cop != 1 )
			gval(0) = pow( R::qchisq(0.975, (double)nc, 1, 0), 0.5 );
		if( cop == 1 && R_IsNA(center(0)) ){
			if( nc == 2 ){
				NumericVector tempd( nv );
				for( int i = 0; i < nv; i++ ){
					tempd(i) = depth( m(i, 0), m(i, 1), Rcpp::as<Rcpp::NumericMatrix>(wrap(m)) );
				}
				double mdep = *std::max_element( tempd.begin(), tempd.end() );
				arma::uvec flag = find( Rcpp::as<arma::vec>(tempd) == mdep );
				if( flag.n_elem == 1 ){
					center = Rcpp::as<Rcpp::NumericMatrix>( wrap(m) )(flag(0), _);
				} else if( flag.n_elem > 1 ){
					center = Rcpp::as<Rcpp::NumericVector>( wrap( mean( m.rows(flag) ) ));
				}
			}
		}
		if( cop == 3 && R_IsNA(center(0)) ){
			center = Rcpp::as<Rcpp::NumericMatrix>(wrap(arma::median( m )))(0,_);
		}
		if( cop == 6 && R_IsNA(center(0)) ){
			center = rmba( Rcpp::as<Rcpp::NumericMatrix>(wrap( m )), 
						   Rcpp::as<Rcpp::IntegerVector>(wrap(5)))(0);
		}
		NumericMatrix m2 = Rcpp::as<Rcpp::NumericMatrix>( wrap(m) );
		NumericVector B(nc), A(nc), BB(nc), temp(nc), temp_idealf(2);
		double bot;
		for( int i=0; i < nv; i++ ){
			B = m2(i,_) - center;
			BB = B*B;
			NumericVector dis(nv);
			bot = sum(BB); 
			if(bot!=0){
				for( int j=0; j<nv; j++ ){
					A = m2(j,_) - center;
					temp = sum( A*B )*B/bot;
					dis(j)=pow(sum( temp*temp ), 0.5);
				}
				temp_idealf=idealf(dis);
				double cu;
				if( !mm(0) )
					cu = median( dis ) + gval( 0 )*( temp_idealf( 1 ) - temp_idealf( 0 ) );
				else if( mm(0) )
					cu = median( dis ) + gval( 0 )*mad( dis );
				for( int k=0; k < nv; k++ ){
					if( dis(k) > cu )
						flag(k) = true;
				}
			}
		}		
		for( int i=0; i < nv; i++ ){
			if(!flag(i))
				keepid.push_back(i);
		}
	}
	return(wrap(keepid));
}


RcppExport SEXP skip(SEXP M, SEXP MM, SEXP OUTPRO_COP, SEXP TR, SEXP Q){
	NumericMatrix m(M);
	NumericVector gval(1, NA_REAL), center(m.ncol(), NA_REAL);
	double tr = Rcpp::as<double>(wrap(TR)), q = Rcpp::as<double>(wrap(Q));
	LogicalVector mm(MM);
	int outpro_cop = Rcpp::as<int>(wrap(outpro_cop));
	IntegerVector keepid = outpro_C(m, gval, center, mm, outpro_cop, tr, q);
	arma::mat val = arma::mean((Rcpp::as<arma::mat>(m)).rows(Rcpp::as<arma::uvec>(keepid)));
	return(wrap(val));
}


NumericVector skip(NumericMatrix m, LogicalVector mm, int outpro_cop, double tr, double q){
	NumericVector gval(1, NA_REAL), center(m.ncol(), NA_REAL);
	IntegerVector keepid = outpro_C(m, gval, center, mm, outpro_cop, tr, q);
	arma::mat val = arma::mean((Rcpp::as<arma::mat>(m)).rows(Rcpp::as<arma::uvec>(keepid)));
	return(wrap(val));
}

NumericVector skip(arma::mat M, LogicalVector mm, int outpro_cop, double tr, double q){
	NumericMatrix m = Rcpp::as<Rcpp::NumericMatrix>(wrap(M));
	NumericVector gval(1, NA_REAL), center(m.ncol(), NA_REAL);
	IntegerVector keepid = outpro_C(m, gval, center, mm, outpro_cop, tr, q);
	arma::mat val = arma::mean((Rcpp::as<arma::mat>(m)).rows(Rcpp::as<arma::uvec>(keepid)));
	return(wrap(val));
}

RcppExport SEXP skip_boot(SEXP M, SEXP MM, SEXP OUTPRO_COP, SEXP BOOTID){
	BEGIN_RCPP
	arma::mat m = Rcpp::as<arma::mat>(wrap(M));
	int outpro_cop = Rcpp::as<int>(wrap(outpro_cop));
	List bootid(BOOTID);
	LogicalVector mm(MM);
	int n = m.n_rows, nc = m.n_cols, nboot = bootid.size();
	NumericVector gval(1, NA_REAL), center(nc, NA_REAL);
	double q = 0.5, tr = 0.2;
	NumericVector output(nboot*nc), temp(nc) ;
	arma::uvec tempid=uvec(n);
	for( int i = 0; i < nboot; i++ ){
		R_CheckUserInterrupt();
		tempid = Rcpp::as<arma::uvec>(bootid(i));
		temp = skip( m.rows(tempid), mm, outpro_cop, tr, q );
		for( int j = 0; j < nc; j++ ){
			output(i*nc + j) = temp(j);
		}
	}
	return(wrap(output));
	END_RCPP
}



NumericVector winval(NumericVector x, double tr=0.2){
	NumericVector xcopy = Rcpp::clone(x);
	std::sort(xcopy.begin(), xcopy.end());
	int n = xcopy.size();
	double ibot = floor(0.2*n);
	double itop = n - ibot - 1;
	double xbot = xcopy(ibot);
	double xtop = xcopy(itop);
	NumericVector output(n);
	for( int i=0; i < n; i++ ){
		if( x(i) > xtop )
			output(i) = xtop;
		else if( x(i) < xbot )
			output(i) = xbot;
		else
			output(i) = x(i);
	}
	return output;
}

double winvar(NumericVector x){
	return arma::as_scalar( arma::var( Rcpp::as<arma::vec>(winval(x)) ) );
}



arma::vec near(arma::vec x, double pt, double fr=1.0){
	//determine which values in x are near pt based on fr * mad
	double m = mad(x);
	NumericVector temp(2);
	arma::vec dflag = zeros(x.n_elem);
	if (m == 0.0){
		temp = idealf(Rcpp::as<Rcpp::NumericVector>(wrap(x)));
		m = ( temp(1) - temp(0) )/(R::qnorm(.75, 0.0, 1.0, 1, 0) - 
								   R::qnorm(.25, 0.0, 1.0, 1, 0));
	}
	if (m == 0.0){
		m = sqrt( winvar( Rcpp::as<Rcpp::NumericVector>(wrap(x))) )/0.4129;
	}
	if (m == 0.0)
		throw std::runtime_error("All measures of dispersion are equal to 0!");
	for( int i=0; i < x.n_elem; i++){
		if( fabs(x(i) - pt) <= (fr*m) ){
			dflag(i) = 1;
		}
	}
	return dflag;
}


RcppExport SEXP near(SEXP X, SEXP PT){
	BEGIN_RCPP
	arma::vec x = Rcpp::as<arma::vec>(wrap(X));
	double pt = Rcpp::as<double>(wrap(PT));
	arma::vec output = near(x, pt);

	return(wrap(output));
	END_RCPP
}



double tmean(NumericVector x, double tr=0.2){
	arma::vec xcopy = Rcpp::as<arma::vec>(x);
	if( tr > 0.5 )
		throw std::runtime_error("tr cannot exceed 0.5");
	else if( tr < 0.0)
		throw std::runtime_error("tr cannot be negative");
	else if( tr == 0.5)
		return median(x);
	else {
		int n = xcopy.n_elem, lo = floor(n * tr), hi = n - lo - 1;
		xcopy = arma::sort(xcopy);
		return arma::as_scalar(arma::mean(xcopy.subvec( span(lo, hi) )));
	}
}

RcppExport SEXP tmean(SEXP X, SEXP TR){
	BEGIN_RCPP
	NumericVector x(X), tr(TR);
	double output;
	output  = tmean(x, tr(0));
	return(wrap(output));
	END_RCPP
}

double tmean(arma::vec x, double tr=0.2){
	arma::vec xcopy(x);
	int n = xcopy.n_elem, lo = floor(n * tr), hi = n - lo - 1;
	xcopy = arma::sort(xcopy);
	return arma::as_scalar(arma::mean(xcopy.subvec( span(lo, hi) )));
}


List aov2depth(List x1, List x2, int nboot=500, double estimator=1, double tr=0.2){
	// 2 by K ANOVA independent groups
	//
	//  Main effect Factor A only
	//
	//Strategy: Use depth of zero based on estimated
	//differences for each column  of the K levels of Factor B
	//That is, testing no main effects for Factor A in 
	//a manner that takes into account the pattern of the
	//measures of location rather then simply averaging 
	//across columns.
	int J = x1.size();
	if( J != x2.size() )
		throw std::runtime_error("x1 and x2 should have same number of groups");
	IntegerVector n1(J), n2(J);
	NumericVector est1(J), est2(J);
	arma::mat difb(nboot +1, J), dif(1, J);
	difb.zeros();

	for( int j=0; j < J; j++ ){
		R_CheckUserInterrupt();
		NumericVector tempx1 = x1(j), tempx2 = x2(j);
		n1(j) = tempx1.size();
		n2(j) = tempx2.size();
		if( estimator == 1){
			est1(j) = tmean(tempx1, tr);
			est2(j) = tmean(tempx2, tr);
		} else if ( estimator == 2){
			est1(j) = median(tempx1);
			est2(j) = median(tempx2);
		} else if ( estimator == 3){
			est1(j) = mean(tempx1);
			est2(j) = mean(tempx2);
		}
		dif(0,j)  = est1(j) - est2(j);

		arma::uvec n1boot(n1(j)), n2boot(n2(j));
		IntegerVector n1boot2(n1(j)), n2boot2(n2(j));

		for( int k=0; k < nboot; k++ ){
			n1boot2 = floor(runif(n1(j), 0, n1(j)));
			n2boot2 = floor(runif(n2(j), 0, n2(j)));
			n1boot  = Rcpp::as<arma::uvec>(n1boot2);
			n2boot  = Rcpp::as<arma::uvec>(n2boot2);

			difb(k, j) = tmean((Rcpp::as<arma::vec>(tempx1)).elem(n1boot), tr) -
					     tmean((Rcpp::as<arma::vec>(tempx2)).elem(n2boot), tr);

			if( estimator == 1){
				difb(k, j) = tmean((Rcpp::as<arma::vec>(tempx1)).elem(n1boot), tr) -
						     tmean((Rcpp::as<arma::vec>(tempx2)).elem(n2boot), tr);
			} else if ( estimator == 2){
				difb(k, j) = arma::median((Rcpp::as<arma::vec>(tempx1)).elem(n1boot)) -
					     	 arma::median((Rcpp::as<arma::vec>(tempx2)).elem(n2boot));
			} else if ( estimator == 3){
				difb(k, j) = arma::mean((Rcpp::as<arma::vec>(tempx1)).elem(n1boot)) -
					     	 arma::mean((Rcpp::as<arma::vec>(tempx2)).elem(n2boot));
			}
		}
	}
	arma::mat m1 = arma::cov(difb.submat(0, 0, nboot-2, J-1));
	arma::mat dis = mahalanobis(difb, dif.t(), m1);
	int bplus = nboot + 1;
	arma::uvec sig_vec = find( dis.col(0) >= dis(bplus-1, 0) );
	double sig = sig_vec.n_elem*1.0/(1.0*bplus);
	return List::create(_["pvalue"] = sig, 
						_["est1"] 	= est1, 
						_["est2"] 	= est2,
						_["dif"] 	= dif,
						_["n1"] 	= n1,
						_["n2"] 	= n2);
}


List ancGLOB_sub3_C(arma::mat x1, arma::vec y1, arma::mat x2, arma::vec y2, bool xout, SEXP PTS,
					int nmin, int nboot, double fr1, double fr2, int estimator, double tr=0.2){
	NumericVector pts;
	List g1, g2;
	IntegerVector n3;
	if(xout){
		arma::mat tempx1 = x1;
		arma::mat tempx2 = x2;
		arma::vec tempy1 = y1;
		arma::vec tempy2 = y2;
		
		NumericVector gval(1, NA_REAL), center(tempx1.n_cols, NA_REAL);
		LogicalVector mm(1, FALSE);
		IntegerVector keepid1 = outpro_C(Rcpp::as<Rcpp::NumericMatrix>(wrap(tempx1)), 
										 gval, center, mm, 3, 0.2, 0.5);
		IntegerVector keepid2 = outpro_C(Rcpp::as<Rcpp::NumericMatrix>(wrap(tempx2)), 
										 gval, center, mm, 3, 0.2, 0.5);
		x1 = tempx1.rows(Rcpp::as<arma::uvec>(keepid1));
		y1 = tempy1(Rcpp::as<arma::uvec>(keepid1));
		x2 = tempx2.rows(Rcpp::as<arma::uvec>(keepid2));
		y2 = tempy2(Rcpp::as<arma::uvec>(keepid2));
	}
	int N1 = y1.n_elem, N2 = y2.n_elem;
	if(Rf_isNull(PTS)){
		int test[5];
		arma::uvec xorder1 = sort_index( x1.col(0) ), isub(5);
		x1 = x1.rows(xorder1);
		y1 = y1.elem(xorder1);
		arma::uvec xorder2 = sort_index( x2.col(0) );
		x2 = x2.rows(xorder2);
		y2 = y2.elem(xorder2);
		arma::uvec n1(N1), n2(N1), vecn(N1);
		arma::uvec tempid1, tempid2;
		IntegerVector sub;
		for( int i=0; i < N1; i++){
			tempid1 = find( near(x1.col(0), x1(i,0), fr1) == 1 );
			tempid2 = find( near(x2.col(0), x1(i,0), fr2) == 1 );
			n1(i) = tempid1.n_elem;
			n2(i) = tempid2.n_elem;
			if ( n1(i) < n2(i) ) 
				vecn(i) = n1(i);
			else 
				vecn(i) = n2(i);
			if ( vecn(i) >= nmin ){
				sub.push_back(i);

			}
		}
		isub(0) = min(sub);
		isub(4) = max(sub);
		isub(2) = floor((isub(0) + isub(4))/2);
		isub(1) = floor((isub(0) + isub(2))/2);
		isub(3) = floor((isub(2) + isub(4))/2);

		for( int i=0; i < 5 && sub.size() > 0; i++ ){
			pts.push_back( arma::as_scalar(x1(isub(i), 0)) );
			arma::vec tempy1 = y1.elem( find(near(x1.col(0), pts(i), fr1) == 1) );
			arma::vec tempy2 = y2.elem( find(near(x2.col(0), pts(i), fr2) == 1) );
			g1.push_back(  tempy1 );
			g2.push_back(  tempy2 );
			n3.push_back(tempy1.n_elem);
		}
	}else{
		pts = Rcpp::as<Rcpp::NumericVector>(wrap(PTS));
		if(pts.size() < 2) 
			throw std::runtime_error("Should have at least two points (use the R function ancova)");
		for( int i=0; i < pts.size(); i++ ){
			arma::vec tempy1 = y1.elem( find(near(x1.col(0), pts(i), fr1) == 1) ); 
			arma::vec tempy2 = y2.elem( find(near(x2.col(0), pts(i), fr2) == 1) );
			//Rprintf("pts: \n", pts(i));
			g1.push_back(  tempy1 );
			g2.push_back(  tempy2 );
			n3.push_back(tempy1.n_elem);
		}
	}

	List output;
	if( g1.size() > 0){
		output = aov2depth(g1, g2, nboot, estimator, tr);
	} else {
		output.push_back(NA_REAL);
	}
	return output;
}



RcppExport SEXP ancGLOB_sub2_C(SEXP BVEC, SEXP N1, SEXP N2, SEXP FR1, SEXP FR2, SEXP PTS, SEXP EST, SEXP TR){
	BEGIN_RCPP
	List bvec(BVEC);
	int niter = bvec.size(), n1 = Rcpp::as<int>(wrap(N1)), n2 = Rcpp::as<int>(wrap(N2));
	int ng1, ng2,  estimator = Rcpp::as<int>(wrap(EST));
	if(n1 > n2){
		ng1 = n2;
		ng2 = n1;
	} else {
		ng1 = n1;
		ng2 = n2;
	}
	double fr1 = Rcpp::as<double>(wrap(FR1)), fr2 = Rcpp::as<double>(wrap(FR2)), tr = Rcpp::as<double>(wrap(TR));
	arma::vec y1(ng1), y2(ng2);
	arma::mat x1(ng1, 1), x2(ng2, 1);
	bool xout=FALSE;
	int nmin = 12, nboot=500;
	arma::vec temp(n1*2+n2*2), pval(niter);
	List templist;
	for( int i = 0; i < niter; i++ ){
		temp = Rcpp::as<arma::vec>(wrap(bvec(i)));
		for( int j = 0; j < ng1; j++ ){
		 	x1(j, 0) = temp(j);
		 	y1(j) = temp(j+ng1);
		}
		for( int j = 0; j < ng2; j++ ){
			x2(j, 0) = temp(j + ng1*2);
		 	y2(j) = temp(j + ng1*2 + ng2);
		}
		templist = ancGLOB_sub3_C(x1, y1, x2, y2, xout, PTS, nmin, nboot, fr1, fr2, estimator, tr);
		pval(i) = Rcpp::as<double>(templist(0));
	}
	return(wrap(pval));
	END_RCPP
}




