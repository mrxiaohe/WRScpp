#include <R.h>
#include <stdio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//IntegerVector order(NumericVector x);

template <class RandomAccessIterator, class StrictWeakOrdering>
void sort(RandomAccessIterator first, RandomAccessIterator last, StrictWeakOrdering comp);

struct val_order{
		int order;
		double value;
};

bool compare(const val_order & a, const val_order & b){return (a.value<b.value);}

IntegerVector order(NumericVector x){
	int n = x.size();
	std::vector<int> output(n);
	std::vector<val_order> index(n);
	for(int i = 0;i < n; i++){
		index[i].value = x(i);
		index[i].order = i;
	}
	std::sort(index.begin(), index.end(), compare); 
	for(int i = 0;i < n; i++){
		output[i] = index[i].order;
	}
	return wrap(output);
}

Environment base("package:base");
Function eigen = base["eigen"]; 


NumericVector eigenval(arma::mat M) {
    List output = eigen(Rcpp::wrap(M));    
    return output(0);
}

double gvar_C(arma::mat M){
	return as_scalar( arma::prod(Rcpp::as<arma::vec>(eigenval( arma::cov(M) ))));
}

RcppExport SEXP mgvar_while(SEXP X, SEXP M){
	BEGIN_RCPP
	arma::mat m = Rcpp::as<arma::mat>(wrap(M));
 	int nr = m.n_rows;	
	arma::vec f = Rcpp::as<arma::vec>(wrap(X));
 	arma::mat m_sub;
 	NumericVector varvec(nr, NA_REAL);
 	int sor, k;
	while( sum(f) > 0 ){
		std::vector<double> chk;
		std::vector<int> remi;
		for(int i=0; i<nr; i++){
 			if(f(i) == 1 ) {
 				m_sub = m.rows( find(f == 0) );
 				m_sub.insert_rows( m_sub.n_rows, m.row(i) );
	 			chk.push_back( gvar_C(m_sub) );
 				remi.push_back( i );
 			}
		}
		sor = sort_index( Rcpp::as<arma::vec>(wrap(chk)) )(0);
 		k = remi[sor];
		varvec(k) = chk[sor];
		f(k) = 0;
		R_CheckUserInterrupt();
	}
	return varvec;
	END_RCPP
}

double median( NumericVector x ) {
    int n = x.size();
  	arma::vec xcopy = arma::sort( Rcpp::as<arma::vec>(wrap(x)) );
    if(n%2==0) {
		return (xcopy(n/2) + xcopy(n/2-1))/2 ;
    } else {
    	return xcopy(n/2);
   }
}

double median( arma::vec x ) {
    int n = x.size();
  	arma::vec xsort = arma::sort( x );
    if(n%2==0) {
	   	return (xsort(n/2) + xsort(n/2-1))/2 ;
    } else {
      	return xsort(n/2);
   }
}

double mad( NumericVector x, double constant=1.4826 ){
	int n = x.size();
	double center = median(x);
	NumericVector diff = x - center;
	NumericVector diff_abs = ifelse(diff<0, -diff, diff);
	return median(diff_abs)*constant;
}

double mad( arma::vec x, double constant=1.4826 ){
	double center = arma::as_scalar(median( x ));
	arma::vec diff = arma::vec( x.n_elem ); 
	for( int i = 0; i < x.n_elem; i++ )
		diff( i ) = fabs( x( i ) - center );
	return arma::as_scalar( median( diff ) ) * constant;
}

NumericVector outer_pos( NumericVector x, NumericVector y ){
	std::vector<double> output;
	int n = x.size();
	double temp;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			temp = x[j] - x[i];
			if(temp > 0){
				output.push_back(( y[j] - y[i]) / temp);
			}
		}
		R_CheckUserInterrupt();
	} 
	return Rcpp::wrap(output);
}

NumericVector outer_pos( arma::vec x, arma::vec y ){
	std::vector<double> output;
	int n = x.size();
	double temp;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			temp = x[j] - x[i];
			if(temp > 0){
				output.push_back(( y[j] - y[i]) / temp);
			}
		}
		R_CheckUserInterrupt();
	} 
	return Rcpp::wrap(output);
}

double hd_C( NumericVector x, double q ){
	int n = x.size();
	double m1 = ( n + 1.0 ) * q, m2=( n + 1.0 ) * ( 1.0 - q );
	double output = 0.0;
	arma::vec xsort = arma::sort( Rcpp::as<arma::vec>(x) );
	for( int i = 1; i <= n; i++ ){
		output += (R::pbeta(i*1.0/n, m1, m2, 1, 0) - 
			       R::pbeta((i-1.0)/n, m1, m2, 1, 0)) * 
		           xsort(i-1);
	}
	return(output);
}

double hd_C( arma::vec x, double q ){
	int n = x.n_elem;
	double m1 = ( n + 1.0 ) * q, m2=( n + 1.0 ) * ( 1.0 - q );
	double output = 0.0;
	arma::vec xsort = arma::sort(x);
	for( int i = 1; i <= n; i++ ){
		output += (R::pbeta(i*1.0/n, m1, m2, 1, 0) - 
			       R::pbeta((i-1.0)/n, m1, m2, 1, 0)) * 
		           xsort(i-1);
	}
	return(output);
}

NumericVector ts_proc( arma::vec x, arma::vec y ){
	arma::uvec ind = arma::sort_index( x );
	return outer_pos( x.elem( ind ), y.elem( ind ) );
}

NumericVector tshd_C( arma::vec x, arma::vec y, bool HD ){
	NumericVector coef(2);
	NumericVector v1v2 = ts_proc( x, y );
	int n_s = v1v2.size();	
	if( !HD ){
		coef( 1 ) = median( v1v2 );
	} else {
		coef( 1 ) = hd_C( v1v2, 0.5 );
	}
	coef( 0 ) = hd_C( y - coef(1) * x, 0.5 );
	return coef;
}


NumericVector tshd_C( NumericVector x, NumericVector y, bool HD ){
	NumericVector coef(2);
	NumericVector v1v2 = ts_proc( Rcpp::as<arma::vec>(x), Rcpp::as<arma::vec>(y) );
	int n_s = v1v2.size();	
	if(!HD){
		coef(1) = median(v1v2);
	}else{
		coef(1) = hd_C(v1v2, 0.5 );
	}
	NumericVector res = y - coef(1)*x;
	coef(0) = hd_C(res, 0.5 );
	return coef;
}

RcppExport SEXP tshd_C( SEXP X, SEXP Y, SEXP hd ){
	BEGIN_RCPP
	return tshd_C( Rcpp::as<arma::vec>( wrap(X) ), 
				   Rcpp::as<arma::vec>( wrap(Y) ), 
				   Rcpp::as<bool>( wrap(hd) ) );
	END_RCPP
}

double pbvar_C( NumericVector x, double beta ){
	int n = x.size();
	double Median = median(x);
	NumericVector w = ifelse( (x-Median)<=0, -(x-Median), (x-Median) );
	double pbvar;
	std::sort( w.begin(), w.end() );
	double omega = w[ floor( (1.0 - beta) * n + 0.5 ) - 1 ];

	if( omega > 0.0 ){
		int len = 0;
		double z_sq_sum = 0.0, div;
		NumericVector z(n);
		for( int i = 0; i < n; i++ ){
			div = ( x(i)-Median )/omega; 
			if( div <= -1.0 )
				z(i) = -1.0;
			else if( div >= 1.0 )
				z(i) = 1.0;
			else {
				z(i) = div;
				len += 1;
			}
			z_sq_sum += pow( z(i), 2.0 );
		}
		pbvar = (double)n*pow(omega, 2)*z_sq_sum/pow((double)len, 2);
		return pbvar;
	} else
		return 0.0;
}

double pbvar_C( arma::vec x, double beta ){
	int n = x.n_elem;
	double Median = arma::as_scalar( arma::median(x) );
	NumericVector w(n);
	for( int i = 0; i < n; i++ ){
		w( i ) = fabs( x(i) - Median );
	}
	double pbvar;
	std::sort( w.begin(), w.end() );
	double omega = w[ floor( (1.0 - beta) * n + 0.5 ) - 1 ];

	if( omega > 0.0 ){
		int len = 0;
		double z_sq_sum = 0.0, div;
		NumericVector z(n);
		for( int i = 0; i < n; i++ ){
			div = ( x(i) - Median )/omega; 
			if( div <= -1.0 )
				z(i) = -1.0;
			else if( div >= 1.0 )
				z(i) = 1.0;
			else {
				z(i) = div;
				len += 1;
			}
			z_sq_sum += pow( z(i), 2.0 );
		}
		pbvar = (double)n*pow(omega, 2)*z_sq_sum/pow((double)len, 2);
		return pbvar;
	} else
		return 0.0;
}


RcppExport SEXP pbvar_C( SEXP X, SEXP BETA ){
	return(wrap(pbvar_C( Rcpp::as<Rcpp::NumericVector>( wrap(X) ),
			             Rcpp::as<double>( wrap(BETA) ))));
}


NumericVector tsp1reg_C( arma::vec x, arma::vec y, bool HD ){
	NumericVector coef(2);
	NumericVector v1v2 = ts_proc( x, y );

	int n_s = v1v2.size();	
	coef(1) = median( v1v2 );
	if(!HD)
		coef(0) = median(y) - coef(1) * median(x);
	else
		coef(0)  = hd_C( y, 0.5 ) - coef( 1 ) * hd_C( x, 0.5 );
	return coef;
}

RcppExport SEXP tsp1reg_C( SEXP X, SEXP Y, SEXP HD ){
	BEGIN_RCPP
	return(wrap(tsp1reg_C( Rcpp::as<arma::vec>( wrap(X) ),
						   Rcpp::as<arma::vec>( wrap(Y) ),
			               Rcpp::as<bool>( wrap(HD) ))));
	END_RCPP
}


List stsregp1_C( arma::vec x, arma::vec y ){
	NumericVector coef(2);
	NumericVector v1v2 = ts_proc( x, y );
	int n_s = v1v2.size();	
	NumericVector allvar(n_s);
	for( int i = 0; i < n_s; i++ ){
		allvar(i) = pbvar_C( y - v1v2(i) * x, 0.2 );
		R_CheckUserInterrupt();
	}
	int b1_id = arma::sort_index(Rcpp::as<arma::vec>(allvar))(0);
	coef(1) = v1v2( b1_id );
	coef(0) = arma::as_scalar( arma::median(y) - coef(1) * arma::median(x) );
	return List::create(_["coef"] = coef, 
						_["res"] = y - coef( 1 ) * x - coef( 0 ));
}


List stsreg_for( arma::mat x, arma::colvec y, int it){
	int ncols = x.n_cols;
	int nrows = x.n_rows;
	NumericVector coef( ncols + 1 ); 
	arma::colvec temp = arma::colvec( ncols );
	for( int i = 0; i < ncols; i++ ){
		temp(i) = tsp1reg_C( arma::conv_to<arma::colvec >::from( x.col(i) ), 
							 arma::conv_to<arma::vec >::from(y), false )(1);
	}
	arma::colvec res = y - x * temp ;
	coef( 0 ) = arma::as_scalar( arma::median( res ) );
	NumericVector temp_coef(2);
	for( int i = 0; i < it; i++ ){
		for( int j = 0; j < ncols; j++ ){
		    temp_coef = stsregp1_C( x.col(j), y - x * temp - coef(0) + temp(j) * x.col(j) )(0);
		    temp( j ) = temp_coef(1);
		}
		coef( 0 ) = arma::as_scalar( arma::median( y - x * temp ) );
		R_CheckUserInterrupt();
	}
	res = y - x * temp - coef( 0 );

	for( int i = 0; i < ncols; i++ )
		coef( i + 1 ) = temp( i ) ;
	return List::create(_["coef"] = coef, _["res"]=res);
}



List stsreg_C( arma::mat x, arma::vec y, int it ){
	if( x.n_cols == 1 ){
		return stsregp1_C( x.col(0), y );
	} else {
		return stsreg_for( x, arma::conv_to<arma::colvec >::from(y), it );
	}
}


RcppExport SEXP stsreg_C( SEXP X, SEXP Y, SEXP IT){
	BEGIN_RCPP
	return stsreg_C( Rcpp::as<arma::mat>(wrap(X)), 
					 Rcpp::as<arma::colvec>(wrap(Y)),
					 Rcpp::as<int>(wrap(IT)));
	END_RCPP
}




List tshdreg_for( arma::mat x, arma::colvec y, int it, double tol, bool HD){
	int ncols = x.n_cols;
	int nrows = x.n_rows;
	arma::colvec temp = arma::colvec( ncols ), tempold = arma::colvec( ncols );
	NumericVector coef( ncols + 1 );

	for( int i = 0; i < ncols; i++ ){
		temp( i ) = tshd_C( x.col(i), y, HD )(1);
		tempold( i ) = temp( i );
	}
	arma::colvec res = y - x * temp;
	coef(0) = hd_C( res,  0.5 );
	arma::colvec diff = arma::colvec( ncols );
	arma::colvec r(nrows);
	NumericVector tempcol(nrows);
	
	for( int i = 0; i < it; i++ ){
		for( int j=0; j < ncols; j++ ){
			r = y - x * temp - coef( 0 ) + temp( j )*x.col( j );
		    temp(j) = tshd_C( x.col( j ), y - x * temp - coef( 0 ) + temp( j )*x.col( j ), HD)(1);
			
			diff(j) = temp( j ) - tempold( j );
			if( diff(j) < 0 )
				diff(j) = fabs( diff(j) );
			tempold( j ) = temp( j );
		}
		R_CheckUserInterrupt();
		if( arma::max( diff ) < tol ) 
			break;
		coef( 0 ) = hd_C( y - x*temp, 0.5 );					
		for( int j = 0; j < ncols; j++ ){
			tempold( j ) = temp( j );
		}
	}
	res = y - x*temp - coef(0);
	for( int i = 0; i < ncols; i++ )
		coef( i + 1 ) = temp( i ) ;
	return List::create( _["coef"] = coef, _["res"] = res );
}



List tshdreg_C( arma::mat x, arma::vec y, int it, double tol, bool HD ){
	if( x.n_cols == 1 ){
		NumericVector res( x.n_rows ), coef = tshd_C( x.col(0), y, HD );
		for( int i = 0; i < x.n_rows; i++ )
			res(i) = y(i) - coef( 1 ) * x(i, 0) - coef(0);
		return List::create( _["coef"] = coef, _["res"] = res );
	} else {
		return tshdreg_for( x, arma::conv_to<arma::colvec >::from(y), it, tol, HD);
	}
}


RcppExport SEXP tshdreg_C( SEXP X, SEXP Y, SEXP IT, SEXP TOL, SEXP HD ){
	BEGIN_RCPP

	return tshdreg_C( Rcpp::as<arma::mat>(wrap(X)), 
					  Rcpp::as<arma::colvec>(wrap(Y)),
					  Rcpp::as<int>(wrap(IT)),
					  Rcpp::as<double>(wrap(TOL)),
					  Rcpp::as<bool>(wrap(HD)));
	END_RCPP
}



List tsreg_C( arma::mat x, arma::colvec y, int it, bool HD){
	int ncols = x.n_cols;
	int nrows = x.n_rows;
	arma::colvec temp = arma::colvec( ncols ), tempold = arma::colvec( ncols );
	NumericVector coef( ncols + 1 );

	for( int i = 0; i < ncols; i++)
		temp(i) = tsp1reg_C( x.col(i), y, false )(1);
	
	arma::colvec res= y - x * temp;
	if(HD)
		coef( 0 ) = median( res );
	else {
		coef( 0 ) = hd_C( res,  0.5 );
	}
	
	arma::colvec r(nrows);
	NumericVector tempcol(nrows);
	for( int i = 0; i < it; i++ ){
		for( int j = 0; j<ncols; j++){
		    temp( j ) = tsp1reg_C( x.col(j), y - x*temp - coef(0) + temp(j)*x.col(j), false)(1);
		}
		R_CheckUserInterrupt();
		if( !HD ) 
			coef( 0 ) = arma::as_scalar( arma::median( y - x*temp ) );
		else 
			coef( 0 ) = hd_C( y - x*temp, 0.5 );
	}
	res = y = x*temp - coef(0); 

	for( int i = 0; i < ncols; i++ )
		coef( i + 1 ) = temp( i ) ;
	return List::create( _["coef"] = coef, _["res"] = res );
}


RcppExport SEXP tsreg_C( SEXP X, SEXP Y, SEXP IT, SEXP HD ){
	BEGIN_RCPP
	return tsreg_C( Rcpp::as<arma::mat>(wrap(X)), 
					Rcpp::as<arma::colvec>(wrap(Y)),
				    Rcpp::as<int>(wrap(IT)),
					Rcpp::as<bool>(wrap(HD)));
	END_RCPP
}


NumericVector unidepth_C(NumericVector x, NumericVector y){
	int nx = x.size();
	int ny = y.size();
	NumericVector output(ny);
	arma::vec p = arma::vec( 2 );
	
	double temp1, temp2;
	for( int i = 0; i < ny; i++ ){
		temp1 = 0;
		temp2 = 0;
		for( int j = 0; j < nx; j++ ){
			if( y(i) <= x(j) )
				temp1 += 1;
			if( y(i) < x(j) )  
				temp2 += 1;
		}
		p(1) = temp1/( (double)nx );
		p(0) = 1 - temp2/( (double)nx );
		output( i ) = arma::as_scalar( arma::min( p ) );
	} 
	return output;
}

NumericVector unidepth_C(arma::vec x, arma::vec y){
	int nx = x.n_elem;
	int ny = y.n_elem;
	NumericVector output(ny);
	arma::vec p = arma::vec( 2 );
	
	double temp1, temp2;
	for( int i = 0; i < ny; i++ ){
		temp1 = 0;
		temp2 = 0;
		for( int j = 0; j < nx; j++ ){
			if( y(i) <= x(j) )
				temp1 += 1;
			if( y(i) < x(j) )  
				temp2 += 1;
		}
		p(1) = temp1/( (double)nx );
		p(0) = 1 - temp2/( (double)nx );
		output( i ) = arma::as_scalar( arma::min( p ) );
	} 
	return output;
}


RcppExport SEXP unidepth_C( SEXP X, SEXP Y ){
	BEGIN_RCPP
	return unidepth_C( Rcpp::as<Rcpp::NumericVector>(wrap(X)), 
					   Rcpp::as<Rcpp::NumericVector>(wrap(Y)) );
	END_RCPP
}

NumericVector subsetVec(NumericVector A, int start, int end) {
	NumericVector B( end - start + 1 ) ;
	std::copy( A.begin() + start, A.begin() + end+1, B.begin() ) ;
	return B;
}


double sign(double x){
    return (((x > 0) ? 1 : ((x == 0)? 0 : -1)));
}


NumericVector fdepthv2_C( NumericMatrix m, int nrows_orig ){
	int nrows = nrows_orig, dis_len = m.nrow(), ic = 0;
	int mdep_nc = ( nrows * nrows - nrows )/2, mdep_nr = nrows;

	if( m.nrow() > nrows_orig ){
		mdep_nr = m.nrow() - nrows_orig;	
	}

	NumericMatrix mdep( mdep_nr , mdep_nc);
	NumericVector dis_temp, dis1_1, dis_2, unidepth_out( mdep_nr );
	NumericVector B(m.ncol()), BB(m.ncol()), A(m.ncol()), temp(m.ncol());
	arma::vec dis = arma::vec( dis_len );
	double bot;

	for( int iall = 0; iall < nrows; iall++ ){
		R_CheckUserInterrupt();
		for( int i = 0; i < nrows; i++ ){
			if( iall < i ){
				ic += 1;
				B = m(i,_) - m(iall,_);
				BB = B * B;
				bot = sum( BB );
				if( bot != 0.0 ){
					if( nrows == m.nrow() ){
						for( int j = 0; j < nrows; j++){
							A = m(j,_) - m(iall,_);
							temp = ( sum( A * B ) * B )/bot;
							dis( j ) = sign( sum( A*B ) )*pow( sum( pow(temp,2) ), 0.5 );
						}
						unidepth_out = unidepth_C( dis, dis );
						for( int k = 0; k < mdep_nr; k++)
							mdep( k, ic - 1 ) = unidepth_out( k );	
					} else {
						for( int j = 0; j < m.nrow(); j++ ){
							A = m( j,_ ) - m( iall, _ );
							temp = ( sum( A * B ) * B )/bot;
							dis( j ) = sign( sum( A*B ) )*pow( sum( pow(temp,2) ), 0.5 );
						}
						unidepth_out = unidepth_C( dis.subvec( 0, nrows-1 ), 
							  				  	   dis.subvec( nrows, nrows + mdep_nr -1 ));
						for( int k = 0; k < mdep_nr; k++ )
							mdep( k, ic - 1 ) = unidepth_out( k );
					}
				} else {
					for( int k = 0; k < mdep_nr; k++)
						mdep( k, ic-1 ) = 0.0;					
				}
			}
		}
	}
	return mdep ;
}
	

RcppExport SEXP fdepthv2_C( SEXP M, SEXP NO ){
	BEGIN_RCPP	
	return fdepthv2_C( Rcpp::as<Rcpp::NumericMatrix>( wrap(M) ),
			 		   Rcpp::as<int>( wrap(NO) )); 
	END_RCPP
}

NumericVector idealf(NumericVector x){
	NumericVector xcopy = Rcpp::clone(x), output(2);
	int n = xcopy.size();
	int j = floor( (double)n/4.0 + 5.0/12.0 ) - 1;
	std::sort( xcopy.begin(), xcopy.end() );
	double g = n/4.0 - j - 1 + 5.0/12.0;
	int k = n - j - 1;
	output(0) = (1.0 - g) * xcopy( j ) + g * xcopy( j + 1 );
	output(1) = (1.0 - g) * xcopy( k ) + g * xcopy( k - 1 );
	return output;
}



 arma::mat mahalanobis(arma::mat x, arma::mat center, arma::mat covs){
	int n = x.n_rows, ngrp = x.n_cols;
	arma::vec mdis(n);
	arma::mat x_i = mat(ngrp, 1);
	arma::mat xt = x.t();
	arma::mat Y;
	for ( int i = 0; i < n; i++ ){
		x_i = xt.col(i) - center;
		Y = arma::solve( covs, x_i );
		mdis(i) = arma::as_scalar(x_i.t() * Y);
	}
	return(mdis);
}

arma::mat mahalanobis(NumericMatrix X, NumericMatrix CENTER, NumericMatrix COV){
	return mahalanobis( Rcpp::as<arma::mat>( wrap(X) ),
						Rcpp::as<arma::mat>( wrap(CENTER) ),
						Rcpp::as<arma::mat>( wrap(COV) ));
}


List rmba(arma::mat x, int csteps){
	//Computes the reweighted MBA estimator
	//Original R code provided by David Olive
	//int p = X.ncol(), n = X.nrow();
	int n 			= x.n_rows, p = x.n_cols;
	arma::mat covs 	= arma::cov(x);
	arma::mat mns 	= mean(x, 0);
	arma::vec md2 	= vec(n);
	double medd2 	= 0.0;

	for(int i = 0; i < csteps; i++){
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

	for( int i = 0; i < csteps; i++ ){
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

NumericVector idealf(arma::vec x){
	NumericVector output(2);
	int n = x.n_elem;
	int j = floor( (double)n/4.0 + 5.0/12.0 ) - 1;
	arma::vec xsort = arma::sort( x );
	double g = n/4.0 - j - 1 + 5.0/12.0;
	int k = n - j - 1;
	output(0) = (1.0 - g) * xsort( j ) + g * xsort( j + 1 );
	output(1) = (1.0 - g) * xsort( k ) + g * xsort( k - 1 );
	return output;
}

LogicalVector outpro_for(NumericMatrix m, NumericVector gval, NumericVector center, bool mm ){
	int nrows = m.nrow();
	NumericVector outid, B(m.ncol()), A(m.ncol()), temp(m.ncol()), BB(m.ncol()), temp_idealf(2);
	LogicalVector flag(nrows, false);
	arma::vec dis = arma::vec( nrows );
	double bot;

	double cu;
	for( int i = 0; i < nrows; i++ ){
		R_CheckUserInterrupt();
		B  = m(i,_) - center;
		BB = B * B;
		bot = sum(BB); 

		if( bot != 0){
			for(int j = 0; j < nrows; j++ ){
				A = m(j,_) - center;
				temp = sum( A*B)*B/bot;
				dis(j) = pow( sum( temp * temp ), 0.5 );
			}
			temp_idealf = idealf(dis);
			if( !mm )
				cu = arma::as_scalar(arma::median( dis )) + gval(0) * ( temp_idealf(1) - temp_idealf(0) );
			else
				cu = arma::as_scalar(arma::median( dis )) + gval(0) * mad(dis);
			for( int k = 0; k < nrows; k++ ){
				if( dis(k) > cu)
					flag(k)=true;
			}
		}
	}
	return(wrap(flag));
}


RcppExport SEXP outpro_for( SEXP M, SEXP GVAL, SEXP CENTER, SEXP MM ){
	BEGIN_RCPP	
	return outpro_for( Rcpp::as<Rcpp::NumericMatrix>(wrap(M)),
					 Rcpp::as<Rcpp::NumericVector>(wrap(GVAL)),
					 Rcpp::as<Rcpp::NumericVector>(wrap(CENTER)),
					 Rcpp::as<bool>(wrap(MM)));
	END_RCPP
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


double depth( double u, double v, NumericMatrix m ){
	std::vector<double> alpha;
	std::vector<int> fv;
	double nums = 0.0, numh = 0.0, sdep = 0.0, p = acos(-1.0), p2 = p*2.0, eps = 0.000001, dv;
	int n = m.nrow(), nt = 0;

	for( int i=0; i < n; i++ ){
		dv = pow( ( m(i, 0) - u )*( m(i, 0) - u ) + ( m(i, 1) - v )*( m(i, 1) - v ), 0.5 );
		if( dv <= eps ){
			nt += 1;
		}else{
			double xu = (m(i, 0) - u)/dv, yu = (m(i, 1) - v)/dv;
			if( fabs(xu) > fabs(yu) ){
				if( m(i,0) >= u ){
					alpha.push_back( asin(yu) );
					if( alpha[ alpha.size()-1 ] < 0.0 ){
						alpha[ alpha.size()-1 ] += p2;
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
			if( alpha[ alpha.size()-1 ] >= (p2 - eps) ){
				alpha[ alpha.size()-1 ] = 0.0;
			}
		}
	}
	int nn = n - nt;
	if( nn <= 1 ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		numh += nt;
		return numh/n ;
	}
	std::sort( alpha.begin(), alpha.end() );
	double angle = alpha[0] - alpha[ alpha.size()-1 ] + p2;
	for( int i=1; i<nn; i++ ){
		double dif = (alpha[i] - alpha[i-1]);  
		if( angle < dif )
			angle = dif;
	}
	if( angle > (p + eps) ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		numh += nt;
		return numh/n ;
	}
	angle = alpha[0];
	int nu = 0;
	for( int i=0; i<nn; i++ ){
		alpha[i] = alpha[i] - angle;
		if( alpha[i] < (p - eps) ){
			nu += 1;
		}
	}
	if( nu >= nn ){
		nums += depths1( (double)nt, 1.0 )*depths1( (double)nn, 2.0 ) +
				depths1( (double)nt, 2.0 )*depths1( (double)nn, 1.0 ) +
				depths1( (double)nt, 3 );
		numh += nt;
		return numh/n ;
	}

//  Mergesort the alpha with their antipodal angles beta and at the same time update I, F(I), and NBAD.
	int ja = 1, jb = 1, nn2 = nn*2, nbad = 0, I = nu, nf = nn, nn1;
	double alphk = alpha[0], betak = alpha[nu] - p, add;
	std::vector<int> fv_id;

	for( int i=0; i<nn2; i++ ){
		add = alphk + eps;
		if ( add < betak ){
			nf += 1;
			if( ja < nn ){
				ja += 1;
				alphk = alpha[ja - 1];
			} else {
				alphk = p2 + 1.0;
			}
		} else {
			I += 1;
			nn1 = nn + 1;
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
					betak = alpha[ jb + nu - 1 ] - p;
				} else {
					betak = alpha[ jb + nu - nn - 1 ] + p;
				}
			} else {
				betak = p2 + 1.0;
			}
		}
	}
	int gi = 0;
	ja = 1;
	arma::vec fv2( fv.size() );

	for( int i=0; i < fv.size(); i++ ){
		fv2( fv_id[i] ) = fv[i];
	}
	//arma::uvec fv_id2 = sort_index( Rcpp::as<arma::vec>(fv_id) );
	angle = alpha[0];
	int dif = nn - fv2(0);
	NumericVector numh_temp(2);
	numh_temp(0) = fv2(0)*1.0;
	numh_temp(1) = dif*1.0;
	numh = min(numh_temp);
	for( int i=1; i < nn; i++ ){
		double aeps = angle + eps;
		if( alpha[i] <= aeps ){
			ja += 1;
		} else {
			gi += ja;
			ja = 1;
			angle = alpha[i];
		}
		int ki = fv2[i] - gi;
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
	return numh/n ;
}


RcppExport SEXP depth(SEXP U, SEXP V, SEXP M){
	BEGIN_RCPP

	return(wrap(depth( Rcpp::as<double>( wrap( U ) ),
			           Rcpp::as<double>( wrap( V ) ),
				       Rcpp::as<NumericMatrix>( wrap( M ) ))));
	END_RCPP
}


IntegerVector outpro_C(NumericMatrix M, NumericVector center, bool mm, 
	                   int cop, double tr, double q, double gval=NA_REAL){
	arma::mat m = Rcpp::as<arma::mat>(M);
	int nv = m.n_rows, nc = m.n_cols;
	LogicalVector flag(nv, false);  
	std::vector<int> outid, keepid; 
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
		if( R_IsNA( gval ) && cop == 1 )
			gval = pow( R::qchisq(0.95, (double)nc, 1, 0), 0.5 );
		if( R_IsNA( gval ) && cop != 1 )
			gval = pow( R::qchisq(0.975, (double)nc, 1, 0), 0.5 );
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
			center = rmba( m , 5)(0);
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
				if( !mm )
					cu = median( dis ) + gval*( temp_idealf( 1 ) - temp_idealf( 0 ) );
				else 
					cu = median( dis ) + gval*mad( dis );
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

arma::uvec outpro_C(arma::mat m, NumericVector center, bool mm, 
	                   int cop, double tr, double q, double gval=NA_REAL){

	int nv = m.n_rows, nc = m.n_cols;
	arma::uvec flag = zeros<uvec>( nv );
	arma::uvec keepid;
	if(nc == 1){
		double Med = median( Rcpp::as<Rcpp::NumericVector>(wrap(m)) );
		double Mad = mad ( Rcpp::as<Rcpp::NumericVector>(wrap(m)) );
		NumericVector dis( nv );
		double crit = pow( R::qchisq( 0.975, 1, 1, 0 ), 0.5 );
		for( int i = 0; i < nv; i++ ){
			dis(i) = ( m(i, 0) - Med )*( m(i, 0) - Med )/( Mad*Mad );
			dis(i) = pow( dis(i), 0.5);
			if( dis(i) > crit ){
				flag( i ) = 1;
			}
		}
	}
	if( nc > 1 ){
		if( R_IsNA( gval ) && cop == 1 )
			gval = pow( R::qchisq(0.95, (double)nc, 1, 0), 0.5 );
		if( R_IsNA( gval ) && cop != 1 )
			gval = pow( R::qchisq(0.975, (double)nc, 1, 0), 0.5 );
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
			center = rmba( m , 5)(0);
		}
		NumericMatrix m2 = Rcpp::as<Rcpp::NumericMatrix>( wrap(m) );
		NumericVector B(nc), A(nc), BB(nc), temp(nc), temp_idealf(2);
		double bot;
		NumericVector dis(nv);
		double cu;
		for( int i=0; i < nv; i++ ){
			B = m2(i,_) - center;
			BB = B*B;
			bot = sum(BB); 
			if(bot!=0){
				for( int j = 0; j < nv; j++ ){
					A = m2(j,_) - center;
					temp = sum( A * B ) * B / bot;
					dis(j) = pow(sum( temp * temp ), 0.5);
				}
				temp_idealf = idealf(dis);
				if( !mm ){
					cu = median( dis ) + gval*( temp_idealf( 1 ) - temp_idealf( 0 ) );

				}else 
					cu = median( dis ) + gval*mad( dis );
				for( int k = 0; k < nv; k++ ){
					if( dis(k) > cu )
						flag( k ) = 1;
				}
			}
		}		
	}
	keepid = find( flag == 0 );
	return keepid;
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
	return(wrap(near( Rcpp::as<arma::vec>(wrap(X)), 
	 				  Rcpp::as<double>(wrap(PT)))));
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


double tmean(arma::vec x, double tr=0.2){
	int n = x.n_elem, lo = floor(n * tr), hi = n - lo - 1;
	arma::vec xcopy = arma::sort( x );
	return arma::as_scalar(arma::mean(xcopy.subvec( span(lo, hi) )));
}

RcppExport SEXP tmean(SEXP X, SEXP TR){
	BEGIN_RCPP
	return(wrap(tmean(Rcpp::as<arma::vec>(wrap(X)), Rcpp::as<double>(wrap(TR)))));
	END_RCPP
}


List aov2depth(List x1, List x2, int nboot=500, double estimator=1, double tr=0.2){
	int J = x1.size();
	if( J != x2.size() )
		throw std::runtime_error("x1 and x2 should have same number of groups");
	NumericVector est1(J), est2(J);
	arma::mat difb(nboot +1, J), dif(1, J);
	difb.zeros();
	arma::vec tempx1, tempx2;
	int n1, n2;

	for( int j=0; j < J; j++ ){
		R_CheckUserInterrupt();
		tempx1 = Rcpp::as<arma::vec>( x1(j) ), tempx2 = Rcpp::as<arma::vec>( x2(j) );
		n1 = tempx1.n_elem;
		n2 = tempx2.n_elem;
		if( estimator == 1){
			est1(j) = tmean(tempx1, tr);
			est2(j) = tmean(tempx2, tr);
		} else if ( estimator == 2){
			est1(j) = median(tempx1);
			est2(j) = median(tempx2);
		} else if ( estimator == 3){
			est1(j) = mean(tempx1);
			est2(j) = mean(tempx2);
		}else if ( estimator == 4){
			est1(j) = hd_C(tempx1, tr);
			est2(j) = hd_C(tempx2, tr);
		}
		dif(0,j)  = est1(j) - est2(j);

		arma::uvec n1boot( n1 ), n2boot( n2 );
		IntegerVector n1boot2( n1 ), n2boot2( n2 );

		for( int k=0; k < nboot; k++ ){
			n1boot2 = floor(runif( n1, 0, n1  ));
			n2boot2 = floor(runif( n2, 0, n2  ));
			n1boot  = Rcpp::as<arma::uvec>( n1boot2 );
			n2boot  = Rcpp::as<arma::uvec>( n2boot2 );
			if( estimator == 1){
				difb(k, j) = tmean( tempx1.elem(n1boot), tr) -
						     tmean( tempx2.elem(n2boot), tr);
			} else if ( estimator == 2){
				difb(k, j) = median( tempx1.elem(n1boot)) -
					     	 median( tempx2.elem(n2boot));
			} else if ( estimator == 3){
				difb(k, j) = mean( tempx1.elem(n1boot)) -
					     	 mean( tempx2.elem(n2boot));
			} else if ( estimator == 4){
				difb(k, j) = hd_C( tempx1.elem(n1boot), tr) -
					     	 hd_C( tempx2.elem(n2boot), tr);
			}
		}
	}
	arma::mat m1 = arma::cov( difb.submat(0, 0, nboot-2, J-1) );
	arma::mat dis = mahalanobis( difb, dif.t(), m1 );
	int bplus = nboot + 1;
	arma::uvec sig_vec = find( dis.col(0) >= dis(bplus-1, 0) );
	double sig = sig_vec.n_elem*1.0/(1.0*bplus);
	return List::create( _["pvalue"] = sig, 
						 _["est1"] 	= est1, 
						 _["est2"] 	= est2,
						 _["dif"] 	= dif,
						 _["n1"] 	= n1,
						 _["n2"] 	= n2);
}


List ancGLOB_sub3_C(arma::mat x1, arma::vec y1, arma::mat x2, arma::vec y2, bool xout, SEXP PTS,
					int nmin, int nboot, double fr1, double fr2, int estimator, double tr=0.2){
	List g1, g2;
	IntegerVector n3;
	if(xout){
		arma::mat tempx1 = x1;
		arma::mat tempx2 = x2;
		arma::vec tempy1 = y1;
		arma::vec tempy2 = y2;
		
		NumericVector center(tempx1.n_cols, NA_REAL);
		double gval = NA_REAL; 
		LogicalVector mm(1, FALSE);
		IntegerVector keepid1 = outpro_C(Rcpp::as<Rcpp::NumericMatrix>(wrap(tempx1)), 
										 center, mm, 3, 0.2, 0.5, gval);
		IntegerVector keepid2 = outpro_C(Rcpp::as<Rcpp::NumericMatrix>(wrap(tempx2)), 
										 center, mm, 3, 0.2, 0.5, gval);
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
		xorder1 = sort_index( x2.col(0) );
		x2 = x2.rows(xorder1);
		y2 = y2.elem(xorder1);
		arma::uvec tempid1, tempid2;
		std::vector<int> n1n2( 2 );
		std::vector<int> sub;
		int vecn;
		for( int i=0; i < N1; i++){
			tempid1 = find( near(x1.col(0), x1(i,0), fr1) == 1 );
			tempid2 = find( near(x2.col(0), x1(i,0), fr2) == 1 );
			n1n2[ 0 ] = tempid1.n_elem;
			n1n2[ 1 ] = tempid2.n_elem;
			vecn = *std::min_element( n1n2.begin(), n1n2.end() );
			if ( vecn >= nmin ){
				sub.push_back(i);
			}
		}
		if(sub.size() > 0){
			std::vector<double> pts;
			isub(0) = *std::min_element( sub.begin(), sub.end() );
			isub(4) = *std::max_element( sub.begin(), sub.end() );
			isub(2) = floor((isub(0) + isub(4))/2);
			isub(1) = floor((isub(0) + isub(2))/2);
			isub(3) = floor((isub(2) + isub(4))/2);

			for( int i=0; i < 5; i++ ){
				pts.push_back( arma::as_scalar(x1(isub(i), 0)) );
				arma::vec tempy1 = y1.elem( find(near(x1.col(0), pts[i], fr1) == 1) );
				arma::vec tempy2 = y2.elem( find(near(x2.col(0), pts[i], fr2) == 1) );
				g1.push_back(  tempy1 );
				g2.push_back(  tempy2 );
				n3.push_back(tempy1.n_elem);
			}
		}
	}else{
		NumericVector pts(PTS);
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
	bool xout = FALSE;
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


IntegerVector outbox_C( NumericVector x, double gval=NA_REAL, bool mbox=false ){
	NumericVector temp = idealf( x );
	double cl, cu;
	int n = x.size();
	if( mbox ){
		if( R_IsNA( gval ) )
			gval = ( 17.63 * n - 23.64 )/( 7.74 * n - 3.71 );
		cl = median( x ) - gval * ( temp( 1 ) - temp( 0 ) );
		cu = median( x ) + gval * ( temp( 1 ) - temp( 0 ) );
	} else {
		if( R_IsNA( gval ) )
			gval = 1.5;
		cl = temp( 0 ) - gval * ( temp( 1 ) - temp( 0 ) );
		cu = temp( 1 ) + gval * ( temp( 1 ) - temp( 0 ) );
	}
	std::vector<int> outid;

	for( int i = 0; i < n; i++ ){
		if( x(i) > cu || x(i) < cl  )
			outid.push_back( i );
	}
	return Rcpp::as<Rcpp::IntegerVector>( wrap( outid ) );
}

RcppExport SEXP outbox_C( SEXP X, SEXP GVAL, SEXP MBOX ){
	BEGIN_RCPP
	return outbox_C( Rcpp::as<NumericVector>( wrap(X) ),
					 Rcpp::as<double>( wrap(GVAL) ),
					 Rcpp::as<bool>( wrap(MBOX) ));

	END_RCPP
}



NumericVector skip(NumericMatrix m, bool mm, int outpro_cop, double tr, double q){
	NumericVector center(m.ncol(), NA_REAL);
	IntegerVector keepid = outpro_C(m, center, mm, outpro_cop, tr, q);
	arma::mat val = arma::mean((Rcpp::as<arma::mat>(m)).rows(Rcpp::as<arma::uvec>(keepid)));
	return(wrap(val));
}

NumericVector skip(arma::mat m, bool mm, int outpro_cop, double tr, double q){
	NumericVector center(m.n_cols, NA_REAL);
	arma::uvec keepid = outpro_C(m, center, mm, outpro_cop, tr, q);
	arma::mat val = arma::mean( m.rows( keepid ) );
	return(wrap(val));
}


RcppExport SEXP skip(SEXP M, SEXP MM, SEXP OUTPRO_COP, SEXP TR, SEXP Q){
	BEGIN_RCPP
	return skip( Rcpp::as<Rcpp::NumericMatrix>( wrap( M )),
				 Rcpp::as<bool>( wrap( MM ) ),
				 Rcpp::as<int>( wrap(OUTPRO_COP) ),
				 Rcpp::as<double>( wrap(TR) ),
				 Rcpp::as<double>( wrap(Q) )); 
	END_RCPP
}

RcppExport SEXP skip_boot(SEXP M, SEXP MM, SEXP OUTPRO_COP, SEXP BOOTID){
	BEGIN_RCPP
	
	arma::mat m = Rcpp::as<arma::mat>( wrap(M) );
	int outpro_cop = Rcpp::as<int>( wrap(OUTPRO_COP) );
	List bootid( BOOTID );
	bool mm = Rcpp::as<bool>( wrap( MM ) );
	int n = m.n_rows, nc = m.n_cols, nboot = bootid.size();
	NumericVector center(nc, NA_REAL);
	double q = 0.5, tr = 0.2;
	NumericVector output( nboot * nc ), temp( nc );
	arma::uvec tempid = uvec( n );
	
	for( int i = 0; i < nboot; i++ ){
		R_CheckUserInterrupt();
		temp = skip( m.rows( Rcpp::as<arma::uvec>( bootid(i) ) ), mm, outpro_cop, tr, q );
		for( int j = 0; j < nc; j++ ){
			output( i * nc + j) = temp(j);
		}
	}
	return output;
	END_RCPP
}

