#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;

 
// [[Rcpp::export]]
NumericMatrix tsHistogram(NumericVector& x,int& tlag, const int& npartitions) {
  int indy, binx, biny;
  NumericMatrix histogram(npartitions,npartitions);
  
  int tslen = x.size();
  for(int indx = 0; indx < (tslen-tlag); indx++) {
    indy = indx + tlag;
		binx = (int) std::min( (int)(x[indx] * npartitions), npartitions - 1);
		biny = (int) std::min( (int)(x[indy] * npartitions), npartitions - 1);
		histogram(binx, biny) += 1.0 / (tslen - tlag);
    
	}
  return histogram;
}


NumericVector marginalHistogram(NumericMatrix& hist){
  int histdim = hist.nrow();
  NumericVector marginalHist(histdim);
  
  for (int i=0;i < histdim; i++){
    marginalHist[i] = sum(hist(i,_));
  }
  
  return marginalHist;
}


double mutInfFromHist(NumericMatrix& hist){
  NumericVector histx = marginalHistogram(hist);
  int histlen = histx.size();
  double mutinf = 0.0;

  //implement sum(hist * log(hist)) - 2 * sum( histx * log(histx) )
  // first term
  for (int i=0; i < histlen; i++){
    for (int j=0; j < histlen; j++){
      if ( hist(i,j) > 0.0){
        mutinf += hist(i,j) * log( hist(i,j) );
      }
    }  
  }
  // second term
  for (int i=0; i < histlen; i++){
    if ( histx[i] > 0.0){
      mutinf -= 2 * histx[i] * log( histx[i] );
    }
  }

  return (mutinf);
}

  
// [[Rcpp::export]]
NumericVector mutualInformation(const NumericVector& tseries,const int& maxlag,
const int& npartitions){
  NumericVector mutinf(maxlag + 1);
  NumericVector x(clone(tseries));
  // reescale time series between 0 and 1
  x = ( x - min(x) ) / ( max(x) - min(x) );
    
  for(int tlag = 0; tlag <= maxlag; tlag++) {
		NumericMatrix hist = tsHistogram(x, tlag, npartitions);
    mutinf[tlag] = mutInfFromHist(hist);
		
	}
	return mutinf;
}
