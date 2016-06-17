#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double calculateFluctuation(NumericVector& yr, int windowSize) {
  int nwindows = floor(yr.size() / windowSize);
  // reuses memory and avoids extra copy
  arma::colvec y(yr.begin(), yr.size(), false);
  
  // Design matrix for regression
  arma::mat X = arma::ones(windowSize, 2);
  X.col(1) = arma::linspace<arma::vec>(1, windowSize, windowSize);
  
  // vector storing the residuals of each window
  arma::colvec resid_vector = arma::vec(nwindows);
  for (int i=0; i < nwindows; i++) {
    arma::colvec y_subvec = y.subvec(i * windowSize, (i + 1) * windowSize - 1);
    arma::colvec coef = arma::solve(X, y_subvec);  // fit model y ~ X
    arma::colvec resid_pow_2 = arma::pow(y_subvec - X * coef, 2);
    resid_vector(i) = std::accumulate(resid_pow_2.begin(),
                 resid_pow_2.end(), 0.0)/ windowSize;
    
  }
  return (sqrt(std::accumulate(resid_vector.begin(),resid_vector.end(),0.0) /
                 nwindows));
}


// [[Rcpp::export]]
NumericVector calculateFluctuationFunction(NumericVector& yr, 
                                           NumericVector& windowSizesVector){
  NumericVector fluctuation_function(windowSizesVector.size());
  for (int i = 0; i < windowSizesVector.size(); i ++) {
    fluctuation_function(i) =  calculateFluctuation(yr, windowSizesVector(i));
  }
  return fluctuation_function;
}

 
