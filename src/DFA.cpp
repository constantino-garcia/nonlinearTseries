#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double calculate_fluctuation(NumericVector& yr, int windowSize) {
  int nWindows = floor((double)yr.size() / windowSize);
  /* reuses memory and avoids extra copy */
  arma::colvec y(yr.begin(), yr.size(), false);
  
  /* Design matrix for regression */
  arma::mat X = arma::ones(windowSize, 2);
  X.col(1) = arma::linspace<arma::vec>(1, windowSize, windowSize);
  
  /* vector storing the residuals of each window */
  arma::colvec windowResiduals = arma::vec(nWindows);
  for (int i=0; i < nWindows; i++) {
    arma::colvec y_subvec = y.subvec(i * windowSize, (i + 1) * windowSize - 1);
    /* fit model y ~ X */
    arma::colvec coef = arma::solve(X, y_subvec);  
    arma::colvec residualsPow2 = arma::pow(y_subvec - X * coef, 2);
    windowResiduals(i) = 
      std::accumulate(residualsPow2.begin(), residualsPow2.end(), 0.0) / windowSize;
    
  }
  return std::sqrt(
    std::accumulate(windowResiduals.begin(), windowResiduals.end(), 0.0) / nWindows
  );
}


// [[Rcpp::export]]
NumericVector calculate_fluctuation_function(NumericVector& yr, 
                                             NumericVector& windowSizesVector){
  int nWindows = windowSizesVector.size();
  NumericVector fluctuationFunction(nWindows);
  for (int i = 0; i < nWindows; i ++) {
    fluctuationFunction[i] =  calculate_fluctuation(yr, windowSizesVector(i));
  }
  return fluctuationFunction;
}

 
