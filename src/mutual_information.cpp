#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;

#include <Rcpp.h>
#include <math.h> 
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix calculate_time_series_histogram(const NumericVector& x, const int& timeLag,
                                              const int& nPartitions) {
  int indy, binx, biny;
  NumericMatrix histogram(nPartitions,nPartitions);
  
  int endPoint = x.size() - timeLag;
  double sampleWeight =  1.0 / static_cast<double>(endPoint);
  for(int indx = 0; indx < endPoint; indx++) {
    indy = indx + timeLag;
    binx = std::min(static_cast<int>(x[indx] * nPartitions), nPartitions - 1);
    biny = std::min(static_cast<int>(x[indy] * nPartitions), nPartitions - 1);
    histogram(binx, biny) += sampleWeight;
  }
  return histogram;
}


NumericVector calculate_marginal_histogram(const NumericMatrix& histogram){
  int histogramDimension = histogram.nrow();
  NumericVector marginalHistogram(histogramDimension);
  NumericMatrix histogramCopy(histogram);
  
  for (int i=0;i < histogramDimension; i++){
    marginalHistogram[i] = sum(histogramCopy(i,_));
  }
  
  return marginalHistogram;
}


double calculate_mutual_information(const NumericMatrix& histogram){
  NumericVector marginalHistogramX = calculate_marginal_histogram(histogram);
  int histogramSize = marginalHistogramX.size();
  double mutualInformation = 0.0;
  
  /* implement sum(histogram * log(histogram)) - 2 * sum( marginalHistogramX * log(marginalHistogramX) ) */
  for (int i=0; i < histogramSize; i++){
    for (int j=0; j < histogramSize; j++){
      if ( histogram(i,j) > 0.0){
        mutualInformation += histogram(i,j) * std::log(histogram(i,j));
      }
    }
    if (marginalHistogramX[i] > 0.0){
      mutualInformation -= 2 * marginalHistogramX[i] * std::log(marginalHistogramX[i]);
    }
  }
  return (mutualInformation);
}

// [[Rcpp::export]]
NumericVector calculate_mutual_information(const NumericVector& tseries,
                                           const int& maxlag, const int& nPartitions){
  NumericVector mutualInformation(maxlag + 1);
  NumericVector x(clone(tseries));
  /* reescale time series between 0 and 1 */
  x = (x - min(x)) / (max(x) - min(x));
  
  for (int timeLag = 0; timeLag <= maxlag; timeLag++) {
    mutualInformation[timeLag] = 
      calculate_mutual_information(calculate_time_series_histogram(x, timeLag, nPartitions));
  }
  
  return mutualInformation;
}


