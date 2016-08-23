#ifndef NONLINEARTSERIES_GENERIC_FUNCTIONS_H
#define NONLINEARTSERIES_GENERIC_FUNCTIONS_H
#include <Rcpp.h>


Rcpp::NumericMatrix build_takens(const Rcpp::NumericVector& timeSeries, 
                                int embeddingDimension, int timeLag);

#endif
