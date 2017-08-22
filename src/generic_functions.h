#ifndef NONLINEARTSERIES_GENERIC_FUNCTIONS_H
#define NONLINEARTSERIES_GENERIC_FUNCTIONS_H
#include <Rcpp.h>

/* Modifies the C mod function in order to return only positives values */
inline int positive_modulo(int x, int m) {
  return (((x % m) + m) % m);
}

Rcpp::NumericMatrix build_takens(const Rcpp::NumericVector& timeSeries, 
                                int embeddingDimension, int timeLag);

#endif
