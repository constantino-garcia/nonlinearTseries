#include <Rcpp.h>
#include "generic_functions.h"
using namespace Rcpp;


NumericMatrix build_takens(const NumericVector& timeSeries,
                           int embeddingDimension, 
                           int timeLag){
  int maxJump = (embeddingDimension - 1) * timeLag;
  IntegerVector jumps(embeddingDimension);
  for (int i=0; i < jumps.size(); i ++){
    jumps[i] = i * timeLag;
  }
  int nTakensVectors = timeSeries.size() - maxJump;
  NumericMatrix takensSpace(nTakensVectors, embeddingDimension);
  
  /* matrix that will store the Takens' vectors. One vector per row */
  for (int i=0; i < nTakensVectors; i++) {
    for (int j=0; j < embeddingDimension; j++){
      takensSpace(i, j) = timeSeries[i + jumps[j]];
    }
  }
  return takensSpace;
}
