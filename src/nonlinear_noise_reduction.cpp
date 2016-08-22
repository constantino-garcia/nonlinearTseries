#include <Rcpp.h>
#include <cmath>
#include "generic_functions.h"
#include "neighbour_search.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector nonlinear_noise_reduction(const NumericVector& timeSeries, 
                               int embeddingDimension,
                               double radius, int nBoxes){
  NumericVector denoisedTimeSeries(clone(timeSeries));
  /* The neighbours of the n-th Takens' vector shall be used to denoise x[n + (m - 1) / 2] */
  int midPosition = static_cast<int>(std::floor( embeddingDimension / 2.0));
  /* A fixed time lag of 1 is used */
  NumericMatrix phaseSpace = build_takens(timeSeries, embeddingDimension, 1);
  neighbour_search neighbourSearcher(phaseSpace, radius, nBoxes);
  int nTakensVectors = phaseSpace.nrow();
  for (int i = 0; i < nTakensVectors; i++){
    IntegerVector neighbours = neighbourSearcher.find_neighbours(i);
    int nNeighbours = neighbours.size();
    if (nNeighbours == 0){
      continue;
    } 
    for (int j=0;j < nNeighbours; j++){
      denoisedTimeSeries[i+midPosition] +=phaseSpace(neighbours[j], midPosition);
    }
    //mean using the neighbours and the proper point
    denoisedTimeSeries[i + midPosition] = denoisedTimeSeries[i + midPosition] / (nNeighbours + 1);
  }
  return denoisedTimeSeries;
}


