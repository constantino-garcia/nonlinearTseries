#include <Rcpp.h>
#include "generic_functions.h"
#include "neighbour_search.h"
using namespace Rcpp;

NumericVector compute_divergence(const NumericVector& timeSeries, 
                                 int embeddingDimension, int timeLag, 
                                 double radius,  int theilerWindow,
                                 int minNumNeighbours, int nRefPoints, 
                                 int maxTimeSteps, int nBoxes){
  NumericMatrix phaseSpace = build_takens(timeSeries,embeddingDimension, timeLag);
  NumericVector Sdn(maxTimeSteps + 1);
  // variables for the box assisted algorithm
  neighbour_search neighbourSearcher(phaseSpace, radius, nBoxes);
  /* number of  Takens' vectors that the algorithm can use as reference vectors.
   * These reference vectors must have another phase space vector after 
   * maxTimeSteps. The algorithm shall use the nTakens first vectors.
   */
  int nTakens = phaseSpace.nrow() - maxTimeSteps;
  int refPointsFound = 0;
  for (int refVectorIndex = 0; (refVectorIndex < nTakens) && (refPointsFound < nRefPoints); refVectorIndex++){
    NumericVector Saux(maxTimeSteps + 1);
    IntegerVector neighboursIndexes = neighbourSearcher.find_neighbours(refVectorIndex);
    int nNeighboursFound = 0;
    for (int j = 0; j < neighboursIndexes.size(); j++) {
      int neighbourIndex = neighboursIndexes[j];
      //check if we can use this possible neighbour to average
      if (neighbourIndex >= nTakens){
        continue;
      }
      if (abs(neighbourIndex - refVectorIndex) > theilerWindow){
        nNeighboursFound++;
        int timeIndex1 = refVectorIndex + (embeddingDimension - 1) * timeLag;
        int timeIndex2 = neighbourIndex + (embeddingDimension - 1) * timeLag;
        // calculate the divergence for every value of timeStep 
        for (int timeStep = 0; timeStep <= maxTimeSteps; timeStep++) {
         Saux[timeStep] +=
            std::abs(timeSeries[timeIndex1 + timeStep] - timeSeries[timeIndex2 + timeStep]);    
        }          
      }
      //check if there are enough neighbours
      if (nNeighboursFound >= minNumNeighbours) {
        // found a reference point with enough neighbours
        refPointsFound++;
        Sdn = Sdn + log(Saux / nNeighboursFound);
      }
    }
  }
  
  // if we have found some reference point, average all refence points
  if (refPointsFound > 0){
    Sdn = Sdn / refPointsFound;
  }
  return Sdn;
}

// [[Rcpp::export]]
NumericMatrix lyapunov_exponent(
                     const NumericVector& timeSeries, 
                     int minEmbeddingDim, int maxEmbeddingDim,
                     int timeLag, double radius, int theilerWindow,
                     int minNumNeighbours, int nRefPoints, int maxTimeSteps,
                     int nBoxes){

  /* TODO: check all arguments... maybe move checking to R code */
  int nrows = maxEmbeddingDim - minEmbeddingDim + 1;
  NumericMatrix divergenceMatrix(nrows,
                                 maxTimeSteps + 1); 
   
  for (int i = 0; i < nrows; i++) {
    divergenceMatrix(i, _) = compute_divergence(timeSeries, minEmbeddingDim + i,
                     timeLag, radius, theilerWindow,
                     minNumNeighbours, nRefPoints, 
                     maxTimeSteps,nBoxes);
  }
  return divergenceMatrix;
}
