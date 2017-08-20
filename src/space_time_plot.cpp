#include <Rcpp.h>
#include "generic_functions.h"
using namespace Rcpp;

/* TODO: until now we assume that vectors have same length */
/*TODO: integrate with neighbour_search */
inline double max_distance_between_rows(const NumericMatrix::Row& v1, const NumericMatrix::Row& v2){
  double distance = -1.0;
  int embeddingDimension = v1.size();
  for (int i=0; i < embeddingDimension; i++){
    distance = std::max(distance,  std::abs(v1[i]- v2[i]));
  }
  return distance;
}


//[[Rcpp::export]]
NumericMatrix space_time_plot(NumericMatrix phaseSpace, NumericVector radii, 
                              int nTimeSteps, int timeStep,
                              int  nPercentages) {
  NumericMatrix spaceTimePlot(nPercentages, nTimeSteps);
  int radiiSize = radii.size();
  int lastPosition = radiiSize - 1;
  double maxRadius = radii[lastPosition];
  int nPhaseSpaceVectors = phaseSpace.nrow();
  /* We shall start at a convenient phase space vector to avoid out of bounds errors */
  int firstPhaseSpaceVector =  nTimeSteps * timeStep;
  /* Hence the number of available phase space vectors is ... */
  int nAvailablePhaseSpaceVectors = (nPhaseSpaceVectors - firstPhaseSpaceVector);

  for (int iTime=0; iTime < nTimeSteps; iTime++) {
    /* Initialice the histogram of distances */
    IntegerVector radiusHistogram(radiiSize, 0);
    for (int iTakensVector = firstPhaseSpaceVector; iTakensVector < nPhaseSpaceVectors; iTakensVector++){
      /* calculate distances for the current timeStep = iTime * timeStep */
      double distance = max_distance_between_rows(phaseSpace(iTakensVector, _), 
                                                  phaseSpace(iTakensVector - iTime * timeStep, _));
      int iHistogram = static_cast<int>(distance / maxRadius *  radiiSize);
      radiusHistogram[std::min(iHistogram,lastPosition)]++;
    }
    for (int iPercentage=0; iPercentage < nPercentages; iPercentage++){
      int nRequiredPoints =  static_cast<int>(
        nAvailablePhaseSpaceVectors * (iPercentage + 1) / static_cast<double>(nPercentages)
      );
      int nPoints=0;
      int i;
      for (i=0; (i < radiiSize) && (nPoints < nRequiredPoints); i++){
        nPoints += radiusHistogram[i];
      }
      /* Store the radius at which we obtain the desired number of neighbours 
       * for a given time delay.
       */
      spaceTimePlot(iPercentage, iTime) = radii[i];
    }
  }
  return spaceTimePlot; 
}
