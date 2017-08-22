#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <numeric>
using namespace Rcpp;

// calculate an aproximation to the point of intersection between the hiperplane
// and the trajectory using linear interpolation
double calculate_crossings(const NumericMatrix& timeSeries, 
                           double sideParam1, double sideParam2,
                           int pos1, int pos2,
                           int dimension, NumericVector& crossPoint){
  // compute the crossing time using linear interpolation
  double crossingTime = (
    pos1 * sideParam2 / (sideParam2 - sideParam1) +
    pos2 * sideParam1 / (sideParam1 - sideParam2)
  );
  // calculate the point of intersection using linear interpolation
 
  for (int i = 0; i < dimension; i++){
    crossPoint[i] = (
      timeSeries(pos1,i) + 
      (timeSeries(pos2, i) - timeSeries(pos1, i)) * (crossingTime - pos1) / static_cast<double>(pos2-pos1)
    );
  }
  return crossingTime;
}

// timeSeries was parsed from a matrix where each rows represented a point. When
// descending by rows, the time is advancing. nPoints is the number of rows
// of the original matrix
// [[Rcpp::export]]
List poincare_map(const NumericMatrix& timeSeries, 
                  const NumericVector& hiperplanePoint, const NumericVector& normalVector){
  int nPoints = timeSeries.nrow();
  int dimension = timeSeries.ncol();
  
  // Create matrices with the maximum posible number of nPoints since NumericMatrix
  // does not grow gracefully. Not sure if this is the best approach, though
  NumericMatrix poincareMapSeries(nPoints, dimension);
  NumericMatrix positivePoincareMapSeries(nPoints, dimension);
  NumericMatrix negativePoincareMapSeries(nPoints, dimension);
  std::vector<double> crossingTime;
  std::vector<double> positiveCrossingTime;
  std::vector<double> negativeCrossingTime;
  int numberCrossings = 0, numberPositiveCrossings = 0, numberNegativeCrossings = 0;
               
  NumericVector incidentVector(dimension);
  NumericVector crossPoint(dimension);
  double referenceSide, currentSide;
  // get the incident vector the Side of the hiperplane in which the point lies
  for (int i = 0; i < dimension; i++) {
    incidentVector[i] = timeSeries(0, i) - hiperplanePoint[i];
  } 
  currentSide = std::inner_product(normalVector.begin(), normalVector.end(),
                                   incidentVector.begin(), 0.0);
  for (int checkingPointIt = 1; checkingPointIt < nPoints; checkingPointIt++){
    referenceSide = currentSide;
    // Let's find the next hiperplane crossing...
    // get the incident vector and calculate the current Side of this point
    for (int i = 0; i < dimension; i++) {
      incidentVector[i] = timeSeries(checkingPointIt,i) - hiperplanePoint[i];  
    }      
    currentSide = std::inner_product(normalVector.begin(), normalVector.end(),
                                     incidentVector.begin(), 0.0);
    // is there a cross?
    if (currentSide * referenceSide < 0.0){
      double currentCrossingTime = calculate_crossings(
        timeSeries, referenceSide, currentSide,
        checkingPointIt - 1, checkingPointIt,
        dimension, crossPoint);
      crossingTime.push_back(currentCrossingTime);
      //store in the poincare map for the positive direction and in the two sided
      //poincare map
      if (referenceSide > 0.0){      
        positiveCrossingTime.push_back(currentCrossingTime);
        for (int i = 0; i < dimension; i++) {
          poincareMapSeries(numberCrossings,i) = crossPoint[i];
          positivePoincareMapSeries(numberPositiveCrossings,i) = crossPoint[i];  
        }
        numberCrossings++;
        numberPositiveCrossings++;
      } else { 
        //store in the poincare map for the negative direction and in the
        // two sided poincare map
        for (int i = 0; i < dimension; i++) {
          poincareMapSeries(numberCrossings, i) = crossPoint[i];
          negativePoincareMapSeries(numberNegativeCrossings,i) = crossPoint[i];  
        }
        negativeCrossingTime.push_back(currentCrossingTime);
        numberCrossings++;
        numberNegativeCrossings++;
      } 
    }
  }
  List ret;
  // TODO; return only the proper submatrices and transform vector to numericVecot
  ret["pm"] = poincareMapSeries(Range(0, numberCrossings -1),
                                Range(0, dimension -1));
  NumericVector pmTime = wrap(crossingTime);
  ret["pm.time"] = pmTime;
  ret["pm.pos"] = positivePoincareMapSeries(Range(0, numberPositiveCrossings -1), 
                                            Range(0, dimension -1));
  NumericVector pmPosTime = wrap(positiveCrossingTime);
  ret["pm.pos.time"] = pmPosTime;
  ret["pm.neg"] =  negativePoincareMapSeries(Range(0, numberNegativeCrossings -1), 
                                             Range(0, dimension -1));
  NumericVector pmNegTime = wrap(negativeCrossingTime);
  ret["pm.neg.time"] = pmNegTime;
  return ret;
}


