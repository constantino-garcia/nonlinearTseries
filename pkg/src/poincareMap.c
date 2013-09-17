#include <R.h>
#include <Rmath.h>
#include <stdlib.h>
#include "genericFunctions.h"

// calculate an aproximation to the point of intersection between the hiperplane
// and the trajectory using linear interpolation
double calculateCrossings(double* timeSeries,int nPoints, double sideParam1,
                                 double sideParam2, int pos1, int pos2,
                                 int dimension, double* crossPoint){
    // compute the crossing time using linear interpolation
    double crossingTime = sideParam2/(sideParam2-sideParam1)*((double)(pos1)) 
                          + sideParam1/(sideParam1-sideParam2)*((double)(pos2));
    // calculate the point of intersection using linear interpolation
    int i;
    for (i=0;i<dimension;i++){
      crossPoint[i] = MAT_ELEM(timeSeries,pos1,i,nPoints) + 
                      (MAT_ELEM(timeSeries,pos2,i,nPoints)-MAT_ELEM(timeSeries,pos1,i,nPoints))*
                      (crossingTime-(double)(pos1))/((double)(pos2-pos1));
    }
    return crossingTime;
                                   
}

// timeSeries was parsed from a matrix where each rows represented a point. When
// descending by rows, the time is advancing. nPoints is the number of rows
// of the original matrix
// dimension is the dimension of each point of the time series
// poincareMaps will store the crossing points
// crossingTime: time at which a cross happens
// numberCrossings that we have found
// poincareMapSeries must be declared in R with the same dimensions as the timeSeries matrix
void poincareMap(double* timeSeries, int* nPoints, int* dimension, 
                double* poincareMapSeries, double* positivePoincareMapSeries,
                double* negativePoincareMapSeries, double* crossingTime, 
                double* positiveCrossingTime, double* negativeCrossingTime,
                int* numberCrossings,  int* numberPositiveCrossings, int* numberNegativeCrossings,
                double* hiperplanePoint, double* normalVector){
  
  int i, checkingPointIt;
  *numberCrossings = 0;  *numberPositiveCrossings = 0; *numberNegativeCrossings = 0;
  double incidentVector[(*dimension)];
  double crossPoint[(*dimension)];
  double referenceSide, currentSide;
  
  // get the incident vector the Side of the hiperplane in which the point lies
  // for the point 
  for (i=0;i<(*dimension);i++){
    incidentVector[i] = MAT_ELEM(timeSeries,0,i,*nPoints) - hiperplanePoint[i];
  } 
  currentSide = scalarProduct (normalVector, incidentVector, *dimension);
  
  for (checkingPointIt=1; checkingPointIt < (*nPoints) ; checkingPointIt++){
     referenceSide = currentSide;
     // Let's find the next hiperplane crossing...
     // get the incident vector a calculate the current Side of this point
      for (i=0;i<(*dimension);i++) {
        incidentVector[i] = MAT_ELEM(timeSeries,checkingPointIt,i,*nPoints) - hiperplanePoint[i];  
      }      
      currentSide = scalarProduct (normalVector, incidentVector, *dimension);
      // is there a cross?
      if (currentSide*referenceSide<0.0){
         crossingTime[*numberCrossings] = calculateCrossings(timeSeries, *nPoints, referenceSide, 
                                          currentSide, checkingPointIt-1,checkingPointIt,
                                          *dimension,crossPoint);
         //store in the poincare map for the positive direction and in the two sided
         //poincare map
         if (referenceSide > 0.0){      
                for (i=0;i<(*dimension);i++) {
                MAT_ELEM(poincareMapSeries,*numberCrossings,i,*nPoints) = crossPoint[i];
                MAT_ELEM(positivePoincareMapSeries,*numberPositiveCrossings,i,*nPoints) = crossPoint[i];  
             }
             positiveCrossingTime[*numberPositiveCrossings] = crossingTime[*numberCrossings];
             (*numberCrossings)++;
             (*numberPositiveCrossings)++;
             
         }else{ //store in the poincare map for the negative direction and in the two sided poincare map
             for (i=0;i<(*dimension);i++) {
               MAT_ELEM(poincareMapSeries,*numberCrossings,i,*nPoints) = crossPoint[i];
               MAT_ELEM(negativePoincareMapSeries,*numberNegativeCrossings,i,*nPoints) = crossPoint[i];  
             }
             negativeCrossingTime[*numberNegativeCrossings] = crossingTime[*numberCrossings];
             (*numberCrossings)++;
             (*numberNegativeCrossings)++;
         } 
      }
    }
}

