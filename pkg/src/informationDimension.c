#include "neighbourSearch.h"
#include "genericFunctions.h"  
#include <math.h>
#include <Rmath.h> //include digamma function
#include <R.h> //include sort functionality

inline double which(double *distVector,int k,int numberValidNeighs){
   R_rsort(distVector,numberValidNeighs);
   return (distVector[k]);
}

/******************************************************************************/
/********************************** d1 ****************************************/
/******************************************************************************/
void d1(double *takens, int *numberTakens, int *embeddingD, double* fixedMass,
        double *eps, double *increasingEpsFactor, int *numberBoxes,int *boxes, 
        int *numberReferenceVectors, int* theilerWindow, int* kMax,
        double* averageLogRadius){
  
  // declare variables
  int error=0;
  int i,iiiii, theilerMargin, k, takensVectorsUsed, remainingReferenceVectors,takensIterator,
      nTakensWithInsufficientNeighs=0, currentReferenceVector, nfound, neigh, neighIt,
      numberValidNeighs;
  double lnFixedMass, currentEps, radius;
  int possibleNeighbours[*numberTakens];
  int neighList[*numberTakens];
  int referenceVectors[*numberReferenceVectors];
  double distVector[*numberTakens];
  /*********************** computing k and N from  fixedMass=k/N **************/
  // when estimating k (from the formula fixedMass=k/N) we have to take into account
  // that there exist a theiler window!! Thus, for every takens' vector i, we cannot
  // use as neighbours those ranging from (i-theilerWindow):(i+theilerWindow). That is,
  // the effective number of takens vector is numberTakens - (2*theilerWindow+1)
  theilerMargin = 2*(*theilerWindow)+1;
  takensVectorsUsed = (*numberTakens);
  k = (int)( (*fixedMass)*( (*numberTakens)-theilerMargin  )) + 1;
  // avoid looking for too many neighbours by imposing a kmax. If k > kmax we will reduce
  // the number of takens vectors considered so that k/N=kMax/N_prime
  if (k>(*kMax) ){
    takensVectorsUsed =  (int)((  (double)( ( (*numberTakens)-theilerMargin  )*(*kMax)  ) )/( (double)(k) ))
                          + theilerMargin;
    k = (*kMax);
  }
  /******************************** estimator for ln(fixedMass) ****************/
  // use grassberger estimator for ln (fixedMass). This estimator takes into account the
  // finite nature of the estimation. See Generalizations of the Hausdorff dimension
  // of fractal measures (Grassberger 1985)
  lnFixedMass = digamma( (double)(k) ) - log10((double)(takensVectorsUsed-theilerMargin)); 
  /****************************** Prepare the iterations **********************/
  // this variable will store the average radious
  *averageLogRadius = 0;
  // vector that will store the number of the reference takens' vector for which
  // we have not yet found enough neighbours. At the beggining they will be all 
  // the first numberReferenceVectors vectors...Thus the remaining reference vectors
  // will be numberReferenceVectors
  for (i=0;i < (*numberReferenceVectors); i++){
    referenceVectors[i]=i;
  }
  
  remainingReferenceVectors = (*numberReferenceVectors);
/*************************** iterate ***************************************/
// Go!: while there is some reference vector with not enough vectors in its neighbourhood
// we will iterate increasing the neighbourhood radius

  for (currentEps =(*eps); remainingReferenceVectors > 0; currentEps*=(*increasingEpsFactor)){
    // use the box assisted algorithm with the current radious
    boxAssistant( takens, numberTakens, embeddingD,
    &currentEps, numberBoxes, boxes, possibleNeighbours);
    // iterate over the vector that stores the number of each reference Takens' vector
    for (takensIterator=0, nTakensWithInsufficientNeighs=0; takensIterator < remainingReferenceVectors; takensIterator++){
      currentReferenceVector=referenceVectors[takensIterator];
      // neighbour search
      neighbourSearchFromBoxes(takens, &currentReferenceVector, numberTakens, embeddingD,
               &currentEps, numberBoxes, boxes, possibleNeighbours, neighList, &nfound);

     /****************************** check neighs **********************************/
     // for each neighbour of the neighList, check if they fulfill the theiler distance
     // and the takensVectorsUsed that we can use
     for (neighIt=0,numberValidNeighs=0;neighIt<nfound;neighIt++){
       neigh=neighList[neighIt];
       // if the current neighbour does not meet the conditions, move to the next
       if ( (abs(currentReferenceVector-neigh)<=(*theilerWindow)) || (neigh > takensVectorsUsed ) ) continue;
       // the current neighbour does meet the conditions (numberValidNeighs++)!! we must compute the distance
       // for using it at the computation of  the radious that encloses k neighbours...
       distVector[numberValidNeighs++]=distance(currentReferenceVector,neigh,takens,*numberTakens,*embeddingD);
     }
     /************* check if we have found enough numberValidNeighs for averaging *****/
     if (numberValidNeighs >= k){
       // which is k-th smallest element of the distVector(which has length=numberValidNeighs)
       radius = which(distVector,k,numberValidNeighs);
       (*averageLogRadius) += log10( radius );
     }else{
       // not enough vectors found... mark for next sweep and increase the number
       // of reference takens with insufficiente neighbours
       referenceVectors[nTakensWithInsufficientNeighs++] = currentReferenceVector;
     }
     
    }
    remainingReferenceVectors = nTakensWithInsufficientNeighs;
  }
  (*averageLogRadius) /= ((double)(*numberReferenceVectors)) ;
}


/******************************************************************************/
/************************* informationDimension *******************************/
/******************************************************************************/
// warning!!!! averageLogRadiusVector and fixedMassVector must have the same length
void informationDimension(double *takens, int *numberTakens, int *embeddingD, double* fixedMassVector,
                          int* fixedMassVectorLength,double *eps, double *increasingEpsFactor,
                          int *numberBoxes,int *boxes, int *numberReferenceVectors, 
                          int* theilerWindow, int* kMax,double* averageLogRadiusVector){
  // declare variables
  int i;
  double fixedMass,averageLogRadius;
  // find the averageLogRadiusVector for each fixedMass value in fixedMassVector
  for (i=0;i<(*fixedMassVectorLength);i++){
    fixedMass = fixedMassVector[i];
    d1(takens, numberTakens, embeddingD, &fixedMass, eps, increasingEpsFactor,
       numberBoxes, boxes, numberReferenceVectors,theilerWindow, kMax, 
       &averageLogRadius);
    //store the averageLogRadius in the averageLogRadiusVector 
    averageLogRadiusVector[i] = averageLogRadius;
  }            
}
