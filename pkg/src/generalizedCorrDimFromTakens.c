#include "neighbourSearch.h"
#include "genericFunctions.h"
#include <math.h>     
#include <stdlib.h>

/*****************************************************************************/
 
 //eps is sorted in decreasing order
// q >1 
void generalizedCorrDimFromTakens(double *takens, int *numberTakens,  int *embedding,
int *q, double *eps,int* numberEps,int *numberBoxes,int *tdist,double* corrVector){
  // auxilar variables  
  int i,j,ep,posNeigh,lastTakens,lastReferenceVector,nfound;
  double denom;
  double epsMax=eps[0];
  // variables for the box assisted algorithm
  int boxes[(*numberBoxes)*(*numberBoxes)+1];
  int possibleNeighbours[*numberTakens];
  int neighList[*numberTakens];
  //box assisted algorithm
  boxAssistant( takens, numberTakens, embedding
  , &epsMax,numberBoxes, boxes, possibleNeighbours);
  //number of takens' vectors in the maximum embedding dimension
  lastTakens= (*numberTakens);
  lastReferenceVector = lastTakens-(*tdist);
  // iterate over the takens' vectors 
  for (i=*tdist;i<lastReferenceVector;i++){
    //find neighbours
    neighbourSearchFromBoxes(takens,&i,numberTakens, embedding,
    &epsMax,numberBoxes,boxes,possibleNeighbours,neighList,&nfound);
    // check neighbours of i
    for (j=0;j<nfound;j++){
      posNeigh=neighList[j];
      //check if we can use this possible neighbour 
      //we can not use it if this t(i) vector do not exist in the following embedded spaces
      //or if we have already count this neighbourhood
      //or they are close in time (do not respect the theiler window)
      if (abs(posNeigh-i)<=(*tdist)) continue;
      //ok: we may use this vector  
      corrVector[0]++;
      // complete the row corresponging to the minimum embedding dimension
      for (ep=1;ep<(*numberEps);ep++){
       if (isNeighbourTakens(i,posNeigh,takens, *numberTakens,*embedding, eps[ep]  )==1)
          corrVector[ep]++;
       else break;       
      }
    }
  }
  
  // denominator
  denom=(  (double)(lastTakens-2*(*tdist)) * pow( (double)((lastTakens-2*(*tdist))-1 ), (*q)-1  ) );
  for (ep=0;ep<(*numberEps);ep++){
     corrVector[ep] = pow( corrVector[ep], (*q)-1  )/denom ;
  }
 
}





