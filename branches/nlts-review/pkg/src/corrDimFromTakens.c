#include "neighbourSearch.h"
#include "genericFunctions.h"
#include <math.h>     
#include <stdlib.h>

//eps is sorted in decreasing order
void corrDimFromTakens(double *takens,int *lenTakens,int* embeddingDim, 
double *eps,int* numberEps,int *numberBoxes,int *tdist,double* corrVector){
  
  // auxilar variables
  int i,j,k,m,ep,posNeigh,numberEmbeddings,lastEpsTocheck,nfound;
  double distance;
  double epsMax=eps[0];
  //correlation Matrix;
  numberEmbeddings=1;
  // variables for the box assisted algorithm
  int boxes[(*numberBoxes)*(*numberBoxes)+1];
  int possibleNeighbours[*lenTakens];
  int neighList[*lenTakens];
  //box assisted algorithm
  boxAssistant( takens, lenTakens, embeddingDim
  , &epsMax,numberBoxes, boxes, possibleNeighbours);
  //total number of distances that we will use in each dimension
  double denom=(  ( (double)((*lenTakens)-(*tdist)) )*( (double)((*lenTakens)-(*tdist)+1) ) )/2.0;
  // iterate over the takens' vectors 
  for (i=0;i<(*lenTakens);i++){
    //find neighbours
    neighbourSearchFromBoxes(takens,&i,lenTakens, embeddingDim,
    &epsMax,numberBoxes,boxes,possibleNeighbours,neighList,&nfound);
    //check neighbours of i
    for (j=0;j<nfound;j++){
      posNeigh=neighList[j];
      //check if we can use this possible neighbour 
      //we can not use it if this t(i) vector do not exist in the following embedded spaces
      //or if we have already count this neighbourhood
      //or they are close in time (do not respect the theiler window)
      if ((posNeigh<i)||(abs(posNeigh-i)<=(*tdist))||(posNeigh>(*lenTakens))) continue;
      //ok: we may use this vector  
      corrVector[0] += 1.0/denom; 
      // complete the row 
      for (ep=1;ep<(*numberEps);ep++){
       if (isNeighbourTakens(i,posNeigh,takens, *lenTakens,*embeddingDim, eps[ep]  )==1)
                corrVector[ep] += 1.0/denom;
       else break;       
      }
    }
  }
}

