#include "neighbourSearch.h"
#include "genericFunctions.h"
#include <math.h>     
#include <stdlib.h>
#include <R.h>

//eps is sorted in decreasing order
void corrDim(double *timeSeries,int *lenTimeSeries,double *takensDimMin, int* tau, int *numberTakens, int *minEmbeddingD
,int *maxEmbeddingD, double *eps,int* numberEps,int *numberBoxes,int *tdist,double* corrMatrix){
  
  // auxilar variables
  int i,j,m,ep,posNeigh,lastTakens,numberEmbeddings,lastEpsTocheck,nfound;
  double distance;
  double epsMax=eps[0];
  //correlation Matrix;
  numberEmbeddings=*maxEmbeddingD-*minEmbeddingD+1;
  //embedding vector
  int embedding[numberEmbeddings];
  embedding[0]=(*minEmbeddingD);
  for (i=1;i<numberEmbeddings;i++)  embedding[i]=embedding[i-1]+1;
  // variables for the box assisted algorithm
  int boxes[(*numberBoxes)*(*numberBoxes)+1];
  int possibleNeighbours[*numberTakens];
  int neighList[*numberTakens];
  //box assisted algorithm
  boxAssistant( takensDimMin, numberTakens, minEmbeddingD
  , &epsMax,numberBoxes, boxes, possibleNeighbours);
  //number of takens' vectors in the maximum embedding dimension
  lastTakens=(*lenTimeSeries)-(*maxEmbeddingD-1)*(*tau);
  //total number of distances that we will use in each dimension
  double denom=(  ( (double)(lastTakens-(*tdist)) )*( (double)(lastTakens-(*tdist)+1) ) )/2.0;
  // iterate over the takens' vectors 
  for (i=0;i<lastTakens;i++){
    //find neighbours
    neighbourSearchFromBoxes(takensDimMin,&i,numberTakens, minEmbeddingD,
    &epsMax,numberBoxes,boxes,possibleNeighbours,neighList,&nfound);
    //check neighbours of i
    for (j=0;j<nfound;j++){
      posNeigh=neighList[j];
      //check if we can use this possible neighbour 
      //we can not use it if this t(i) vector do not exist in the following embedded spaces
      //or if we have already count this neighbourhood
      //or they are close in time (do not respect the theiler window)
      if ((posNeigh<i)||(abs(posNeigh-i)<(*tdist))||(posNeigh>lastTakens)) continue;
      //ok: we may use this vector  
      MAT_ELEM(corrMatrix,0,0,numberEmbeddings)=MAT_ELEM(corrMatrix,0,0,numberEmbeddings)+1.0/denom; 
      // complete the row corresponging to the minimum embedding dimension
      for (ep=1;ep<(*numberEps);ep++){
       if (isNeighbourTakens(i,posNeigh,takensDimMin, *numberTakens,*minEmbeddingD, eps[ep]  )==1)
          MAT_ELEM(corrMatrix,0,ep,numberEmbeddings)=MAT_ELEM(corrMatrix,0,ep,numberEmbeddings)+1.0/denom;
       else break;       
      }
      // the rest of the embedding dimensions
      lastEpsTocheck=ep;
      //iterate over the whole matrix
      for (m=1;m<numberEmbeddings;m++){
        for (ep=0;ep<lastEpsTocheck;ep++){
          //we just have to chek the last position of the vector in the new embedding dimension (as they were neighbours in the
          //embedding[m-1] dimensional space)
          distance=(double)(fabs(timeSeries[i+(embedding[m]-1)*(*tau)]-timeSeries[posNeigh+(embedding[m]-1)*(*tau)]));
          if ( distance < eps[ep] ) {MAT_ELEM(corrMatrix,m,ep,numberEmbeddings)=MAT_ELEM(corrMatrix,m,ep,numberEmbeddings)+1.0/denom;}
          //there is  no need to check other eps in this dimension (because they are smaller)
          //Also, if i and posNeigh aren't neighbours in a m-dimensional space for a given eps,
          //nor will in a bigger dimension or smaller eps
          else {
            lastEpsTocheck=ep;
            break;
          }
        }  
      }
    }
  }
}

/*****************************************************************************/
 
 //eps is sorted in decreasing order
// q >1 
void generalizedCorrDim(double *timeSeries,int *lenTimeSeries,double *takensDimMin, int* tau, int *numberTakens, int *minEmbeddingD
,int *maxEmbeddingD, int *q, double *eps,int* numberEps,int *numberBoxes,int *tdist,double* corrMatrix){
  
  
  // auxilar variables  
  int i,j,m,ep,posNeigh,lastTakens,lastReferenceVector,numberEmbeddings,lastEpsTocheck,nfound;
  double distance, denom;
  double epsMax=eps[0];
  //correlation Matrix;
  numberEmbeddings=*maxEmbeddingD-*minEmbeddingD+1;
  //embedding vector
  int embedding[numberEmbeddings];
  embedding[0]=(*minEmbeddingD);
  for (i=1;i<numberEmbeddings;i++)  embedding[i]=embedding[i-1]+1;
  // variables for the box assisted algorithm
  int boxes[(*numberBoxes)*(*numberBoxes)+1];
  int possibleNeighbours[*numberTakens];
  int neighList[*numberTakens];
  //box assisted algorithm
  boxAssistant( takensDimMin, numberTakens, minEmbeddingD
  , &epsMax,numberBoxes, boxes, possibleNeighbours);
  //number of takens' vectors in the maximum embedding dimension
  lastTakens=(*lenTimeSeries)-(*maxEmbeddingD-1)*(*tau);
  lastReferenceVector = lastTakens-(*tdist);
  // iterate over the takens' vectors 
  for (i=*tdist;i<lastReferenceVector;i++){
    //find neighbours
    neighbourSearchFromBoxes(takensDimMin,&i,numberTakens, minEmbeddingD,
    &epsMax,numberBoxes,boxes,possibleNeighbours,neighList,&nfound);
    // check neighbours of i
    for (j=0;j<nfound;j++){
      posNeigh=neighList[j];
      //check if we can use this possible neighbour 
      //we can not use it if this t(i) vector do not exist in the following embedded spaces
      //or if we have already count this neighbourhood
      //or they are close in time (do not respect the theiler window)
      if ((abs(posNeigh-i)<(*tdist))||(posNeigh>lastTakens)) continue;
      //ok: we may use this vector  
      MAT_ELEM(corrMatrix,0,0,numberEmbeddings)=MAT_ELEM(corrMatrix,0,0,numberEmbeddings)+1.0; 
      // complete the row corresponging to the minimum embedding dimension
      for (ep=1;ep<(*numberEps);ep++){
       if (isNeighbourTakens(i,posNeigh,takensDimMin, *numberTakens,*minEmbeddingD, eps[ep]  )==1)
          MAT_ELEM(corrMatrix,0,ep,numberEmbeddings)=MAT_ELEM(corrMatrix,0,ep,numberEmbeddings)+1.0;
       else break;       
      }
      // the rest of the embedding dimensions
      lastEpsTocheck=ep;
      //iterate over the whole matrix
      for (m=1;m<numberEmbeddings;m++){
        for (ep=0;ep<lastEpsTocheck;ep++){
          //we just have to chek the last position of the vector in the new embedding dimension (as they were neighbours in the
          //embedding[m-1] dimensional space)
          distance=(double)(fabs(timeSeries[i+(embedding[m]-1)*(*tau)]-timeSeries[posNeigh+(embedding[m]-1)*(*tau)]));
          if ( distance < eps[ep] ) {MAT_ELEM(corrMatrix,m,ep,numberEmbeddings)=MAT_ELEM(corrMatrix,m,ep,numberEmbeddings)+1.0;}
          //there is  no need to check other eps in this dimension (because they are smaller)
          //Also, if i and posNeigh aren't neighbours in a m-dimensional space for a given eps,
          //nor will in a bigger dimension or smaller eps
          else {
            lastEpsTocheck=ep;
            break;
          }
        }  
      }
    }
  }
  
  // denominator
  denom=(  (double)(lastTakens-2*(*tdist)) * pow( (double)((lastTakens-2*(*tdist))-1 ), (*q)-1  ) );
  for (m=0;m<numberEmbeddings;m++){
    for (ep=0;ep<(*numberEps);ep++){
      MAT_ELEM(corrMatrix,m,ep,numberEmbeddings) = pow( (double)MAT_ELEM(corrMatrix,m,ep,numberEmbeddings), (*q)-1  )/denom ;
    }
  }
}






