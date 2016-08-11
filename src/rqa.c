#include "neighbourSearch.h"
#include "genericFunctions.h"  
#include <math.h>
#include <R.h> //include sort functionality


// neighs is a matrix with ntakens rows. The n-th row has the neighbours of the 
// n-th Takens' vector. n is also present as neighbour. The n-th position of the 
// nneighs vector has the number of neighs of the n-th takens' vector.
// verticalHistogram: vector with ntakens elements
void getVerticalHistogram(int* neighs, int* nneighs, int* ntakens, int *vmin, 
int* verticalHistogram){
  int i, j, count, lne;
  for (i = 0; i < (*ntakens); i++){
    verticalHistogram[i] = 0;
  }
  // find vertical lines in every column/row
  for (i = 0; i < (*ntakens); i++){
      // count number of neighbours + itself
      lne = nneighs[i];
      j=1;
      while (j < lne){
        count=1;
        // while there is a vertical line, update the length (count)
        while( (j < lne) && (MAT_ELEM(neighs,i,j,*ntakens) == (MAT_ELEM(neighs,i,(j-1),*ntakens)+1) ) ){
          j++;
          count++;
        }
        // update the histogram if the current length > vmin
        if (count >= (*vmin)){
          verticalHistogram[count-1]++;
        }
        j++;
      }
  } 
}  

int isContainedInNeighbourhood(int possibleNeigh, int i,int* neighs,int ntakens,int nneighs){
  int iter;
  int found = 0;
  for (iter = 0; iter < nneighs; iter++ ){
    if (possibleNeigh == MAT_ELEM(neighs,i,iter,ntakens)){
      found = 1;
      break;
    }
  }
  return found;
  
}



void updateLengthHistogram(int i,int j,int* neighs, int* nneighs, int ntakens,int* diagonalHistogram,int lmin){
  //update diagonalHistogram
  int actualLength=1;
  // advance in diagonal direction
  i++;
  j++;
  // advance while there exist diagonal 
  while ( ( j< ntakens ) &&  (isContainedInNeighbourhood(j,i,neighs,ntakens,nneighs[i]) == 1)  ){
      // update the length of the diagonal
      actualLength++;
      // advance
      j++;
      i++;
  }
  // update the Histogram
  if (actualLength >= lmin){
    diagonalHistogram[actualLength-1]++;
  }
}

//diagonalHistogram and recurrenceHistogram have ntakens elements

// compute the diagonal histogram and the recurrenceHistogram
void getDiagonalHistogramRecurrenceHistogram(int* neighs, int* nneighs, int* ntakens, int *lmin,
      int* diagonalHistogram,int* recurrenceHistogram){
  int i,j, currentNeigh, lastPosition;
  for (i=0;i <  (*ntakens) ;i++){
    diagonalHistogram[i]=0;
    recurrenceHistogram[i]=0;
  }
  // Treat the first row separately: It will be easier if we find diagonals that
  // "originate" at the first row separately
  for (j=0;j < nneighs[0]; j++){
    currentNeigh = MAT_ELEM(neighs,0,j,*ntakens);
    // update recurrenceHistogram: the recurrence histogram takes into accoutn
    // the distance between neigbours. In this case:currentNeigh-0, and we store
    //it in tre previous position (position 0 stores distance 1)
    recurrenceHistogram[currentNeigh-1]++;
    // update diagonal histogram
    updateLengthHistogram(0,currentNeigh,neighs,nneighs, *ntakens,diagonalHistogram,*lmin);
  }
  // find diagonals that "originate" at other rows
  lastPosition = *ntakens-*lmin;
  for (i=1;i < lastPosition; i++){
    for (j=0;j < nneighs[i]; j++){
       currentNeigh = MAT_ELEM(neighs,i,j,*ntakens);
      //check if we have already took into acount the pair i,j
      if (currentNeigh <= i) continue;
      // update recurrenceHistogram
      recurrenceHistogram[currentNeigh-i-1]++;
      //check if this diagonal has its origin in any row before
      if (isContainedInNeighbourhood(currentNeigh-1,i-1,neighs,*ntakens,nneighs[i-1]) == 1){
        continue;
      } else{
        //update the diagonalHistogram and the recurrenceHistogram
        updateLengthHistogram(i,currentNeigh,neighs,nneighs, *ntakens,diagonalHistogram,*lmin);  
      }      
    }
  }
  //complete the recurrence vector if needed
  if (lastPosition < (*ntakens)){
    for (i=lastPosition;i < (*ntakens); i++){
      for (j=0;j < nneighs[i]; j++){
        currentNeigh = MAT_ELEM(neighs,i,j,*ntakens);
         if (currentNeigh <= i) continue;
        // update recurrenceHistogram
        recurrenceHistogram[currentNeigh-i-1]++;
      }
    }
  }  
  
  for (i=0;i <  (*ntakens) ;i++){
    diagonalHistogram[i]=2*diagonalHistogram[i];
    recurrenceHistogram[i]=2*recurrenceHistogram[i];
  }
  diagonalHistogram[*ntakens-1] = 1;
}

/****************** getHistograms *******************************************/ 
void getHistograms(int* neighs, int* nneighs, int* ntakens, int *vmin, int *lmin,
int* verticalHistogram, int* diagonalHistogram,int* recurrenceHistogram){
  // auxiliar variables
  getVerticalHistogram(neighs, nneighs, ntakens, vmin, verticalHistogram);
  getDiagonalHistogramRecurrenceHistogram(neighs, nneighs, ntakens, lmin,
      diagonalHistogram, recurrenceHistogram);
      
} 


