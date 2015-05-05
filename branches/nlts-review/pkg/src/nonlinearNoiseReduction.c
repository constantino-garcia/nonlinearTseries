#include <R.h>
#include <math.h> 
#include "genericFunctions.h"
#include "neighbourSearch.h"

//This function computes the accumulated histogram of the "takens" takens'
//vectors using a wrapped grid of "numberBoxes" per dimension. Each box
//has a size of "eps".
//takens: array of takens' vectors build from the time series using tau=1
//numberTakens: number of takens' vectors stored in takens
//embeddingD: dimension of the embedded space
//eps: size of each box of the grid.
//numberBoxes: number of boxes used for the grid in each direction
//boxes: each box of the array will point to all the vectors that fell on it
//possibleNeighbours: the numbers of each takens vector are stored here
// deppending on which box they fell. 
void nonlinearNoiseReduction(double *timeSeries, double *takens, int *numberTakens, int *embeddingD
, double *eps,int *numberBoxes){
  int i,j, posNeigh, midPos, nfound;
  // variables for the box assisted algorithm
  int boxes[(*numberBoxes)*(*numberBoxes)+1];
  int possibleNeighbours[*numberTakens];
  int neighList[*numberTakens];
  //box assisted algorithm
  boxAssistant( takens, numberTakens, embeddingD
  , eps,numberBoxes, boxes, possibleNeighbours);
  
  //the neighbours of the Sn takens' vector will be used to denoise x(n+midPos)...
  midPos = (int) floor(   ((double)(*embeddingD))/2.0);
  for (i = 0; i < (*numberTakens);i++){
    neighbourSearchFromBoxes(takens,&i,numberTakens, embeddingD,
    eps,numberBoxes,boxes,possibleNeighbours,neighList,&nfound);
    //no neigbours found...
    if (nfound==0) continue;
    //there exist a neighbourhood...denoise using this neihgbourhood
    for (j=0;j < nfound;j++){
      posNeigh=neighList[j];
      timeSeries[i+midPos] = timeSeries[i+midPos] + MAT_ELEM(takens,posNeigh,midPos,*numberTakens);
    }
    //mean using the neighbours and the proper point
    timeSeries[i+midPos] = timeSeries[i+midPos]/(nfound+1);
  }

}
