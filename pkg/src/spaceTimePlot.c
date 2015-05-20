#include "genericFunctions.h"
#include <math.h>


//takens: array of takens' vectors
//numberTakens: number of takens' vectors stored in takens
//embeddingD: dimension of the embedded space
//eps: size of each box of the grid.
void spaceTimePlot( double *takens, int *numberTakens, int *embeddingD
,double *eps,int* leps, int* numberTimeSteps,int* timeStep,int*  numberPercentages,double* spaceTimePlot){
  //auxiliar variables
  int time,i,ntakens, positionEps, lastPosition, percentageIterator,npoints,neededPoints,firstPosition;
  int histEps[*leps];
  double dis;
  lastPosition = (*leps)-1;
  double maxEps=eps[lastPosition];
  
  for (time=0; time < (*numberTimeSteps);time++){
    // inicialice the histogram of distances
    for(i=0;i<(*leps);i++) histEps[i]=0;
    // calculate distances for the current timeStep = time*timeStep
    firstPosition = ( (*numberTimeSteps)*(*timeStep) ); // first position we can use
    for (ntakens=firstPosition;ntakens<(*numberTakens);ntakens++){
      dis = distance(ntakens,ntakens-time*(*timeStep),takens,*numberTakens,*embeddingD);
      // update the histogram of distances
      positionEps = (int)( (dis/maxEps)*((double)(*leps)) );
      histEps[MIN(positionEps,lastPosition)]++;
    }
    for ( percentageIterator=0; percentageIterator<(* numberPercentages); percentageIterator++){
      neededPoints =  (*numberTakens-firstPosition)*(percentageIterator+1)/(*numberPercentages);
      for (npoints=0,i=0;(i<(*leps))&&(npoints<neededPoints);i++)
        npoints += histEps[i];
      //store the resulting eps at which we obtain the desired number of neighbours for a given 
      //time delay
      MAT_ELEM(spaceTimePlot,percentageIterator,time,*numberPercentages) = eps[i];
    }
  }
  
}
