#include "neighbourSearch.h"
#include <math.h>
#include <stdlib.h>
/*this function computes the function S(N) ("Sdn"") for a given array of takens' 
* vectors "takens" (length of "takens" = "numberTakens") in a 
* "embeddingD"-dimensional space, and
* neighbourhood size "eps". This function will compute "Sdn" for n=0,1,...,"nmax".
* This functions truncates the computation of
* "Sdn" if "nminRP" reference points are found, each of which have at least 
* "neighMin" neighbours. Only such neighbous are considered which at least "tdist"
* time steps apart. This is done to avoid temporal correlations. 
* To calculate the neighbourhoods, a box assisted algorithm
* is used. The wrapped grid for the box assisted algorithm uses "numberBoxes" 
* per dimension.
*/
void maxLyapunov(double *timeSeries,double *takens, int* tau, int *numberTakens, int *embeddingD
, double *eps,double *Sdn,int *nmax, int *nminRP,int *neighMin,int *numberBoxes,int *tdist){
  // auxilar variables
  int i,ii,j,k,lastTakens,posNeigh,nfound,nf,rpfound;
  double distance0;
  double Saux[(*nmax)+1];
  // variables for the box assisted algorithm
  int boxes[(*numberBoxes)*(*numberBoxes)+1];
  int possibleNeighbours[*numberTakens];
  int neighList[*numberTakens];
  //box assisted algorithm
  boxAssistant( takens, numberTakens, embeddingD
  , eps,numberBoxes, boxes, possibleNeighbours);
  // initializing Sdn
  for (i=0;i<=(*nmax);i++) Sdn[i]=0;
  //last takens' vector that you can use
  lastTakens=(*numberTakens)-1-(*nmax);
  // iterate over the takens' vectors
  for (i=0,rpfound=0;(i<lastTakens)&&(rpfound<(*nminRP));i++){
    //inicializate Saux for averaging the divergence S[k] for the neighbours of i
    for (ii=0;ii<=(*nmax);ii++) Saux[ii]=0;
    //find neighbours
    neighbourSearchFromBoxes(takens,&i,numberTakens, embeddingD,
    eps,numberBoxes,boxes,possibleNeighbours,neighList,&nfound);
    //check neighbours of i
    for (nf=0,j=0;j<nfound;j++){
      posNeigh=neighList[j];
      //check if we can use this possible neighbour to average
      if (posNeigh>=lastTakens) continue;
      //avoid temporal correlations
      if ( abs(posNeigh-i)>=(*tdist) ){
        nf++;
        // calculate the divergence for every value of k 
        for (k=0;k<=(*nmax);k++) {
          distance0=fabs( (double)(timeSeries[ i+((*embeddingD)-1)*(*tau)]-timeSeries[posNeigh+((*embeddingD)-1)*(*tau)] ));
          if (distance0!=0)
            Saux[k]+=(fabs( (double)(timeSeries[ i+((*embeddingD)-1)*(*tau)+k]-
              timeSeries[posNeigh+((*embeddingD)-1)*(*tau)+k] ))/distance0);    
            
        }        
      }
      //check if there are enough neighbours
      if (nf>=(*neighMin)) {
        //found a reference point with enough neighbours
        rpfound++;
        for (k=0;k<=(*nmax);k++) Sdn[k]+=log(Saux[k]/nf);
      }
    }
  }
 
  // if we have found some reference point, average all refence points
  if (rpfound>0){
    for (k=0;k<=(*nmax);k++) Sdn[k]=Sdn[k]/rpfound;
  }
  
}
