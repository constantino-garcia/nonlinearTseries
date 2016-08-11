
#ifndef NEIGHS_H
#define NEIGHS_H

#include "genericFunctions.h"

#define OFFSET 0//auxiliar parameter

//computes the position of a takens' vector that falls into the (a,b)
//unwrapped box in the wrapped 2d box grid formed with "numberBoxes" in 
//each dimension
static inline int position(int a, int b,int numberBoxes){
  //a specifies the file, b the colum in the imaginary unwrapped 2D cover with boxes
  return ( 
             numberBoxes * ( MOD(a+OFFSET,numberBoxes) )
             +  MOD(b+OFFSET,numberBoxes)
    
         );
}


void boxAssistant( double *takens, int *numberTakens, int *embeddingD
, double *eps,int *numberBoxes,int *boxes,int *possibleNeighbours);

void neighbourSearchFromBoxes(double *takens,int *positionTakens, int *numberTakens, int *embeddingD
, double *eps,int *numberBoxes,int *boxes,int *possibleNeighbours,int *neighList,int *nfound);


#endif //NEIGHS_H








