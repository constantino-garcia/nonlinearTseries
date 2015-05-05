#include <R.h>
#include <math.h> 
#include "genericFunctions.h"
#include "neighbourSearch.h"

//This function computes the accumulated histogram of the "takens" takens'
//vectors using a wrapped grid of "numberBoxes" per dimension. Each box
//has a size of "eps".
//takens: array of takens' vectors
//numberTakens: number of takens' vectors stored in takens
//embeddingD: dimension of the embedded space
//eps: size of each box of the grid.
//numberBoxes: number of boxes used for the grid in each direction
//boxes: each box of the array will point to all the vectors that fell on it
//possibleNeighbours: the numbers of each takens vector are stored here
// deppending on which box they fell. 
void boxAssistant( double *takens, int *numberTakens, int *embeddingD
, double *eps,int *numberBoxes,int *boxes,int *possibleNeighbours){
  //auxiliar variables
  int i,j;
  int lastPosition=(*embeddingD)-1;
  //the matrix has numberBoxes*numberBoxes + the last position to store the total length
  int lengthBoxes=(*numberBoxes)*(*numberBoxes)+1;
  int wrappedBoxPosition,xBoxPosition,yBoxPosition;
  //initialize the boxes 
  for (i=0;i<lengthBoxes;i++) boxes[i]=0;
  //count number of taken vectors in each box. We use the first and last coordinates
  //since they will probably be the least correlated ones
  for (i=0;i<(*numberTakens);i++) {
      xBoxPosition=(int)(MAT_ELEM(takens,i,0,*numberTakens)/(*eps));    
      yBoxPosition=(int)(MAT_ELEM(takens,i,lastPosition,*numberTakens)/(*eps));
      wrappedBoxPosition=position(xBoxPosition,yBoxPosition,*numberBoxes);
      boxes[wrappedBoxPosition]++;
  
  }
  //accumulate histogram
  for (i=1;i<lengthBoxes;i++) boxes[i]+=boxes[i-1];
  //fill list of pointers to possible neighbours
  for (i=0;i<(*numberTakens);i++){
    xBoxPosition=(int)(MAT_ELEM(takens,i,0,*numberTakens)/(*eps));    
    yBoxPosition=(int)(MAT_ELEM(takens,i,lastPosition,*numberTakens)/(*eps));
    wrappedBoxPosition=position(xBoxPosition,yBoxPosition,*numberBoxes);
    possibleNeighbours[--boxes[wrappedBoxPosition]]=i;
  }
}

//find neighbours of the "positionTakens"-th vector from the "takens" takens'
//vector array. The neighbours are found using a box assisted algorithm that
//creates a wrapped grid of "numberBoxes" per dimension. Each box
//has a size of "eps".
//takens: array of takens' vectors
//positionTakens: number of the vector which you want to calculate all of its neighbors
//numberTakens: number of takens' vectors stored in takens
//embeddingD: dimension of the embedded space
//eps: size of each box of the grid.
//numberBoxes: number of boxes used for the grid in each direction
//boxes: each box of the array will point to all the vectors that fell on it
//possibleNeighbours: the numbers of each takens vector are stored here
// deppending on which box they fell.
//neighList: neighbours of the positionTakens-th vector.
//nfound: number of neighbours of the positionTakens-th vector.
void neighbourSearchFromBoxes(double *takens,int *positionTakens, int *numberTakens, int *embeddingD
, double *eps,int *numberBoxes,int *boxes,int *possibleNeighbours,int *neighList,int *nfound){
  //auxiliar variable
  int xBoxPosition,yBoxPosition,lastPosition,i,j,auxiliarBoxPos,boxptr,possibleNeigh;
  lastPosition=(*embeddingD)-1;  
  xBoxPosition=(int)(MAT_ELEM(takens,*positionTakens,0,*numberTakens)/(*eps));    
  yBoxPosition=(int)(MAT_ELEM(takens,*positionTakens,lastPosition,*numberTakens)/(*eps));
  //printf("caigo en %d\n",position(xBoxPosition,yBoxPosition,*numberBoxes));
  //look for vector neighbours in the 9 neighbours of the (iBoxPos,jBoxPos) box
  for ((*nfound)=0,i=xBoxPosition-1;i<=xBoxPosition+1;i++){
    for (j=yBoxPosition-1;j<=yBoxPosition+1;j++){
       auxiliarBoxPos=position(i,j,*numberBoxes);
       //printf("boxPos %d\n",auxiliarBoxPos);
       //printf("voy de  %d a %d\n",boxes[auxiliarBoxPos+1]-1,boxes[auxiliarBoxPos]);
       //avoid empty boxes
       for (boxptr=boxes[auxiliarBoxPos+1]-1;boxptr>=boxes[auxiliarBoxPos];boxptr--){
        possibleNeigh=possibleNeighbours[boxptr];
        if (possibleNeigh==(*positionTakens)) continue;
        if (isNeighbourTakens(*positionTakens,possibleNeigh,takens,*numberTakens,*embeddingD,*eps)==1){
          neighList[(*nfound)++]=possibleNeigh;
        }
      }
    }
  }
}
