#ifndef GENERIC_FUNCTIONS_H
#define GENERIC_FUNCTIONS_H

#include <math.h>


#define MAT_ELEM(mat,i,j,nfiles) (mat[(j)*(nfiles)+(i)])
#define MAX(a,b) (a)>(b) ? (a) : (b)
#define MIN(a,b) (a)<(b) ? (a) : (b)


//this function modifies the C mod function in order to return only positives
//values
inline int MOD(int x,int m){
 return ( ((x%m) + m) %m);
}

//computes the max distance between the n1-th and n2-th takens vector.
//the neighbourhood is defined in a embeddingD-dimension space
inline double distance(int n1,int n2,double *takens, int numberTakens,int embeddingD){
  int i;
  double dis=0.0;
  double a,b;
  for (i=0;i<embeddingD;i++){
    a=MAT_ELEM(takens,n1,i ,numberTakens);
    b=MAT_ELEM(takens,n2,i,numberTakens);
    dis = MAX(dis , ( (double)(fabs(a-b)) ) );
  }
  return (dis);
}

//computes if vectors n1 and n2 in the takens array are neighbours.
//the neighbourhood is defined in a embeddingD-dimension space with
//the max norm and eps radious.
inline int isNeighbourTakens(int n1,int n2,double *takens, int numberTakens,int embeddingD,double eps){
  int i;
  double a,b;
  for (i=0;i<embeddingD;i++){
    a=MAT_ELEM(takens,n1,i ,numberTakens);
    b=MAT_ELEM(takens,n2,i,numberTakens);
    if ( ((double)(fabs(a-b)))>=eps) return 0;
  }
  return (1);
}


// scalar product
inline double scalarProduct (double* vector1, double* vector2, int len){
  double result = 0.0;
  int i;
  
  for (i=0;i<len;i++){
    result += vector1[i]*vector2[i];
  }
  return (result);
    
}

#endif //GENERIC_FUNCTIONS_H
