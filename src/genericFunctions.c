#include <R.h>
#include <math.h> 
#include "genericFunctions.h"


double distance(int n1,int n2,double *takens, int numberTakens,int embeddingD){
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


int isNeighbourTakens(int n1,int n2,double *takens, int numberTakens,int embeddingD,double eps){
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
double scalarProduct (double* vector1, double* vector2, int len){
  double result = 0.0;
  int i;
  
  for (i=0;i<len;i++){
    result += vector1[i]*vector2[i];
  }
  return (result);
    
}
