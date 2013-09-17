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
double distance(int n1,int n2,double *takens, int numberTakens,int embeddingD);

//computes if vectors n1 and n2 in the takens array are neighbours.
//the neighbourhood is defined in a embeddingD-dimension space with
//the max norm and eps radious.
 int isNeighbourTakens(int n1,int n2,double *takens, int numberTakens,int embeddingD,double eps);

// scalar product
double scalarProduct (double* vector1, double* vector2, int len);

#endif //GENERIC_FUNCTIONS_H
