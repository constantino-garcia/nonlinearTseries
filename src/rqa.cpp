#include <Rcpp.h>
using namespace Rcpp;
#include "generic_functions.h"  


// neighs is a matrix with ntakens rows. The n-th row has the neighbours of the 
// n-th Takens' vector. n is also present as neighbour. The n-th position of the 
// nneighs vector has the number of neighs of the n-th takens' vector.
// verticalHistogram: vector with ntakens elements
void get_vertical_histogram(IntegerMatrix& neighs, IntegerVector& nneighs, 
                          int ntakens, int vmin, 
                          IntegerVector& verticalHistogram){
  for (int i = 0; i < ntakens; i++){
    verticalHistogram[i] = 0;
  }
  // find vertical lines in every column/row
  for (int i = 0; i < ntakens; i++){
    // count number of neighbours + itself
    int lne = nneighs[i];
    int j=1;
    while (j < lne){
      int count = 1;
      // while there is a vertical line, update the length (count)
      while( (j < lne) && (neighs(i, j) == (neighs(i, j - 1) + 1)) ){
        j++;
        count++;
      }
      // update the histogram if the current length > vmin
      if (count >= vmin){
        // lengths of 1 in position 0, 2 in position 1, etc.
        verticalHistogram[count - 1]++;
      }
      j++;
    }
  } 
}  

bool is_contained_in_neighbourhood(int possibleNeigh, int i, IntegerMatrix& neighs,
                                int ntakens, int nneighs){
  bool found = false;
  for (int iter = 0; iter < nneighs; iter++ ){
    if (possibleNeigh == neighs(i,iter)){
      found = true;
      break;
    }
  }
  return found;
  
}

void update_length_histogram(int i, int j, IntegerMatrix& neighs,
                           IntegerVector& nneighs, int ntakens,
                           IntegerVector& diagonalHistogram, int lmin){
  // update diagonalHistogram
  int actualLength = 1;
  // advance in diagonal direction
  i++;
  j++;
  // advance while there exist diagonal
  while ( (j < ntakens) &&  (is_contained_in_neighbourhood(j, i, neighs, ntakens, nneighs[i])) ){
    // update the length of the diagonal
    actualLength++;
    // advance
    j++;
    i++;
  }
  // update the Histogram
  if (actualLength >= lmin){
    // diagonal lengths of 1 in position 0, 2 in position 1, etc.
    diagonalHistogram[actualLength - 1]++;
  }
}

//diagonalHistogram and recurrenceHistogram have ntakens elements
// compute the diagonal histogram and the recurrenceHistogram
void get_diagonal_recurrence_histogram(
    IntegerMatrix& neighs, IntegerVector& nneighs, int ntakens, int lmin,
    IntegerVector& diagonalHistogram, IntegerVector& recurrenceHistogram){
  for (int i = 0; i < ntakens ;i++){
    diagonalHistogram[i]=0;
    recurrenceHistogram[i]=0;
  }
  // Treat the first row separately: It will be easier if we find diagonals that
  // "originate" at the first row separately
  for (int j = 0; j < nneighs[0]; j++){
    int currentNeigh = neighs(0,j);
    // update recurrenceHistogram: the recurrence histogram takes into account
    // the distance between neigbours. In this case: currentNeigh - 0, and we store
    // it in tre previous position (position 0 stores distance 1)
    recurrenceHistogram[currentNeigh - 1]++;
    // update diagonal histogram
    update_length_histogram(0, currentNeigh, neighs, nneighs, 
                          ntakens, diagonalHistogram, lmin);
  }
  // find diagonals that "originate" at other rows
  int lastPosition = ntakens - lmin;
  for (int i = 1; i < lastPosition; i++){
    for (int j = 0; j < nneighs[i]; j++){
      int currentNeigh = neighs(i,j);
      //check if we have already took into acount the pair i,j
      if (currentNeigh <= i) continue;
      // update recurrenceHistogram
      recurrenceHistogram[currentNeigh - i - 1]++;
      //check if this diagonal has its origin in any row before
      if (is_contained_in_neighbourhood(currentNeigh - 1, i - 1, neighs,
                                     ntakens, nneighs[i-1])) {
        continue;
      } else{
        //update the diagonalHistogram and the recurrenceHistogram
        update_length_histogram(i,currentNeigh,neighs,nneighs, ntakens,
                              diagonalHistogram, lmin);
      }
    }
  }
  //complete the recurrence vector if needed
  if (lastPosition < ntakens){
    for (int i = lastPosition; i < ntakens; i++){
      for (int j = 0; j < nneighs[i]; j++){
        int currentNeigh = neighs(i,j);
        if (currentNeigh <= i) continue;
        // update recurrenceHistogram
        recurrenceHistogram[currentNeigh - i - 1]++;
      }
    }
  }
  
  // take into account that the matrix is symmetric
  for (int i = 0; i <  ntakens; i++){
    diagonalHistogram[i] = 2 * diagonalHistogram[i];
    recurrenceHistogram[i] = 2 * recurrenceHistogram[i];
  }
  // Add the main diagonal, of length ntakens
  diagonalHistogram[ntakens - 1] = 1;
}

// [[Rcpp::export]]
List get_rqa_histograms(IntegerMatrix& neighs, IntegerVector& nneighs, 
                         int ntakens, int vmin, int lmin){
  // auxiliar variables
  IntegerVector verticalHistogram(ntakens, 0.0);
  IntegerVector diagonalHistogram(ntakens, 0.0);
  IntegerVector recurrenceHistogram(ntakens, 0.0);
  
  get_vertical_histogram(neighs, nneighs, ntakens, vmin, verticalHistogram);
  get_diagonal_recurrence_histogram(neighs, nneighs, ntakens, lmin,
                                    diagonalHistogram, recurrenceHistogram);
 List ret;
 ret["diagonalHist"] = diagonalHistogram;
 ret["recurrenceHist"] = recurrenceHistogram;
 ret["verticalHist"] = verticalHistogram;
 return ret; 
}  

