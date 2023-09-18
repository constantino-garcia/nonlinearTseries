#include <Rcpp.h>
using namespace Rcpp;
#include "generic_functions.h"  


// neighs is a matrix with ntakens rows. The n-th row has the neighbours of the 
// n-th Takens' vector. n is also present as neighbour. The n-th position of the 
// nneighs vector has the number of neighs of the n-th takens' vector.
// verticalHistogram: vector with ntakens elements
void get_vertical_histogram(List& neighs,  int ntakens, int vmin, 
                          IntegerVector& verticalHistogram){
  // find vertical lines in every column/row
  for (int i = 0; i < ntakens; i++){
    // count number of neighbours + itself
    IntegerVector ith_takens_neighs = as<IntegerVector>(neighs[i]);
    int lne = ith_takens_neighs.length();
    int j = 1;
    while (j < lne){
      int count = 1;
      // while there is a vertical line, update the length (count)
      while( (j < lne) && (ith_takens_neighs[j] == (ith_takens_neighs[j - 1] + 1)) ) {
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

bool is_contained_in_neighbourhood(int possibleNeigh, int i, List& neighs){
  bool found = false;
  IntegerVector ith_takens_neighs = as<IntegerVector>(neighs[i]);
  int nneighs = ith_takens_neighs.length();
  for (int iter = 0; iter < nneighs; iter++ ){
    if (possibleNeigh == ith_takens_neighs[iter]){
      found = true;
      break;
    }
  }
  return found;
}

void update_length_histogram(int i, int j, List& neighs, int ntakens,
                           IntegerVector& diagonalHistogram, int lmin){
  // update diagonalHistogram
  int actualLength = 1;
  // advance in diagonal direction
  i++;
  j++;
  // advance while there exist diagonal
  while ( (j < ntakens) &&  (is_contained_in_neighbourhood(j, i, neighs)) ){
    // update the length of the diagonal
    actualLength++;
    // advance
    j++;
    i++;
  }
  // update the Histogram
  if (actualLength >= lmin){
    // diagonal lengths of 1 in position 0, 2 in position 1, etc.
    // +=2 to take into account that the matrix is symmetric
    diagonalHistogram[actualLength - 1] += 2;
  }
}

//diagonalHistogram and recurrenceHistogram have ntakens elements
// compute the diagonal histogram and the recurrenceHistogram
void get_diagonal_recurrence_histogram(
    List& neighs, int ntakens, int lmin,
    IntegerVector& diagonalHistogram, IntegerVector& recurrenceHistogram){
  // Treat the first row separately: It will be easier if we find diagonals that
  // "originate" at the first row separately
  IntegerVector first_takens_neighs = as<IntegerVector>(neighs[0]);
  int nneighs = first_takens_neighs.length();
  for (int j = 0; j < nneighs; j++){
    int currentNeigh = first_takens_neighs[j];
    // update recurrenceHistogram: the recurrence histogram takes into account
    // the distance between neigbours. In this case: currentNeigh - 0, and we store
    // it in the previous position (position 0 stores distance 1).
    // +=2 to take into account that the RQA matrix is symmetric
    recurrenceHistogram[currentNeigh - 1] += 2;
    // update diagonal histogram
    update_length_histogram(0, currentNeigh, neighs,
                            ntakens, diagonalHistogram, lmin);
  }
  // find diagonals that "originate" at other rows
  int lastPosition = ntakens - lmin;
  for (int i = 1; i < lastPosition; i++){
    IntegerVector ith_takens_neighs = as<IntegerVector>(neighs[i]);
    nneighs = ith_takens_neighs.length();
    for (int j = 0; j < nneighs; j++){
      int currentNeigh = ith_takens_neighs[j];
      //check if we have already took into account the pair i,j
      if (currentNeigh <= i) continue;
      // +=2 to take into account that the RQA matrix is symmetric
      recurrenceHistogram[currentNeigh - i - 1] += 2;
      //check if this diagonal has its origin in any row before
      if (is_contained_in_neighbourhood(currentNeigh - 1, i - 1, neighs)) {
        continue;
      } else{
        //update the diagonalHistogram and the recurrenceHistogram
        update_length_histogram(i,currentNeigh,neighs, ntakens,
                              diagonalHistogram, lmin);
      }
    }
  }
  //complete the recurrence vector if needed
  if (lastPosition < ntakens){
    for (int i = lastPosition; i < ntakens; i++){
      IntegerVector ith_takens_neighs = as<IntegerVector>(neighs[i]);
      nneighs = ith_takens_neighs.length();
      for (int j = 0; j < nneighs; j++){
        int currentNeigh = ith_takens_neighs[j];
        if (currentNeigh <= i) continue;
        // +=2 to take into account that the matrix is symmetric
        recurrenceHistogram[currentNeigh - i - 1] += 2;
      }
    }
  }
  // Add the main diagonal, of length ntakens (there is only one!)
  diagonalHistogram[ntakens - 1] = 1;
}

// [[Rcpp::export]]
List get_rqa_histograms(List& neighs, int ntakens, int vmin, int lmin){
  // auxiliar variables
  IntegerVector verticalHistogram(ntakens, 0);
  IntegerVector diagonalHistogram(ntakens, 0);
  IntegerVector recurrenceHistogram(ntakens, 0);
  
  get_vertical_histogram(neighs, ntakens, vmin, verticalHistogram);
  get_diagonal_recurrence_histogram(neighs, ntakens, lmin,
                                    diagonalHistogram, recurrenceHistogram);
 List ret;
 ret["diagonalHist"] = diagonalHistogram;
 ret["recurrenceHist"] = recurrenceHistogram;
 ret["verticalHist"] = verticalHistogram;
 return ret; 
}  

