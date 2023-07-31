#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void neighsList2Sparse(List&  neighs,NumericMatrix& neighs_matrix) {
  NumericVector x;
  int neighs_matrix_nrow = neighs_matrix.nrow();
  for (int i = 0; i < neighs_matrix_nrow; i++){
    x = Rcpp::as<NumericVector>(neighs[i]); 
    x.sort();
    // add diagonal elements
    if ( x.size() > 0){
      int j = 0;
      int final_size = (x.size() + 1);
      // insert the diagonal neighbour in order
      while ( ( j < x.size() ) && ((x(j) - 1) < i)  ){
        // substract one to get C indexes
        neighs_matrix( i, j ) = x(j) -1;
        j++;
      }
      
      neighs_matrix( i, j++ ) = i;
      
      for(; j < final_size; j++){
        /* substract one to the position to take into account that we
        we have introduced the diagonal element */
        neighs_matrix( i, j ) = x(j-1) -1;
      }
    }else{
      neighs_matrix(i, 0) = i;  
    }
  }
}




// [[Rcpp::export]]
void neighsList2SparseRCreator(const List&  neighs,const int& ntakens,
NumericMatrix& neighs_matrix) {
  NumericVector x;
  int irow = 0;
  for (int i = 0;i < ntakens; i++){
    x = Rcpp::as<NumericVector>(neighs[i]); 
    // i should start in 1 and ends in Takens as required by the R sparseMatrix
    // object constructor
    neighs_matrix(irow,0) = i + 1;
    neighs_matrix(irow++,1) = i + 1;
    for (int j = 0; j < x.size(); j++){
      if ((i + 1) < x(j)) {
        neighs_matrix(irow,0) = i + 1;
        neighs_matrix(irow++,1) = x(j);
      }
    }
  }
}
