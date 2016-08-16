#include <Rcpp.h>
#include "NeighbourSearchAlgorithm.h"
using namespace Rcpp;

void transformCppIndexes(IntegerVector& indexes) {
  if (!Rf_isNull(indexes) && indexes.size() > 0) {
    std::transform(indexes.begin(), indexes.end(), indexes.begin(),
                   std::bind2nd(std::plus<int>(), 1));  
  }
}

/*
 * Search for all the neighbours of the vectorIndex-th phase space vector. This
 * function is only a wrapper around the methods provide by the 
 * NeighbourSearchAlgorihm class, which handles the "translation" from the 
 * R-like indexes (1:N) to the C++ indexes (0:N-1).
 */
// [[Rcpp::export]]
IntegerVector getVectorNeighbours(const NumericMatrix& phaseSpace,    
                                  int vectorIndex, double radius, int numberBoxes){
  NeighbourSearchAlgorithm neighbourSearcher(phaseSpace, 
                                             radius, numberBoxes);
  IntegerVector neighbours = neighbourSearcher.getVectorNeighbours(vectorIndex - 1);
  
  // transform C++ indexes to R indexes
  transformCppIndexes(neighbours);
  return neighbours;
}


/*
 * Search for the neighbours of each of the vectors in phase space. This
 * function is only a wrapper around the methods provide by the 
 * NeighbourSearchAlgorihm class, which handles the "translation" from the 
 * R-like indexes (1:N) to the C++ indexes (0:N-1).
 */
// [[Rcpp::export]]
List getAllNeighbours(const NumericMatrix& phaseSpace,    
                      double radius, int numberBoxes){
  NeighbourSearchAlgorithm neighbourSearcher(phaseSpace, 
                                             radius, numberBoxes);
  List allNeighbours = neighbourSearcher.getAllNeighbours();
  // transform C++ indexes to R indexes
  for (int i = 0; i < allNeighbours.size(); i++){
    IntegerVector iNeighbours = as<IntegerVector>(allNeighbours[i]);
    transformCppIndexes(iNeighbours);
    allNeighbours[i] = iNeighbours;
  }
  return  allNeighbours;
}

/*** R
ntakens = 500
nrepeat = 3
radius = 2
embeddingD = sample(2:10, 1)
itak = sample(1:ntakens, 1)
takens = matrix(rnorm(ntakens * embeddingD), nrow = ntakens)
number.boxes = 4

R = neighbourSearch(takens, itak, radius, number.boxes)
C = rcppNeighbourSearch(takens, itak, radius, number.boxes)
testthat::expect_equal(R,C)

library(microbenchmark)
print(
  microbenchmark(
    R = neighbourSearch(takens, itak, radius, number.boxes),
    C = rcppNeighbourSearch(takens, itak, radius, number.boxes),
    times = 100)
)

radius = 5.4
R = findAllNeighbours(takens, radius, number.boxes)
C = rcppFindAllNeighbours(takens, radius, number.boxes)
testthat::expect_equal(R,C)
*/
