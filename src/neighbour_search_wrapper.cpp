#include <Rcpp.h>
#include "neighbour_search.h"
using namespace Rcpp;

void transformCppIndexes(IntegerVector& indexes) {
  if (!Rf_isNull(indexes) && indexes.size() > 0) {
    std::transform(
      indexes.begin(), indexes.end(), indexes.begin(),
      std::bind(std::plus<int>(), std::placeholders::_1, 1)
    );  
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
  neighbour_search neighbourSearcher(phaseSpace, 
                                             radius, numberBoxes);
  IntegerVector neighbours = neighbourSearcher.find_neighbours(vectorIndex - 1);
  
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
  neighbour_search neighbourSearcher(phaseSpace, 
                                             radius, numberBoxes);
  List allNeighbours = neighbourSearcher.find_all_neighbours();
  // transform C++ indexes to R indexes
  for (int i = 0; i < allNeighbours.size(); i++){
    IntegerVector iNeighbours = as<IntegerVector>(allNeighbours[i]);
    transformCppIndexes(iNeighbours);
    allNeighbours[i] = iNeighbours;
  }
  return  allNeighbours;
}

