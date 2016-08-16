#ifndef RCPP_NEIGHBOUR_SEARCH_H
#define RCPP_NEIGHBOUR_SEARCH_H
#include <Rcpp.h>

class neighbour_search {
public:
  neighbour_search();
  neighbour_search(const Rcpp::NumericMatrix& phaseSpace, double radius, int numberBoxes);
  Rcpp::NumericMatrix get_phase_space() const;
  int get_dimension() const;
  int get_number_vectors() const;
  Rcpp::List find_all_neighbours() const;
  Rcpp::IntegerVector find_neighbours(int vectorIndex) const;
  bool are_neighbours(int vectorIndex1, int vectorIndex2, 
                     double neighbourhoodRadius) const;
  
private:
  Rcpp::NumericMatrix mPhaseSpace;
  int mEmbeddingDim;
  int mNumberVectors;
  double mRadius;
  Rcpp::IntegerVector mBoxes;
  Rcpp::IntegerVector mPossibleNeighbours;
  
  int get_wrapped_position(int row, int col) const;
  Rcpp::IntegerVector box_assisted_search(int vectorIndex,
                                          Rcpp::IntegerVector& neighbourWorkspace) const;
  
};

#endif
