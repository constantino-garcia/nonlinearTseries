#ifndef RCPP_NEIGHBOUR_SEARCH_H
#define RCPP_NEIGHBOUR_SEARCH_H
#include <Rcpp.h>
#include <vector>

class neighbour_search {
public:
  neighbour_search();
  neighbour_search(const Rcpp::NumericMatrix& phaseSpace, double radius,
                   int numberBoxes);
  void set_radius(double radius);
  Rcpp::NumericMatrix get_phase_space() const;
  int get_dimension() const;
  int get_number_vectors() const;
  Rcpp::List find_all_neighbours();
  Rcpp::List find_all_neighbours(int theilerWindow);
  Rcpp::IntegerVector find_neighbours(int vectorIndex);
  Rcpp::IntegerVector find_neighbours(int vectorIndex, int theilerWindow);
  double calculate_max_distance(int vectorIndex1, int vectorIndex2) const;
  bool are_neighbours(int vectorIndex1, int vectorIndex2, 
                      double neighbourhoodRadius) const;
  static bool comply_theiler_window(int vectorIndex1, int vectorIndex2, int theilerWindow);
  
  
private:
  Rcpp::NumericMatrix mPhaseSpace;
  int mEmbeddingDim;
  int mNumberVectors;
  double mRadius;
  Rcpp::IntegerVector mSearchingWorkspace;
  Rcpp::IntegerVector mBoxes;
  Rcpp::IntegerVector mPossibleNeighbours;
  
  int get_wrapped_position(int row, int col) const;
  void prepare_box_assisted_search();
  Rcpp::IntegerVector box_assisted_search(int vectorIndex,
                                          int theilerWindow);
  
};

#endif
