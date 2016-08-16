#ifndef RCPP_NEIGHBOUR_SEARCH_H
#define RCPP_NEIGHBOUR_SEARCH_H
#include <Rcpp.h>

class NeighbourSearchAlgorithm {
private:
  Rcpp::NumericMatrix mPhaseSpace;
  int mEmbeddingDim;
  int mNumberVectors;
  double mRadius;
  Rcpp::IntegerVector mBoxes;
  Rcpp::IntegerVector mPossibleNeighbours;
  
  int getWrappedBoxPosition(int row, int col) const;
  Rcpp::IntegerVector boxAssistedNeighbourSearch(int vectorIndex,
                                           Rcpp::IntegerVector& neighbourWorkspace) const;
public:
  NeighbourSearchAlgorithm();
  NeighbourSearchAlgorithm(const Rcpp::NumericMatrix& phaseSpace, double radius, int numberBoxes);
  Rcpp::NumericMatrix getPhaseSpace() const;
  int getEmbeddingDim() const;
  int getNumberVectors() const;
  Rcpp::IntegerVector getBoxes() const;
  Rcpp::List getAllNeighbours() const;
  Rcpp::IntegerVector getVectorNeighbours(int vectorIndex) const;
  bool areNeighbours(int vectorIndex1, int vectorIndex2, double neighbourhoodRadius) const;
};





#endif
