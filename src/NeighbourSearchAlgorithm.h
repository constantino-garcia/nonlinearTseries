#ifndef NEIGHBOUR_SEARCH_H
#define NEIGHBOUR_SEARCH_H
#include <Rcpp.h>

using namespace Rcpp;
class NeighbourSearchAlgorithm {
private:
  const Rcpp::NumericMatrix& mPhaseSpace;
  int mEmbeddingDim;
  int mNumberVectors;
  double mRadius;
  Rcpp::IntegerVector mBoxes;
  Rcpp::IntegerVector mPossibleNeighbours;
  
  int getWrappedBoxPosition(int row, int col) const;
  bool areNeighbours(const int vectorIndex1, const int vectorIndex2) const;
  IntegerVector boxAssistedNeighbourSearch(int vectorIndex,
                                           IntegerVector& neighbourWorkspace) const;
public:
  NeighbourSearchAlgorithm(const Rcpp::NumericMatrix& phaseSpace, double radius, int numberBoxes);
  Rcpp::IntegerVector getBoxes() const;
  Rcpp::List getAllNeighbours() const;
  Rcpp::IntegerVector getVectorNeighbours(int vectorIndex) const;
};





#endif
