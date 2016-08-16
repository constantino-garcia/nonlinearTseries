#ifndef RCPP_CORR_DIM_H
#define RCPP_CORR_DIM_H

#include <Rcpp.h>
#include <vector>
#include "NeighbourSearchAlgorithm.h"


class CorrDimAlgorithm {
private:
  NeighbourSearchAlgorithm mNeighbourSearcher;
  Rcpp::NumericVector mTimeSeries;
  int mTimeLag;
  int mTheilerDistance;
  Rcpp::NumericVector mRadiusVector;
  int mMinEmbeddingDim;
  int mMaxEmbeddingDim;
  /* Maximum number of Takens vector that can be used in ALL dimensions */
  int mNumberTakens;
  /* Index of the first Takens vector that can be used as reference vector */
  int mFirstReferenceVector;
  /* number of neighbours depending on the radius and the embedding dimension
   * for each reference vector in phase space
   */
  std::vector<Rcpp::NumericMatrix> mNumberOfNeighbors;

  void calculateCorrMatrix();
  void accumulateWeightedProbability(Rcpp::NumericMatrix& accumulated,
                                     const Rcpp::NumericMatrix& x, double exponent);
  void updateNeighboursMatrix(Rcpp::NumericMatrix& currentNeighbourMatrix, 
                                                int refVectorIndex);
public:
  CorrDimAlgorithm(const Rcpp::NumericVector& timeSeries,
                   int timeLag,
                   int theilerDistance,
                   Rcpp::NumericVector& radiusVector,
                   int minEmbeddingDim, int maxEmbeddingDim,
                   int numberBoxes);
  Rcpp::NumericMatrix getCorrSum(int order);
};

#endif
