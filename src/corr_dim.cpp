#include <Rcpp.h>
#include "neighbour_search.h"
#include "generic_functions.h"
using namespace Rcpp;

struct correlation_sum_information {
  const NumericVector& mTimeSeries; 
  int mTimeLag;
  int mTheilerWindow;
  NumericVector& mRadii;
  int mMinEmbeddingDim;
  int mMaxEmbeddingDim;
  int mCorrSumOrder;
  int mNumTakens;
  int mNumRefVectors;
  
  correlation_sum_information(const NumericVector& timeSeries,  int timeLag,
                              int theilerWindow, NumericVector& radii,
                              int minEmbeddingDim, int maxEmbeddingDim,
                              int corrSumOrder, int nTakens, int nRefVectors):
    mTimeSeries(timeSeries),
    mTimeLag(timeLag),
    mTheilerWindow(theilerWindow),
    mRadii(radii),
    mMinEmbeddingDim(minEmbeddingDim),
    mMaxEmbeddingDim(maxEmbeddingDim),
    mCorrSumOrder(corrSumOrder),
    mNumTakens(nTakens),
    mNumRefVectors(nRefVectors){}
};

void count_neighbours(NumericMatrix& currentNeighbourCount,
                      neighbour_search& neighbourSearcher,
                      int refVectorIndex,
                      correlation_sum_information& corrSumInfo){
  int currentNeighbourCountNumberCol = currentNeighbourCount.ncol();
  int currentNeighbourCountNumberRow = currentNeighbourCount.nrow();
  IntegerVector neighboursIndexes = 
    neighbourSearcher.find_neighbours(refVectorIndex);
  int neighboursIndexesSize = neighboursIndexes.size();
  /* check each of the neighbours of refVectorIndex */
  for (int j = 0; j < neighboursIndexesSize; j++){
    int neighbourIndex = neighboursIndexes[j];
    /* Check if we can use this possible neighbour
     * we can not use it if this t(i) vector do not exist in the following 
     * embedded spaces or if we have already count this neighbourhood
     * or they are close in time (do not respect the theiler window)
     */
    bool isCloserTheilerWindow = 
      std::abs(neighbourIndex - refVectorIndex) <= corrSumInfo.mTheilerWindow;
    bool isInvalidTakensVector = neighbourIndex >= corrSumInfo.mNumTakens;
    if (isCloserTheilerWindow || isInvalidTakensVector) {
      continue;
    }
    /* OK: we may use this vector */
    currentNeighbourCount(0,0)++;
    /* Complete the row corresponging to the minimum embedding dimension. */
    int ep;
    for (ep = 1; ep < currentNeighbourCountNumberCol; ep++){
      if (neighbourSearcher.are_neighbours(refVectorIndex, neighbourIndex, 
                                           corrSumInfo.mRadii[ep])){
        currentNeighbourCount(0,ep)++;
      } else {
        break;
      }
    }
    int lastRadiusTocheck = ep;
    for (int m=1; m < currentNeighbourCountNumberRow; m++){
      int embeddingDim = corrSumInfo.mMinEmbeddingDim + m;
      for (ep = 0; ep < lastRadiusTocheck; ep++){
        /* We just have to chek the last position of the vector in the new
         * embedding dimension (as they were neighbours in the
         * mEmbeddingVector[m-1] dimensional space)
         */
        double distance =  std::abs(
          corrSumInfo.mTimeSeries[refVectorIndex + (embeddingDim - 1) * corrSumInfo.mTimeLag] -
            corrSumInfo.mTimeSeries[neighbourIndex + (embeddingDim - 1) * corrSumInfo.mTimeLag]
        );
        if (distance < corrSumInfo.mRadii[ep]) {
          currentNeighbourCount(m,ep)++;
        }else {
          /* There is  no need to check other radius in this dimension since 
           * they are sorted in descending order. Also, if refVectorIndex and
           * neighbourIndex  aren't neighbours in a m-dimensional space for a 
           * given radius, nor will in a bigger dimension or smaller radius.
           */
          lastRadiusTocheck = ep;
          break;
        }
      } /* end looping different radius */
    } /* end looping different embeddings */
  } /* end looping neighbours */
}



void update_correlation_sum(NumericMatrix& corrSum,
                            NumericMatrix& currentNeighbourCount,
                            double exponent){
  int nrow = corrSum.nrow();
  int ncol = corrSum.ncol();
  
  for (int i=0; i < nrow; i++) {
    for (int j=0; j < ncol; j++){
      corrSum(i,j) += std::pow(currentNeighbourCount(i,j), exponent);
    }
  }
}



void  calculate_weighted_neighbour_count(NumericMatrix& corrSum, 
                                         neighbour_search& neighbourSearcher,
                                         correlation_sum_information& corrSumInfo){
  int nDimensions = corrSumInfo.mMaxEmbeddingDim - corrSumInfo.mMinEmbeddingDim + 1;
  int nRadii = corrSumInfo.mRadii.size();
  double exponent = static_cast<double>(corrSumInfo.mCorrSumOrder - 1);
  
  /* The case corrSumOder == 2 permits to save computation time since we can
   * accumulate directly the neighbour count of each reference vector, without
   *  the use of temporal matrices (we can accumulate directly on the final
   *  correlation sum).
   */
  if (corrSumInfo.mCorrSumOrder == 2) {
    for (int refVectorNumber = 0, refVectorIndex = corrSumInfo.mTheilerWindow;
         refVectorNumber < corrSumInfo.mNumRefVectors;
         refVectorNumber++, refVectorIndex++) {
      count_neighbours(corrSum, neighbourSearcher,
                       refVectorIndex, corrSumInfo);
      
    }
  } else {
    for (int refVectorNumber = 0, refVectorIndex = corrSumInfo.mTheilerWindow;
         refVectorNumber < corrSumInfo.mNumRefVectors;
         refVectorNumber++, refVectorIndex++) {
      NumericMatrix currentNeighbourCount(nDimensions, nRadii);
      count_neighbours(currentNeighbourCount, neighbourSearcher,
                       refVectorIndex, corrSumInfo);
      update_correlation_sum(corrSum, currentNeighbourCount, 
                             exponent);
    }
  }
}


//[[Rcpp::export]]
NumericMatrix generalized_correlation_sum(const NumericVector& timeSeries, 
                                          int timeLag, int theilerWindow,
                                          NumericVector& radii,
                                          int minEmbeddingDim,
                                          int maxEmbeddingDim,
                                          int corrSumOrder, int numberBoxes){
  /* check embedding dimensions correctness */
  if (minEmbeddingDim > maxEmbeddingDim) {
    throw std::invalid_argument("minEmbeddingDim > maxEmbeddingDim");
  }
  /*
  * Check if in the maximum embedding dimension (worst case) there will be at
  * least 2 phase space vectors to use for computing the correlation sum.
  */ 
  int minimutimeSeriesSize = (2 + (maxEmbeddingDim - 1) * timeLag - 2 * theilerWindow);
  if (timeSeries.size() < minimutimeSeriesSize) {
     throw std::invalid_argument("There aren't enough phase space vectors");
  }
  
  /* order the radii (decreasing order) */
  std::sort(radii.begin(), radii.end(), std::greater<double>());
  neighbour_search neighbourSearcher(build_takens(timeSeries, minEmbeddingDim, timeLag),
                                             radii[0], numberBoxes);
  int nTakens = timeSeries.size() - (maxEmbeddingDim - 1) * timeLag;
  /* First reference vector is theilerWindow. The last reference vector (not
   * included) would be  endReferenceVector = endTakens - theilerWindow. These
   * vectors are selected since they have the same number of possible neighbours,
   * which eases the normalization of the correlation sum.
   */
  int nRefVectors = nTakens - 2 * theilerWindow;
  /* store all the information required to compute the correlation sum in a 
   * handy structure
   */
  correlation_sum_information corrSumInfo(timeSeries, timeLag,
                                          theilerWindow, radii,
                                          minEmbeddingDim, 
                                          maxEmbeddingDim, corrSumOrder,
                                          nTakens, nRefVectors);
  NumericMatrix corrSum(maxEmbeddingDim - minEmbeddingDim + 1, radii.size());

  calculate_weighted_neighbour_count(corrSum, 
                                     neighbourSearcher,
                                     corrSumInfo);
  
  /* Normalize results */
  double denominator = 
  static_cast<double>(nRefVectors * std::pow(nRefVectors - 1, corrSumOrder - 1));
  int nrow = corrSum.nrow();
  int ncol = corrSum.ncol();
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      corrSum(i,j) = corrSum(i, j) / denominator ;
    }
  }
  
  return corrSum;
}

