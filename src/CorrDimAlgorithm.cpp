#include "CorrDimAlgorithm.h"
#include "genericFunctions.h"
#include <cmath>
#include <stdexcept> 
using namespace Rcpp;

//TODO: move to genericFunctions
NumericMatrix buildTakens(const NumericVector& timeSeries, int embeddingDim, int timeLag){
  int maxJump = (embeddingDim - 1) * timeLag;
  IntegerVector jumpsVector(embeddingDim);
  for (int i=0; i < jumpsVector.size(); i ++){
    jumpsVector[i] = i * timeLag;
  }
  NumericMatrix takensSpace(timeSeries.size() - maxJump, embeddingDim);

  // matrix that will store the takens' vectos. One vector per row
  for (int i=0; i < takensSpace.nrow(); i++) {
    for (int j=0; j < takensSpace.ncol(); j++){
      takensSpace(i, j) = timeSeries[i + jumpsVector[j]];
    }
  }
  return takensSpace;
}

CorrDimAlgorithm::CorrDimAlgorithm(const NumericVector& timeSeries,
                                   int timeLag,
                                   int theilerDistance,
                                   NumericVector& radiusVector,
                                   int minEmbeddingDim,
                                   int maxEmbeddingDim,
                                   int numberBoxes):
  mTimeSeries(timeSeries),
  mTimeLag(timeLag),
  mTheilerDistance(theilerDistance),
  mRadiusVector(radiusVector),
  mMinEmbeddingDim(minEmbeddingDim),
  mMaxEmbeddingDim(maxEmbeddingDim),
  mNumberTakens(timeSeries.size() - (maxEmbeddingDim - 1) * timeLag),
  mFirstReferenceVector(mTheilerDistance) {
  /* check embedding dimensions correctness */
  if (mMinEmbeddingDim > mMaxEmbeddingDim) {
    throw std::invalid_argument("minEmbeddingDim > maxEmbeddingDim");
  }
  /*
   * Check if in the maximum embedding dimension (worst case) there will be at
   * least 2 phase space vectors to use for computing the correlation sum.
   */ 
  int minimumTimeSeriesSize = 
        (2 + (mMaxEmbeddingDim - 1) * mTimeLag - 2 * mTheilerDistance);
  if (mTimeSeries.size() >= minimumTimeSeriesSize) {
    mNeighbourSearcher =
      NeighbourSearchAlgorithm(buildTakens(mTimeSeries, mMinEmbeddingDim, timeLag),
                               max(radiusVector), numberBoxes);
  } else {
    throw std::invalid_argument("There aren't enough phase space vectors");
  }
   /* First reference vector is mTheilerDistance. The last reference vector (not
   * included) would be  endReferenceVector = endTakens - mTheilerDistance. These
   * vectors are selected since they have the same number of possible neighbours,
   * which eases the normalization of the correlation sum.
   */
  int nRefVectors = mNumberTakens - 2 * mTheilerDistance;
  mNumberOfNeighbors = std::vector<NumericMatrix>(nRefVectors);
  for (int i = 0; i < mNumberOfNeighbors.size(); i++){
    mNumberOfNeighbors[i] = NumericMatrix(mMaxEmbeddingDim - mMinEmbeddingDim + 1, mRadiusVector.size());
  }
  // order the radiusVector (decreasing order) and the embeddings vector (ascending)
  std::sort(mRadiusVector.begin(), mRadiusVector.end(), std::greater<double>());
  calculateCorrMatrix();
}

void CorrDimAlgorithm::calculateCorrMatrix(){
  double radiusMax = mRadiusVector[0]; 
  for (int refVectorNumber = 0, refVectorIndex = mFirstReferenceVector; refVectorNumber < mNumberOfNeighbors.size(); 
  refVectorNumber++, refVectorIndex++) {
    updateNeighboursMatrix(mNumberOfNeighbors[refVectorNumber], refVectorIndex);
  } 
} 

void CorrDimAlgorithm::updateNeighboursMatrix(NumericMatrix& currentNeighbourMatrix, int refVectorIndex){
    //find neighbours
    IntegerVector neighboursIndexes = mNeighbourSearcher.getVectorNeighbours(refVectorIndex);
    // check each of the neighbours of refVectorIndex
    for (int j = 0; j < neighboursIndexes.size(); j++){
      int neighbourIndex = neighboursIndexes[j];
      /* Check if we can use this possible neighbour
       * we can not use it if this t(i) vector do not exist in the following embedded spaces
       * or if we have already count this neighbourhood
       * or they are close in time (do not respect the theiler window)
       */
      if ((std::abs(neighbourIndex - refVectorIndex) <= mTheilerDistance) || (neighbourIndex >= mNumberTakens)) {
        continue;
      }
      //ok: we may use this vector
      currentNeighbourMatrix(0,0)++; 
      // complete the row corresponging to the minimum embedding dimension
      int ep;
      for (ep = 1; ep < currentNeighbourMatrix.ncol(); ep++){
        if (mNeighbourSearcher.areNeighbours(refVectorIndex, neighbourIndex, mRadiusVector[ep])){
          currentNeighbourMatrix(0,ep)++;
        } else {
          break;
        }
      }
      // the rest of the embedding dimensions
      int lastRadiusTocheck = ep;
      //iterate over the whole matrix
      for (int m=1; m < currentNeighbourMatrix.nrow(); m++){
        int embeddingDim = mMinEmbeddingDim + m;
        for (ep = 0; ep < lastRadiusTocheck; ep++){ 
          /* We just have to chek the last position of the vector in the new 
           * embedding dimension (as they were neighbours in the mEmbeddingVector[m-1] 
           * dimensional space)
           */ 
          double distance = std::abs(mTimeSeries[refVectorIndex + (embeddingDim - 1) * mTimeLag] - 
                                     mTimeSeries[neighbourIndex + (embeddingDim - 1) * mTimeLag]);
          if (distance < mRadiusVector[ep]) {
            currentNeighbourMatrix(m,ep)++;
          }else {
            /* There is  no need to check other radius in this dimension since they
             * aret sorted in descending order. Also, if i and neighbourIndex aren't 
             * neighbours in a m-dimensional space for a given radius, nor will in
             * a bigger dimension or smaller radius.
             */
            lastRadiusTocheck = ep;
            break;
          }
        } // end looping different radius
      } // end looping different embeddings
    } // end looping neighbours of the i-th reference vector 
}
  
void CorrDimAlgorithm::accumulateWeightedProbability(NumericMatrix& accumulated, const NumericMatrix& x, double exponent){
  for (int i=0; i < accumulated.nrow(); i++) {
    for (int j=0; j < accumulated.ncol(); j++){
      accumulated(i,j) += std::pow(x(i,j), exponent);
    } 
  }
}

NumericMatrix CorrDimAlgorithm::getCorrSum(int order){
  NumericMatrix corrSum(mMaxEmbeddingDim - mMinEmbeddingDim + 1, mRadiusVector.size());
  double exponent = static_cast<double>(order - 1);

  for (int i=0; i < mNumberOfNeighbors.size(); i++){
    accumulateWeightedProbability(corrSum, mNumberOfNeighbors[i], exponent);
  }
 
  /* number of Takens vectors used for searching for neighbours */
  int nReferenceVectors = static_cast<double>(mNumberOfNeighbors.size());
  /* maximum number of neighbours for each of the Takens vectors */
  int nNeighbours = nReferenceVectors - 1;
  double denominator = nReferenceVectors * std::pow(nNeighbours, exponent);
  
  for (int i = 0; i < corrSum.nrow(); i++){
    for (int j = 0; j < corrSum.ncol(); j++){
      corrSum(i,j) = corrSum(i, j) / denominator ;
    }
  }
  return corrSum;
}

