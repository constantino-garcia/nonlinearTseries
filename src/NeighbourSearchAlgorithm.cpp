#include "NeighbourSearchAlgorithm.h"
#include "genericFunctions.h"
using namespace Rcpp;

NeighbourSearchAlgorithm::NeighbourSearchAlgorithm(const NumericMatrix& phaseSpace, 
                                                   double radius, int numberBoxes) :
  mPhaseSpace(phaseSpace), mRadius(radius), mBoxes(numberBoxes * numberBoxes + 1),
  mPossibleNeighbours(phaseSpace.nrow()) {
  // The grid that we shall construct  has numberBoxes*numberBoxes boxes.
  // We add other position (the last position) to store the total length
  
  int numberTakens = mPhaseSpace.nrow();
  int lastPosition = mPhaseSpace.ncol() - 1; //last position of a given mPhaseSpace vector
  int xBoxPosition, yBoxPosition, wrappedBoxPosition;
  //count number of taken vectors in each box. We use the first and last coordinates
  //since they will probably be the least correlated ones
  for (int i = 0; i < numberTakens; i++) {
    xBoxPosition = (int) (mPhaseSpace(i, 0) / radius);
    yBoxPosition = (int) (mPhaseSpace(i, lastPosition) / radius);
    wrappedBoxPosition = getWrappedBoxPosition(xBoxPosition, yBoxPosition);
    mBoxes[wrappedBoxPosition]++;
  }
  
  //accumulate histogram
  //std::partial_sum(mBoxes.begin(), mBoxes.end(), cumBoxes.begin());
  for (int i=1;i<mBoxes.size();i++) mBoxes[i]+=mBoxes[i-1];  
  //fill list of pointers to possible neighbours
  for (int i = 0; i < (numberTakens); i++) {
    xBoxPosition = (int) (mPhaseSpace(i, 0) / radius);
    yBoxPosition = (int) (mPhaseSpace(i, lastPosition) / radius);
    wrappedBoxPosition = getWrappedBoxPosition(xBoxPosition, yBoxPosition);
    mPossibleNeighbours[--mBoxes[wrappedBoxPosition]] = i;
  }
}


//computes the position of a takens' vector that falls into the (row, col)
//unwrapped box in the wrapped 2d box grid.
inline int NeighbourSearchAlgorithm::getWrappedBoxPosition(int row, int col) const{
  int numberBoxes = (int) std::sqrt(mBoxes.size() - 1);
  return (numberBoxes * positiveModulo(row, numberBoxes) + positiveModulo(col, numberBoxes));
}

bool NeighbourSearchAlgorithm::areNeighbours(const int vectorIndex1, const int vectorIndex2) const{
  // use max metric
  for (int i = 0; i < mPhaseSpace.ncol(); i++) {
    if (std::abs(mPhaseSpace(vectorIndex1, i) - mPhaseSpace(vectorIndex2, i)) >= mRadius){
      return false;
    }
  }
  return true;
}


//find neighbours of the "vectorIndex"-th vector from the "takens" takens'
//vector array. The neighbours are found using a box assisted algorithm that
//creates a wrapped grid of "numberBoxes" per dimension. Each box
//has a size of "eps".
//takens: array of takens' vectors
//vectorIndex: number of the vector which you want to calculate all of its neighbors
//numberTakens: number of takens' vectors stored in takens
//embeddingD: dimension of the embedded space
//eps: size of each box of the grid.
//numberBoxes: number of boxes used for the grid in each direction
//boxes: each box of the array will point to all the vectors that fell on it
//possibleNeighbours: the numbers of each takens vector are stored here
// deppending on which box they fell.
//neighbourWorkspace: neighbours of the vectorIndex-th vector.
//nfound: number of neighbours of the vectorIndex-th vector.
// length(neihgList) should be equal to takens.nrow()
IntegerVector NeighbourSearchAlgorithm::getVectorNeighbours(int vectorIndex) const{
  IntegerVector neighbourWorkspace(mPhaseSpace.nrow());
  return boxAssistedNeighbourSearch(vectorIndex, neighbourWorkspace);
}

List NeighbourSearchAlgorithm::getAllNeighbours() const{
  int nVectors = mPhaseSpace.nrow();
  Rcpp::List neighbourList(nVectors);
  Rcpp::IntegerVector neighbourWorkspace(nVectors);
  for (int i = 0; i < nVectors; i++) {
    neighbourList[i] = boxAssistedNeighbourSearch(i, neighbourWorkspace);
  }
  return neighbourList;
}

IntegerVector NeighbourSearchAlgorithm::boxAssistedNeighbourSearch(int vectorIndex, 
                                                                   IntegerVector& neighbourWorkspace) const{
  int nfound = 0;
  int embeddingDim = mPhaseSpace.ncol();
  int lastPosition = embeddingDim - 1;
  int xBoxPosition = (int) (mPhaseSpace(vectorIndex, 0) / mRadius);
  int yBoxPosition = (int) (mPhaseSpace(vectorIndex, lastPosition) / mRadius);
  //look for vector neighbours in the 9 neighbours of the (iBoxPos,jBoxPos) box
  for (int i = xBoxPosition - 1; i <= xBoxPosition + 1; i++) {
    for (int j = yBoxPosition - 1; j <= yBoxPosition + 1; j++) {
      int auxiliarBoxPos = getWrappedBoxPosition(i, j);
      //avoid empty boxes
      for (int boxptr = mBoxes[auxiliarBoxPos + 1] - 1; boxptr >= mBoxes[auxiliarBoxPos]; boxptr--) {
        int possibleNeigh = mPossibleNeighbours[boxptr];
        if (possibleNeigh == vectorIndex) continue;
        if (areNeighbours(vectorIndex, possibleNeigh)) {
          neighbourWorkspace[nfound++] = possibleNeigh;
        }
      }
    }
  }
  return (IntegerVector(neighbourWorkspace.begin(), neighbourWorkspace.begin() + nfound));
}

