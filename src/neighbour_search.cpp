#include "neighbour_search.h"
#include "generic_functions.h"
using namespace Rcpp;

neighbour_search::neighbour_search():
  mPhaseSpace(), mEmbeddingDim(0),
  mNumberVectors(0), mRadius(0.0),
  mBoxes(0), mPossibleNeighbours() {
}

neighbour_search::neighbour_search(const NumericMatrix& phaseSpace, 
                                   double radius, int nBoxes) :
  mPhaseSpace(phaseSpace), mEmbeddingDim(phaseSpace.ncol()),
  mNumberVectors(phaseSpace.nrow()), mRadius(radius),
  mSearchingWorkspace(phaseSpace.nrow()),
  /* The grid used to aid the search has nBoxes * nBoxes boxes. We add another 
  * /position (the last position) to store the total length.
  */ 
  mBoxes(nBoxes * nBoxes + 1), mPossibleNeighbours(phaseSpace.nrow()) {
  prepare_box_assisted_search();
}

void neighbour_search::prepare_box_assisted_search(){
  std::fill(mBoxes.begin(), mBoxes.end(), 0);
  std::fill(mPossibleNeighbours.begin(), mPossibleNeighbours.end(), 0);
  int nTakens = mPhaseSpace.nrow();
  int lastPosition = mEmbeddingDim - 1; 
  int xBoxPosition, yBoxPosition, wrappedBoxPosition;
  /* count number of taken vectors in each box. We use the first and last coordinates
  * since they will probably be the least correlated ones.
  */
  for (int i = 0; i < nTakens; i++) {
    xBoxPosition = static_cast<int>(mPhaseSpace(i, 0) / mRadius);
    yBoxPosition = static_cast<int>(mPhaseSpace(i, lastPosition) / mRadius);
    wrappedBoxPosition = get_wrapped_position(xBoxPosition, yBoxPosition);
    mBoxes[wrappedBoxPosition]++;
  }
  
  std::partial_sum(mBoxes.begin(), mBoxes.end(), mBoxes.begin());
  /* fill list of pointers to possible neighbours. */
  for (int i = 0; i < (nTakens); i++) {
    xBoxPosition = static_cast<int>(mPhaseSpace(i, 0) / mRadius);
    yBoxPosition = static_cast<int>(mPhaseSpace(i, lastPosition) / mRadius);
    wrappedBoxPosition = get_wrapped_position(xBoxPosition, yBoxPosition);
    mPossibleNeighbours[--mBoxes[wrappedBoxPosition]] = i;
  }
}



inline bool neighbour_search::comply_theiler_window(int vectorIndex1, int vectorIndex2, 
                                                    int theilerWindow) {
  if (theilerWindow < 0) {
    return true;
  } else {
    return (std::abs(vectorIndex1 - vectorIndex2) > theilerWindow);
  }
}





void neighbour_search::set_radius(double radius){
  mRadius = radius;
  prepare_box_assisted_search();
}



NumericMatrix neighbour_search::get_phase_space() const{
  return mPhaseSpace;  
}



int neighbour_search::get_dimension() const {
  return mEmbeddingDim;
}



int neighbour_search::get_number_vectors() const{
  return mNumberVectors;
}

/*
* Computes the position of a takens' vector that falls into the (row, col)
* unwrapped 2D box in the wrapped box grid.
*/
inline int neighbour_search::get_wrapped_position(int row, int col) const{
  int nBoxes = static_cast<int>(std::sqrt(mBoxes.size() - 1));
  return (nBoxes * positive_modulo(row, nBoxes) + positive_modulo(col, nBoxes));
}

double neighbour_search::calculate_max_distance(int vectorIndex1, int vectorIndex2) const {
  double distance = -1.0;
  for (int i=0; i < mEmbeddingDim; i++){
    distance = std::max(distance, 
                        std::abs(mPhaseSpace(vectorIndex1, i)- mPhaseSpace(vectorIndex2, i))
    );
  }
  return distance;
}  



/* check if two Phase Space vectors are neighbours using the max metric */
bool neighbour_search::are_neighbours(int vectorIndex1, int vectorIndex2,
                                      double neighbourhoodRadius) const{
  for (int i = 0; i < mEmbeddingDim; i++) {
    if (std::abs(mPhaseSpace(vectorIndex1, i) - mPhaseSpace(vectorIndex2, i)) >=
        neighbourhoodRadius){
      return false;
    }
  }
  return true;
}

IntegerVector neighbour_search::find_neighbours(int vectorIndex, int theilerWindow) {
  return box_assisted_search(vectorIndex, theilerWindow);
}

IntegerVector neighbour_search::find_neighbours(int vectorIndex) {
  /* A negative value for the Theiler window indicates that it can be ignored */
  return find_neighbours(vectorIndex, -1);
}

List neighbour_search::find_all_neighbours(int theilerWindow) {
  int nVectors = mPhaseSpace.nrow();
  Rcpp::List neighbourList(nVectors);
  for (int i = 0; i < nVectors; i++) {
    neighbourList[i] = box_assisted_search(i, theilerWindow);
  }
  return neighbourList;
}

List neighbour_search::find_all_neighbours() {
  /* A negative value for the Theiler window indicates that it can be ignored */
  return find_all_neighbours(-1);
}

IntegerVector 
  neighbour_search::box_assisted_search(int vectorIndex, 
                                        int theilerWindow) {
    int nfound = 0;
    int embeddingDim = mEmbeddingDim;
    int lastPosition = embeddingDim - 1;
    int xBoxPosition = static_cast<int>(mPhaseSpace(vectorIndex, 0) / mRadius);
    int yBoxPosition = static_cast<int>(mPhaseSpace(vectorIndex, lastPosition) / mRadius);
    /* Look for vector neighbours in the 9 neighbours of the (iBoxPos,jBoxPos) box. */
    for (int i = xBoxPosition - 1; i <= xBoxPosition + 1; i++) {
      for (int j = yBoxPosition - 1; j <= yBoxPosition + 1; j++) {
        int auxiliarBoxPos = get_wrapped_position(i, j);
        /* avoid empty boxes */
        for (int boxptr = mBoxes[auxiliarBoxPos + 1] - 1; 
             boxptr >= mBoxes[auxiliarBoxPos]; boxptr--) {
          int possibleNeigh = mPossibleNeighbours[boxptr];
          if (possibleNeigh == vectorIndex){
            continue;
          }
          if (are_neighbours(vectorIndex, possibleNeigh, mRadius) && 
              comply_theiler_window(vectorIndex, possibleNeigh, theilerWindow)) {
            mSearchingWorkspace[nfound++] = possibleNeigh;
          }
        }
      }
    }
    return (IntegerVector(mSearchingWorkspace.begin(),
                          mSearchingWorkspace.begin() + nfound));
  }

