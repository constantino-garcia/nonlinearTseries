#include <Rcpp.h>
#include <vector>
#include "generic_functions.h"
#include "neighbour_search.h"
// TODO: ELIMINATE
//#include <gperftools/profiler.h>
using namespace Rcpp;


double which_is_k_smallest(std::vector<double>& v, int k) {
  std::nth_element(v.begin(), v.begin() + k - 1, v.end());
  return v[k - 1]; 
}

double calculate_average_log_radius_with_fixed_mass(NumericMatrix& phaseSpace,  
                                                    double fixedMass,
                                                    double radius, double increasingRadiusFactor,
                                                    int nBoxes, int nReferenceVectors,
                                                    int theilerWindow, int nNeighbours) {
  // Rcout << "1" << std::endl;
  int nPhaseSpaceVectors = phaseSpace.nrow();
  // int embeddingDimension = phaseSpace.ncol();  // Unused
  // Rcout << "le fuck" << std::endl;
  
  neighbour_search neighbourSearcher(phaseSpace, radius, nBoxes);
  
  // Rcout << "2" << std::endl;
  IntegerVector referenceVectors(nReferenceVectors);
  std::vector<double> distVector;
  /*********************** computing k and N from  fixedMass=k/N **************/
  /* when estimating k (from the formula fixedMass=k/N) we have to take into account
   * that there exist a theiler window!! Thus, for every phaseSpace' vector i, we cannot
   * use as neighbours those ranging from (i-theilerWindow):(i+theilerWindow). That is,
   * the effective number of phaseSpace vectors is
   * nPhaseSpaceVectors - 2 * theilerWindow - 1
   */
  int nUsablePhaseSpaceVectors = nPhaseSpaceVectors - 2 * theilerWindow - 1;
  int k = static_cast<int>(fixedMass * nUsablePhaseSpaceVectors) + 1;
  // avoid looking for too many neighbours by imposing a max number nNeighbours. If
  // k > nNeighbours we will reduce the number of phaseSpace vectors
  // considered so that k/N=nNeighbours/N_prime
  int takensVectorsUsed = nPhaseSpaceVectors;
  
  if (k > nNeighbours){
    takensVectorsUsed = 
      static_cast<int>(nUsablePhaseSpaceVectors * nNeighbours / static_cast<double>(k) +
      2 * theilerWindow + 1);
    k = nNeighbours;
  }
  //Rprintf("Original: %d, final: %d ",nPhaseSpaceVectors,takensVectorsUsed);
  // Rcout << "3" << std::endl;
  /******************************** estimator for ln(fixedMass) ****************/
  // use grassberger estimator for ln (fixedMass). This estimator takes into account the
  // finite nature of the estimation. See Generalizations of the Hausdorff dimension
  // of fractal measures (Grassberger 1985). It is currently unused...
  // double lnFixedMass = R::digamma(k) - std::log(takensVectorsUsed - theilerMargin); // Unused
  // double log10FixedMass = lnFixedMass / std::log(10.0); // Unused
  // /****************************** Prepare the iterations **********************/
  // this variable will store the average radious
  double averageLogRadius = 0.0;
  // vector that will store the number of the reference phaseSpace' vector for which
  // we have not yet found enough neighbours. At the beggining they will be all
  // the first numberReferenceVectors vectors...Thus the remaining reference vectors
  // will be numberReferenceVectors
  
  // Rcout << "4" << std::endl;
  for (int i = 0; i < referenceVectors.size(); i++){
    referenceVectors[i] = i;
  }
  
  int remainingReferenceVectors = nReferenceVectors;
  /*************************** iterate ***************************************/
  // Go!: while there is some reference vector with not enough vectors in its neighbourhood
  // we will iterate increasing the neighbourhood radius
  // Rcout << "5" << std::endl;
  for (double currentRadius = radius; remainingReferenceVectors > 0;
  currentRadius *= increasingRadiusFactor){
    // Rcout << "otra vuelta" << std::endl;
    // use the box assisted algorithm with the current radious
    //neighbour_search neighbourSearcher(phaseSpace, currentRadius, nBoxes);
    neighbourSearcher.set_radius(currentRadius);
    int nTakensWithInsufficientNeighs=0;
  
  //Rcout << takensVectorsUsed << std::endl;
    // iterate over the vector that stores the number of each reference phaseSpace' vector
    for (int takensIterator=0; takensIterator < remainingReferenceVectors; takensIterator++){
      // Rcout << "tak_it "<< takensIterator << std::endl;
      int currentReferenceVector=referenceVectors[takensIterator];
      IntegerVector possibleNeighbours = 
        neighbourSearcher.find_neighbours(currentReferenceVector, theilerWindow);
      
      // Rcout << "end search" << std::endl;
      /****************************** check neighs **********************************/
      // for each neighbour of the possibleNeighbours, check if they fulfill the theiler distance
      // and the takensVectorsUsed that we can use
      int nValidNeighbours=0;
      distVector.clear();
      for (int neighIt=0; neighIt < possibleNeighbours.size(); neighIt++){
        int currentNeighbour = possibleNeighbours[neighIt];
        // if the current neighbour does not meet the conditions, move to the next
        if (currentNeighbour > takensVectorsUsed) {
          continue;
        }
        /* The current neighbour does meet the conditions (nValidNeighbours++)!
         * We must compute the distance for using it at the computation of
         *  the radius that encloses k neighbours...
         */
        nValidNeighbours++;
        
        // Rcout << "otro a la saca" << std::endl;
        // if (currentReferenceVector < 2){
        //   
        // Rprintf("%d %d %f\n", currentReferenceVector,currentNeighbour,
        //         neighbourSearcher.calculate_max_distance(currentReferenceVector,
        //                                        currentNeighbour));
        // 
        // 
        // }
        distVector.push_back(neighbourSearcher.calculate_max_distance(currentReferenceVector,
                                                            currentNeighbour)
        );
      }
      /* check if we have found enough nValidNeighbours for averaging */
      if (nValidNeighbours >= k){
        
        // Rcout << ">=k" << std::endl;
        /* which is k-th smallest element of the distVector */
        radius = which_is_k_smallest(distVector, k);
        if (currentReferenceVector < 3){
          //Rprintf("%d: %f\n",currentReferenceVector, radius);    
        }
       
        averageLogRadius += std::log10(radius);
      }else{
        // not enough vectors found... mark for next sweep and increase the number
        // of reference phaseSpace with insufficiente neighbours
        
        // Rcout << "a fregar" << std::endl;
        referenceVectors[nTakensWithInsufficientNeighs++] = currentReferenceVector;
      }
      
    }
    remainingReferenceVectors = nTakensWithInsufficientNeighs;
  }
  // Rcout << "end" <<std::endl;
  return (averageLogRadius / static_cast<double>(nReferenceVectors));
}

// [[Rcpp::export]]
NumericMatrix rcpp_information_dimension(const NumericVector& timeSeries,
                                         const IntegerVector& embeddingDimensions,
                                         int timeLag, const NumericVector& fixedMasses,
                                         double radius, double increasingRadiusFactor,
                                         int nBoxes, int nReferenceVectors,
                                         int theilerWindow, int nNeighbours){
  //ProfilerStart("/tmp/speedup_infDimProf_2.txt");
  NumericMatrix averageLogRadii(embeddingDimensions.size(), fixedMasses.size());
  
  // find the averageLogRadius for each embeddingDimension and fixed mass 
  // Rcout << "vamos!" << std::endl;
  for (int i = 0; i < embeddingDimensions.size(); i++) {
    NumericMatrix phaseSpace = build_takens(timeSeries, embeddingDimensions[i], timeLag);
    for(int j = 0; j < fixedMasses.size(); j++){
      // Rcout << "(,) = " << i << ", " << j << std::endl;
      averageLogRadii(i,j) =
        calculate_average_log_radius_with_fixed_mass(phaseSpace,
                                                    fixedMasses[j],
                                                    radius, increasingRadiusFactor,
                                                    nBoxes, nReferenceVectors,
                                                    theilerWindow, nNeighbours);
      //TODO: store the corrections for lnFixedMass
      //lnFixedMassVector[i] = lnFixedMass;
    }
  }
 // ProfilerStop();
  return averageLogRadii;
}
