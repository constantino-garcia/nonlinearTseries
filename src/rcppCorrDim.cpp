#include <Rcpp.h>
#include "CorrDimAlgorithm.h"
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix rcppGeneralizedCorrDim(NumericVector timeSeries, 
                                     int timeLag, int theilerDistance,
                                     NumericVector radiusVector,
                                     int minEmbeddingDim, int maxEmbeddingDim,
                                     int corrSumOrder, int numberBoxes){
  CorrDimAlgorithm corrDimCalculator(timeSeries, timeLag,
                                     theilerDistance, radiusVector,
                                     minEmbeddingDim, maxEmbeddingDim,
                                     numberBoxes);
  return corrDimCalculator.getCorrSum(corrSumOrder);
  
}

/*** R
q = 3
# Henon
h = henon(n.sample = 5000, n.transient = 100,
          a = 1.4, b = 0.3, start = c(0.5702737, 0.3953898),
          do.plot = FALSE)
  ts = h$x
  mmin = 2
mmax = 5
time.lag = 1
rmin = 10 ^ -2.6
rmax = 0.01
np = 20
theiler.window = 5
# message("\nComputing Renyi entropy for the Henon attractor (1,4,0.3)\n")
x = corrDim(time.series = ts, min.embedding.dim = mmin,
            max.embedding.dim = mmax, corr.order = q, 
            time.lag = time.lag, min.radius = rmin, max.radius = rmax,
            n.points.radius = np, do.plot = FALSE, theiler.window = theiler.window,
            number.boxes = 100)
  
  rcppx = rcppCorrDim(time.series = ts, min.embedding.dim = mmin,
                      max.embedding.dim = mmax, corr.order = q, 
                      time.lag = time.lag, min.radius = rmin, max.radius = rmax,
                      n.points.radius = np, do.plot = FALSE, theiler.window = theiler.window,
                      number.boxes = 100)
  
  testthat::expect_equal(x, rcppx,tolerance= 1e-3)

library(microbenchmark)
microbenchmark(R=corrDim(time.series = ts, min.embedding.dim = mmin,
                         max.embedding.dim = mmax, corr.order = q,
                         time.lag = time.lag, min.radius = rmin, max.radius = rmax,
                         n.points.radius = np, do.plot = FALSE, theiler.window = theiler.window,
                         number.boxes = 100),
               C = rcppCorrDim(time.series = ts, min.embedding.dim = mmin,
                               max.embedding.dim = mmax, corr.order = q,
                               time.lag = time.lag, min.radius = rmin, max.radius = rmax,
                               n.points.radius = np, do.plot = FALSE, theiler.window = theiler.window,
                               number.boxes = 100),
               times=50)
  
*/
