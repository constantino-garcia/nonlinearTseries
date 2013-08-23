########################### Determining embedding dimension ####################
#' Estimate the embedding dimension
#' @description
#' This function determines the minimum embedding dimension from a scalar time 
#' series using the algorithm proposed by L. Cao (see references).
#' @details
#' The Cao's algorithm uses 2 functions in order to estimate the embedding dimension
#' from a time series: the E1(d) and the E2(d) functions, where d denotes the dimension.
#' 
#' E1(d) stops changing when d is greater than or equal to the embedding dimension, staying close to 1.
#' On the other hand, E2(d) is used to distinguish deterministic signals from stochastic signals. For 
#' deterministic signals, there exist some d such that E2(d)!=1. For stochastic signals,
#' E2(d) is approximately 1 for all the values. 
#' @note
#' The current implementation of this function is fully written in R (as a prototype).
#' Thus it requires heavy computations and may be quite slow. Future versions of the package
#' will solve this issue.
#' 
#' In the current version of the package, the automatic detection of stochastic 
#' signals has not been implemented yet.
#' @param time.series The original time series.
#' @param number.points Number of points from the time series that will be used to estimate
#' the embedding dimension. By default, all the points in the time series are used.
#' @param time.lag Time lag used to build the Takens' vectors needed to estimate the
#' embedding dimension (see \link{buildTakens}). Default: 1.
#' @param max.embedding.dim Maximum possible embedding dimension for the time series. Default: 15.
#' @param threshold Numerical value between 0 and 1. The embedding dimension is estimated
#' using the E1(d) function. E1(d) stops changing when d is greater than or equal to
#' embedding dimension, staying close to 1. This value establishes a threshold for 
#' considering that E1(d) has stopped to change. Default: 0.95
#' @param do.plot Logical value. If TRUE (default value), a plot of E1(d) and E2(d) is shown.
#' @param theiler.window Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#'  theiler.window time steps in order to be considered neighbours. By using a Theiler window, we exclude temporally correlated 
#'  vectors from our estimations. Default: 1.
#' @references 
#' Cao, L. Practical method for determining the minimum embedding dimension of a scalar time series. Physica D: Nonlinear Phenomena,
#' 110,1, pp. 43-50 (1997).
#' @author Constantino A. Garcia
#' @examples 
#' \dontrun{
#' h = henon(do.plot=FALSE) 
#' dimension = estimateEmbeddingDim(h$x, time.lag=1, max.embedding.dim=6,
#'              theiler.window=10, threshold=0.9, do.plot=TRUE)
#'              }
#' @export estimateEmbeddingDim
estimateEmbeddingDim = function(time.series,  number.points = length(time.series), 
                                          time.lag = 1,  max.embedding.dim = 15,  threshold = 0.95, 
                                          do.plot = TRUE,  theiler.window = 1){
  if (max.embedding.dim < 3) stop("max.embedding.dim should be greater that 2...\n")
  time.series.len = length(time.series)
  data = time.series[(time.series.len/2-number.points/2+1):(time.series.len/2+number.points/2)]
  #if no d verifies E1(d) >= threshold,  then we shall return 0
  embedding.dim = 0 
  # First iteration: get E1(1) and E2(1)
  E.parameters = getCaoParameters(data,  1,  time.lag,  theiler.window)
  E.parameters.next.dim = getCaoParameters(data,  2,  time.lag,  theiler.window)
  E.vector = c(E.parameters$E,  E.parameters.next.dim$E)
  E.star.vector = c(E.parameters$E.star,  E.parameters.next.dim$E.star)
  E1.vector = c(E.vector[[2]]/E.vector[[1]])
  E2.vector = c(E.star.vector[[2]]/E.star.vector[[1]])
  # compute from d = 3 to d = max.embedding.dim
  for (dimension in 3:max.embedding.dim){
    #compute E parameters, E1 and E2
    E.parameters = getCaoParameters(data, dimension, time.lag, theiler.window)
    E.vector[[dimension]] = E.parameters$E
    E.star.vector[[dimension]] = E.parameters$E.star
    E1.vector[[dimension-1]] = E.vector[[dimension]]/E.vector[[dimension-1]]
    E2.vector[[dimension-1]] = E.star.vector[[dimension]]/E.star.vector[[dimension-1]]
    #compute if E1(d)>=threshold...If it is the first time it happens(embedding.dim==0), store
    # the dimension
    if ((embedding.dim==0)&&(E1.vector[[dimension-1]]>=threshold)){
      embedding.dim = dimension-1
    }
  }
  #plot graphics
  if (do.plot){
    plot(1:length(E1.vector), E1.vector, 'b', main = "Computing the embedding dimension", xlab="dimension (d)", ylab="E1(d) & E2(d)", cex = 0.1, ylim = c(0, max(E2.vector)))
    lines(1:length(E2.vector), E2.vector, 'b', col = 4, cex = 0.1)
    abline(h = c(1, threshold), col = 2, lty = c(3, 3))
    legend("bottomright", col = c(1, 4, 2), lty = c(1, 1, 3), lwd = c(2.5, 2.5, 2.5), legend = c("E1(d)", "E2(d)", "limits for E1(d)"))
    
  }
  return (embedding.dim)
  
}


# private function
# auxiliar function to compute E,  E1 and E2 based on the 
# L.Cao article: Practical method for determining the minimum embedding dimension of a scalar time series.
getCaoParameters = function(data, m, time.lag, theiler.window){
  if (theiler.window>length(data)) stop("invalid theiler window\n")
  #construct takens vectors of dimensions m and m+1
  takens = buildTakens(data,  m,  time.lag)
  takens.next.dimension = buildTakens(data,  m+1,  time.lag)
  #get distance in the m dimension to compute nearest neighbour
  if (m==1){
    mutual.distance = as.matrix(dist(data,  method = "maximum"))
  }else{
    mutual.distance = as.matrix(dist(takens,  method = "maximum"))
  } 
  # avoid problems eliminating same vectors from the Takens' matrix
  mutual.distance[mutual.distance==0] = Inf
  #number of iterations needed 
  max.iter = nrow(takens.next.dimension)
  # the a(i, d) parameter from the Cao's article (Equation 1) will be call here
  # min.dist.ratio. On the other hand,  the expression inside the summatory in 
  # equation 4 will be called stochastic.parameter
  min.dist.ratio = c()
  stochastic.parameter = c()
  #computing...
  for (takens.position in 1:max.iter){
    # getClosest neighbour respecting the theiler window. We do not use the 
    # findAllNeighbours algorithm since we need to ensure that every point has a neighbour! 
    # We shall use a matrix of distances. We eliminate vectors that
    # do not respect the Theiler window by putting Infinites
    if ((takens.position-theiler.window-1)<1) left_index = c() else left_index = 1:(takens.position-theiler.window-1)
    if ((takens.position+theiler.window+1)>max.iter) right_index = c() else right_index=(takens.position+theiler.window+1):max.iter
    present.distances = mutual.distance[takens.position, 1:max.iter]
    present.distances[-c(left_index, right_index)]=Inf
    
    # get closest neighbour
    closest.neigh = as.numeric(which.min(present.distances))
    numerator = as.numeric(dist(rbind(takens.next.dimension[takens.position, ], takens.next.dimension[closest.neigh, ]),  method = "maximum"))
    min.dist.ratio[[takens.position]] = numerator/mutual.distance[takens.position, closest.neigh]
    stochastic.parameter[[takens.position]] = abs(data[[takens.position+m*time.lag]]-data[[closest.neigh+m*time.lag]])
  }
  return (list(E = mean(min.dist.ratio), E.star = mean(stochastic.parameter)))
}
