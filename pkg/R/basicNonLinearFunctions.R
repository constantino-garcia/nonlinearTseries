# private function
# Estimates an appropiate number of boxes for using in the box assisted algorithm
# using the time series being analyzed and the radius of the grid.
# @details
# The estimation takes the difference between the maximum and the minimum of the time
# series and divides it by the radius of the grid. If the number of boxes is too large
# (over a kMaximumNumberBoxes), the number of boxes is fixed to kMaximumNumberBoxes.
# If the number of boxes is too small (below kMinimumNumberBoxes), the number of boxes
# is fixed to kMinimumNumberBoxes.
# @param time.series the original time.series from which the surrogate data is generated
# @param radius width of each of the  rectangles that will be used to build the grid
# in the box assisted algorithm
# @return the number of boxes to used in the box assisted algorithm
# @author Constantino A. Garcia
estimateNumberBoxes = function(time.series, radius){
   kMinimumNumberBoxes = 10
   kMaximumNumberBoxes = 500
   number.boxes = as.integer( ( max(time.series)-min(time.series) )/radius  )
   if (number.boxes > kMaximumNumberBoxes) number.boxes = kMaximumNumberBoxes
   if (number.boxes < kMinimumNumberBoxes) number.boxes = kMinimumNumberBoxes
   
   return (number.boxes)
}

# private function
# checks if takens' vectors v1 and v2 are neighbours
# the neighbourhood is determined using the max norm and an radius radious.
# embedding.dim is not strictly necessary but it is used to accelerate the computations
isNeighbour=function(v1,v2,embedding.dim,radius){
  for (i in 1:embedding.dim){
    if ( abs(v1[[i]]-v2[[i]]) >=radius) return (FALSE);
  }
  return (TRUE);
}


################################################################################
#' Build the Takens' vectors 
#' @description
#' This function builds the Takens' vectors from a given time series. The set
#' of Takens' vector is the result of embedding the time series in a m-dimensional
#' space. That is, the \eqn{n^{th}} Takens' vector is defined as 
#' \deqn{T[n]=\{time.series[n], time.series[n+ timeLag],...time.series[n+m*timeLag]\}.}
#' Taken's theorem states that we can then reconstruct an equivalent dynamical 
#' system to the original one (the 
#' dynamical system that generated the observed time series) by using the Takens' vectors.

#' @param time.series The original time series.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the time.series.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @return A matrix containing the Takens' vectors (one per row).
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @examples 
#'\dontrun{
#'# Build the Takens vector for the Henon map using the x-coordinate time series
#' h = henon(n.sample=  3000,n.transient= 100, a = 1.4, b = 0.3, 
#' start = c(0.73954883, 0.04772637), do.plot = FALSE)
#' takens = buildTakens(h$x,embedding.dim=2,time.lag=1)
#' # using the x-coordinate time series we are able to reconstruct
#' # the state space of the Henon map
#' plot(takens)}
#' @export buildTakens
buildTakens=function (time.series, embedding.dim, time.lag) 
{
  N = length(time.series)
  #vector of increments
  maxjump = (embedding.dim - 1) * time.lag
  jumpsvect = seq(0, maxjump, time.lag)
  #matrix that will store the takens' vectos. One vector per row
  numelem = N - maxjump
  takens = matrix(nrow = numelem, ncol = embedding.dim)
  #build takens' vectors
  for (i in 1:numelem) {
    takens[i, 1:embedding.dim] = time.series[jumpsvect + i]
  }
  return(takens)
}


#' Estimate an appropiate time lag for the Takens' vectors
#' @description
#' Given a time series (time.series), an embedding dimension (m) and a 
#' time lag (timeLag), the \eqn{n^{th}} 
#' Takens' vector is defined as 
#' \deqn{T[n]=\{time.series[n], time.series[n+ timeLag],...time.series[n+m*timeLag]\}.}
#' This function estimates an appropiate time lag by using the autocorrelation function.
#' @details 
#' A basic criteria for estimating a proper time lag is based on the following reasoning:
#' if the time lag used to build the Takens' vectors is too small, the coordinates will
#' be too highly temporally correlated and the embedding will tend to cluster around 
#' the diagonal in the phase space. If the time lag is chosen too large, the resulting coordinates may be almost uncorrelated
#' and the resulting embedding will be very complicated.  Thus, there is a wide variety of methods
#' for estimating an appropiate time lag  based on the study of the autocorrelation function
#' of a given time series:
#'  \itemize{
#'    \item Select the time lag where the autocorrelation function decays to 0 (first.zero method).
#'    \item Select the time lag where the autocorrelation function decays to 1/e (first.e.decay method).
#'    \item Select the time lag where the autocorrelation function reaches its first minimum (first.minimum method).
#'    \item Select the time lag where the autocorrelation function decays to the value specified by the user (first.value method and value parameter).
#' }   
#' @param time.series The original time series.
#' @param method The method that we shall use to estimate the time lag (see the Details section). Available methods
#' are \emph{"first.zero"}, \emph{"first.e.decay"} (default), \emph{"first.minimum"} and \emph{"first.value"}. 
#' @param value Numeric value indicating the value that the autocorrelation function must cross in order to
#' select the time lag. It is used only with the "first.value" method.
#' @param lag.max Maximum lag at which to calculate the acf. By default, the length of the time.series
#' is used.
#' @param do.plot Logical value. If TRUE (default value), a plot of the autocorrelation function is shown.
#' @return The estimated time lag.
#' @note If the autocorrelation function does not cross the specifiged value, an error is thrown. This may be solved
#' by increasing the lag.max or selecting a higher value to which the autocorrelation function must decay.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @export timeLag
timeLag = function (time.series, method="first.e.decay",value = 1/exp(1),lag.max=NULL,do.plot=TRUE){
  ##############################################################################
  ########################### internal functions definitions ###################
  #############################################################################
  # internal function. It will calculate the first time lag where
  # the autocorrelation function decays to a given value
  position.acf.decays=function(time.series,value,lag.max=NULL,do.plot){
      if (is.null(lag.max)){
        l.time.series = length(time.series)
        lag.max = floor(l.time.series/2)
      }
      ac = stats::acf(time.series,lag.max=lag.max,plot=do.plot,type="correlation")
      ac = ac$acf
      cross.position <- which(ac < value)
      # if the autocorrelation function does not cross the value specified
      # by the user, a warning is given
      if (length(cross.position) == 0){
        stop(paste("The autocorrelation function does not cross ", value,
                   ". Choose another \"cross\" value!\n",sep=""))
      }else{   
        cross.position = cross.position[[1]]
        # select closest point to given value
        if ((cross.position > 1) && (abs(ac[[cross.position]] - value) > abs(ac[[cross.position - 1]] - value)))
          cross.position = cross.position - 1
      }
      # If possible: convert from positions to time lags by substracting 1
      if ((cross.position) > 1){cross.position = cross.position-1}
      return (cross.position)
  }
  
  # internal function. It will calculate the first time lag where
  # the autocorrelation function has a minimum
  first.acf.minimum=function(time.series,lag.max=NULL,do.plot){
    if (is.null(lag.max)){
      l.time.series = length(time.series)
      lag.max = l.time.series
    }
    ac = stats::acf(time.series,lag.max=lag.max,plot=do.plot,type="correlation")
    ac = ac$acf
    first.minimum = which(diff(ac) >= 0)
    # if the autocorrelation function is monotonically decreasing, an error is given
    if (length(first.minimum) == 0) {
      stop(paste("The autocorrelation function does not have a minimum\n",sep=""))
    }else{
      first.minimum = first.minimum[[1]]
    } 
    # convert from positions to time lags by substracting 1
    return (first.minimum-1)
  }
  #############################################################################
  #############################################################################
  time.lag=0
  switch(method,
         "first.zero"={time.lag=position.acf.decays(time.series,0,lag.max,do.plot=do.plot)},
         "first.e.decay"={time.lag=position.acf.decays(time.series,1/exp(1),lag.max,do.plot=do.plot)},
         "first.value"={
                        if ((value < 0)||(value>1))
                          stop(paste("Incorrect \"cross\" value", value,
                                      ". \"Cross\" value must be between 0 and 1.
                                     Choose another \"cross\" value!\n",sep=""))
                        time.lag=position.acf.decays(time.series,value,lag.max,do.plot=do.plot)},
         "first.minimum"={time.lag=first.acf.minimum(time.series,lag.max,do.plot=do.plot)},
         {stop("Incorrect method. Available methods are:\"fist.zero\",
               \"first.e.decay\", \"fist.value\" and \"fist.minimum\".\n")})
  return (time.lag)
  
}

