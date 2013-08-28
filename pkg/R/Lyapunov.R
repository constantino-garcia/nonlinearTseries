################################################################################
#' Maximum lyapunov exponent
#' @description
#' Functions for estimating the maximal Lyapunov exponent of a dynamical system from 1-dimensional time series
#' using Takens' vectors.
#' @details It is a well-known fact that close trajectories diverge exponentially fast in a chaotic system. The 
#' averaged exponent that determines the divergence rate is called the Lyapunov exponent (usually denoted with \eqn{\lambda}{lambda}). 
#' If \eqn{\delta(0)}{delta(0)} is the distance between two Takens' vectors in the embedding.dim-dimensional space, we expect that the distance
#' after a time \eqn{t} between the two trajectories arising from this two vectors fulfills:
#' \deqn{\delta (n) \sim \delta (0)\cdot exp(\lambda \cdot t)}{\delta (n) is.approximately \delta (0) exp(\lambda *t).}
#' The lyapunov exponent is estimated using the slope obtained by performing a linear regression of 
#' \eqn{S(t)=\lambda \cdot t \sim log(\delta (t)/\delta (0))}{S(t)=\lambda *t is.approximately log(\delta (t)/\delta (0))} 
#' on  \eqn{t}. \eqn{S(t)} will be estimated by averaging the divergence of several reference points.
#' 
#' The user should plot \eqn{S(t) Vs t} when looking for the maximal lyapunov exponent and, if for some temporal range
#' \eqn{S(t)} shows a linear behaviour, its slope is an estimate of the maximal Lyapunov exponent per unit of time. The estimate
#'  routine allows the user to get always an estimate of the maximal Lyapunov exponent, but the user must check that there is a linear region in the  
#' \eqn{S(t) Vs t}. If such a region does not exist, the estimation should be discarded.
#' @param takens A matrix containing the Takens' vectors (one per row) that will be used to estimate the maximal Lyapunov
#' exponent (see \link{buildTakens}). If the Takens' vectors are not specified, the user must specify the time series
#' (time.series), the embedding dimension (embedding.dim) and the time lag (time.lag) that
#' shall be used to construct the Takens' vectors.
#' @param time.series The original time series from which the maximal Lyapunov exponent will be estimated
#' @param embedding.dim Integer denoting the dimension in which we shall embed the time series (see \link{buildTakens}).
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors (see \link{buildTakens}).
#' @param radius Maximum distance in which will look for nearby trajectories.
#' @param theiler.window Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#'  \emph{theiler.window} time steps in order to be considered neighbours. By using a Theiler window, we exclude temporally correlated 
#'  vectors from our estimations. 
#' @param min.neighs Minimum number of neighbours that a Takens' vector must have to be considered
#' a reference point.
#' @param min.ref.points Number of reference points that the routine will try to use. The routine stops when it finds 
#' \emph{min.ref.points} reference points, saving computation time.
#' @param min.time.steps Integer denoting the number of time steps marking the start of the linear region.
#' @param max.time.steps Integer denoting the number of time steps marking the end of the linear region.
#' @param number.boxes Number of boxes that will be used in the box assisted algorithm (see \link{neighbourSearch}).
#' @param sampling.period Sampling period of the time series. When dealing with a discrete
#' system, the \emph{sampling.period} should be set to 1.
#' @param do.plot Logical value. If TRUE (default value), a plot of \eqn{S(t)} Vs  \eqn{t} is shown.
#' @return A list with three components named  \eqn{time} and \eqn{s.function}. \eqn{time}
#' is a vector containing the temporal interval where the system evolves. It ranges from 0 to \emph{\eqn{max.time.steps \cdot sampling.period}{max.time.steps * sampling.period}}.
#' \eqn{s.function} is a vector containing the 
#' values of the \eqn{S(t)} for each t in the time vector.
#' @references  
#' Eckmann, Jean-Pierre and Kamphorst, S Oliffson and Ruelle, David and Ciliberto, S and others. Liapunov exponents from time series.
#' Physical Review A, 34-6, 4971--4979, (1986).
#' 
#' Rosenstein, Michael T and Collins, James J and De Luca, Carlo J.A practical method for calculating largest Lyapunov exponents from small data sets.
#' Physica D: Nonlinear Phenomena, 65-1, 117--134, (1993).
#' @author Constantino A. Garcia
#' @examples 
#'  \dontrun{
#'  # Estimating the  maximal Lyapunov exponent of the Henon attractor
#'  h=henon(start = c(0.63954883, 0.04772637), do.plot = FALSE)
#'  estimation = maxLyapunov(time.series = h$x, embedding.dim=2,time.lag=1,
#'  radius=0.001,theiler.window=4, min.neighs=2, min.ref.points=500 ,
#'  min.time.steps=5, max.time.steps=10)
#'  cat("expected: ",0.41," calculated: ",estimate(estimation),"\n")
#'  }
#' @rdname maxLyapunov
#' @export maxLyapunov
#' @exportClass maxLyapunov
#' @useDynLib nonlinearTseries
maxLyapunov=function(time.series,takens=NULL,embedding.dim=2,time.lag=1,radius,theiler.window=1,min.neighs=5,
                     min.ref.points=500,min.time.steps=0,max.time.steps=10,number.boxes=NULL,sampling.period=1,do.plot=TRUE){
  #C parameters
  if (is.null(takens)) takens=buildTakens(time.series,embedding.dim=embedding.dim,time.lag=time.lag)
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(time.series, radius)
  numberTakens=nrow(takens)
  Sdn=rep(0,max.time.steps+1)
  
  sol=.C("maxLyapunov",timeSeries=as.double(time.series),takens=as.double(takens),
                   tau=as.integer(time.lag),numberTakens=as.integer(numberTakens),
                   embeddingD=as.integer(embedding.dim),eps=as.double(radius),
                   Sdn=as.double(Sdn),nmax=as.integer(max.time.steps),
                   nminRP=as.integer(min.ref.points),neighMin=as.integer(min.neighs),
                   numberBoxes=as.integer(number.boxes),tdist=as.integer(theiler.window),
                   PACKAGE="nonlinearTseries")

  # create lyapunov structure
  time=(0:max.time.steps)*sampling.period
  max.lyapunov.structure = list(time=time,s.function=sol$Sdn)
  class(max.lyapunov.structure) = "maxLyapunov"
  #plot
  if (do.plot){
    plot(max.lyapunov.structure)
  }
  #return results
  return (max.lyapunov.structure)
}

#' @return The \emph{getTime} function returns the 
#' time in which the divergence of close trajectories was computed
#' @rdname maxLyapunov
#' @export getTime
getTime = function(x){
  return (x$time)
}

#' @return The \emph{getDivergence} function returns the 
#' rate of divergence of close trajectories needed for the maximum Lyapunov
#' exponent estimation
#' @rdname maxLyapunov
#' @export getDivergence
getDivergence = function(x){
  return (x$s.function)
}

#' @rdname maxLyapunov
#' @method plot maxLyapunov  
#' @param ... Additional parameters.
#' @method plot maxLyapunov
#' @S3method plot maxLyapunov 
#' @method plot
plot.maxLyapunov= function(x, ...){
  plot(x$time,x$s.function,xlab="t",ylab=expression("S(t)"),main="Estimating maximal Lyapunov exponent")
}

#' @return In order to obtain an estimation of the Lyapunov exponent the user can use the
#' \emph{estimate} function. The  \eqn{estimate} function allows the user to obtain
#' the maximal Lyapunov exponent obtained by performing a linear regression of \eqn{S(t)} on  \eqn{t}
#' in the region especified with the \emph{regression.range}.
#' @param x A \emph{maxLyapunov} object.
#' @param regression.range Vector with 2 components denoting the range where the function will perform linear regression.
#' @rdname maxLyapunov
#' @S3method estimate maxLyapunov 
#' @method estimate maxLyapunov
#' @method estimate
estimate.maxLyapunov= function(x,regression.range = NULL,do.plot=FALSE,...){
  if (is.null(regression.range)){
    min.time = x$time[[2]] # the first position is always 0
    max.time = tail(x$time,1)
  }else{
    min.time = regression.range[[1]]
    max.time = regression.range[[2]]
  }
  
  
  indx = which(x$time >= min.time & x$time <= max.time )
  y.values = x$s.function[indx]
  x.values = x$time[indx]
  fit=lm(y.values~x.values)
  lyapunov.estimate=fit$coefficients[[2]]
  if (do.plot){
    plot(x)
    lines(x.values,fit$fitted.values,col=4)
    legend("bottomright",col=c(1,4),lty=c(1,1),lwd=c(2.5,2.5), legend=c("S(t) Vs. t", 
                                                                        paste("lyapunov estimate =", sprintf("%.2f",lyapunov.estimate))))
  }
  return (lyapunov.estimate)
}



