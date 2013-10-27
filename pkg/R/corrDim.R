################################################################################
#' Correlation sum, correlation dimension and generalized correlation dimension 
#' (order q >1).
#' @description
#' Functions for estimating the correlation sum and the correlation dimension of a dynamical 
#' system from 1-dimensional time series using Takens' vectors.
#' @details 
#' The correlation dimension is the most common measure of the fractal dimensionality
#' of a geometrical object embedded in a phase space. In order to estimate the correlation
#' dimension, the correlation sum is defined over the points from the phase space:
#' \deqn{C(r) = \{(number\;of\;points\;(x_i,x_j)\;verifying\;that\;distance\;(x_i,x_j)<r\})/N^2}{C(r) = {number of points(xi,xj)  verifying distance(xi,xj)<r}/N^2}
#' However, this estimator is biased when the pairs in the sum are not statistically independent. For example,
#' Taken's vectors that are close in time, are usually close in the phase space due to the non-zero autocorrelation
#' of the original time series. This is solved by using the so-called Theiler window: two Takens' vectors must be
#' separated by, at least, the time steps specified with this window in order to be considered neighbours. By using a Theiler window,
#' we exclude temporally correlated vectors from our estimations. 
#' 
#' The correlation dimension is estimated using the slope obtained by performing a linear regression of 
#' \eqn{\log10(C(r))\;Vs.\;\log10(r)}{log10(C(r)) Vs. log10(r)}. Since this dimension is supposed to be an invariant of the system, it should not
#' depend on the dimension of the Taken's vectors used to estimate it. Thus, the user should plot \eqn{\log10(C(r))\;Vs.\;\log10(r)}{log10(C(r)) Vs. log10(r)} for several embedding
#' dimensions when looking for the correlation 
#' dimension and, if for some range \eqn{\log10(C(r))}{log10(C(r))} shows a similar linear behaviour in different embedding dimensions (i.e. parallel
#' slopes), these slopes are an estimate of the
#' correlation dimension. The \emph{estimate} routine allows the user to get always an estimate of the correlation dimension,
#' but the user must check that there is a linear region in the correlation sum over different dimensions. 
#' If such a region does not exist, the estimation should be discarded.
#' 
#' Note that the correlation sum  C(r) may be interpreted as:
#' \eqn{C(r) = <p(r)>,}
#' that is: the mean probability of finding a neighbour in a ball of radius r surrounding
#' a point in the phase space. Thus, it is possible to define a generalization of the correlation dimension by writing:
#' \deqn{C_q(r) = <p(r)^{(q-1)}>}{Cq(r) = <p(r)^(q-1)>.}
#' Note that the correlation sum \deqn{C(r) = C_2(r)}{C(r) = C2(r).}
#' 
#' It is possible to determine generalized dimensions Dq using the slope obtained by performing a linear regression of 
#' \eqn{log10(Cq(r))\;Vs.\;(q-1)log10(r)}. The case q=1 leads to the information dimension, that is treated separately
#' in this package (\link{infDim}). The considerations discussed for the correlation dimension estimate
#' are also valid for these generalized dimensions. 
#' @param time.series The original time series from which the correlation sum will be estimated.
#' @param min.embedding.dim Integer denoting the minimum dimension in which we shall embed the time.series (see \link{buildTakens}). 
#' @param max.embedding.dim Integer denoting the maximum dimension in which we shall embed the time.series (see \link{buildTakens}).Thus,
#' we shall estimate the correlation dimension between \emph{min.embedding.dim} and \emph{max.embedding.dim}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors (see \link{buildTakens}).
#' @param min.radius Minimum distance used to compute the correlation sum C(r).
#' @param max.radius Maximum distance used to compute the correlation sum C(r).
#' @param corr.order Order of the generalized correlation Dimension q. It must be greater than 1 (corr.order>1). Default, corr.order=2.
#' @param n.points.radius The number of different radius where we shall estimate.
#' C(r). Thus,  we will estimate C(r) in n.points.radius between min.radius and max.radius.
#' @param theiler.window Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#'  theiler.window time steps in order to be considered neighbours. By using a Theiler window, we exclude temporally correlated 
#'  vectors from our estimations. 
#' @param do.plot Logical value. If TRUE (default value), a plot of the correlation sum is shown.
#' @param number.boxes Number of boxes that will be used in the box assisted algorithm (see \link{neighbourSearch}). If the user does not specify it, the function
#' uses a proper number of boxes.
#' @return  A \emph{corrDim} object that consist of a list with four components named \emph{radius}, \emph{embedding.dims}, \emph{order} and \emph{corr.matrix}.
#' \emph{radius} is a vector containing the different radius where we have evaluated C(r). \emph{embedding.dims} is a vector containing
#' all the embedding dimensions in which we have estimated C(r). \emph{order} stores the order of the generalized correlation dimension
#' that has been used. Finally, \emph{corr.matrix} stores all the correlation
#' sums that have been computed. Each row stores the correlation sum for a concrete embedding dimension whereas each colum
#' stores the correlation sum for a specific radius. 
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @rdname corrDim
#' @export corrDim
#' @exportClass corrDim
#' @useDynLib nonlinearTseries
corrDim = function ( time.series, min.embedding.dim=2, max.embedding.dim = 5, time.lag=1, 
                     min.radius,max.radius,corr.order=2,n.points.radius=5,theiler.window=100,do.plot=TRUE,number.boxes=NULL){
 
  #estimate number of boxes for the box assisted algorithm
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(time.series, min.radius)
  # getTakens vectors for the  minimium embedding dimension. This is done for simplicity
  # in the C code
  takensDimMin=buildTakens(time.series=time.series,embedding.dim=min.embedding.dim,time.lag=time.lag)
  numberTakens=nrow(takensDimMin)
  #other C params
  lenTimeSeries=length(time.series)
  numberEmbeddings=max.embedding.dim-min.embedding.dim+1
  corr.matrix=matrix(0,nrow=numberEmbeddings,ncol=n.points.radius)
  #radius vector. It is prepared to have equally spaced points in the log10(radius) axis
  log.radius=seq(log10(max.radius),log10(min.radius),len=n.points.radius)
  radius=10^log.radius
  
  #call the C program
  if (corr.order==2){
    sol=.C("corrDim", timeSeries=as.double(time.series),lenTimeSeries = as.integer(lenTimeSeries),
           takensDimMin=as.double(takensDimMin),tau=as.integer(time.lag), numberTakens = as.integer(numberTakens),
           minEmbeddingD=as.integer(min.embedding.dim),maxEmbeddingD=as.integer(max.embedding.dim),
           eps=as.double(radius),numberEps=as.integer(n.points.radius),
           numberBoxes=as.integer(number.boxes),tdist=as.integer(theiler.window),corrMatrix=as.double(corr.matrix),
           PACKAGE="nonlinearTseries")
  }else{
    sol=.C("generalizedCorrDim", time.series=as.double(time.series),lenTimeSeries = as.integer(lenTimeSeries),
           takensDimMin=as.double(takensDimMin),tau=as.integer(time.lag), numberTakens = as.integer(numberTakens),
           minEmbeddingD=as.integer(min.embedding.dim),maxEmbeddingD=as.integer(max.embedding.dim), q = as.integer(corr.order),
           eps=as.double(radius),numberEps=as.integer(n.points.radius),
           numberBoxes=as.integer(number.boxes),tdist=as.integer(theiler.window),corrMatrix=as.double(corr.matrix),
           PACKAGE="nonlinearTseries")
  }

  #get the correlation sum matrix
  corr.matrix=matrix(sol$corrMatrix,nrow=numberEmbeddings,ncol=n.points.radius)
  dimnames(corr.matrix) = list(min.embedding.dim:max.embedding.dim,radius)
  #eliminate columns with at least one 0
  wh=which(corr.matrix==0,arr.ind=TRUE)
  wh=unique(wh[,'col'])
  if (length(wh>0)){
    corr.matrix=corr.matrix[,-wh,drop=FALSE]
    #eliminate the corresponding radius values in the radius vector
    radius=radius[-wh]
  }
  # create the corrDim object
  corr.dim = list(corr.matrix = corr.matrix,embedding.dims = min.embedding.dim:max.embedding.dim,radius=radius,corr.order=corr.order)
  class(corr.dim) = "corrDim"
  # plot if necessary
  if (do.plot){
    plot(corr.dim)
  }
  
  return(corr.dim)
  
}

#' @return The \emph{nlOrder} function returns the order of the correlation sum.
#' @rdname corrDim
#' @method nlOrder corrDim
#' @S3method nlOrder corrDim
nlOrder.corrDim = function(x){
  return (x$corr.order)
}


#' Returns the correlation sums stored in the \emph{corrDim} object
#' @param x A \emph{corrDim} object.
#' @return The \emph{corrMatrix} function returns the correlations matrix 
#' storing the correlation sums that have been computed for all the embedding 
#' dimensions.
#' @seealso \code{\link{corrDim}}
#' @export corrMatrix
corrMatrix = function(x){
  UseMethod("corrMatrix")
}

#' @return The \emph{corrMatrix} function returns the correlations matrix 
#' storing the correlation sums that have been computed for all the embedding 
#' dimensions.
#' @rdname corrDim
#' @method corrMatrix corrDim
#' @S3method corrMatrix corrDim
corrMatrix.corrDim = function(x){
  return (x$corr.matrix)
}

#' @return The \emph{radius} function returns the radius on which the correlation sum function has been evaluated.
#' @rdname corrDim
#' @method radius corrDim
#' @S3method radius corrDim
radius.corrDim = function(x){
  return (radius.default(x))
}

#' @return The \emph{embeddingDims} function returns the embedding dimensions on which 
#' the correlation sum function has been evaluated.
#' @rdname corrDim
#' @method embeddingDims corrDim
#' @S3method embeddingDims corrDim
embeddingDims.corrDim = function(x){
  return (embeddingDims.default(x))
}


#' @rdname corrDim
#' @param ... Additional parameters.
#' @method print corrDim
#' @S3method print corrDim
print.corrDim = function(x, ...){
  print(x$corr.matrix)
}


#' @return The \emph{plot} function shows two graphics of the correlation integral:
#' a log-log plot of the correlation sum Vs the radius and the local slopes of 
#' \eqn{log10(C(r))\;Vs\;log10(C(r)).}{log10(C(r)) Vs log10(C(r)).}
#' @rdname corrDim
#' @method plot corrDim
#' @S3method plot corrDim
#' @method plot
plot.corrDim = function(x, ...){
   current.par = par() 
   # Check if it is possible to compute local slopes and
   # set layout depending on it
   if ( ncol(x$corr.matrix) > 1) {
      # it is possible... 3 regions
      layout(rbind(1,2,3), heights=c(4,4,2))
    }else{
      # not possible... just 2 regions (not local slopes)
      layout(rbind(1,2), heights=c(8,2))
    }
    number.embeddings=nrow(x$corr.matrix)
    ### log-log plot
    xlab = ifelse(x$corr.order==2,{"Radius r"},{paste("Radius r^",x$corr.order-1,"",sep="")})
    plot(x$radius^(x$corr.order-1),x$corr.matrix[1,],log="xy",'b',col=1,cex=0.3,ylim=range(x$corr.matrix),
         xlab=xlab,ylab="Correlation Sum C(r)",main="Correlation Sum Vs radius")
    i=2
    while(i<=number.embeddings){
      lines(x$radius^(x$corr.order-1),x$corr.matrix[i,],'b',col=i,cex=0.3)
      i = i + 1
    }
    #### local slopes
    # Check if it is possible to compute local slopes
    if ( ncol(x$corr.matrix) > 1) {
      lcm = log10(x$corr.matrix)
      dlcm= matrix(
                    t(apply(lcm,MARGIN=1,differentiate, 
                          h = (x$corr.order-1)*(log10(x$radius[[2]])-log10(x$radius[[1]]))
                          )),
                    nrow = number.embeddings)
      dlcm=10^dlcm
      radius.axis = differentiateAxis(x$radius)
      xlab = ifelse(x$corr.order==2,{"Radius r"},{paste("Radius r^",x$corr.order-1,"",sep="")})
      plot(radius.axis^(x$corr.order-1),dlcm[1,],'b',log="xy",cex=0.3,col=1,ylim=range(dlcm),
           xlab=xlab,ylab="Local slopes of the Correlation Sum",main="Local slopes of the Correlation Sum Vs Radius")
      i=2
      while(i <= number.embeddings){
        lines(radius.axis^(x$corr.order-1),dlcm[i,],'b',cex=0.3,col=i)
        i=i+1
      }
    }
    ### add legend
    par(mar=c(0, 0, 0, 0))
    plot.new()
    legend("center","groups",ncol=ceiling(number.embeddings/2), 
           col=1:number.embeddings,lty=rep(1,number.embeddings),
           lwd=rep(2.5,number.embeddings),
           legend=x$embedding.dims, title="Embedding dimension")
    par(mar=current.par$mar)
    par(mfrow=c(1,1))
    
}

#' @return The \emph{estimate} function estimates the correlation dimension of the 
#' \emph{corr.dim} object by averaging the slopes of the embedding dimensions specified in
#' the \emph{use.embeddings} parameter. The slopes are determined by performing a linear regression
#' over the radius' range specified in \emph{regression.range}.If \emph{do.plot} is TRUE,
#' a graphic of the regression over the data is shown.
#' @param use.embeddings A numeric vector specifying which embedding dimensions should the \emph{estimate} function use to compute
#' the correlation dimension.
#' @param x A \emph{corrDim} object.
#' @param regression.range Vector with 2 components denoting the range where the function will perform linear regression.
#' @rdname corrDim
#' @S3method estimate corrDim
#' @method estimate corrDim
#' 
estimate.corrDim=function(x, regression.range = NULL, do.plot=FALSE,
                          use.embeddings = NULL,...){
  corr.matrix = corrMatrix(x)
  if (!is.null(use.embeddings)){
    corr.matrix = corr.matrix[as.character(use.embeddings),]
  }
  average=0
  #x axis
  q = nlOrder(x)
  radius=radius(x)
  numberEmbeddings=nrow(corr.matrix)
  log.radius = log10(radius)
  if (is.null(regression.range)){
    r.min = min(radius)
    r.max = max(radius)
  }else{
    # transform the regression range in the corresponding radius
    r.min = (regression.range[[1]])^(1/(q-1))
    r.max = (regression.range[[2]])^(1/(q-1))
  }
  lcm = log10(corr.matrix)
  if (do.plot){
    plot(radius^(q-1),lcm[1,],'b',log="x",col=1,
         cex=0.3,ylim=c(lcm[numberEmbeddings,
         ncol(corr.matrix)],lcm[1,1]),...)
    i=2
    while(i<=numberEmbeddings){
      lines(radius^(q-1),lcm[i,],'b',col=i,cex=0.3)
      i=i+1
    } 
  }
  #average over differents embedding dimensions
  for (i in 1:numberEmbeddings){
    new.corr = eliminateDuplicates(corr.matrix[i,] , radius)
    indx = which(new.corr$radius >= r.min & new.corr$radius <= r.max)
    y.values = log10(new.corr$correlation[indx])
    x.values = (q-1)*log10(new.corr$radius[indx])
    fit=lm(y.values ~ x.values)
    if (do.plot){
      lines(new.corr$radius[indx]^(q-1),fit$fitted.values,
            col="blue",type="b",cex=0.3)
    }
    #print(fit$coefficients[[2]])
    average=average + fit$coefficients[[2]]
  }
  average=average/numberEmbeddings
  #return the correlation dimension estimate
  return(average)
}

#private function
#eliminate duplicate correlation.sums with different radius
eliminateDuplicates = function(correlation.sum , radius){
  len.correlation.sum  = length(correlation.sum)
  unique.correlation.sum = unique(correlation.sum)
  unique.radius = c()
  len.unique.correlation.sum = length(unique.correlation.sum)
  if (len.unique.correlation.sum < len.correlation.sum){
    radius.position = 1
    for (corr in unique.correlation.sum){
      index = which(correlation.sum == corr)
      unique.radius[[radius.position]] = median(radius[index])
      radius.position = radius.position + 1
    }
  }else{
    unique.radius = radius
  }
  return (list(correlation = unique.correlation.sum,radius = unique.radius))
}

