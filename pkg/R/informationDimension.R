################################################################################
#' Information dimension
#' @description
#' Functions for estimating the information dimension of a dynamical 
#' system from 1-dimensional time series using Takens' vectors
#' @details 
#' The information dimension is a particular case of the generalized correlation dimension
#' when setting the order q = 1. It is possible to demonstrate that the information dimension
#' \eqn{D_1}{D1} may be defined as:
#' \eqn{D_1=lim_{r \rightarrow 0} <\log p(r)>/\log(r)}{D1=lim{r->0} <ln p(r)>/ln(r)}.
#' Here, \eqn{p(r)} is the probability of finding a neighbour in a neighbourhood of size \eqn{r} and 
#' <> is the mean value. Thus, the information dimension specifies how the average
#' Shannon information scales with the radius \eqn{r}. The user should compute
#' the information dimension for different embedding dimensions for checking 
#' if \eqn{D_1}{D1} saturates.
#' 
#' In order to estimate \eqn{D_1}{D1}, the algorithm looks for the scaling behaviour of the the average
#' radius that contains a given portion (a "fixed-mass") of the total points in the phase space. By performing
#' a linear regression of \eqn{\log(p)\;Vs.\;\log(<r>)}{ln p Vs ln <r>} (being \eqn{p} the fixed-mass of the total points), an estimate
#' of \eqn{D_1}{D1} is obtained. 
#' 
#' The algorithm also introduces a variation of \eqn{p} for achieving a better performance: 
#' for small values of \eqn{p}, all the points in the time series (\eqn{N}) are considered for obtaining
#' \eqn{p=n/N}. Above a maximum number of neighbours \eqn{kMax}, the algorithm obtains \eqn{p} by decreasing the number
#' of points considerd  from the time series  \eqn{M<N}. Thus \eqn{p = kMax/M}.
#' 
#' Even with these improvements, the calculations for the information dimension are heavier than
#' those needed for the correlation dimension. 
#' @param time.series The original time series from which the information dimension will be estimated.
#' @param min.embedding.dim Integer denoting the minimum dimension in which we shall embed the time.series (see \link{buildTakens}). 
#' @param max.embedding.dim Integer denoting the maximum dimension in which we shall embed the time.series (see \link{buildTakens}).Thus,
#' we shall estimate the information dimension between \emph{min.embedding.dim} and \emph{max.embedding.dim}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors (see \code{\link{buildTakens}}).
#' @param min.fixed.mass Minimum percentage of the total points that the algorithm shall use for the estimation.
#' @param max.fixed.mass Maximum percentage of the total points that the algorithm shall use for the estimation.
#' @param number.fixed.mass.points The number of different \emph{fixed mass} fractions between \emph{min.fixed.mass}
#' and \emph{max.fixed.mass} that the algorithm will use for estimation.
#' @param radius Initial radius for searching neighbour points in the phase space. Ideally, it should be small
#' enough so that the fixed mass contained in this radius is slightly greater than the \emph{min.fixed.mass}. However,
#' whereas the radius is not too large (so that the performance decreases) the choice is not critical.
#' @param increasing.radius.factor Numeric value. If no enough neighbours are found within \emph{radius}, the radius
#' is increased by a factor \emph{increasing.radius.factor} until succesful. Default: sqrt(2) = 1.414214.
#' @param number.boxes Number of boxes that will be used in the box assisted algorithm (see \link{neighbourSearch}).
#' @param number.reference.vectors Number of reference points that the routine will try to use, saving computation time.
#' @param theiler.window Integer denoting the Theiler window:  Two Takens' vectors must be separated by more than
#'  theiler.window time steps in order to be considered neighbours. By using a Theiler window, we exclude temporally correlated 
#'  vectors from our estimations. 
#' @param kMax Maximum number of neighbours used for achieving p with all the points from the time series (see Details). Default: 100.
#' @param do.plot Logical value. If TRUE (default value), a plot of the correlation sum is shown.
#' @return A \emph{infDim} object that consist of a list with two components: \emph{log.radius} and \emph{fixed.mass}. \emph{log.radius} contains
#' the average log10(radius) in which the \emph{fixed.mass} can be found.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#  @examples
#  \dontrun{
#  s = sinaiMap(a=0.3,n.sample=5000,start=c(0.23489,0.8923),do.plot=FALSE)
#  inf.dim = infDim(time.series = s$x, embedding.dim = 2, time.lag = 1,
#                         min.fixed.mass=0.01, max.fixed.mass=0.03,
#                         number.fixed.mass.points=1000, radius =0.1, 
#                         increasing.radius.factor = sqrt(2), number.boxes=100, 
#                         number.reference.vectors=100, theiler.window = 200, 
#                         kMax = 100,do.plot=FALSE)
#  estimate(inf.dim)}
#' @rdname infDim
#' @export infDim
#' @exportClass infDim
#' @useDynLib nonlinearTseries
#' @seealso \code{\link{corrDim}}.
infDim <- 
  function(time.series, min.embedding.dim=2, 
           max.embedding.dim = min.embedding.dim,time.lag=1,
           min.fixed.mass, max.fixed.mass, number.fixed.mass.points = 10,
           radius, increasing.radius.factor = sqrt(2), number.boxes=NULL,
           number.reference.vectors, theiler.window = 1,
           kMax = 100,do.plot=TRUE){
   
    embeddings =  min.embedding.dim:max.embedding.dim
    n.embeddings = length(embeddings)
    fixed.mass.vector = 10^(seq(log10(min.fixed.mass),log10(max.fixed.mass),
                                length.out=number.fixed.mass.points))
    infDim.matrix = matrix(0, ncol = number.fixed.mass.points, nrow = n.embeddings)
    dimnames(infDim.matrix)  = list(embeddings,fixed.mass.vector)
    for (m in embeddings){
      infDim.matrix[as.character(m),] = infDimSingleDimension(
                 time.series, m, time.lag, min.fixed.mass,
                 max.fixed.mass, number.fixed.mass.points, radius, 
                 increasing.radius.factor, number.boxes,
                 number.reference.vectors, theiler.window,
                 kMax)
    }
   
    information.dimension.structure = list( fixed.mass = fixed.mass.vector, 
                                            log.radius =  infDim.matrix,
                                            embedding.dims = embeddings)
    class(information.dimension.structure) = "infDim"
    
    if (do.plot){
      plot(information.dimension.structure)
    }

    return(information.dimension.structure)
}


# Private function... computes the average log radius for a given fixed mass
# vector and a single embedding dimension.
infDimSingleDimension <- 
  function(time.series, embedding.dim, time.lag, min.fixed.mass,
                max.fixed.mass, number.fixed.mass.points, radius, 
                increasing.radius.factor, number.boxes,
                number.reference.vectors, theiler.window,
                kMax){
  #estimates number.boxes if it has not been specified
  takens = buildTakens(time.series,embedding.dim,time.lag)
  numberTakens = nrow(takens)
  embedding.dim = ncol(takens)
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(takens, radius)
  fixedMassVector = 10^(seq(log10(min.fixed.mass),log10(max.fixed.mass),length.out=number.fixed.mass.points))
  boxes=rep(0,number.boxes*number.boxes+1)
  averageLogRadiusVector = rep(0,number.fixed.mass.points)
  
  c.result=.C("informationDimension", takens = as.double(takens), 
              numberTakens = as.integer(numberTakens), embeddingD = as.integer(embedding.dim),
              fixedMassVector = as.double(fixedMassVector), 
              fixedMassVectorLength = as.integer(number.fixed.mass.points),
              eps = as.double(radius),
              increasingEpsFactor = as.double(increasing.radius.factor),
              numberBoxes = as.integer(number.boxes), boxes = as.integer(boxes), 
              numberReferenceVectors = as.integer( number.reference.vectors), 
              theilerWindow = as.integer(theiler.window), kMax = as.integer(kMax),
              averageLogRadiusVector = as.double(averageLogRadiusVector),
              PACKAGE="nonlinearTseries")
 
  
  return (c.result$averageLogRadiusVector)  
}

#' Obtain the fixed mass vector used in the information dimension algorithm.
#' @param x A \emph{infDim} object.
#' @return A numeric vector representing the fixed mass vector used
#' in the information dimension algorithm represented by the \emph{infDim} object.
#' @seealso \code{\link{infDim}}
#' @export fixedMass
fixedMass = function(x){
  UseMethod("fixedMass")
}

#' @return The \emph{fixedMass} function returns the fixed mass vector used
#' in the information dimension algorithm.
#' @rdname infDim
#' @method fixedMass infDim
#' @S3method fixedMass infDim
fixedMass.infDim = function(x){
  return (x$fixed.mass)
}


#' Obtain the the average log(radius) computed
#' on the information dimension algorithm.
#' @param x A \emph{infDim} object.
#' @return A numeric vector representing the average log(radius) computed
#' on the information dimension algorithm represented by the \emph{infDim} object.
#' @seealso \code{\link{infDim}}
#' @export logRadius
logRadius = function(x){
  UseMethod("logRadius")
}


#' @return The \emph{logRadius} function returns the average log(radius) computed
#' on the information dimension algorithm.
#' @rdname infDim
#' @method fixedMass infDim
#' @S3method fixedMass infDim
logRadius.infDim = function(x){
  return (x$log.radius)
}


#' @return The \emph{embeddingDims} function returns the 
#' embeddings in which the information dimension was computed
#' @rdname infDim
#' @method embeddingDims infDim
#' @S3method embeddingDims infDim
embeddingDims.infDim = function(x){
  return (embeddingDims.default(x))
}

#' @return The 'estimate' function estimates the information dimension of the 
#' 'infDim' object by by averaging the slopes of the
#' embedding dimensions specified in the \emph{use.embeddings} parameter. The
#' slopes are determined  by performing a linear regression
#' over the fixed mass' range specified in 'regression.range'. If do.plot is TRUE,
#' a graphic of the regression over the data is shown.
#' @param x A \emph{infDim} object.
#' @param regression.range Vector with 2 components denoting the range where the function will perform linear regression.
#' @param use.embeddings A numeric vector specifying which embedding dimensions should 
#' the \emph{estimate} function use to compute the information dimension.
#' @rdname infDim
#' @S3method estimate infDim
#' @method estimate infDim
#' @method estimate
estimate.infDim = function(x, regression.range=NULL, do.plot=TRUE,
                           use.embeddings = NULL,...){
  if (is.null(regression.range)){
    min.fixed.mass = min(x$fixed.mass) # the first position is always 0
    max.fixed.mass = max(x$fixed.mass)
  }else{
    min.fixed.mass = regression.range[[1]]
    max.fixed.mass = regression.range[[2]]
  }
  if (!is.null(use.embeddings)){
    log.radius = x$log.radius[as.character(use.embeddings),]
  }else{
    use.embeddings = x$embedding.dims
    log.radius=x$log.radius
  }
  
  n.embeddings = length(use.embeddings)
  indx = which(x$fixed.mass>=min.fixed.mass & x$fixed.mass<=max.fixed.mass)
  x.values = log10(x$fixed.mass[indx])
  information.dimension=c()
  for(i in 1:n.embeddings){
    y.values = log.radius[as.character(use.embeddings[[i]]),indx]
    reg = lm(y.values~x.values)
    information.dimension = c(information.dimension,1/reg$coefficients[[2]])
    # plotting
    if (do.plot){
      if (i!=1){
        lines(x$fixed.mass,x$log.radius[as.character(use.embeddings[[i]]),],
              type = "p", col = i)
      }else{
        plot(x$fixed.mass,x$log.radius[as.character(use.embeddings[[i]]),],
             log="x",main="Information Dimension",ylab="<log10(radius)>",
             xlab="fixed mass (p)",ylim=range(log.radius))
      }
      lines(x$fixed.mass[indx],reg$fitted.values,col="blue")
    }  
  }
  legend("bottomright",col=1:n.embeddings,lty=rep(1,n.embeddings),
         lwd=rep(2.5,n.embeddings), 
         legend=use.embeddings,title="Embedding dimension")
  information.dimension = mean(information.dimension)
  
  return(information.dimension)
  
}


#' @return The 'plot' function shows two graphics of the information dimension estimate:
#' a graphic of <log10(radius)> Vs fixed mass and a graphic of the local slopes of the information
#' dimension Vs the fixed mass, both in a semi-log scale.
#' @param ... Additional parameters.
#' @rdname infDim
#' @S3method plot infDim
#' @method plot infDim
#' @method plot
plot.infDim = function(x, ...){
  current.par = par() 
  # Check if it is possible to compute local slopes and
  # set layout depending on it
  if ( length(x$fixed.mass) > 1) {
    # it is possible... 3 regions
    layout(rbind(1,2,3), heights=c(4,4,2))
  }else{
    # not possible... just 2 regions (not local slopes)
    layout(rbind(1,2), heights=c(8,2))
  }
  n.embeddings = length(x$embedding.dims)
  for (i in 1:n.embeddings){
    if (i!=1){
      lines(x$fixed.mass,x$log.radius[as.character(x$embedding.dims[[i]]),],
            col = i, type = "p")
    }else{
      plot(x$fixed.mass,x$log.radius[as.character(x$embedding.dims[[i]]),],
           log="x",main="Information Dimension",
           ylab="<log10(radius)>",xlab="fixed mass (p)",
           ylim=range(x$log.radius))  
    } 
    
  }
  #local slopes
  if (length(x$fixed.mass) > 1){
    lfm = log10(x$fixed.mass)
    derivative = t(apply(x$log.radius,MARGIN=1,differentiate,h=lfm[[2]] - lfm[[1]]))
    fixed.mass.axis = differentiateAxis(x$fixed.mass)
    for (i in 1:n.embeddings){
      if (i!=1){
        lines(fixed.mass.axis,derivative[i,],'b',cex=0.3,
              col = i, type = "p")
      }else{
        plot(fixed.mass.axis,derivative[i,],log="x",'b',cex=0.3,col=1,
             xlab="fixed mass p",ylab="local slope d1(p)", 
             main="Local slopes for the Information dimension estimate",
             ylim=range(derivative))  
      } 
      
    }
   
  }
  ### add legend
  par(mar=c(0, 0, 0, 0))
  plot.new()
  legend("center","groups",ncol=ceiling(n.embeddings/2), 
         col=1:n.embeddings,lty=rep(1,n.embeddings),
         lwd=rep(2.5,n.embeddings),
         legend=x$embedding.dims, title="Embedding dimension")
  par(mar=current.par$mar)
  par(mfrow=c(1,1))
}
