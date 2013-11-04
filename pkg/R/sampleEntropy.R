################################################################################
#' Sample Entropy (also known as Kolgomorov-Sinai Entropy)
#' @description
#' The Sample Entropy measures the complexity of a time series. Large values of 
#' the Sample Entropy indicate high complexity whereas that smaller values characterize
#' more regular signals.
#' @details  The sample entropy is computed using:
#' \deqn{h_q(m,r) = log(C_q(m,r)/C_{q}(m+1,r))}{hq(m,r) = log(Cq(m,r)/Cq(m+1,r)),}
#' where \emph{m} is the embedding dimension and \emph{r} is the radius of the neighbourhood. When 
#' computing the correlation dimensions we use the linear regions from the correlation
#' sums in order to do the estimates. Similarly, the sample entropy \eqn{h_q(m,r)}{hq(m,r)} 
#' should not change for both various \emph{m} and \emph{r}.
#' @param corrDim.object A \emph{corrDim} object from which the Sample Entropy
#' of the time series characterized by \emph{corrDim} shall be estimated.
#' @param do.plot do.plot Logical value. If TRUE (default value), a plot of the sample entropy is shown.
#' @return A \emph{sampleEntropy} object that contains a list storing the sample entropy (\emph{sample.entropy}),
#' the embedding dimensions ( \emph{embedding.dims}) and radius (\emph{radius}) for which the sample entropy has 
#' been computed, and the order of the sample entropy (\emph{entr.order}). The sample entropy
#' is stored as a matrix in which each row contains the computations for a given embedding dimension and 
#' each column stores the computations for a given radius.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @examples
#' \dontrun{
#' h=henon(n.sample = 15000, n.transient = 100, a = 1.4, b = 0.3,
#' start = c(0.78,0.8165), do.plot = FALSE)
#' gen.corr.dim=corrDim(time.series=h$x,min.embedding.dim=2,max.embedding.dim=9,
#'                      corr.order=2, time.lag=1,min.radius=0.025,
#'                      max.radius=0.01,n.points.radius=20, do.plot=FALSE,
#'                      theiler.window=10,number.boxes=100)
#' se=sampleEntropy(gen.corr.dim, do.plot=FALSE)
#' estimate(se)}
#' @author Constantino A. Garcia
#' @rdname sampleEntropy
#' @export sampleEntropy
#' @exportClass sampleEntropy
sampleEntropy = function (corrDim.object, do.plot=TRUE){ 
  radius = radius(corrDim.object)
  corr.matrix = corrMatrix(corrDim.object)
  embeddings = embeddingDims(corrDim.object)
  number.embeddings = length(embeddings) - 1
  entropy = matrix(0,nrow= number.embeddings,ncol=length(radius))
  for (i in 1:number.embeddings){
    entropy[i,] = log(corr.matrix[i,]/corr.matrix[i+1,])
  }
  dimnames(entropy)=list(head(embeddings,-1),radius)
  sample.entropy = list(sample.entropy = entropy,embedding.dims = head(embeddings,-1),
                        entr.order=nlOrder(corrDim.object), radius=radius)
  class(sample.entropy)="sampleEntropy"
  
  if (do.plot){
   plot(sample.entropy) 
  }
  
  return (sample.entropy)

}

#' Returns the sample entropy function \eqn{h_q(m,r)} used for the computations
#' of the sample entropy.
#' @param x A \emph{sampleEntropy} object.
#' @return A numeric matrix representing the sample entropy function
#' \eqn{h_q(m,r)} obtained in #' the sample entropy computations represented
#' by the \emph{sampleEntropy} object.
#' @seealso \code{\link{sampleEntropy}}
#' @export sampleEntropyFunction
sampleEntropyFunction = function(x){
  UseMethod("sampleEntropyFunction")
}

#' @return The \emph{sampleEntropyFunction} returns the sample entropy function
#' \eqn{h_q(m,r)} used for the computations. The sample
#' entropy function is represented by a matrix. Each row represents a given
#' embedding dimension whereas that each column representes a different radius.
#' @rdname sampleEntropy
#' @method sampleEntropyFunction sampleEntropy
#' @S3method sampleEntropyFunction sampleEntropy
sampleEntropyFunction.sampleEntropy = function(x){
  return (x$sample.entropy)
}

#' @return The \emph{nlOrder} function returns the order of the sample entropy.
#' @rdname sampleEntropy
#' @method nlOrder sampleEntropy
#' @S3method nlOrder sampleEntropy
nlOrder.sampleEntropy = function(x){
  return(x$entr.order)
}

#' @return The \emph{radius} function returns the radius on which the sample entropy
#'  function has been evaluated.
#' @rdname sampleEntropy
#' @method radius sampleEntropy
#' @S3method radius sampleEntropy
radius.sampleEntropy = function(x){
  return (radius.default(x))
}

#' @return The \emph{embeddingDims} function returns the embedding dimensions 
#' on which the sample entropy function has been evaluated.
#' @rdname sampleEntropy
#' @method embeddingDims sampleEntropy
#' @S3method embeddingDims sampleEntropy
embeddingDims.sampleEntropy = function(x){
  return (embeddingDims.default(x))
}


#' @return The \emph{plot} function shows the graphics for the sample entropy.
#' @rdname sampleEntropy
#' @method plot sampleEntropy
#' @param ... Additional parameters.
#' @S3method plot sampleEntropy
plot.sampleEntropy = function(x,...){
  xlab = expression("ln("*epsilon*")")
  ylab = expression(h[q]*"("*epsilon*")")
  main = expression("Sample entropy (q = "*x$entr.order*")"*h[x$entr.order]*"("*epsilon*")")
  
  number.embeddings = length(x$embedding.dims)
  
  current.par = par()
  layout(rbind(1,2), heights=c(8,2))
  for (i in 1:number.embeddings){
     if (i == 1) {
       plot(x$radius,x$sample.entropy[1,],xlab = xlab,ylab = ylab,main=main,'l',col=1,ylim=range(x$sample.entropy))
     }else{
       lines(x$radius,x$sample.entropy[i,],col=i)
     }
  }
  par(mar=c(0, 0, 0, 0))
  # c(bottom, left, top, right)
  plot.new()
  legend("center","groups",ncol=ceiling(number.embeddings/2),col=1:number.embeddings,lty=rep(1,number.embeddings),
         lwd=rep(2.5,number.embeddings),
         legend=x$embedding.dims,title="Embedding dimension")
  par(mar=current.par$mar)
  par(mfrow=c(1,1))
  
}


#' @details For each embedding dimension the sample
#' entropy is estimated by averaging  \deqn{h_q(m,r) = log(C_q(m,r)/C_{q}(m+1,r))}{hq(m,r) = log(Cq(m,r)/Cq(m+1,r))}
#' over the range specified by \emph{regression range} in the \emph{estimate} function.
#' @param x A \emph{sampleEntropy} object.
#' @param regression.range Vector with 2 components denoting the range where the function will perform linear regression.
#' @param use.embeddings A numeric vector specifying which embedding dimensions should the \emph{estimate} function use to compute
#' the sample entropy.
#' @return The  \emph{estimate} function returns a vector storing the sample entropy estimate for each embedding dimension.
#' @rdname sampleEntropy
#' @method estimate sampleEntropy
#' @S3method estimate sampleEntropy
estimate.sampleEntropy = function(x,regression.range=NULL,do.plot=TRUE,use.embeddings=NULL,...){
  if (is.null(regression.range)){
    regression.range = range(x$radius)
  }  
  if (is.null(use.embeddings)){
    use.embeddings = x$embedding.dims
  }
  #select only embeddings used in the object
  use.embeddings = intersect(as.numeric(use.embeddings),x$embedding.dims)
  len.use.embeddings = length(use.embeddings)
  if(len.use.embeddings<1) {stop("The embeddings specified are not present!\n")}
  #range
  indx = which(x$radius >= regression.range[[1]] & x$radius <=regression.range[[2]])
    
  if (len.use.embeddings == 1){
    sample.entropy.estimate = mean(x$sample.entropy[as.character(use.embeddings),indx])
  }else{
    sample.entropy.estimate = apply(x$sample.entropy[as.character(use.embeddings),indx],MARGIN=1,FUN=mean)  
  }  
  if (do.plot){
    #create new sample entropy object with the embedding.dims and radius range used
    x = list(sample.entropy = x$sample.entropy[as.character(use.embeddings),],
             embedding.dims = use.embeddings,entr.order=x$entr.order, radius=x$radius)
    class(x)="sampleEntropy"
    plotSampleEntropyEstimate(x,sample.entropy.estimate)
  }
  return(sample.entropy.estimate)
}

plotSampleEntropyEstimate = function(sampleEntropy.object,sample.entropy.estimate){
  xlab = expression("ln("*epsilon*")")
  ylab = expression(h[q]*"("*epsilon*")")
  main = expression("Sample entropy (q = "*sampleEntropy.object$entr.order*")"*h[sampleEntropy.object$entr.order]*"("*epsilon*")")
  main=paste("Sample entropy (q = ",sampleEntropy.object$entr.order,")")
  number.embeddings = length(sampleEntropy.object$embedding.dims)
  
  current.par = par()
  layout(rbind(1,2), heights=c(8,2))
  for (i in 1:number.embeddings){
    if (i == 1) {
      plot(sampleEntropy.object$radius,sampleEntropy.object$sample.entropy[1,],xlab = xlab,ylab = ylab,main=main,'l',col=1,ylim=range(sampleEntropy.object$sample.entropy))
      abline(h=sample.entropy.estimate[[i]],lty=3)
    }else{
      lines(sampleEntropy.object$radius,sampleEntropy.object$sample.entropy[i,],col=i)
      abline(h=sample.entropy.estimate[[i]],lty=3)
    }
  }
  par(mar=c(0, 0, 0, 0))
  # c(bottom, left, top, right)
  plot.new()
  legend("center","groups",ncol=ceiling(number.embeddings/2),col=1:number.embeddings,lty=rep(1,number.embeddings),
         lwd=rep(2.5,number.embeddings),
         legend=sampleEntropy.object$embedding.dims,title="Embedding dimension")
  par(mar=current.par$mar)
  
}
