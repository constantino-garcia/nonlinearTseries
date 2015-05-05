################################################################################
#' Recurrence Quantification Analysis (RQA)
#' @description
#' The Recurrence Quantification Analysis (RQA) is an advanced technique for the nonlinear
#' analysis that allows to quantify the number and duration of the recurrences in the 
#' phase space.
#' @param time.series The original time series from which the phase-space reconstruction is performed.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the \emph{time.series}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param takens Instead of specifying the \emph{time.series}, the \emph{embedding.dim} and the \emph{time.lag}, the user
#' may specify directly the Takens' vectors. 
#' @param radius Maximum distance between two phase-space points to be considered a recurrence.
#' @param lmin Minimal length of a diagonal line to be considered in the RQA. Default \emph{lmin} = 2.
#' @param vmin Minimal length of a vertical line to be considered in the RQA. Default \emph{vmin} = 2.
#' @param do.plot Logical. If TRUE, the recurrence plot is shown. However, plotting the recurrence matrix is computationally 
#'  expensive. Use with caution.
#' @param distanceToBorder In order to avoid border effects, the \emph{distanceToBorder} points near the 
#' border of the recurrence matrix are ignored when computing the RQA parameters. Default, \emph{distanceToBorder} = 2.
#' @return A \emph{rqa}  object that consist of a list with the most important RQA parameters:
#' \itemize{
#'  \item \emph{REC}: Recurrence. Percentage of recurrence points in a Recurrence Plot.
#'  \item \emph{DET}: Determinism. Percentage of recurrence points that form diagonal lines.
#'  \item \emph{LAM}: Percentage of recurrent points that form vertical lines.
#'  \item \emph{RATIO}: Ratio between \emph{DET} and \emph{RR}.
#'  \item \emph{Lmax}: Length of the longest diagonal line.
#'  \item \emph{Lmean}: Mean length of the diagonal lines. The main diagonal is not taken into account.
#'  \item \emph{DIV}: Inverse of \emph{Lmax}.
#'  \item \emph{Vmax}: Longest vertical line.
#'  \item \emph{Vmean}: Average length of the vertical lines. This parameter is also referred to as the Trapping time.
#'  \item \emph{ENTR}: Shannon entropy of the diagonal line lengths distribution
#'  \item \emph{TREND}: Trend of the number of recurrent points depending on the distance to the main diagonal
#'  \item \emph{diagonalHistogram}: Histogram of the length of the diagonals.
#'  \item \emph{recurrenceRate}: Number of recurrent points depending on the distance to the main diagonal.
#' }
#' 
#' @references Zbilut, J. P. and C. L. Webber. Recurrence quantification analysis. Wiley Encyclopedia of Biomedical Engineering  (2006).
#' @examples
#' \dontrun{
#' rossler.ts =  rossler(time=seq(0, 10, by = 0.01),do.plot=FALSE)$x
#' rqa.params=rqa(time.series = rossler.ts, embedding.dim=2, time.lag=1,
#'                radius=1.2,lmin=2,do.plot=FALSE,distanceToBorder=2)
#'                }
#' @author Constantino A. Garcia
#' @rdname rqa
#' @export rqa
#' 
rqa=function(takens = NULL, time.series=NULL, embedding.dim=2, time.lag = 1,radius,lmin = 2,vmin = 2,do.plot=FALSE,distanceToBorder=2){
  if(is.null(takens)){
    takens = buildTakens( time.series, embedding.dim = embedding.dim, time.lag = time.lag)  
  } 
  ntakens = nrow(takens)
  # distance to the border of the matrix to use in the linear regression that estimates
  #the trend
  maxDistanceMD=ntakens-distanceToBorder
  if (maxDistanceMD <=1) maxDistanceMD=2 # this should not happen
  
  neighs=findAllNeighbours(takens,radius)
  if (do.plot) {recurrencePlotAux(neighs)}
  hist=getHistograms(neighs,ntakens,lmin,vmin)
  # calculate the number of recurrence points from the recurrence rate. The recurrence
  # rate counts the number of points at every distance in a concrete side of the main diagonal.
  # Thus, sum all points for all distances, multiply by 2 (count both sides) and add the main
  # diagonal
  numberRecurrencePoints=sum(hist$recurrenceHist)+ntakens
  # calculate the recurrence rate dividing the number of recurrent points at a given
  # distance by all points that could be at that distance
  recurrence_rate_vector=hist$recurrenceHist[1:(ntakens-1)]/((ntakens-1):1)
  #percentage of recurrent points
  REC=(numberRecurrencePoints)/ntakens^2
  diagP=calculateDiagonalParameters(ntakens,numberRecurrencePoints,lmin,hist$diagonalHist,recurrence_rate_vector,maxDistanceMD)  
  #paramenters dealing with vertical lines
  vertP=calculateVerticalParameters(ntakens,numberRecurrencePoints,vmin,hist$verticalHist)
  #join all computations
  rqa.parameters=c(REC=REC,RATIO=diagP$DET/REC,diagP,vertP,list(diagonalHistogram=hist$diagonalHist,recurrenceRate=recurrence_rate_vector))
  class(rqa.parameters) = "rqa"
  return(rqa.parameters)
}


################################################################################
#' Recurrence Plot 
#' @description
#' Plot the recurrence matrix.
#' @details
#' WARNING: This function is computationally very expensive. Use with caution.
#' @param time.series The original time series from which the phase-space reconstruction is performed.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the \emph{time.series}.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param takens Instead of specifying the \emph{time.series}, the \emph{embedding.dim} and the \emph{time.lag}, the user
#' may specify directly the Takens' vectors. 
#' @param radius Maximum distance between two phase-space points to be considered a recurrence.
#' @references Zbilut, J. P. and C. L. Webber. Recurrence quantification analysis. Wiley Encyclopedia of Biomedical Engineering  (2006).
#' @author Constantino A. Garcia
#' @export recurrencePlot
#' @import Matrix
#' @useDynLib nonlinearTseries
recurrencePlot=function(takens = NULL, time.series, embedding.dim, time.lag,radius){
  if(is.null(takens)){
    takens = buildTakens( time.series, embedding.dim = embedding.dim, time.lag = time.lag)  
  } 
  neighs=findAllNeighbours(takens,radius)
  recurrencePlotAux(neighs)
}

#private 
recurrencePlotAux=function(neighs){
  ntakens=length(neighs)
  neighs.matrix = neighbourListSparseNeighbourMatrix(neighs,ntakens)
  # need a print because it is a trellis object!!
  print(image(neighs.matrix,xlab="Number of Takens' vector", ylab="Number of Takens' vector"))
    
}

neighbourListSparseNeighbourMatrix  = function(neighs,ntakens){
  neighs.matrix = Diagonal(ntakens)
  for (i in 1:ntakens){
    if (length(neighs[[i]])>0){
      for (j in neighs[[i]]){
        neighs.matrix[i,j] = 1
      }
    }
  }
  return (neighs.matrix)
}


neighbourListToCsparseNeighbourMatrix = function(neighs,ntakens){
  max.number.neighs = -1
  # store the number of neighs of each takens' vector
  number.neighs = rep(0,ntakens)
  for (i in 1:ntakens){
    number.neighs[[i]] = length(neighs[[i]]) + 1 # add the proper vector the its neighbourhood
    if (number.neighs[[i]] > max.number.neighs){
      max.number.neighs = number.neighs[[i]]
    }
  }
  # add one because i is neighbour of i (itself)
  neighs.matrix = matrix(-1,ncol=max.number.neighs,nrow=ntakens)
  for (i in 1:ntakens){
    # substract 1 to convert to C notation!!!
    if (number.neighs[[i]] > 1) {
      neighs.matrix[i,(1:number.neighs[[i]])] = c(i,neighs[[i]]) -1
    } else{
      neighs.matrix[i,1] = i-1
    } 
    neighs.matrix[i,(1:number.neighs[[i]])] = sort( neighs.matrix[i,1:(number.neighs[[i]])])
  }
  
  return (list(neighs = neighs.matrix, nneighs = number.neighs ))
}


countRecurrencePoints=function(neighs,ntakens){
  count=0
  for (i in 1:ntakens){
    count = count + length(neighs[[i]])
  }
  return (count)
}


calculateVerticalParameters=function(ntakens,numberRecurrencePoints,vmin=2,verticalLinesHistogram){
  #begin parameter computations
  num=sum((vmin:ntakens)*verticalLinesHistogram[vmin:ntakens])
  LAM=num/numberRecurrencePoints
  Vmean=num/sum(verticalLinesHistogram[vmin:ntakens])
  if (is.nan(Vmean)) Vmean=0
  #pick the penultimate
  histogramWithoutZeros=which(verticalLinesHistogram>0)
  if (length(histogramWithoutZeros)>0) Vmax=tail(histogramWithoutZeros,1) else Vmax=0
  
  #results
  params=list(LAM=LAM,Vmax=Vmax,Vmean=Vmean)
  return(params)
}

calculateDiagonalParameters=function(ntakens,numberRecurrencePoints,lmin=2,lDiagonalHistogram,recurrence_rate_vector,maxDistanceMD){
  #begin parameter computations
  num=sum((lmin:ntakens)*lDiagonalHistogram[lmin:ntakens]);
  DET=num/numberRecurrencePoints
  Lmean=num/sum(lDiagonalHistogram[lmin:ntakens])
  aux.index=lmin:(ntakens-1)
  LmeanWithoutMain=(sum((aux.index)*lDiagonalHistogram[aux.index]))/(sum(lDiagonalHistogram[aux.index]))
  #pick the penultimate
  Lmax=tail(which(lDiagonalHistogram>0),2)[1]
  if (Lmax==ntakens) Lmax=0
  DIV=1/Lmax
  pl=lDiagonalHistogram/sum(lDiagonalHistogram)
  diff_0=which(pl>0)
  ENTR=-sum(pl[diff_0]*log(pl[diff_0]));
  
  # use only recurrent points with a distance to the main diagonal < maxDistance
  recurrence_rate_vector=recurrence_rate_vector[1:maxDistanceMD]
  mrrv=mean(recurrence_rate_vector)
  #auxiliar vector for the linear regresion: It is related to the general regression
  #formula xi-mean(x)
  auxiliarVector=(1:maxDistanceMD-(maxDistanceMD+1)/2);auxiliarVector2=auxiliarVector*auxiliarVector
  num=sum(auxiliarVector*((recurrence_rate_vector-mrrv)/2) ) # divide by two because we are having into account just one side of the main diag
  den=sum(auxiliarVector2)
  TREND=num/den
  #results
  params=list(DET=DET,DIV=DIV,Lmax=Lmax,Lmean=Lmean,LmeanWithoutMain=LmeanWithoutMain,ENTR=ENTR,TREND=TREND)
  return(params)
}

getHistograms=function(neighs,ntakens,lmin,vmin){
  
  # the neighbours are labeled from 0 to ntakens-1
  c.matrix = neighbourListToCsparseNeighbourMatrix(neighs,ntakens)
  verticalHistogram = rep(0,ntakens)
  diagonalHistogram = rep(0,ntakens)
  recurrenceHistogram = rep(0,ntakens)
  # auxiliar variables
  hist = .C("getHistograms", neighs = as.integer(c.matrix$neighs),
            nneighs = as.integer(c.matrix$nneighs), ntakens = as.integer(ntakens), 
            vmin = as.integer(vmin), lmin = as.integer(lmin),
            verticalHistogram = as.integer(verticalHistogram),
            diagonalHistogram = as.integer(diagonalHistogram),
            recurrenceHistogram = as.integer(recurrenceHistogram),
            PACKAGE="nonlinearTseries" )
  
  
  return(list(diagonalHist=hist$diagonalHistogram,recurrenceHist=hist$recurrenceHistogram,
       verticalHist=hist$verticalHistogram))
  
}  





