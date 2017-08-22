#' neighbour search
#' @description
#' This function finds all the neighbours of a given Takens' vector. The 
#' neighbours are found using a box assisted algorithm that creates a wrapped 
#' grid with a given number of  boxes per dimension. 
#' @param takens The matrix containing all the Takens' vectors (see 
#' \link{buildTakens}).
#' @param positionTakens Integer denoting the Takens' vector whose neighbours 
#' will be searched.
#' @param radius Distance in which the algorithm will search for neighbours.
#' @param number.boxes Integer denoting the number of boxes per dimension that 
#' will be used to construct a wrapped grid (see Schreiber). If the user does 
#' not specify a number of boxes, this function estimates a proper number.
#' @return A containing all the neighbours of the 
#' \emph{positionTakens-th} Takens' vector. If the list is empty, that means 
#' that there is no neighbour of the \emph{positionTakens-th} Takens' vector 
#' in the given radius.
#' @seealso \code{\link{findAllNeighbours}}.
#' @references  Schreiber, T. Efficient neighbor searching in nonlinear time 
#' series analysis. Int. J. Bifurcation and Chaos, 5, p. 349, (1995).
#' @author Constantino A. Garcia
#' @export neighbourSearch
#' @useDynLib nonlinearTseries
neighbourSearch = function(takens,positionTakens,radius,number.boxes=NULL) {
  #estimates number.boxes if it has not been specified
  if (is.null(number.boxes)) {
    number.boxes = estimateNumberBoxes(takens, radius)
  }
  
  cneighs = .Call('_nonlinearTseries_getVectorNeighbours', 
                  PACKAGE = 'nonlinearTseries', 
                  takens, positionTakens, radius, number.boxes)
  # TODO: remove list conversion
  nfound = length(cneighs)
  if (nfound == 0) {
    finalNeighs = list(nfound = 0, neighList = list())
  } else {
    finalNeighs = list(nfound = nfound, neighList = cneighs)
  }
  # end TODO 
  
  finalNeighs = propagateTakensAttr(finalNeighs, takens)
  attr(finalNeighs,"takens.index") = positionTakens
  attr(finalNeighs,"radius") = radius
  finalNeighs
}

#' neighbour search
#' @description
#' This function finds all the neighbours of all the vectors from Takens'
#' vector array. The neighbours are found using a box assisted algorithm that
#' creates a wrapped grid of a given number of boxes per dimension. 
#' @param takens The matrix containing all the Takens' vectors 
#' (see \link{buildTakens}).
#' @param radius Distance in which the algorithm will search for neighbours.
#' @param number.boxes Integer denoting the number of boxes per dimension that 
#' will be used to construct a wrapped grid (see Schreiber). If the user does 
#' not specify a number of boxes, this function estimates a proper number.
#' @return A list in which the n-th position contains another list with all 
#' the neighbours of the n-th Takens' vector. If the list is empty, that means 
#' that there is no neighbour of the n-th Takens' vector in the given radius.
#' @references Schreiber, T. Efficient neighbor searching in nonlinear time 
#' series analysis. Int. J. Bifurcation and Chaos, 5, p. 349, (1995).
#' @author Constantino A. Garcia
#' @examples 
#' \dontrun{
#' # Find all the neighbours Takens' vectors build from the Henon time
#' # series. The size of the neighbourhood is set to 0.1.
#' h=henon(start = c(0.63954883, 0.04772637), do.plot = FALSE)
#' takens = buildTakens(h$x,embedding.dim=2,time.lag=1)
#' neighbours=findAllNeighbours(takens,0.1)
#' }
#' @seealso \code{\link{neighbourSearch}}.
#' @export findAllNeighbours
findAllNeighbours = function(takens, radius, number.boxes = NULL) {
  #estimates number.boxes if it has not been specified
  if (is.null(number.boxes)) {
    number.boxes = estimateNumberBoxes(takens, radius)
  }
  allneighs = .Call('_nonlinearTseries_getAllNeighbours',
                    PACKAGE = 'nonlinearTseries',
                    takens, radius, number.boxes)
  # TODO: eliminate in future conversion to list
  allneighs = lapply(allneighs, function(x){
    if (length(x) == 0) {
      return(list())
    } else {
      return(x)
    }
  })
  
  allneighs = propagateTakensAttr(allneighs, takens)
  attr(allneighs, "radius") = radius
  allneighs
}

#brute force algorithm that computes all distances  to find
#each neighbourhood. This function will be used to test the
#box assisted algorithm
findAllNeighboursDist = function(data, radius) {
  dst = as.matrix(dist(data, method = "maximum"))
  dst[col(dst) == row(dst)] = Inf
  sol = list()
  l = nrow(data)
  for (i in 1:l) {
    aux = as.vector(which(dst[i, ] < radius))
    if (length(aux) == 0) {
      sol[[i]] = list()
    } else {
      sol[[i]] = aux
    }
  }
  sol
}
