#' Poincare map
#' @description
#' Computes the Poincare map of the reconstructed trajectories in the 
#' phase-space.
#' @description
#' The Poincare map is a classical dynamical system technique that replaces the 
#' n-th dimensional trajectory in the phase space with an (n-1)-th order 
#' discrete-time called the Poincare map. The points of the Poincare map are 
#' the intersection of the trajectories in the phase-space with a certain 
#' Hyper-plane.
#' @details
#' This function computes the Poincare map taking the Takens' vectors as the 
#' continuous trajectory in the phase space. The  \emph{takens} param has been 
#' included so that the user may specify the real phase-space instead of using 
#' the phase-space reconstruction (see examples).
#' @param time.series The original time series from which the phase-space 
#' reconstruction is done.
#' @param embedding.dim Integer denoting the dimension in which we shall 
#' embed the \emph{time.series}.
#' @param time.lag Integer denoting the number of time steps that will be 
#' use to construct the  Takens' vectors.
#' @param takens Instead of specifying the \emph{time.series}, the 
#' \emph{embedding.dim} and the \emph{time.lag}, the user
#' may specify directly the Takens' vectors. 
#' @param normal.hiperplane.vector The normal vector of the hyperplane that 
#' will be used to compute the Poincare map. If the vector is not specifyed 
#' the program choses the vector (0,0,...,1).
#' @param hiperplane.point  A point on the hyperplane (an hyperplane is 
#' defined with a point and a normal vector).
#' @return Since there are three different Poincare maps, an R list is returned 
#' storing all the information related which all of these maps:
#' \itemize{
#'    \item The positive Poincare map is formed by all the intersections with
#'     the hyperplane in positive direction (defined by the normal vector). The 
#'     \emph{pm.pos} returns the points of the map whereas that 
#'     \emph{pm.pos.time} returns the number of time steps since the beginning
#'      where the intersections occurred. Note that these time steps probably 
#'      won't be integers since the algorithm uses an interpolation procedure 
#'      for calculating the intersection with the hyperplane.
#'    \item Similarly we define a negative Poincare map (\emph{pm.neg} and 
#'    \emph{pm.neg.time}).
#'    \item  Finally, we may define a two-side Poincare map that stores all the
#'     intersections (no matter the direction of the intersection) (\emph{pm}
#'     and \emph{pm.time}).
#' }
#' @examples
#' \dontrun{
#' r=rossler(a = 0.2, b = 0.2, w = 5.7, start=c(-2, -10, 0.2),
#' time=seq(0,300,by = 0.01), do.plot=FALSE)
#' takens=cbind(r$x,r$y,r$z)
#' # calculate poincare sections
#' pm=poincareMap(takens = takens,normal.hiperplane.vector = c(0,1,0), 
#'  hiperplane.point=c(0,0,0) )
#' if (requireNamespace("rgl", quietly = TRUE)) {
#'    rgl::plot3d(takens,size=0.7)
#'    rgl::points3d(pm$pm,col="red")
#' }
#' }
#' @references Parker, T. S., L. O. Chua, and T. S. Parker (1989). Practical
#' numerical algorithms for chaotic systems. Springer New York
#' @author Constantino A. Garcia
#' @rdname poincareMap
#' @export poincareMap
#' @useDynLib nonlinearTseries
poincareMap = function(time.series = NULL, embedding.dim = 2,
                       time.lag = 1, takens = NULL,
                       normal.hiperplane.vector = NULL,
                       hiperplane.point) {
  
  if (is.null(takens)) {
    takens = buildTakens(time.series,embedding.dim,time.lag)
  }
  dimension = ncol(takens)
  n.points = nrow(takens)
  if (is.null(normal.hiperplane.vector)) {
    normal.hiperplane.vector = rep(0, dimension)
    normal.hiperplane.vector[[dimension]] = 1
    hiperplane.point = rep(0, dimension)
  } 
  if (length(normal.hiperplane.vector) != dimension) {
    stop("The hiperplane was defined in a wrong dimensional space\n")
  }
  # this variables will store information about the poincare map
  p.map = .Call("_nonlinearTseries_poincare_map", timeSeries = takens, 
                hiperplanePoint = hiperplane.point,
                normalVector = normal.hiperplane.vector,
                PACKAGE = "nonlinearTseries")
  
  # add attributes
  p.map = propagateTakensAttr(p.map, takens)
  attr(p.map, "normal.hiperplane.vector") = normal.hiperplane.vector
  attr(p.map, "hiperplane.point") = hiperplane.point 
  
  p.map
}
