################################################################################
#' Estimate several chaotic invariants using linear regression
#' @description
#' Several chaotic invariants are estimated by using linear regression. This function
#' provides a common interface for the estimate of all these parameters (see \code{\link{corrDim}}, 
#' \code{\link{dfa}} and \code{\link{maxLyapunov}} for examples).
#' @param x An object containing all the information needed for the estimate.
#' @param regression.range Range of values on the x-axis on which the regression is performed.
#' @param do.plot Logical value. If TRUE (default value), a plot of the regression is shown.
#' @param ... Additional parameters.
#' @return An estimate of the proper chaotic invariant.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @export estimate
estimate = function(x, regression.range, do.plot,...){
  UseMethod("estimate")
}
  
