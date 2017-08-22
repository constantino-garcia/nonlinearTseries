# References: http://www.mpipks-dresden.mpg.de/~tisean/TISEAN_2.1/docs/chaospaper/node23.html#SECTION00061000000000000000
#' Nonlinear noise reduction
#' @description
#' Function for denoising a given time series using nonlinear analysis 
#' techniques. 
#' @details
#' This function takes a given time series and denoises it. The denoising
#' is achieved by averaging each Takens' vector in an m-dimensional space
#' with his neighbours (time lag=1). Each neighbourhood is specified with balls 
#' of a given radius
#' (max norm is used).
#' @param time.series The original time series to denoise.
#' @param embedding.dim Integer denoting the dimension in which we shall embed 
#' the \emph{time.series}.
#' @param radius The radius used to looking for neighbours in the phase space 
#' (see details).
#' @return A vector containing the denoised time series.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis 
#' (Cambridge university press)
#' @author Constantino A. Garcia
#' @rdname nonLinearNoiseReduction
#' @export nonLinearNoiseReduction
#' @useDynLib nonlinearTseries
nonLinearNoiseReduction = function(time.series, 
                                   embedding.dim, radius){
  # TODO: provide a better calculation of n.boxes
  n.boxes = 400
  .Call('_nonlinearTseries_nonlinear_noise_reduction',
        PACKAGE = 'nonlinearTseries', 
        as.numeric(time.series), embedding.dim,
        radius, n.boxes)
}
