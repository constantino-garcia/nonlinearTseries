#' Time Reversibility statistic
#' @details  The time simmetry statistic  measures the asymmetry of a time 
#' series under time reversal by implementing the third order statistic:
#' \deqn{E[s_n - s_{n-\tau})^3] }{E [s(n) - s(n-tau))^3]}.
#' Since linear stochastic series are symmetric under time reversal, this 
#' statistic may be used for testing the assertion that the data was generated
#' from a stationary linear stochastic process or not.
#' @param time.series The time series used to compute the statistic
#' @param tau Time delay used to compute the third order statistic.
#' @return The time simmetry statistic for the delays specified with
#' \emph{tau}.
#' @seealso \code{\link{timeAsymmetry}}
#' @author Constantino A. Garcia
#' @rdname timeAsymmetry2
#' @export timeAsymmetry2
timeAsymmetry2 = function(time.series, tau) {
  len.ts = length(time.series)
  statistic = numeric(length(tau))
  for (i in seq_along(tau)) {
    statistic[i] = mean(
      (time.series[(tau[[i]] + 1):len.ts] - time.series[1:(len.ts - tau[[i]])]) ^ 3
    )  
  }
  statistic
}

#' Time Reversibility statistic
#' @details  The time simmetry statistic  measures the asymmetry of a time series 
#' under time reversal by calculating:
#' \deqn{E[s_n\cdot s_{n+1}^2] - E[s_n^2\cdot s_{n+1}] }{E[s_n * s_{n+1}^2] - E[s_n^2 * s_{n+1}] }.
#' Since linear stochastic series are symmetric under time reversal, this 
#' statistic may be used for testing the assertion that the data was generated
#' from a stationary linear stochastic process or not.
#' @param time.series The time series used to compute the statistic.
#' @return The time simmetry statistic.
#' @author Constantino A. Garcia
#' @rdname timeAsymmetry
#' @export timeAsymmetry
timeAsymmetry = function(time.series){
  len = length(time.series)
  mean(time.series[1:(len - 1)] *
         time.series[2:len] ^ 2 -  time.series[1:(len - 1)] ^ 2 * time.series[2:len])
}
