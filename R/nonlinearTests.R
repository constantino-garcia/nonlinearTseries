#' Nonlinearity test
#' @details
#' This function runs a set of nonlinearity tests implemented in other R packages including:
#' \itemize{
#'    \item Teraesvirta's neural metwork test for nonlinearity (\code{\link[tseries]{terasvirta.test}}).
#'    \item White neural metwork test for nonlinearity (\code{\link[tseries]{white.test}}).
#'    \item Keenan's one-degree test for nonlinearity (\code{\link[TSA]{Keenan.test}}).
#'    \item Perform the McLeod-Li test for conditional heteroscedascity (ARCH). (\code{\link[TSA]{McLeod.Li.test}}).
#'    \item Perform the Tsay's test for quadratic nonlinearity in a time series. (\code{\link[TSA]{Tsay.test}}).
#'    \item Perform the Likelihood ratio test for threshold nonlinearity. (\code{\link[TSA]{tlrt}}).
#' }
#' @param time.series The original time.series from which the surrogate data is generated.
#' @param verbose Logical value. If TRUE, a summary of each of the tests is shown.
#' @return A list containing the results of each of the tests.
#' @export nonlinearityTest
#' @import tseries
#' @import TSA
#' @seealso \code{\link{keenanTest}}
nonlinearityTest = function(time.series, verbose = TRUE) {
  nltests = list()
  # apply all tests
  # testing linearity in mean
  terasvirta = tseries::terasvirta.test(x = ts(time.series), type = "Chisq")
  nltests$Terasvirta = terasvirta
  if (verbose) {
    cat("\t\t ** Teraesvirta's neural network test  **\n")  
    cat("\t\t Null hypothesis: Linearity in \"mean\" \n")
    cat("\t\t X-squared = ",terasvirta$statistic," df = ",terasvirta$parameter,
        " p-value = ", terasvirta$p.value,"\n\n")  
  }  
  white = tseries::white.test(ts(time.series))
  nltests$White = white
  if (verbose) {
    cat("\t\t ** White neural network test  **\n")  
    cat("\t\t Null hypothesis: Linearity in \"mean\" \n")
    cat("\t\t X-squared = ",white$statistic," df = ",white$parameter,
        " p-value = ", white$p.value,"\n\n")    
  }
  # nonlinearity against the null hypothesis that the time series follows 
  # some AR process.
  keenan = TSA::Keenan.test(time.series)
  nltests$Keenan = keenan
  if (verbose) {
    cat("\t\t ** Keenan's one-degree test for nonlinearity  **\n")  
    cat("\t\t Null hypothesis: The time series follows some AR process\n")
    cat("\t\t F-stat = ",keenan$test.stat," p-value = ", keenan$p.value,"\n\n")    
  }
  # test for conditional heteroscedascity
  mcleodLi = TSA::McLeod.Li.test(y=time.series,plot=FALSE) 
  nltests$McLeodLi = mcleodLi
  if (verbose) {
    cat("\t\t ** McLeod-Li test  **\n")  
    cat("\t\t Null hypothesis: The time series follows some ARIMA process\n")
    cat("\t\t Maximum p-value = ",  max( unlist(mcleodLi) ),"\n\n")    
  }
  # Tsay test for nonlinearity
  tsay = TSA::Tsay.test(time.series)
  nltests$Tsay = tsay
  if (verbose) {
    cat("\t\t ** Tsay's Test for nonlinearity **\n")  
    cat("\t\t Null hypothesis: The time series follows some AR process\n")
    cat("\t\t F-stat = ",tsay$test.stat," p-value = ", tsay$p.value,"\n\n")    
  }
  # Likelihood ratio test for threshold nonlinearity
  tarTest = TSA::tlrt(time.series)
  nltests$TarTest = tarTest
  if (verbose) {
    cat("\t\t ** Likelihood ratio test for threshold nonlinearity **\n")  
    cat("\t\t Null hypothesis: The time series follows some AR process\n")
    cat("\t\t Alternativce hypothesis: The time series follows some TAR process\n")
    cat("\t\t X-squared = ",tarTest$test.stat," p-value = ", tarTest$p.value,"\n\n")    
  }
  nltests
}


#' Keenan's test
#' @details
#' Keenan's test: test for nonlinearity against the null hypothesis that the time
#' series follows some AR process.
#' @param time.series The original time.series.
#' @param ... Additional arguments for the (\code{\link[stats]{ar}}) function.
#' @return A list containing the results of the test, including:
#' \itemize{
#'    \item test.stat: the F-squared test statistic
#'    \item df1 and df2: the degrees of freedom of the test statistic.
#'    \item p.value.
#'    \item order: order of the AR process used for testing.
#' }
#' @references Keenan, D. M. (1985), A Tukey nonadditivity-type test for 
#' time series Nonlinearity, Biometrika, 72, 39-44.
#' @export keenanTest
#' @seealso \code{\link{nonlinearityTest}}
keenanTest = function(time.series, ...) {
  # In the following commentaries, y = time.series to follow the paper's original
  # notation
  
  # Step 1
  # Regress y[t] on {1, y[t-1], ..., y[t-order]}.
  # We implement it using the ar function, since this also permits to estimate
  # a proper order for the AR process.
  y_model = ar(time.series, ...)
  y_predictions = time.series - y_model$resid
  # The first predictions and residuals are NA, since there is no enough 
  # data to compute them
  y_predictions = y_predictions[!is.na(y_predictions)]
  y_residuals = y_model$resid[!is.na(y_model$resid)]
  
  # Step 2
  # We have to regress (\hat{y}[t])^2 on {1, y[t-1], ..., y[t-order]}. We do
  # it computing the regression matrix X:
  X = buildTakens(head(time.series, -1), y_model$order, 1)
  X = cbind(X)
  y2_model = lm.fit(x = X, y = y_predictions ^ 2)
  y2_residuals = y2_model$residuals
  
  # Step 3
  residuals_model = lm(y_residuals ~ y2_residuals + 0)
  
  df1 = 1
  n = length(time.series)
  df2 = n - 2 * y_model$order - 2
  f_stat = summary(residuals_model)$fstatistic[1] * df2 / (n - y_model$order - 1)
  names(f_stat) = NULL
  p_value = pf(f_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
  
  list(
    'test.stat' = f_stat, 
    'df1' = df1,
    'df2' = df2, 
    'p.value' = p_value,
    'order' = y_model$order
  )
}

