#' Nonlinearity test
#' @details
#' This function runs a set of nonlinearity tests implemented by this and other 
#' R packages, including:
#' \itemize{
#'    \item Teraesvirta's neural metwork test for nonlinearity (\code{\link[tseries]{terasvirta.test}}).
#'    \item White neural metwork test for nonlinearity (\code{\link[tseries]{white.test}}).
#'    \item Keenan's one-degree test for nonlinearity (\code{\link{keenanTest}}).
#'    \item Perform the McLeod-Li test for conditional heteroscedascity (ARCH). (\code{\link{mcleodLiTest}}).
#'    \item Perform the Tsay's test for quadratic nonlinearity in a time series. (\code{\link{tsayTest}}).
#'    \item Perform the Likelihood ratio test for threshold nonlinearity. (\code{\link{thresholdTest}}).
#' }
#' @param time.series The original time.series from which the surrogate data is generated.
#' @param verbose Logical value. If TRUE, a summary of each of the tests is shown.
#' @return A list containing the results of each of the tests.
#' @export nonlinearityTest
#' @import tseries
#' @seealso \code{\link{keenanTest}}, \code{\link{tsayTest}}, 
#' \code{\link{mcleodLiTest}},\code{\link{thresholdTest}}
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
  keenan = keenanTest(time.series)
  nltests$Keenan = keenan
  if (verbose) {
    cat("\t\t ** Keenan's one-degree test for nonlinearity  **\n")  
    cat("\t\t Null hypothesis: The time series follows some AR process\n")
    cat("\t\t F-stat = ",keenan$test.stat," p-value = ", keenan$p.value,"\n\n")    
  }
  # test for conditional heteroscedascity
  mcleodLi = mcleodLiTest(time.series) 
  nltests$McLeodLi = mcleodLi
  if (verbose) {
    cat("\t\t ** McLeod-Li test  **\n")  
    cat("\t\t Null hypothesis: The time series follows some ARIMA process\n")
    cat("\t\t Maximum p-value = ",  max( unlist(mcleodLi) ),"\n\n")    
  }
  # Tsay test for nonlinearity
  tsay = tsayTest(time.series)
  nltests$Tsay = tsay
  if (verbose) {
    cat("\t\t ** Tsay's Test for nonlinearity **\n")  
    cat("\t\t Null hypothesis: The time series follows some AR process\n")
    cat("\t\t F-stat = ",tsay$test.stat," p-value = ", tsay$p.value,"\n\n")    
  }
  # Likelihood ratio test for threshold nonlinearity
  tarTest = thresholdTest(time.series)
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
#' @description
#' Keenan's test: test for nonlinearity against the null hypothesis that the time
#' series follows some AR process.
#' @param time.series The original time.series.
#' @param ... Additional arguments for the \code{\link[stats]{ar}} function.
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
#' @seealso \code{\link{nonlinearityTest}}, \code{\link{tsayTest}}, 
#' \code{\link{mcleodLiTest}},\code{\link{thresholdTest}}
keenanTest = function(time.series, ...) {
  # In the following commentaries, y = time.series to follow the paper's original
  # notation
  y = as.numeric(time.series)
  # Step 1
  # Regress y[t] on {1, y[t-1], ..., y[t-order]}.
  # We implement it using the ar function, since this also permits to estimate
  # a proper order for the AR process.
  y_model = ar(y, ...)
  y_predictions = y - y_model$resid
  # The first predictions and residuals are NA, since there is no enough 
  # data to compute them
  y_predictions = y_predictions[!is.na(y_predictions)]
  y_residuals = y_model$resid[!is.na(y_model$resid)]
  
  # Step 2
  # We have to regress (\hat{y}[t])^2 on {1, y[t-1], ..., y[t-order]}. We do
  # it computing the regression matrix X:
  if (y_model$order > 0) {
    X = buildTakens(head(y, -1), y_model$order, 1)
    X = cbind(1, X) 
  } else {
    X = matrix(rep(1, length(y_predictions)), ncol = 1)
  }
  y2_model = lm.fit(x = X, y = y_predictions ^ 2)
  y2_residuals = y2_model$residuals
  
  # Step 3
  residuals_model = lm(y_residuals ~ y2_residuals + 0)
  
  df1 = 1
  n = length(y)
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


#' Tsay's test
#' @description
#' Tsay's test: test for nonlinearity against the null hypothesis that the time
#' series follows some AR process. This is a generalization of Keenan's test.
#' @param time.series The original time.series.
#' @param order Order used for the AR model.
#' @return A list containing the results of the test, including:
#' \itemize{
#'    \item test.stat: the F-squared test statistic
#'    \item p.value.
#'    \item order: order of the AR process used for testing.
#' }
#' @references Tsay, R. S. (1986), Nonlinearity test for time series, Biometrika,
#'  73, 461-466.
#' @export tsayTest
#' @seealso \code{\link{nonlinearityTest}}, \code{\link{keenanTest}},
#' \code{\link{mcleodLiTest}},\code{\link{thresholdTest}}
tsayTest = function(time.series, order) {
  y = as.numeric(time.series)
  N = length(y)
  # Step 1
  # Regress y[t] on {1, y[t-1], ..., y[t-order]}.
  # If the order is missing, we estimate a proper one using the ar function
  if (missing(order)) {
    order = max(1, ar(y)$order)
  }
  X = buildTakens(head(y, -1), order, 1)
  linear_model = lm(tail(y, -order) ~ X)
  
  # Step 2
  # Regress y[t] on {1, y[t-1], ..., y[t-order]} AND some extra nonlinear terms:
  # y[t -1]^2 , y[t-1] * y[t -2], y[t-1] * y[t-3], ... y[t-order]^2
  lagpad = function(x, k, order, len) {
    c(rep(NA, k), x)[(order + 1):len]
  }
  nonlinear_X = matrix(NA, nrow = nrow(X), 
                       ncol = order * (order + 1) / 2)
  index = 1
  for (lag_1 in 1:order) {
    for (lag_2 in lag_1:order) {
      nonlinear_X[, index] = lagpad(y, lag_1, order, N) * lagpad(y, lag_2, order, N)
      index = index + 1
    }
  }
  nonlinear_model = lm(tail(y, -order) ~ cbind(X, nonlinear_X))
  
  # Step 3
  # Test if the nonlinear terms are significant or not
  compare_models = anova(linear_model, nonlinear_model)
  list('test.stat' = compare_models$F[2],
       'p.value' = compare_models$`Pr(>F)`[2],
       order = order)
}


#' McLeod-Li test
#' @description
#' McLeod-Li test for conditional heteroscedascity (ARCH).
#' @param time.series The original time.series.
#' @param lag.max Maximum number of lags for which to test for 
#' conditional heteroscedascity.
#' @return A list containing the \emph{p.values} for each of the 
#' Ljung-Box tests computed using lags ranging from 1 to \emph{lag.max}.
#' @references Tsay, Ruey S., and Rong Chen. Nonlinear time series analysis. 
#' Vol. 891. John Wiley & Sons, 2018.
#' @export mcleodLiTest
#' @seealso \code{\link{nonlinearityTest}}, \code{\link{keenanTest}},
#'  \code{\link{tsayTest}}, \code{\link{thresholdTest}}
mcleodLiTest = function(time.series, lag.max) {
  y = as.numeric(time.series)
  if (missing(lag.max)) {
    lag.max = 10 * log10(length(y))
  }
  p_values = sapply(
    1:lag.max, 
    function(lag) Box.test(y ^ 2, lag = lag, type='Ljung')$p.value
  )
  list(p.values = p_values)
}


#' Threshold nonlinearity test
#' @description Computes the likelihood ratio test for threshold nonlinearity 
#' with H0 being an AR process and H1 a TAR model.
#' @param time.series The original time.series.
#' @param p The order of the AR process.
#' @param d Delay used for the threshold value in the TAR process.
#' @param lower.percent The threshold value is searched over an interval 
#' defined by \emph{lower.percent} and \emph{upper.percent} of the time series
#' values. This defines the lower percent.
#' @param upper.percent The threshold value is searched over an interval 
#' defined by \emph{lower.percent} and \emph{upper.percent} of the time series
#' values. This defines the upper percent.
#' @return A list containing
#' \itemize{
#'    \item p.value: p-value of the test
#'    \item test.statistic: Likelihood ratio test statistic.
#'    \item percentiles: Since the search for the threshold parameter may 
#'    occur in a narrower interval than the one specified by the user, the
#'    effective lower and upper percents are returned here.
#' }
#' @note Adapted from the \emph{tlrt} function of the \emph{TSA} package.
#' @seealso \code{\link{nonlinearityTest}}, \code{\link{keenanTest}},
#'  \code{\link{tsayTest}}, \code{\link{mcleodLiTest}}
#' @export thresholdTest
#' @author Kung-Sik Chan
#' @references Chan, K.S. (1990). Percentage points of likelihood ratio tests 
#' for threshold autoregression. Journal of Royal Statistical Society,
#' B 53, 3, 691-696.
thresholdTest = function(time.series, p, d = 1, 
                         lower.percent = 0.25, upper.percent = 0.75) {
  y = as.numeric(time.series)
  if (missing(p)) { 
    p = ar(y)$order 
  }
  XY = tarRegressionMatrix(y, p, d)
  
  effective_len = dim(XY)[1]
  nb_coefficients = p + 1  # + 1 to account for the intercept
  
  rss1 = rss2 = rep(Inf, effective_len)
  repeated_min = sum(y[!is.na(NA)] == min(y, na.rm = TRUE))
  repeated_max = sum(y[!is.na(NA)] == max(y, na.rm = TRUE))
  start = max(c(nb_coefficients, repeated_min + 1, d))
  start = max(start, floor(lower.percent * effective_len))
  
  end = effective_len - max(c(nb_coefficients, repeated_max + 1, d)) - 1
  end = min(end, ceiling(upper.percent * effective_len))
  
  effective.lower.percentile = start / effective_len
  effective.upper.percentile = end / effective_len
  
  # We iterate adding points to the regression problem to compute the best
  # TAR threshold. Take into account that XY is sorted by amplitude of the
  # d-delayed time series
  R = qr.R(qr(XY[1:start, ]))
  # y is at the last column of the XY matrix
  y_index = dim(R)[2]
  rss1[start] = R[y_index, y_index] ^ 2
  for (i in (start + 1):end) {
    R = qr.R(qr(rbind(R, XY[i, ])))
    rss1[i] = R[y_index, y_index] ^ 2
  }
  
  R = qr.R(qr(XY[(end + 1):effective_len, ]))
  rss2[end] = R[y_index, y_index] ^ 2
  for (i in rev(start:(end - 1))) {
    R = qr.R(qr(rbind(R, XY[i + 1, ])))
    rss2[i] = R[y_index, y_index] ^ 2
  }
  rss = rss1 + rss2
  rss.H1 = min(rss)
  
  R = qr.R(qr(XY))
  rss.H0 = R[y_index, y_index] ^ 2
  
  test.stat = effective_len * (rss.H0 - rss.H1) / rss.H1
  
  list(
    'percentiles' = c(effective.lower.percentile, effective.upper.percentile) * 100,
    'test.statistic' = test.stat, 
    'p.value' = tlrtPValue(test.stat, 
                           effective.lower.percentile, 
                           effective.upper.percentile, 
                           p = p)
  )
}


tarRegressionMatrix = function(x, p, d) {
  # The matrix will contain all the possible delays of x required for building 
  # an AR(p) process and the 'target' time series
  
  n = length(x)
  start = max(p, d) + 1
  XY = matrix(NA, ncol = p + 2, nrow = n - start + 1)
  
  XY[, 1] = 1  # the intercept
  for (i in 2:(p + 1)) {
    XY[, i] =  x[(start - i + 1):(n - i + 1)]
  }
  XY[, p + 2] = x[start:n]
  # We now sort depending on the d-delayed time series (for the threshold condition).
  # This is used letter to find the best possible threshold using a sequential
  # least squares
  sorting_indices = order(x[(start - d):(n - d)])
  XY[sorting_indices, ]
}


# This function computes approximate p-values for the threshold test. 
# Adapted from the p.value.tlrt of the TSA package
# @author Kung-Sik Chan
tlrtPValue = function(y, a = 0.25, b = 0.75, p = 0) {
  t1 = function(x) {
    Fval = pnorm(x)
    0.5 * log(Fval/(1 - Fval))
  }
  lower = qnorm(a)
  upper = qnorm(b)
  if (p == 0) {
    temp = t1(upper) - t1(lower)
  }
  if (p > 0) {
    tp1 = function(x) {
      Fval = pnorm(x)
      f = dnorm(x)
      b = 2 * Fval - x * f
      c = Fval * (Fval - x * f) - f * f
      root = 0.5 * (b + sqrt(b * b - 4 * c))
      0.5 * log(root/(1 - root))
    }
    tp2 = function(x) {
      Fval = pnorm(x)
      f = dnorm(x)
      b = 2 * Fval - x * f
      c = Fval * (Fval - x * f) - f * f
      root = 0.5 * (b - sqrt(b * b - 4 * c))
      0.5 * log(root/(1 - root))
    }
    temp = (p - 1) * (t1(upper) - t1(lower)) + tp1(upper) - 
      tp1(lower) + tp2(upper) - tp2(lower)
  }
  if (p > 0) {
    return(1 - exp(-2 * dchisq(y, df = p + 1) * (y/(p + 1) - 1) * temp))
  } else {
    z = sqrt(y)
    return(sqrt(2/pi) * exp(-y/2) * (temp * (z - 1/z) +  1/z))
  }
}
