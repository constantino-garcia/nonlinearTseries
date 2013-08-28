\name{dfa}
\alias{dfa}
\alias{estimate.dfa}
\alias{getFluctuationFunction}
\alias{getWindowSizes}
\alias{plot.dfa}
\title{Detrended Fluctuation Analysis}
\usage{
  dfa(time.series, window.size.range = c(10, 300),
    npoints = 20, do.plot = TRUE)

  getWindowSizes(x)

  getFluctuationFunction(x)

  \method{plot}{dfa} (x, ...)

  \method{estimate}{dfa} (x, regression.range = NULL,
    do.plot = FALSE, ...)
}
\arguments{
  \item{time.series}{The original time series to be
  analyzed.}

  \item{npoints}{The number of different window sizes that
  will be used to estimate the Fluctuation function in each
  zone.}

  \item{window.size.range}{Range of values for the windows
  size that will be used to estimate the fluctuation
  function. Default: c(10,300).}

  \item{do.plot}{logical value. If TRUE (default value), a
  plot of the Fluctuation function is shown.}

  \item{...}{Additional parameters.}

  \item{x}{A \emph{dfa} object.}

  \item{regression.range}{Vector with 2 components denoting
  the range where the function will perform linear
  regression.}
}
\value{
  A \emph{dfa} object.

  The \emph{getWindowSizes} function returns the windows
  sizes used to detrend the time series

  The \emph{getFluctuationFunction} function returns the
  Fluctuation function of the DFA
}
\description{
  Functions for performing Detrended Fluctuation Analysis
  (DFA), a widely used technique for detecting long range
  correlations in time series. These functions are able to
  estimate several scaling exponents from the time series
  being analyzed. These scaling exponents characterize
  short or long-term fluctuations, depending of the range
  used for regression (see details).
}
\details{
  The Detrended Fluctuation Analysis (DFA) has become a
  widely used technique for detecting long range
  correlations in time series. The DFA procedure may be
  summarized as follows: \enumerate{ \item Integrate the
  time series to be analyzed. The time series resulting
  from the integration will be referred to as the profile.
  \item Divide the profile into N non-overlapping segments.
  \item Calculate the local trend for each of the segments
  using least-square regression.  Compute the total error
  for each of the segments.  \item Compute the average of
  the total error over all segments and take its root
  square. By repeating the previous steps for several
  segment sizes (let's denote it by t), we obtain the
  so-called Fluctuation function \eqn{F(t)}.  \item If the
  data presents long-range power law correlations:
  \eqn{F(t) \sim t^\alpha}{F(t) proportional t^alpha} and
  we may estimate using regression.  \item Usually, when
  plotting \eqn{\log(F(t))\;Vs\;log(t)}{log(F(t)) Vs
  log(t)} we may distinguish two linear regions.  By
  regression them separately, we obtain two scaling
  exponents, \emph{\eqn{\alpha_1}{alpha1}} (characterizing
  short-term fluctuations) and
  \emph{\eqn{\alpha_2}{alpha2}} (characterizing long-term
  fluctuations). } Steps 1-4 are performed using the
  \emph{dfa} function. In order to obtain a estimate of
  some scaling exponent, the user must use the
  \emph{estimate} function specifying the regression range
  (window sizes used to detrend the series).
}
\examples{
\dontrun{
 noise = rnorm(5000)
 dfa.analysis = dfa(time.series = noise, npoints = 10,
               window.size.range=c(10,1000), do.plot=FALSE)
 cat("Theorical: 0.5---Estimated: ", estimate(dfa.analysis),"\\n")
}
}
\author{
  Constantino A. Garcia
}

