\name{estimateEmbeddingDim}
\alias{estimateEmbeddingDim}
\title{Estimate the embedding dimension}
\usage{
  estimateEmbeddingDim(time.series,
    number.points = length(time.series), time.lag = 1,
    max.embedding.dim = 15, threshold = 0.95,
    do.plot = TRUE, theiler.window = 1)
}
\arguments{
  \item{time.series}{The original time series.}

  \item{number.points}{Number of points from the time
  series that will be used to estimate the embedding
  dimension. By default, all the points in the time series
  are used.}

  \item{time.lag}{Time lag used to build the Takens'
  vectors needed to estimate the embedding dimension (see
  \link{buildTakens}). Default: 1.}

  \item{max.embedding.dim}{Maximum possible embedding
  dimension for the time series. Default: 15.}

  \item{threshold}{Numerical value between 0 and 1. The
  embedding dimension is estimated using the E1(d)
  function. E1(d) stops changing when d is greater than or
  equal to embedding dimension, staying close to 1. This
  value establishes a threshold for considering that E1(d)
  has stopped to change. Default: 0.95}

  \item{do.plot}{Logical value. If TRUE (default value), a
  plot of E1(d) and E2(d) is shown.}

  \item{theiler.window}{Integer denoting the Theiler
  window: Two Takens' vectors must be separated by more
  than theiler.window time steps in order to be considered
  neighbours. By using a Theiler window, we exclude
  temporally correlated vectors from our estimations.
  Default: 1.}
}
\description{
  This function determines the minimum embedding dimension
  from a scalar time series using the algorithm proposed by
  L. Cao (see references).
}
\details{
  The Cao's algorithm uses 2 functions in order to estimate
  the embedding dimension from a time series: the E1(d) and
  the E2(d) functions, where d denotes the dimension.

  E1(d) stops changing when d is greater than or equal to
  the embedding dimension, staying close to 1. On the other
  hand, E2(d) is used to distinguish deterministic signals
  from stochastic signals. For deterministic signals, there
  exist some d such that E2(d)!=1. For stochastic signals,
  E2(d) is approximately 1 for all the values.
}
\note{
  The current implementation of this function is fully
  written in R (as a prototype). Thus it requires heavy
  computations and may be quite slow. Future versions of
  the package will solve this issue.

  In the current version of the package, the automatic
  detection of stochastic signals has not been implemented
  yet.
}
\examples{
\dontrun{
h = henon(do.plot=FALSE)
dimension = estimateEmbeddingDim(h$x, time.lag=1, max.embedding.dim=6,
             theiler.window=10, threshold=0.9, do.plot=TRUE)
             }
}
\author{
  Constantino A. Garcia
}
\references{
  Cao, L. Practical method for determining the minimum
  embedding dimension of a scalar time series. Physica D:
  Nonlinear Phenomena, 110,1, pp. 43-50 (1997).
}

