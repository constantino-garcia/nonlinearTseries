\name{corrDim}
\alias{corrDim}
\alias{corrMatrix.corrDim}
\alias{embeddingDims.corrDim}
\alias{estimate.corrDim}
\alias{nlOrder.corrDim}
\alias{plot.corrDim}
\alias{print.corrDim}
\alias{radius.corrDim}
\title{Correlation sum, correlation dimension and generalized correlation dimension
(order q >1).}
\usage{
  corrDim(time.series, min.embedding.dim = 2,
    max.embedding.dim = 5, time.lag = 1, min.radius,
    max.radius, corr.order = 2, n.points.radius = 5,
    theiler.window = 100, do.plot = TRUE,
    number.boxes = NULL)

  \method{nlOrder}{corrDim} (x)

  \method{corrMatrix}{corrDim} (x)

  \method{radius}{corrDim} (x)

  \method{embeddingDims}{corrDim} (x)

  \method{print}{corrDim} (x, ...)

  \method{plot}{corrDim} (x, ...)

  \method{estimate}{corrDim} (x, regression.range = NULL,
    do.plot = FALSE, use.embeddings = NULL, ...)
}
\arguments{
  \item{time.series}{The original time series from which
  the correlation sum will be estimated.}

  \item{min.embedding.dim}{Integer denoting the minimum
  dimension in which we shall embed the time.series (see
  \link{buildTakens}).}

  \item{max.embedding.dim}{Integer denoting the maximum
  dimension in which we shall embed the time.series (see
  \link{buildTakens}).Thus, we shall estimate the
  correlation dimension between \emph{min.embedding.dim}
  and \emph{max.embedding.dim}.}

  \item{time.lag}{Integer denoting the number of time steps
  that will be use to construct the Takens' vectors (see
  \link{buildTakens}).}

  \item{min.radius}{Minimum distance used to compute the
  correlation sum C(r).}

  \item{max.radius}{Maximum distance used to compute the
  correlation sum C(r).}

  \item{corr.order}{Order of the generalized correlation
  Dimension q. It must be greater than 1 (corr.order>1).
  Default, corr.order=2.}

  \item{n.points.radius}{The number of different radius
  where we shall estimate. C(r). Thus, we will estimate
  C(r) in n.points.radius between min.radius and
  max.radius.}

  \item{theiler.window}{Integer denoting the Theiler
  window: Two Takens' vectors must be separated by more
  than theiler.window time steps in order to be considered
  neighbours. By using a Theiler window, we exclude
  temporally correlated vectors from our estimations.}

  \item{do.plot}{Logical value. If TRUE (default value), a
  plot of the correlation sum is shown.}

  \item{number.boxes}{Number of boxes that will be used in
  the box assisted algorithm (see \link{neighbourSearch}).
  If the user does not specify it, the function uses a
  proper number of boxes.}

  \item{...}{Additional parameters.}

  \item{use.embeddings}{A numeric vector specifying which
  embedding dimensions should the \emph{estimate} function
  use to compute the correlation dimension.}

  \item{x}{A \emph{corrDim} object.}

  \item{regression.range}{Vector with 2 components denoting
  the range where the function will perform linear
  regression.}
}
\value{
  A \emph{corrDim} object that consist of a list with four
  components named \emph{radius}, \emph{embedding.dims},
  \emph{order} and \emph{corr.matrix}. \emph{radius} is a
  vector containing the different radius where we have
  evaluated C(r). \emph{embedding.dims} is a vector
  containing all the embedding dimensions in which we have
  estimated C(r). \emph{order} stores the order of the
  generalized correlation dimension that has been used.
  Finally, \emph{corr.matrix} stores all the correlation
  sums that have been computed. Each row stores the
  correlation sum for a concrete embedding dimension
  whereas each colum stores the correlation sum for a
  specific radius.

  The \emph{nlOrder} function returns the order of the
  correlation sum.

  The \emph{corrMatrix} function returns the correlations
  matrix storing the correlation sums that have been
  computed for all the embedding dimensions.

  The \emph{radius} function returns the radius on which
  the correlation sum function has been evaluated.

  The \emph{embeddingDims} function returns the embedding
  dimensions on which the correlation sum function has been
  evaluated.

  The \emph{plot} function shows two graphics of the
  correlation integral: a log-log plot of the correlation
  sum Vs the radius and the local slopes of
  \eqn{log10(C(r))\;Vs\;log10(C(r)).}{log10(C(r)) Vs
  log10(C(r)).}

  The \emph{estimate} function estimates the correlation
  dimension of the \emph{corr.dim} object by averaging the
  slopes of the embedding dimensions specified in the
  \emph{use.embeddings} parameter. The slopes are
  determined by performing a linear regression over the
  radius' range specified in \emph{regression.range}.If
  \emph{do.plot} is TRUE, a graphic of the regression over
  the data is shown.
}
\description{
  Functions for estimating the correlation sum and the
  correlation dimension of a dynamical system from
  1-dimensional time series using Takens' vectors.
}
\details{
  The correlation dimension is the most common measure of
  the fractal dimensionality of a geometrical object
  embedded in a phase space. In order to estimate the
  correlation dimension, the correlation sum is defined
  over the points from the phase space: \deqn{C(r) =
  \{(number\;of\;points\;(x_i,x_j)\;verifying\;that\;distance\;(x_i,x_j)<r\})/N^2}{C(r)
  = {number of points(xi,xj) verifying
  distance(xi,xj)<r}/N^2} However, this estimator is biased
  when the pairs in the sum are not statistically
  independent. For example, Taken's vectors that are close
  in time, are usually close in the phase space due to the
  non-zero autocorrelation of the original time series.
  This is solved by using the so-called Theiler window: two
  Takens' vectors must be separated by, at least, the time
  steps specified with this window in order to be
  considered neighbours. By using a Theiler window, we
  exclude temporally correlated vectors from our
  estimations.

  The correlation dimension is estimated using the slope
  obtained by performing a linear regression of
  \eqn{\log10(C(r))\;Vs.\;\log10(r)}{log10(C(r)) Vs.
  log10(r)}. Since this dimension is supposed to be an
  invariant of the system, it should not depend on the
  dimension of the Taken's vectors used to estimate it.
  Thus, the user should plot
  \eqn{\log10(C(r))\;Vs.\;\log10(r)}{log10(C(r)) Vs.
  log10(r)} for several embedding dimensions when looking
  for the correlation dimension and, if for some range
  \eqn{\log10(C(r))}{log10(C(r))} shows a similar linear
  behaviour in different embedding dimensions (i.e.
  parallel slopes), these slopes are an estimate of the
  correlation dimension. The \emph{estimate} routine allows
  the user to get always an estimate of the correlation
  dimension, but the user must check that there is a linear
  region in the correlation sum over different dimensions.
  If such a region does not exist, the estimation should be
  discarded.

  Note that the correlation sum C(r) may be interpreted as:
  \eqn{C(r) = <p(r)>,} that is: the mean probability of
  finding a neighbour in a ball of radius r surrounding a
  point in the phase space. Thus, it is possible to define
  a generalization of the correlation dimension by writing:
  \deqn{C_q(r) = <p(r)^{(q-1)}>}{Cq(r) = <p(r)^(q-1)>.}
  Note that the correlation sum \deqn{C(r) = C_2(r)}{C(r) =
  C2(r).}

  It is possible to determine generalized dimensions Dq
  using the slope obtained by performing a linear
  regression of \eqn{log10(Cq(r))\;Vs.\;(q-1)log10(r)}. The
  case q=1 leads to the information dimension, that is
  treated separately in this package (\link{infDim}). The
  considerations discussed for the correlation dimension
  estimate are also valid for these generalized dimensions.
}
\author{
  Constantino A. Garcia
}
\references{
  H. Kantz and T. Schreiber: Nonlinear Time series Analysis
  (Cambridge university press)
}

