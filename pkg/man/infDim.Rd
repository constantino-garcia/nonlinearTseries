\name{infDim}
\alias{embeddingDims.infDim}
\alias{estimate.infDim}
\alias{fixedMass.infDim}
\alias{infDim}
\alias{logRadius.infDim}
\alias{plot.infDim}
\title{Information dimension}
\usage{
  infDim(time.series, min.embedding.dim = 2,
    max.embedding.dim = min.embedding.dim, time.lag = 1,
    min.fixed.mass, max.fixed.mass,
    number.fixed.mass.points = 10, radius,
    increasing.radius.factor = sqrt(2),
    number.boxes = NULL, number.reference.vectors,
    theiler.window = 1, kMax = 100, do.plot = TRUE)

  \method{fixedMass}{infDim} (x)

  \method{fixedMass}{infDim} (x)

  \method{embeddingDims}{infDim} (x)

  \method{estimate}{infDim} (x, regression.range = NULL,
    do.plot = TRUE, use.embeddings = NULL, ...)

  \method{plot}{infDim} (x, ...)
}
\arguments{
  \item{time.series}{The original time series from which
  the information dimension will be estimated.}

  \item{min.embedding.dim}{Integer denoting the minimum
  dimension in which we shall embed the time.series (see
  \link{buildTakens}).}

  \item{max.embedding.dim}{Integer denoting the maximum
  dimension in which we shall embed the time.series (see
  \link{buildTakens}).Thus, we shall estimate the
  information dimension between \emph{min.embedding.dim}
  and \emph{max.embedding.dim}.}

  \item{time.lag}{Integer denoting the number of time steps
  that will be use to construct the Takens' vectors (see
  \code{\link{buildTakens}}).}

  \item{min.fixed.mass}{Minimum percentage of the total
  points that the algorithm shall use for the estimation.}

  \item{max.fixed.mass}{Maximum percentage of the total
  points that the algorithm shall use for the estimation.}

  \item{number.fixed.mass.points}{The number of different
  \emph{fixed mass} fractions between \emph{min.fixed.mass}
  and \emph{max.fixed.mass} that the algorithm will use for
  estimation.}

  \item{radius}{Initial radius for searching neighbour
  points in the phase space. Ideally, it should be small
  enough so that the fixed mass contained in this radius is
  slightly greater than the \emph{min.fixed.mass}. However,
  whereas the radius is not too large (so that the
  performance decreases) the choice is not critical.}

  \item{increasing.radius.factor}{Numeric value. If no
  enough neighbours are found within \emph{radius}, the
  radius is increased by a factor
  \emph{increasing.radius.factor} until succesful. Default:
  sqrt(2) = 1.414214.}

  \item{number.boxes}{Number of boxes that will be used in
  the box assisted algorithm (see \link{neighbourSearch}).}

  \item{number.reference.vectors}{Number of reference
  points that the routine will try to use, saving
  computation time.}

  \item{theiler.window}{Integer denoting the Theiler
  window: Two Takens' vectors must be separated by more
  than theiler.window time steps in order to be considered
  neighbours. By using a Theiler window, we exclude
  temporally correlated vectors from our estimations.}

  \item{kMax}{Maximum number of neighbours used for
  achieving p with all the points from the time series (see
  Details). Default: 100.}

  \item{do.plot}{Logical value. If TRUE (default value), a
  plot of the correlation sum is shown.}

  \item{x}{A \emph{infDim} object.}

  \item{regression.range}{Vector with 2 components denoting
  the range where the function will perform linear
  regression.}

  \item{use.embeddings}{A numeric vector specifying which
  embedding dimensions should the \emph{estimate} function
  use to compute the information dimension.}

  \item{...}{Additional parameters.}
}
\value{
  A \emph{infDim} object that consist of a list with two
  components: \emph{log.radius} and \emph{fixed.mass}.
  \emph{log.radius} contains the average log10(radius) in
  which the \emph{fixed.mass} can be found.

  The \emph{fixedMass} function returns the fixed mass
  vector used in the information dimension algorithm.

  The \emph{logRadius} function returns the average
  log(radius) computed on the information dimension
  algorithm.

  The \emph{embeddingDims} function returns the embeddings
  in which the information dimension was computed

  The 'estimate' function estimates the information
  dimension of the 'infDim' object by by averaging the
  slopes of the embedding dimensions specified in the
  \emph{use.embeddings} parameter. The slopes are
  determined by performing a linear regression over the
  fixed mass' range specified in 'regression.range'. If
  do.plot is TRUE, a graphic of the regression over the
  data is shown.

  The 'plot' function shows two graphics of the information
  dimension estimate: a graphic of <log10(radius)> Vs fixed
  mass and a graphic of the local slopes of the information
  dimension Vs the fixed mass, both in a semi-log scale.
}
\description{
  Functions for estimating the information dimension of a
  dynamical system from 1-dimensional time series using
  Takens' vectors
}
\details{
  The information dimension is a particular case of the
  generalized correlation dimension when setting the order
  q = 1. It is possible to demonstrate that the information
  dimension \eqn{D_1}{D1} may be defined as:
  \eqn{D_1=lim_{r \rightarrow 0} <\log
  p(r)>/\log(r)}{D1=lim{r->0} <ln p(r)>/ln(r)}. Here,
  \eqn{p(r)} is the probability of finding a neighbour in a
  neighbourhood of size \eqn{r} and <> is the mean value.
  Thus, the information dimension specifies how the average
  Shannon information scales with the radius \eqn{r}. The
  user should compute the information dimension for
  different embedding dimensions for checking if
  \eqn{D_1}{D1} saturates.

  In order to estimate \eqn{D_1}{D1}, the algorithm looks
  for the scaling behaviour of the the average radius that
  contains a given portion (a "fixed-mass") of the total
  points in the phase space. By performing a linear
  regression of \eqn{\log(p)\;Vs.\;\log(<r>)}{ln p Vs ln
  <r>} (being \eqn{p} the fixed-mass of the total points),
  an estimate of \eqn{D_1}{D1} is obtained.

  The algorithm also introduces a variation of \eqn{p} for
  achieving a better performance: for small values of
  \eqn{p}, all the points in the time series (\eqn{N}) are
  considered for obtaining \eqn{p=n/N}. Above a maximum
  number of neighbours \eqn{kMax}, the algorithm obtains
  \eqn{p} by decreasing the number of points considerd from
  the time series \eqn{M<N}. Thus \eqn{p = kMax/M}.

  Even with these improvements, the calculations for the
  information dimension are heavier than those needed for
  the correlation dimension.
}
\author{
  Constantino A. Garcia
}
\references{
  H. Kantz and T. Schreiber: Nonlinear Time series Analysis
  (Cambridge university press)
}
\seealso{
  \code{\link{corrDim}}.
}

