\name{spaceTimePlot}
\alias{contourLines.spaceTimePlot}
\alias{plot.spaceTimePlot}
\alias{spaceTimePlot}
\title{Space Time plot}
\usage{
  spaceTimePlot(takens = NULL, time.series = NULL,
    embedding.dim = 2, time.lag = 1, max.radius = NULL,
    time.step = 1, number.time.steps = NULL,
    numberPercentages = 10, do.plot = TRUE)

  \method{contourLines}{spaceTimePlot} (x)

  \method{plot}{spaceTimePlot} (x, ...)
}
\arguments{
  \item{time.series}{The original time series being
  analyzed.}

  \item{embedding.dim}{Integer denoting the dimension in
  which we shall embed the time series.}

  \item{time.lag}{Integer denoting the number of time steps
  that will be use to construct the Takens' vectors.}

  \item{takens}{Instead of specifying the
  \emph{time.series}, the \emph{embedding.dim} and the
  \emph{time.lag}, the user may specify directly the
  Takens' vectors.}

  \item{max.radius}{Maximum neighbourhood radius in which
  the algorithm will look for finding neighbours. This
  parameter may be used to avoid heavy computations. If the
  user does not specify a radius, the algorithm estimates
  it.}

  \item{time.step}{Integer denoting the number of discrete
  steps between two calculations of the space time plot.}

  \item{number.time.steps}{Integer denoting the number of
  temporal jumps in steps of \emph{time.step} in which we
  want to compute the space time plot.}

  \item{numberPercentages}{Number of contour lines to be
  computed. Each contour line represent a concrete
  percentage of points (see Details).}

  \item{do.plot}{Logical. If TRUE, the time space plot is
  shown.}

  \item{x}{A \emph{spaceTimePlot} object.}

  \item{...}{Additional parameters.}
}
\value{
  A \emph{timeSpacePlot} object that consist, essentially,
  of a matrix storing the values for each contour line.
  Each row stores the value for a given percentage of
  points. Each column stores the value of the radius you
  have to go to find a given fraction of neighbour pairs
  (the rows), depending on their temporal separation (the
  colums). This matrix can be accessed by using the
  \emph{contourlines} method.

  The \emph{contourLines} function returns the contour
  lines of the space time plot.
}
\description{
  The space time separation is a broadly-used method of
  detecting non-stationarity and temporal correlations in
  the time series being analyzed. The space time separation
  plot is also used to select a proper Theiler window by
  selecting a temporal separation enough to saturate the
  contour lines.
}
\details{
  Each contour line of the space time plot indicate the
  distance you have to go (y-axis) to find a given fraction
  of neighbour pairs, depending on their temporal
  separation (x-axis).

  WARNING: The parameter \emph{number.time.steps} should be
  used with caution since this method performs heavy
  computations.
}
\examples{
\dontrun{
 tak = buildTakens(sin(2*pi*0.005*(0:5000)),2,1)
 stp.test = spaceTimePlot(takens=tak,number.time.steps=400,do.plot=TRUE)
 }
}
\author{
  Constantino A. Garcia
}
\references{
  H. Kantz and T. Schreiber: Nonlinear Time series Analysis
  (Cambridge university press)
}

