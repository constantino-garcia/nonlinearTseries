\name{recurrencePlot}
\alias{recurrencePlot}
\title{Recurrence Plot}
\usage{
  recurrencePlot(takens = NULL, time.series, embedding.dim,
    time.lag, radius)
}
\arguments{
  \item{time.series}{The original time series from which
  the phase-space reconstruction is performed.}

  \item{embedding.dim}{Integer denoting the dimension in
  which we shall embed the \emph{time.series}.}

  \item{time.lag}{Integer denoting the number of time steps
  that will be use to construct the Takens' vectors.}

  \item{takens}{Instead of specifying the
  \emph{time.series}, the \emph{embedding.dim} and the
  \emph{time.lag}, the user may specify directly the
  Takens' vectors.}

  \item{radius}{Maximum distance between two phase-space
  points to be considered a recurrence.}
}
\description{
  Plot the recurrence matrix.
}
\details{
  WARNING: This function is computationally very expensive.
  Use with caution.
}
\author{
  Constantino A. Garcia
}
\references{
  Zbilut, J. P. and C. L. Webber. Recurrence quantification
  analysis. Wiley Encyclopedia of Biomedical Engineering
  (2006).
}

