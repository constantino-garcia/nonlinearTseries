\name{rqa}
\alias{rqa}
\title{Recurrence Quantification Analysis (RQA)}
\usage{
  rqa(takens = NULL, time.series = NULL, embedding.dim = 2,
    time.lag = 1, radius, lmin = 2, vmin = 2,
    do.plot = FALSE, distanceToBorder = 2)
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

  \item{lmin}{Minimal length of a diagonal line to be
  considered in the RQA. Default \emph{lmin} = 2.}

  \item{vmin}{Minimal length of a vertical line to be
  considered in the RQA. Default \emph{vmin} = 2.}

  \item{do.plot}{Logical. If TRUE, the recurrence plot is
  shown. However, plotting the recurrence matrix is
  computationally expensive. Use with caution.}

  \item{distanceToBorder}{In order to avoid border effects,
  the \emph{distanceToBorder} points near the border of the
  recurrence matrix are ignored when computing the RQA
  parameters. Default, \emph{distanceToBorder} = 2.}
}
\value{
  A \emph{rqa} object that consist of a list with the most
  important RQA parameters: \itemize{ \item \emph{REC}:
  Recurrence. Percentage of recurrence points in a
  Recurrence Plot.  \item \emph{DET}: Determinism.
  Percentage of recurrence points that form diagonal lines.
  \item \emph{LAM}: Percentage of recurrent points that
  form vertical lines.  \item \emph{RATIO}: Ratio between
  \emph{DET} and \emph{RR}.  \item \emph{Lmax}: Length of
  the longest diagonal line.  \item \emph{Lmean}: Mean
  length of the diagonal lines. The main diagonal is not
  taken into account.  \item \emph{DIV}: Inverse of
  \emph{Lmax}.  \item \emph{Vmax}: Longest vertical line.
  \item \emph{Vmean}: Average length of the vertical lines.
  This parameter is also referred to as the Trapping time.
  \item \emph{ENTR}: Shannon entropy of the diagonal line
  lengths distribution \item \emph{TREND}: Trend of the
  number of recurrent points depending on the distance to
  the main diagonal \item \emph{diagonalHistogram}:
  Histogram of the length of the diagonals.  \item
  \emph{recurrenceRate}: Number of recurrent points
  depending on the distance to the main diagonal. }
}
\description{
  The Recurrence Quantification Analysis (RQA) is an
  advanced technique for the nonlinear analysis that allows
  to quantify the number and duration of the recurrences in
  the phase space.
}
\examples{
\dontrun{
rossler.ts =  rossler(time=seq(0, 10, by = 0.01),do.plot=FALSE)$x
rqa.params=rqa(time.series = rossler.ts, embedding.dim=2, time.lag=1,
               radius=1.2,lmin=2,do.plot=FALSE,distanceToBorder=2)
               }
}
\author{
  Constantino A. Garcia
}
\references{
  Zbilut, J. P. and C. L. Webber. Recurrence quantification
  analysis. Wiley Encyclopedia of Biomedical Engineering
  (2006).
}

