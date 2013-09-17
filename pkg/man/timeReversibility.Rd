\name{timeReversibility}
\alias{timeReversibility}
\title{Time Reversibility statistic}
\usage{
  timeReversibility(time.series, tau)
}
\arguments{
  \item{time.series}{The time series used to compute the
  statistic}

  \item{tau}{Time delay used to compute the third order
  statistic.}
}
\value{
  The time simmetry statistic for the delays specified with
  \emph{tau}.
}
\description{
  Time Reversibility statistic
}
\details{
  The time simmetry statistic measures the asymmetry of a
  time series under time reversal by implementing the third
  order statistic: \deqn{E[s_n - s_{n-\tau})^3] }{E [s(n) -
  s(n-tau))^3]}. Since linear stochastic series are
  symmetric under time reversal, this statistic may be used
  for testing the assertion that the data was generated
  from a stationary linear stochastic process or not.
}
\author{
  Constantino A. Garcia
}

