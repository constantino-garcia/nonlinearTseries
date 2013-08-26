\name{rossler}
\alias{rossler}
\title{Rossler system}
\usage{
  rossler(a = 0.2, b = 0.2, c = 5.7,
    start = c(-2, -10, 0.2), time = seq(0, 50, by = 0.01),
    do.plot = TRUE)
}
\arguments{
  \item{start}{A 3-dimensional numeric vector indicating
  the starting point for the time series. Default: c(-2,
  -10, 0.2).}

  \item{a}{The \emph{a} parameter. Default:0.2.}

  \item{b}{The \emph{b} parameter. Default: 0.2.}

  \item{c}{The \emph{c} parameter. Default: 5.7.}

  \item{time}{The temporal interval at which the system
  will be generated. Default: time=seq(0,50,by = 0.01).}

  \item{do.plot}{Logical value. If TRUE (default value), a
  plot of the generated Lorenz system is shown.}
}
\value{
  A list with four vectors named \emph{time}, \emph{x},
  \emph{y} and \emph{z} containing the time, the
  x-components, the y-components and the z-components of
  the Rossler system, respectively.
}
\description{
  Generates a 3-dimensional time series using the Rossler
  equations.
}
\details{
  The Rossler system is a system of ordinary differential
  equations defined as: \deqn{\dot{x} = -(y + z)}{dx/dt =
  -(y + z)} \deqn{\dot{y} = x+a \cdot y}{dy/dt = x + a*y}
  \deqn{\dot{z} = b + z*(x-c)}{dz/dt = b + z*(x-c)} The
  default selection for the system parameters (\emph{a} =
  0.2, \emph{b} = 0.2, \emph{c} = 5.7) is known to produce
  a deterministic chaotic time series.
}
\note{
  Some initial values may lead to an unstable system that
  will tend to infinity.
}
\examples{
\dontrun{
r.ts = rossler(time=seq(0,30,by = 0.01))
}
}
\author{
  Constantino A. Garcia
}
\references{
  Strogatz, S.: Nonlinear dynamics and chaos: with
  applications to physics, biology, chemistry and
  engineering (Studies in Nonlinearity)
}
\seealso{
  \code{\link{henon}, \link{logisticMap}, \link{rossler},
  \link{ikedaMap}, \link{cliffordMap}, \link{sinaiMap},
  \link{gaussMap}}
}

