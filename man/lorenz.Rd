% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonLinearSystems.R
\name{lorenz}
\alias{lorenz}
\title{Lorenz system}
\usage{
lorenz(
  sigma = 10,
  beta = 8/3,
  rho = 28,
  start = c(-13, -14, 47),
  time = seq(0, 50, by = 0.01),
  do.plot = deprecated()
)
}
\arguments{
\item{sigma}{The \eqn{\sigma}{sigma} parameter. Default: 10.}

\item{beta}{The \eqn{\beta}{beta} parameter. Default: 8/3.}

\item{rho}{The \eqn{\rho}{rho} parameter. Default: 28.}

\item{start}{A 3-dimensional numeric vector indicating the starting point for the time series.
Default: c(-13, -14, 47).}

\item{time}{The temporal interval at which the system will be generated. 
Default: time=seq(0,50,by = 0.01).}

\item{do.plot}{Logical value. If TRUE, a plot of the 
generated Lorenz system is shown. Before version 0.2.11, default value was
TRUE; versions 0.2.11 and later use FALSE as default.}
}
\value{
A list with four vectors named \emph{time}, \emph{x}, \emph{y} 
and \emph{z} containing the time, the x-components, the 
y-components and the z-components of the Lorenz system, respectively.
}
\description{
Generates a 3-dimensional time series using the Lorenz equations.
}
\details{
The Lorenz system is a system of ordinary differential equations defined as:
\deqn{\dot{x} = \sigma(y-x)}{dx/dt = sigma*( y - x )}
\deqn{\dot{y} = \rho x-y-xz}{dy/dt = rho*x - y - xz}
\deqn{\dot{z} = -\beta z + xy}{dz/dt = -beta*z + xy}
The default selection for the system parameters (\eqn{\sigma=10, \rho=28, \beta=8/3}{sigma=10, rho=28, beta=8/3}) is known to
produce a deterministic chaotic time series.
}
\note{
Some initial values may lead to an unstable system that will tend to 
infinity.
}
\examples{
\dontrun{
lor=lorenz(time=seq(0,30,by = 0.01))
# plotting the x-component 
plot(lor$time,lor$x,type="l")
}
}
\references{
Strogatz, S.: Nonlinear dynamics and chaos: with applications 
to physics, biology, chemistry and engineering (Studies in Nonlinearity)
}
\seealso{
\code{\link{henon}, \link{logisticMap}, \link{rossler}, 
\link{ikedaMap}, \link{cliffordMap}, \link{sinaiMap}, \link{gaussMap}}
}
\author{
Constantino A. Garcia
}
