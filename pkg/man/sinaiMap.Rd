\name{sinaiMap}
\alias{sinaiMap}
\title{Sinai map}
\usage{
  sinaiMap(a = 0.1, start = runif(2), n.sample = 5000,
    n.transient = 500, do.plot = TRUE)
}
\arguments{
  \item{start}{A 2-dimensional vector indicating the
  starting values for the x and y Sinai coordinates. If the
  starting point is not specified, it is generated
  randomly.}

  \item{a}{The \emph{a} parameter. Default: 0.1}

  \item{n.sample}{Length of the generated time series.
  Default: 5000 samples.}

  \item{n.transient}{Number of transient samples that will
  be discarded. Default: 500 samples.}

  \item{do.plot}{Logical value. If TRUE (default value), a
  plot of the generated Sinai system is shown.}
}
\value{
  A list with two vectors named x and y containing the
  x-components and the y-components of the Sinai map,
  respectively.
}
\description{
  Generates a 2-dimensional time series using the Sinai
  map.
}
\details{
  The Sinai map is defined as follows: \deqn{x_{n+1} =
  (x_{n} + y_{n} + a \cdot cos(2 \cdot pi \cdot y_{n}) )mod
  1}{x[n+1] = (x[n] + y[n] + a*cos(2*pi*y[n]) )mod 1}
  \deqn{y_{n+1} = (x_{n} + 2 \cdot y_n)mod 1}{y[n+1] =
  (x[n] + 2*y[n])mod 1} The default selection for the
  \emph{a} parameter is known to produce a deterministic
  chaotic time series.
}
\note{
  Some initial values may lead to an unstable system that
  will tend to infinity.
}
\examples{
\dontrun{
sinai.map = sinaiMap(n.sample = 1000, n.transient=10,do.plot=TRUE)
# accessing the x coordinate and plotting it
plot(ts(sinai.map$x))
}
}
\author{
  Constantino A. Garcia
}
\references{
  Mcsharry, P. E. and P. R. Ruffino (2003). Asymptotic
  angular stability in nonlinear systems: rotation numbers
  and winding numbers. Dynamical Systems 18(3), 191-200.
}
\seealso{
  \code{\link{henon}, \link{logisticMap}, \link{lorenz},
  \link{rossler}, \link{ikedaMap}, \link{cliffordMap},
  \link{gaussMap}}
}

