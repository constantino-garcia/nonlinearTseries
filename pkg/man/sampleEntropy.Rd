\name{sampleEntropy}
\alias{estimate.sampleEntropy}
\alias{getSampleEntropy}
\alias{plot.sampleEntropy}
\alias{sampleEntropy}
\title{Sample Entropy (also known as Kolgomorov-Sinai Entropy)}
\usage{
  sampleEntropy(corrDim.object, do.plot = TRUE)

  getSampleEntropy(x)

  \method{plot}{sampleEntropy} (x, ...)

  \method{estimate}{sampleEntropy} (x,
    regression.range = NULL, do.plot = TRUE, ...)
}
\arguments{
  \item{corrDim.object}{A \emph{corrDim} object from which
  the Sample Entropy of the time series characterized by
  \emph{corrDim} shall be estimated.}

  \item{do.plot}{do.plot Logical value. If TRUE (default
  value), a plot of the sample entropy is shown.}

  \item{...}{Additional parameters.}

  \item{x}{A \emph{sampleEntropy} object.}

  \item{regression.range}{Vector with 2 components denoting
  the range where the function will perform linear
  regression.}
}
\value{
  A \emph{sampleEntropy} object that contains a list
  storing the sample entropy (\emph{sample.entropy}), the
  embedding dimensions ( \emph{embedding.dims}) and radius
  (\emph{radius}) for which the sample entropy has been
  computed, and the order of the sample entropy
  (\emph{order}). The sample entropy is stored as a matrix
  in which each row contains the computations for a given
  embedding dimension and each column stores the
  computations for a given radius.

  The \emph{getSampleEntropy} returns the sample entropy
  function depending of the radius used for the
  computations.

  The \emph{plot} function shows the graphics for the
  sample entropy.

  The \emph{estimate} function returns a vector storing the
  sample entropy estimate for each embedding dimension.
}
\description{
  The Sample Entropy measures the complexity of a time
  series. Large values of the Sample Entropy indicate high
  complexity whereas that smaller values characterize more
  regular signals.
}
\details{
  The sample entropy is computed using: \deqn{h_q(m,r) =
  log(C_q(m,r)/C_{q}(m+1,r))}{hq(m,r) =
  log(Cq(m,r)/Cq(m+1,r)),} where \emph{m} is the embedding
  dimension and \emph{r} is the radius of the
  neighbourhood. When computing the correlation dimensions
  we use the linear regions from the correlation sums in
  order to do the estimates. Similarly, the sample entropy
  \eqn{h_q(m,r)}{hq(m,r)} should not change for both
  various \emph{m} and \emph{r}.

  For each embedding dimension the sample entropy is
  estimated by averaging \deqn{h_q(m,r) =
  log(C_q(m,r)/C_{q}(m+1,r))}{hq(m,r) =
  log(Cq(m,r)/Cq(m+1,r))} over the range specified by
  \emph{regression range} in the \emph{estimate} function.
}
\examples{
\dontrun{
h=henon(n.sample = 15000, n.transient = 100, a = 1.4, b = 0.3,
start = c(0.78,0.8165), do.plot = FALSE)
gen.corr.dim=corrDim(time.series=h$x,min.embedding.dim=2,max.embedding.dim=9,
                     corr.order=2, time.lag=1,min.radius=0.025,
                     max.radius=0.01,n.points.radius=20, do.plot=FALSE,
                     theiler.window=10,number.boxes=100)
se=sampleEntropy(gen.corr.dim, do.plot=FALSE)
estimate(se)}
}
\author{
  Constantino A. Garcia
}
\references{
  H. Kantz and T. Schreiber: Nonlinear Time series Analysis
  (Cambridge university press)
}

