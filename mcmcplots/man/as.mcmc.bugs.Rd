\name{as.mcmc.bugs}
\alias{as.mcmc.bugs}

\title{Convert a bugs Object to an mcmc or mcmc.list Object}

\description{
  Converts a \code{bugs} object to an \code{mcmc} object.
}

\usage{
\method{as.mcmc}{bugs}(x)
}

\arguments{
  \item{x}{\code{bugs} object}
}

\value{
  If \code{x} contains multiple chains, the function returns an \code{mcmc.list} object.  Otherwise, the function returns an \code{mcmc} object.
}

\author{ S. McKay Curtis }

\seealso{
  \code{bugs} in \pkg{R2WinBUGS}
}

\examples{
\dontrun{
## Data object "schools.sim" generated from the examples
## in the bugs function of the R2WinBUGS package.
outmcmc <- as.mcmc(schools.sim)

## Gelman Rubin diagnostics
coda:::gelman.diag(outmcmc)
coda:::mcmc.plot(outmcmc)
}
}

\keyword{ manip }
