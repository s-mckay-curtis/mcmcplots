\name{as.mcmc.bugs}
\alias{as.mcmc.bugs}

\title{Convert a bugs Object to an mcmc or mcmc.list Object}

\description{
This function converts a \code{bugs} object to an \code{mcmc} object so that the functions of the \code{coda} package may be used to analyze the posterior draws in the \code{bugs} object.
}

\usage{
\method{as.mcmc}{bugs}(x)
}

\arguments{
  \item{x}{\code{bugs} object}
}

\details{
  None.
}

\value{
  If \code{x} contains multiple chains, the function returns an \code{mcmc.list} object.  Otherwise, the function returns an \code{mcmc} object.
}

\references{
  None.
}

\author{ S. McKay Curtis }

\note{
  None.
}

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
