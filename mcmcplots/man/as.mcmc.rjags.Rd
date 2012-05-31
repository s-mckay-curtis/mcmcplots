\name{as.mcmc.rjags}
\alias{as.mcmc.rjags}

\title{Convert an rjags Object to an mcmc or mcmc.list Object.}

\description{
  This is a function to convert and object produced by the \code{jags} function in the \pkg{R2jags} package.   If the \code{rjags} object contains more than one chain, then this function will return an \code{mcmc.list} object. Otherwise, the function returns an \code{mcmc} object.
}

\usage{
\method{as.mcmc}{rjags}(x)
}

\arguments{
  \item{x}{ \code{rjags} object.}
}

\details{
  No details.
}

\value{
  An \code{mcmc.list} or \code{mcmc} object.
}

\references{
  None.
}

\author{ S. McKay Curtis }

\note{
  None.
}

\seealso{
  \pkg{coda} package.
}

\examples{
## None ##
}

\keyword{manip}
