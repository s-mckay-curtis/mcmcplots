\name{caterpoints}
\alias{caterpoints}
\title{
  Points on a "caterplot"
}
\description{
  Adds points to a caterplot.
}

\usage{
caterpoints(x, parnames, horizontal = TRUE, ...)
}

\arguments{
  \item{x}{ vector of points to add to a \code{caterplot}. }
  \item{parnames}{ an optional vector of parameter names.  If specified, \code{x} must have a \code{names} attribute.   The argument \code{parnames} will be used to subset the vector \code{x} for plotting, as in \code{x[parnames]}. }
  \item{horizontal}{ logical value that should match the argument of the same name in the original call to \code{caterplot}.}
  \item{\dots}{ further arguments passed to the function \code{points}. }
}

\details{
  None.
}

\value{
  None.
}

\references{
  None.
}

\author{
  S. McKay Curtis
}

\note{
  None.
}

\seealso{
  \code{\link{caterplot}}
}

\examples{
## Create fake MCMC output
nc <- 10; nr <- 1000
pnames <- c(paste("alpha[", 1:5, "]", sep=""), paste("gamma[", 1:5, "]", sep=""))
means <- rpois(10, 20)
fakemcmc <- as.mcmc.list(lapply(1:3, function(i) mcmc(matrix(rnorm(nc*nr, rep(means, each=nr)), nrow=nr, dimnames=list(NULL,pnames)))))

## caterplot plots of the fake MCMC output
par(mfrow=c(2,2))
caterplot(fakemcmc, "alpha", collapse=FALSE)
caterpoints(runif(5, -0.1, 0.1), pch="x", col="red")
caterplot(fakemcmc, "alpha", horizontal=FALSE)
caterpoints(runif(5, -0.1, 0.1), horizontal=FALSE, pch="x", col="red")
}

\keyword{aplot}
