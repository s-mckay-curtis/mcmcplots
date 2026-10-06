\name{convert.mcmc.list}
\alias{convert.mcmc.list}
\title{
  Convert an object to mcmc.list object
}

\description{
  Attempts to convert any object to an \code{mcmc.list} object. Intended for internal use only.
}

\usage{
convert.mcmc.list(x)
}

\arguments{
  \item{x}{
    an object.
  }
}

\details{
  Handles objects of class \code{mcmc}, \code{mcmc.list}, lists of matrices, 3D arrays with dimensions (iterations x chains x variables), and objects from the \pkg{posterior} package (\code{draws_array}, \code{draws_df}, etc.), as well as Stan and brms model outputs.
}

\value{
  An object of class \code{mcmc.list}.
}

\author{
  S. McKay Curtis
}

\keyword{ utilities }
