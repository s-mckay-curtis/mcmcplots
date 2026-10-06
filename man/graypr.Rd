\name{.graypr}
\alias{.graypr}
\alias{.cleanpr}
\title{Create a Styled Plotting Region}
\description{
  Functions to draw styled plotting regions: \code{.cleanpr} creates a modern clean canvas with subtle gridlines, while \code{.graypr} colors the plotting region gray. These functions are intended for internal use only.
}
\usage{
.graypr(x.axis = TRUE, y.axis = TRUE,
x.major = TRUE, y.major = TRUE,
x.minor = TRUE, y.minor=TRUE,
x.malty = 1, y.malty = 1,
x.milty = 1, y.milty = 1)
}
\arguments{
  \item{x.axis}{ if \code{TRUE} an x axis will be plotted.}
  \item{y.axis}{ if \code{TRUE} a y axis will be plotted.}
  \item{x.major}{ if \code{TRUE}, major gridlines on the x axis will be plotted.}
  \item{y.major}{ if \code{TRUE}, major gridlines on the y axis will be plotted.}
  \item{x.minor}{ if \code{TRUE}, minor gridlines on the x axis will be plotted.}
  \item{y.minor}{ if \code{TRUE}, minor gridlines on the y axis will be plotted.}
  \item{x.malty}{ line type to be used on the x-axis, major gridlines.}
  \item{y.malty}{ line type to be used on the y-axis, major gridlines.}
  \item{x.milty}{ line type to be used on the x-axis, minor gridlines.}
  \item{y.milty}{ line type to be used on the y-axis, minor gridlines.}
}
\seealso{
  \code{rect}
}
\author{S. McKay Curtis}
\keyword{aplot}
