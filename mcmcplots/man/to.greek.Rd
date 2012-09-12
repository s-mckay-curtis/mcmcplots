\name{.to.greek}
\alias{.to.greek}
\title{
Parse strings containing names of greek letters.
}
\description{
If a character vector contains names of greek letters, apply the
\code{link{parse}} function.
}
\usage{
.to.greek(instr)
}
\arguments{
  \item{instr}{
	a character vector
}
}
\details{
  None.
}
\value{
  Parsed representation of the character vector
}
\author{
  Ilya Goldin
}
\note{
  None.
}

\seealso{
\code{\link{parse}}
}
\examples{
\dontrun{.to.greek(c("mybeta", "YourDelta"))}
}
\keyword{ utilities }
