\name{initsFS}
\alias{initsFS}
\title{create initial values}
\description{
  Initial values for nodes are created after the suggestion in Sylvia
  Fruehwirth-Schnatter's book.
}
\usage{
initsFS(x, k, restrict, initialValues = list())
}
\arguments{
  \item{x}{a numeric vector.}
  \item{k}{number of segments.}
  \item{initialValues}{additional initial values specifications.}
  \item{restrict}{one of "none", "mu", "tau".}
}
\details{
  The initial values for \eqn{\mu} are determined by the quantiles of
  the data, those for \eqn{\eta} give equal weight on each segment and
  those for \eqn{\tau} are equal for all segments and estimated by the
  inverse of the IQR of the data divided by 1.34 and squared.
}
\value{
  A list with initial values for the parameter indicated by the
  name of the respective list element is returned. 
}
\author{Bettina Gruen}
\keyword{utilities}% at least one, from doc/KEYWORDS

