\name{JAGScall}
\alias{JAGScall}
\title{Calls jags}
\description{
  Using functionality from package rjags the JAGS model is created and
  posterior draws are sampled.
}
\usage{
JAGScall(model, y, prefix, control, ...)
}
\arguments{
  \item{model}{\code{JAGSmodel} object or output from \code{BMMmodel}.}
  \item{y}{a numeric vector.}
  \item{prefix}{character: prefix for .bug-file.}
  \item{control}{named list or \code{JAGScontrol} object.}
  \item{\dots}{additional parameters handed over to \code{BMMmodel}.}
}
\value{
  Returns a \code{"jags"} object.
}
\author{Bettina Gruen}
\seealso{\code{\link{JAGSrun}}}
\keyword{utilities}

