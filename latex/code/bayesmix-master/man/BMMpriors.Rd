\name{BMMpriors}
\alias{BMMpriors}
\title{Creates a 'BMMpriors' object}
\description{This function enables a comfortable creation of
  \code{BMMpriors} objects, which can be used for specifying the priors 
  of a \code{BMMmodel}.
}
\usage{
BMMpriors(specification, y, eps = 10^-16)
}
\arguments{
  \item{specification}{named list including \code{kind}, \code{parameter},
    \code{hierarchical} and \code{mod}.}
  \item{y}{a numeric vector.}
  \item{eps}{a numeric value indicating the smallest value for flat priors.}
}
\details{In \code{specification} \code{kind} can be used for specifying
  if an "independent" or a "conditionallyconjugate"
  prior shall be used. \code{parameter} can be a named list of values for
  the prior variables or a function name (e.g., "priorsUncertain",
  "priorsFish", "priorsRaftery"). \code{hierarchical} can be \code{NULL}
  or "tau" if a hierarchical prior shall be taken for
  \eqn{\tau}. \code{mod} is a named list which provides the possibility
  to override the values from \code{parameter}. 
}
\value{
  There is an object \code{BMMpriors} returned with components
  \item{name}{vector indicating which kind of prior is specified and if
    it is an hierarchical prior and if appropriate with respect to which
    variable.}
  \item{var}{list of variables in the prior and their specified values.} 
}
\author{Bettina Gruen}
\examples{
  data("fish", package = "bayesmix")
  priors <- BMMpriors(y = fish)
}
\keyword{utilities}

