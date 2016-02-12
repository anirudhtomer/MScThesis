\name{priors}
\alias{priorsFish}
\alias{priorsRaftery}
\alias{priorsUncertain}
\title{Creates list of prior specifications}
\description{Given the data values for the priors are determined.
}
\usage{
priorsFish(y, eps = 10^-16)
priorsRaftery(y)
priorsUncertain(y, eps = 10^-16)

}
\arguments{
  \item{y}{a numeric vector.}
  \item{eps}{a numeric value indicating the smallest value for flat priors.}
}
\details{Values for the prior parameter \code{b0}, \code{B0}, \code{nu0}
  and \code{S0} are determined.
}
\value{
  There is a list returned with named components of the prior parameters.
}
\author{Bettina Gruen}
\keyword{utilities}

