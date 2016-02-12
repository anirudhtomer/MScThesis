\name{JAGScontrol}
\alias{JAGScontrol}
\title{Control parameters for the sampling.}
\description{
  The information on the number of burn-in draws and monitored draws is
  specified. Furthermore, it includes the information which variables
  shall be monitored and possibly a seed and a random number generator (RNG).
}
\usage{
JAGScontrol(variables, n.iter = 1000, thin = 1, burn.in = 0, seed,
            rng = c("base::Wichmann-Hill",
                    "base::Marsaglia-Multicarry",
                    "base::Super-Duper",
                    "base::Mersenne-Twister"))
}
\arguments{
  \item{variables}{names of variables which shall be monitored.}
  \item{n.iter}{number of monitored draws.}
  \item{thin}{thinning interval.}
  \item{burn.in}{number of discarded burn-in draws.}
  \item{seed}{integer setting the seed for the RNG.}
  \item{rng}{specification of RNG.}
}
\value{
  An object of class \code{JAGScontrol} is returned which is a list
  containing the specified information.
}
\author{Bettina Gruen}
\seealso{\code{\link{JAGSrun}}}
\examples{
control <- JAGScontrol(variables = "mu")
control
}
\keyword{utilities}

