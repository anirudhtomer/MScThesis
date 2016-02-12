\name{BMMmodel}
\alias{BMMmodel}
\alias{print.JAGSmodel}
\title{Creates bugs model.} 
\description{
  Creates the BUGS-model specification and return the values for the
  initialization, the prior specification and the observations.
}
\usage{
BMMmodel(y, k, priors, inits = "initsFS", aprioriWeights = 1,
         no.empty.classes = FALSE, restrict = "none", ...)
}
\arguments{
  \item{y}{a numeric vector.}
  \item{k}{integer indicating the number of segments.}
  \item{priors}{specification of priors by a named list or a
    \code{BMMpriors} object.}
  \item{inits}{specification of initial values by a named list or string
    indicating the function to be called.}
  \item{aprioriWeights}{specification of prior of the a-priori weights. If
    \code{aprioriWeights} does not have length = k, there is an equal
    prior for the a-priori weights assumed.}
  \item{no.empty.classes}{logical: should it be prevented that empty
    classes arise during sampling.}
  \item{restrict}{one of "none", "mu", "tau".}
  \item{...}{further parameters for the function specified in
    \code{inits}.} 
}
\details{By default the function \code{initsFS} is called for generating
  initial values. Any other function specified by \code{inits} is
  assumed to have at least \code{x}, \code{k} and \code{restrict} as
  input parameters.

  The parameter \code{restrict} indicates if a
  location-shift model ("tau"), a scale contaminated model
  ("mu") or a model where both variables vary over components shall
  be fitted.

  If the logical \code{no.empty.classes} is \code{TRUE} there
  are observations added to the model that the classes are not empty. This
  signifies that the likelihood when sampling the class affiliations is
  changed thus that any data point which is sampled and is the last one
  in its class stays there.  
}
\value{
  If \code{y} is specified an object of class \code{BMMmodel} is
  returned with components:
  \item{inits}{named list.}
  \item{data}{named list.}
  \item{bugs}{text for .bug-file with prefix missing.}

  If \code{y} is missing an object of class \code{BMMsetup} is
  returned containing the parameter specifications. When \code{JAGScall}
  is called with this object as model argument, \code{BMMmodel} is
  called with \code{y} and the other parameters as input arguments
  before calling JAGS.
}
\author{Bettina Gruen}
\seealso{\code{\link{JAGSrun}}, \code{\link{initsFS}}}
\examples{
data("fish", package = "bayesmix")
model <- BMMmodel(fish, k = 4, priors = list(kind = "independence",
                  parameter = "priorsFish", hierarchical = "tau"),
                  initialValues = list(S0 = 2))
model
}
\keyword{utilities}

