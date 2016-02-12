\name{BMMposteriori}
\alias{BMMposteriori}
\title{Plots aposteriori probabilities of data points}
\description{
  Given a \code{jags} object with model of class \code{BMMmodel} the
  aposteriori probabilities are determined. If \code{plot=TRUE}, the
  resulting object of class \code{BMMposteriori} is plotted.
}
\usage{
BMMposteriori(object, class, caption = NULL, plot = TRUE,
              auto.layout = TRUE, ...)
}
\arguments{
  \item{object}{a \code{jags} object with model of class \code{BMMmodel}.}
  \item{class}{a vector of integers indicating for which classes the
    posterior probabilities shall be plotted. The default is all.}
  \item{caption}{captions to appear above the plots.}
  \item{plot}{logical indicating if a plot shall be made.}
  \item{auto.layout}{logical: if \code{TRUE} puts all classes in the
    posterior probabilities plot on the same figure.}
  \item{\dots}{further graphical parameters may also be supplied as
    arguments.} 
}
\details{
  Given a \code{jags} object with model of class \code{BMMmodel} the a
  posteriori probabilities are calculated for the unique data points
  with respect to the components specified by \code{class}.
}
\value{
  There is a \code{BMMposteriori} object returned which is a list
  including the following components
  \item{data}{vector of unique data points.}
  \item{post}{a matrix including the posteriori probability of the data
    points for each class.} 
}
\author{Bettina Gruen}

\seealso{\code{\link{plot.BMMposteriori}}}
\keyword{hplot}

