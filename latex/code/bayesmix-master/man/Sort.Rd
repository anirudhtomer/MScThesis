\name{Sort}
\alias{Sort}
\title{Sorts MCMC chains according to certain variables}
\description{Ascending sorting of results of \code{jags} object with
  model of class \code{BMMmodel} with respect to a given variable. 
}
\usage{
Sort(x, by = NULL)
}
\arguments{
  \item{x}{a \code{jags} object with model of class \code{BMMmodel}.}
  \item{by}{variable name according to which the segments shall be
    ordered.}
}
\details{If \code{by} is not specified, the first variable in the
  corresponding vector of the \code{jags} object is taken.
}
\value{
  The input object with results sorted in ascending order according
  to the variable given in \code{by} is returned.  
}
\author{Bettina Gruen}
\section{Warning}{If there arise problems, the original object is
  returned with a warning.}
\keyword{utilities}

