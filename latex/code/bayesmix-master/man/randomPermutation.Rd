\name{randomPermutation}
\alias{randomPermutation}
\title{Randomly permute segments for MCMC draws}
\description{
  Random permutation of segment labels for each draw in order to get a
  better estimate of the unrestricted likelihood.
}
\usage{
randomPermutation(x)
}
\arguments{
  \item{x}{a \code{jags} object with model of class \code{BMMmodel}.}
}
\details{
  The draws are permuted with respect to the different classes
  \code{k}. 
}
\value{
  The input object with permuted results for each draw is returned.
}
\author{Bettina Gruen}
\section{Warning}{Any variables where there are neither \code{k}
  different chains nor only one chain observed are dropped.} 
\keyword{utilities}


