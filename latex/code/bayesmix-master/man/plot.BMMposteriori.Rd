\name{plot.BMMposteriori}
\alias{plot.BMMposteriori}
\title{Plots a posteriori probabilities of data points}
\description{
  Plot method for object of class \code{BMMposteriori}, typically
  called by \code{BMMposteriori}.
}
\usage{
\method{plot}{BMMposteriori}(x, caption, main = "", ...)
}
\arguments{
  \item{x}{a \code{BMMposteriori} object.}
  \item{caption}{captions to appear above the plots.}
  \item{main}{title to each plot-in addition to the above \code{caption}.}
  \item{\dots}{further graphical parameters may also be supplied as
    arguments.} 
}
\details{
  This function is called by \code{BMMposteriori} if \code{plot =
    TRUE}. 
}
\author{Bettina Gruen}

\seealso{\code{\link{BMMposteriori}}}
\keyword{hplot}

