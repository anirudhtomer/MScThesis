\name{darwin}
\alias{darwin}
\docType{data}
\title{Differences in heights between plants}
\description{
  A numeric vector containing 15 observations of differences in heights
  between pairs of self-fertilized and cross-fertilized plants grown in
  the same condition.
}
\usage{data(darwin)}
\format{A numeric vector of length 15.
}
\details{
  Darwin's data set contains two extremely small values. Therefore, this
  data set can be used for outlier modelling.

}
\source{
  Abraham, B. and G. Box (1978)
  Linear models and spurious observations.
  \emph{Applied Statistics}, \bold{27}, 131--8.  
}
\examples{
data("darwin", package = "bayesmix")
## Estimated sample density
plot(density(darwin[[1]]), ylim = c(0, 0.02), main = "Outlier modelling")
ss <- seq(-100, 100, by = 1)
## Normal density with estimated mean and sd of whole sample
lines(ss, dnorm(ss, mean = mean(darwin[[1]]), sd = sd(darwin[[1]])), col = "red")
## Normal density with estimated mean and sd of sample, where the 2
## extremely small values are removed
lines(ss, dnorm(ss, mean = mean(darwin[-c(1:2),1]),
      sd = sd(darwin[-c(1:2),1])), col = "green")
}
\keyword{datasets}

