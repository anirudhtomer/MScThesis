\name{fish}
\alias{fish}
\docType{data}
\title{Fish length data}
\description{
  A numeric vector containing 256 observations of fish lengths.
}
\usage{data(fish)}
\format{A numeric vector of length 256.
}
\details{This data set can be used for modeling unobserved
  heterogeneity, as it can be assumed that underlying categories
  present in the data are the age groups to which the fish belong.  
}
\source{
  D. M. Titterington, A. F. M. Smith and U.E. Makov (1985)
  \emph{Statistical Analysis of Finite Mixture Distributions.}
  Wiley.
}
\examples{
data("fish", package = "bayesmix")
ss <- seq(-3, 13, by = 0.01)
hist(fish[[1]], 20, freq = FALSE, main = "Fish data")
lines(ss, dnorm(ss, mean(fish[[1]]), sd(fish[[1]])), col = "red")
}
\keyword{datasets}

