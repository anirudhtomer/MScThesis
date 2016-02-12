\name{BMMdiag}
\alias{BMMdiag}
\title{Plot Identifiability Diagnostics for jags Object}
\description{
  Two different plots are currently provided: a plot of different
  variables against each other and a plot of the same variable against
  its values in the other classes.
}
\usage{
BMMdiag(object, which = 1:2, variables, ask = interactive(), fct1,
        fct2, xlim, ylim, auto.layout = TRUE, caption = NULL,
        main = "", ...)
}

\arguments{
  \item{object}{a \code{jags} object with model of class \code{BMMmodel}.}
  \item{which}{if only one of the plots is required, specify its number.}
  \item{variables}{if \code{variables} is missing, the names are taken
    from the \code{jags} object.}
  \item{ask}{prompt user before each page of plots}
  \item{fct1}{string: name of transformation function for variable on x-axis.}
  \item{fct2}{string: name of transformation function for variable on y-axis.}
  \item{xlim}{if no range for \code{xlim} is specified, a sensible range is
    taken.}
  \item{ylim}{if no range for \code{ylim} is specified, a sensible range is
    taken.}
  \item{auto.layout}{logical: if \code{TRUE} puts each of the two
    different plots on one figure.}
  \item{caption}{captions to appear above the plots.}
  \item{main}{title to each plot (in addition to the above 'caption').}
  \item{...}{further graphical parameters (see 'plot.xy' and 'par') may
    also be supplied as arguments.}
}
      
\details{The plots help determining which variable will induce a unique
  labelling when taken for ordering of the segments and indicate if the
  model is overfitted by specifying too  many segments.
}
\author{Bettina Gruen}
\keyword{hplot}
