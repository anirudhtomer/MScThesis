\name{NEWS}
\title{News for Package 'bayesmix'}
\section{Changes in bayesmix version 0.7-4}{
  \itemize{
    \item Functions from packages \pkg{stats} and \pkg{graphics} are now
      correctly imported.
  }
}
\section{Changes in bayesmix version 0.7-3}{
  \itemize{
    \item Packages \pkg{rjags} and \pkg{coda} moved from Depends to Imports. 
  }
}
\section{Changes in bayesmix version 0.7-2}{
  \itemize{
    \item Avoid to use methods \code{mean(<data.frame>)} and \code{sd(<data.frame>)}
      (deprecated in \R-2.14)
  }
}
\section{Changes in bayesmix version 0.7-1}{
  \itemize{
    \item \pkg{bayesmix} was modified to work with \pkg{rjags} >= 2.1.0
    with changes if no empty classes during sampling.
  }
}
\section{Changes in bayesmix version 0.7}{
  \itemize{
    \item \pkg{bayesmix} now uses functionality from package \pkg{rjags}
    for the MCMC simulations instead of calling the stand-alone software
    JAGS. 
  }
}
\section{Changes in bayesmix version 0.6-0}{
  \itemize{
    \item \pkg{bayesmix} modified to work with jags-1.0.2
  }
}
\section{Changes in bayesmix version 0.5-4}{
  \itemize{
    \item Dumping the data modified due to changes in \R.
  }
}
\section{Changes in \pkg{bayesmix} version 0.5-3}{
  \itemize{
    \item Patch from Martyn Plummer for jags-0.75 and \pkg{coda}-0.9-1. 
    \item First version released on CRAN: 0.5-2.
  }
}
