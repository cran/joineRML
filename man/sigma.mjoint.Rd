% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigma.mjoint.R
\name{sigma.mjoint}
\alias{sigma.mjoint}
\title{Extract residual standard deviation(s) from an \code{mjoint} object}
\usage{
\method{sigma}{mjoint}(object, ...)
}
\arguments{
\item{object}{an object inheriting from class \code{mjoint} for a joint model
of time-to-event and multivariate longitudinal data.}

\item{...}{additional arguments; currently none are used.}
}
\value{
a number (standard deviation) if \eqn{K = 1} (univariate model), or a
  vector if \eqn{K>1} (multivariate model).
}
\description{
Extract residual standard deviation(s) from an \code{mjoint}
  object.
}
\references{
Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
Springer Verlag; 2000.
}
\seealso{
\code{\link[lme4]{sigma}} in the \strong{lme4} package.
}
\author{
Graeme L. Hickey (\email{graemeleehickey@gmail.com})
}
\keyword{methods}
