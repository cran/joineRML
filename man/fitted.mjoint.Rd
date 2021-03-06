% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted.mjoint.R
\name{fitted.mjoint}
\alias{fitted.mjoint}
\title{Extract \code{mjoint} fitted values}
\usage{
\method{fitted}{mjoint}(object, level = 0, ...)
}
\arguments{
\item{object}{an object inheriting from class \code{mjoint} for a joint model
of time-to-event and multivariate longitudinal data.}

\item{level}{an optional integer giving the level of grouping to be used in
extracting the fitted values from object. Level values increase from outermost
to innermost grouping, with level 0 corresponding to the population
fitted values and level 1 corresponding to subject-specific fitted values Defaults
to level 0.}

\item{...}{additional arguments; currently none are used.}
}
\value{
A \code{list} of length \emph{K} with each element a vector of
  fitted values for the \emph{k}-th longitudinal outcome.
}
\description{
The fitted values at level \emph{i} are obtained by adding
  together the population fitted values (based only on the fixed effects
  estimates) and the estimated contributions of the random effects.
}
\references{
Pinheiro JC, Bates DM. \emph{Mixed-Effects Models in S and S-PLUS.} New York:
Springer Verlag; 2000.
}
\seealso{
\code{\link{mjoint}}, \code{\link{residuals.mjoint}}
}
\author{
Graeme L. Hickey (\email{graemeleehickey@gmail.com})
}
\keyword{methods}
