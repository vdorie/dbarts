\name{dbarts-deprecated}
\alias{dbarts-deprecated}
\alias{bart2}
\alias{rbart_vi}
\alias{predict.rbart}
\alias{extract.rbart}
\alias{fitted.rbart}
\alias{residuals.rbart}
\title{Deprecated and Removed Functions in dbarts}
\description{
  Names kept reachable for one release past their removal or rename, each pointing at its successor. They are removed in dbarts 1.1-0.

  \code{bart2} is the former name of the modern front door, which is now \code{\link{bart}}; a \code{bart2} call is forwarded unchanged, with a warning the first time in a session.

  \code{rbart_vi} and its methods fit a grouped random-intercept BART model. That model now lives in the \pkg{stan4bart} package (\code{stan4bart::stan4bart}), whose prior on the group spread differs from the one dbarts drew from, so a refit moves the results rather than reproducing them. These stubs error.
}
\usage{
rbart_vi(\dots)

\method{predict}{rbart}(object, \dots)
\method{extract}{rbart}(object, \dots)
\method{fitted}{rbart}(object, \dots)
\method{residuals}{rbart}(object, \dots)
}
\arguments{
  \item{object}{
    An \code{rbart} fit object saved by dbarts 0.9-x.
  }
  \item{\dots}{
    Ignored; the stubs do not run. \code{bart2}'s own arguments are \code{\link{bart}}'s, copied rather than restated - \code{formals(bart2)} is \code{formals(bart)} - so a caller reading the defaults programmatically still finds them.
  }
}
\value{
  \code{bart2} returns whatever \code{\link{bart}} returns. The \code{rbart_vi} stubs do not return.
}
\seealso{
  \code{\link{bart}}, \code{\link{bartBT}}
}
\keyword{internal}
