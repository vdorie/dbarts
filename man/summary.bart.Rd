\name{summary.bart}
\alias{summary.bart}
\alias{summary.rbart}
\alias{print.summary.bart}
\alias{as_draws_array.bart}
\alias{as_draws_array.rbart}
\alias{as_draws_df.bart}
\alias{as_draws_df.rbart}
\title{Convergence Diagnostics and Posterior-Package Draws for BART Fits}
\description{
  \code{summary} reports a per-variable posterior summary of the scalar
  parameters (\code{sigma}, \code{k}, and, for \code{\link{rbart_vi}}
  fits, \code{tau}) of a \code{bart}/\code{bart2}/\code{rbart_vi} fit,
  along with split-\eqn{\hat{R}}{Rhat} and effective sample size when
  the \pkg{posterior} package is installed.

  \code{as_draws_array} and \code{as_draws_df} convert a fit's
  chain-dimensioned draws to \pkg{posterior}'s array/data frame
  conventions. Both are methods for \pkg{posterior}'s generics and so
  require it to be loaded (\code{library(posterior)} or a qualified
  \code{posterior::as_draws_array} call); dbarts registers them only
  when \pkg{posterior} is available, as it is a \code{Suggests}, not a
  hard dependency.
}
\usage{
\method{summary}{bart}(object, vars = c("sigma", "k", "tau"), \dots)
\method{summary}{rbart}(object, vars = c("sigma", "k", "tau"), \dots)
\method{print}{summary.bart}(x, \dots)

\method{as_draws_array}{bart}(x, vars = c("sigma", "k", "tau"), \dots)
\method{as_draws_array}{rbart}(x, vars = c("sigma", "k", "tau"), \dots)
\method{as_draws_df}{bart}(x, vars = c("sigma", "k", "tau"), \dots)
\method{as_draws_df}{rbart}(x, vars = c("sigma", "k", "tau"), \dots)
}
\arguments{
  \item{object, x}{
    An object of class \code{bart} or \code{rbart}, as returned by
    \code{\link{bart}}, \code{\link{bart2}}, or \code{\link{rbart_vi}}.
  }
  \item{vars}{
    Character vector of fields to gather. Requested fields absent from
    \code{object} (e.g. \code{k} when unmodeled, \code{tau} for a plain
    \code{bart}/\code{bart2} fit) are silently dropped. \code{sigma},
    \code{k}, and \code{tau} contribute one draws variable each; any
    other field (\code{varcount}, \code{varprobs}, \code{yhat.train},
    \code{yhat.test}, \code{ranef}) contributes one variable per
    column, named \code{"field[column]"}. \eqn{f(x)} draws
    (\code{yhat.*}) are reachable this way but are not summarized
    automatically, as they carry one variable per observation.
  }
  \item{\dots}{Unused.}
}
\value{
  \code{summary} returns an object of class \code{summary.bart} with
  elements \code{call}, \code{stats} (a data frame with one row per
  requested scalar, or \code{NULL} if none of \code{vars} are present
  on the fit; \code{rhat}, \code{ess_bulk}, and \code{ess_tail} columns
  are added when \pkg{posterior} is installed), and \code{posterior}
  (whether those columns were computed). Its print method notes when
  \pkg{posterior} is missing, and separately when any R-hat exceeds
  1.01; neither withholds the rest of the summary or errors.

  \code{as_draws_array} and \code{as_draws_df} return a \pkg{posterior}
  \code{draws_array}/\code{draws_df}: the fit's native (chain, sample[,
  variable]) storage, transposed to \pkg{posterior}'s (iteration,
  chain, variable) convention.
}
\seealso{
\code{\link{bart}}, \code{\link{rbart_vi}}
}
\examples{
\donttest{
fit <- bart2(y ~ x, data.frame(y = rnorm(100), x = rnorm(100)), n.chains = 2L,
             n.samples = 20L, n.burn = 20L, n.trees = 5L, n.threads = 1L,
             verbose = FALSE)
summary(fit)

if (requireNamespace("posterior", quietly = TRUE)) {
  library(posterior)
  as_draws_array(fit)
}
}
}
\keyword{htest}
