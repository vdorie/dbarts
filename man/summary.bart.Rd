\name{summary.bart}
\alias{summary.bart}
\alias{print.summary.bart}
\title{Convergence Diagnostics for BART Fits}
\description{
  Reports a per-variable posterior summary of the scalar parameters
  (\code{sigma} and \code{k}) of a \code{bart}/\code{bart2} fit, along
  with split-\eqn{\hat{R}}{Rhat} and bulk/tail effective sample size,
  computed by dbarts itself (no \pkg{posterior} dependency). A
  heteroscedastic fit (\code{\link{dbarts}}'s \code{variance}) has no
  scalar residual scale and reports \code{mean.s} in place of
  \code{sigma}; see \code{vars}.

  \code{bart2}'s four own-class fits (\code{"bartMultinomial"},
  \code{"bartOrdinal"}, \code{"bartNegbin"}, \code{"bartHurdle"}; see
  \code{\link{bart2}}) summarize through this same method, each exposing
  the scalar posterior parameters that family carries rather than its
  per-observation channels - never \code{yhat.train} itself.

  \code{\link{draws}} returns the same fit's chain-dimensioned draws as
  a plain array, the shape this summary is computed from.
}
\usage{
\method{summary}{bart}(object, vars = c("sigma", "k", "tau"), \dots)
\method{print}{summary.bart}(x, \dots)
}
\arguments{
  \item{object, x}{
    An object of class \code{bart}, as returned by
    \code{\link{bart}} or \code{\link{bart2}}; for the four own-class
    fits, a \code{bart2} fit of the matching class.
  }
  \item{vars}{
    Character vector of fields to summarize. Requested fields absent
    from \code{object} (e.g. \code{k} when unmodeled, or \code{tau},
    which no shipped family carries) are silently dropped. \code{sigma},
    \code{k}, and \code{tau} contribute one summarized variable each; any
    other field (\code{varcount}, \code{varprobs}, \code{yhat.train},
    \code{yhat.test}) contributes one variable per column, named
    \code{"field[column]"}.

    On a heteroscedastic fit the \code{sigma} token resolves to
    \code{mean.s}: the mean over the training observations of that
    draw's variance surface \eqn{s(x)}, one value per draw; see
    \code{\link{draws}} for the full rule.

    A \code{"bartMultinomial"} fit has no \code{vars} argument: its only
    scalar posterior parameter is the per-category mean predicted
    probability, reported as \code{meanProb[<level>]}.
  }
  \item{\dots}{Unused.}
}
\value{
  An object of class \code{summary.bart} with elements \code{call},
  \code{stats} (a data frame with one row per requested scalar, or
  \code{NULL} if none of \code{vars} are present on the fit), and
  \code{vars}. \code{stats} has columns \code{mean}, \code{median},
  \code{sd}, \code{mad}, \code{q5}, \code{q95}, \code{rhat},
  \code{ess_bulk}, and \code{ess_tail} - the same columns
  \code{posterior::summarise_draws} reports for a plain array, computed
  without that package. Its print method notes when any R-hat exceeds
  1.01; this does not withhold the rest of the summary or error, as
  dbarts does not refuse to summarize a non-converged fit.
}
\seealso{
\code{\link{bart}}, \code{\link{bart2}}, \code{\link{draws}}
}
\examples{
\donttest{
fit <- bart2(y ~ x, data.frame(y = rnorm(100), x = rnorm(100)), n.chains = 2L,
             n.samples = 20L, n.burn = 20L, n.trees = 5L, n.threads = 1L,
             verbose = FALSE)
summary(fit)
}
}
\keyword{htest}
