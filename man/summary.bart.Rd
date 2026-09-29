\name{summary.bart}
\alias{summary.bart}
\alias{print.summary.bart}
\title{Convergence Diagnostics for BART Fits}
\description{
  Reports a per-variable posterior summary of the scalar parameters
  (\code{sigma} and \code{k}) of a \code{bartBT}/\code{bart} fit, along
  with split-\eqn{\hat{R}}{Rhat} and bulk/tail effective sample size,
  computed by dbarts itself (no \pkg{posterior} dependency). A
  heteroscedastic fit (\code{\link{dbarts}}'s \code{variance}) has no
  scalar residual scale and reports \code{mean.s} in place of
  \code{sigma}; see \code{vars}.

  \code{bart}'s four own-class fits (\code{"bartMultinomial"},
  \code{"bartOrdinal"}, \code{"bartNegbin"}, \code{"bartHurdle"}; see
  \code{\link{bart}}) summarize through this same method, each exposing
  the scalar posterior parameters that family carries rather than its
  per-observation channels - never \code{yhat.train} itself. A chain-
  separated (rather than pooled) view of these same draws is
  \code{\link{extract}}, whose \code{"sigma"}, \code{"k"}, \code{"sd"}, and other
  scalar types return a vector or, with \code{combineChains = FALSE}, a
  chains-by-samples matrix.
}
\usage{
\method{summary}{bart}(object, vars = c("sigma", "k", "sd"), \dots)
\method{print}{summary.bart}(x, \dots)
}
\arguments{
  \item{object, x}{
    An object of class \code{bart}, as returned by
    \code{\link{bartBT}} or \code{\link{bart}}; for the four own-class
    fits, a \code{bart} fit of the matching class.
  }
  \item{vars}{
    Character vector of fields to summarize. Requested fields absent
    from \code{object} (e.g. \code{k} when unmodeled, or on a fit whose
    leaf prior was named by \code{sd}, which carries \code{sd} draws
    instead) are silently dropped. \code{sigma}, \code{k} and \code{sd}
    contribute one summarized variable each; any other field (\code{varcount}, \code{varprobs},
    \code{yhat.train}, \code{yhat.test}) contributes one variable per
    column, named \code{"field[column]"}.

    On a heteroscedastic fit the \code{sigma} token resolves to
    \code{mean.s}: the mean over the training observations of that
    draw's variance surface \eqn{s(x)}, one value per draw. The
    \code{sigma} such a fit stores is the variance forest
    parameterization's fixed unit residual times the range of the
    response, a constant with no posterior content, so it is not
    reported as a parameter; \code{s.train} itself is reachable by
    name, contributing one variable per observation.

    For the four own-class fits, \code{vars} is scoped to that
    family's own vocabulary (see \code{\link{bart}}): a
    \code{"bartOrdinal"} fit's \code{"thresholds"} contributes
    \code{threshold[1]} (pinned at 0) through \code{threshold[K - 1]};
    a \code{"bartNegbin"} fit's \code{"dispersion"} contributes the
    per-draw dispersion \eqn{r}; a \code{"bartMultinomial"} fit has no
    \code{vars} argument at all - its only scalar posterior parameter
    is the per-category mean predicted probability, reported as
    \code{prob[<level>]}; a \code{"bartHurdle"} fit applies \code{vars}
    to both components separately.
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
\code{\link{bart}}, \code{\link{bartBT}}, \code{\link{extract}}
}
\examples{
\donttest{
fit <- bart(y ~ x, data.frame(y = rnorm(100), x = rnorm(100)), n.chains = 2L,
             n.samples = 20L, n.burn = 20L, n.trees = 5L, n.threads = 1L,
             verbose = FALSE)
summary(fit)
}
}
\keyword{htest}
