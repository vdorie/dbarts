\name{summary.bart}
\alias{summary.bart}
\alias{summary.rbart}
\alias{print.summary.bart}
\alias{as_draws_array.bart}
\alias{as_draws_array.rbart}
\alias{as_draws_df.bart}
\alias{as_draws_df.rbart}
\alias{as_draws_array.bartMultinomial}
\alias{as_draws_array.bartOrdinal}
\alias{as_draws_array.bartNegbin}
\alias{as_draws_array.bartHurdle}
\alias{as_draws_df.bartMultinomial}
\alias{as_draws_df.bartOrdinal}
\alias{as_draws_df.bartNegbin}
\alias{as_draws_df.bartHurdle}
\title{Convergence Diagnostics and Posterior-Package Draws for BART Fits}
\description{
  \code{summary} reports a per-variable posterior summary of the scalar
  parameters (\code{sigma}, \code{k}, and, for \code{\link{rbart_vi}}
  fits, \code{tau}) of a \code{bart}/\code{bart2}/\code{rbart_vi} fit,
  along with split-\eqn{\hat{R}}{Rhat} and effective sample size when
  the \pkg{posterior} package is installed. A heteroscedastic fit
  (\code{\link{dbarts}}'s \code{variance}) has no scalar residual scale
  and reports \code{mean.s} in place of \code{sigma}; see \code{vars}.

  \code{as_draws_array} and \code{as_draws_df} convert a fit's
  chain-dimensioned draws to \pkg{posterior}'s array/data frame
  conventions. Both are methods for \pkg{posterior}'s generics and so
  require it to be loaded (\code{library(posterior)} or a qualified
  \code{posterior::as_draws_array} call); dbarts registers them only
  when \pkg{posterior} is available, as it is a \code{Suggests}, not a
  hard dependency.

  \code{bart2}'s four own-class fits (\code{"bartMultinomial"},
  \code{"bartOrdinal"}, \code{"bartNegbin"}, \code{"bartHurdle"}; see
  \code{\link{bart2}}) are draws objects too, through their own
  \code{summary}/\code{as_draws_array}/\code{as_draws_df} methods
  (documented under \code{\link{bart2}}'s \sQuote{Value}), each exposing
  the scalar posterior parameters that family carries rather than its
  per-observation channels - never \code{yhat.train} itself.
}
\usage{
\method{summary}{bart}(object, vars = c("sigma", "k", "tau"), \dots)
\method{summary}{rbart}(object, vars = c("sigma", "k", "tau"), \dots)
\method{print}{summary.bart}(x, \dots)

\method{as_draws_array}{bart}(x, vars = c("sigma", "k", "tau"), \dots)
\method{as_draws_array}{rbart}(x, vars = c("sigma", "k", "tau"), \dots)
\method{as_draws_df}{bart}(x, vars = c("sigma", "k", "tau"), \dots)
\method{as_draws_df}{rbart}(x, vars = c("sigma", "k", "tau"), \dots)

\method{as_draws_array}{bartMultinomial}(x, vars = "meanProb", \dots)
\method{as_draws_df}{bartMultinomial}(x, vars = "meanProb", \dots)
\method{as_draws_array}{bartOrdinal}(
  x, vars = c("cutpoints", "sigma", "k", "tau"), \dots)
\method{as_draws_df}{bartOrdinal}(
  x, vars = c("cutpoints", "sigma", "k", "tau"), \dots)
\method{as_draws_array}{bartNegbin}(
  x, vars = c("dispersion", "sigma", "k", "tau"), \dots)
\method{as_draws_df}{bartNegbin}(
  x, vars = c("dispersion", "sigma", "k", "tau"), \dots)
\method{as_draws_array}{bartHurdle}(x, vars = c("sigma", "k", "tau"), \dots)
\method{as_draws_df}{bartHurdle}(x, vars = c("sigma", "k", "tau"), \dots)
}
\arguments{
  \item{object, x}{
    An object of class \code{bart} or \code{rbart}, as returned by
    \code{\link{bart}}, \code{\link{bart2}}, or \code{\link{rbart_vi}};
    for the four own-class methods, a \code{bart2} fit of the matching
    class.
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

    On a heteroscedastic fit the \code{sigma} token resolves to
    \code{mean.s}: the mean over the training observations of that
    draw's variance surface \eqn{s(x)}, one value per draw. The
    \code{sigma} such a fit stores is the variance forest
    parameterization's fixed unit residual times the range of the
    response, a constant with no posterior content, so it is not
    reported as a parameter. \code{s.train} itself is reachable by
    name, contributing one variable per observation.

    For the four own-class methods, \code{vars} is scoped to that
    family's own vocabulary (see \code{\link{bart2}}): a
    \code{"bartOrdinal"} fit's \code{"cutpoints"} contributes
    \code{cutpoint[1]} (pinned at 0) through \code{cutpoint[K - 1]}; a
    \code{"bartNegbin"} fit's \code{"dispersion"} contributes the
    per-draw dispersion \eqn{r}; a \code{"bartMultinomial"} fit has a
    single channel, so its \code{vars} only ever means \code{"meanProb"}
    and it always reports \code{meanProb[<level>]}, its only scalar
    posterior parameter (\code{summary} on such a fit carries no
    \code{vars} at all and refuses one by name); a \code{"bartHurdle"} fit
    applies \code{vars} to both components and labels the result
    \code{occupancy.<field>}/\code{positive.<field>} (a dot, not a
    bracket, since \pkg{posterior} parses a bracket as an index).
  }
  \item{\dots}{Unused.}
}
\value{
  \code{summary} returns an object of class \code{summary.bart} with
  elements \code{call}, \code{stats} (a data frame with one row per
  requested scalar, or \code{NULL} if none of \code{vars} are present
  on the fit), and \code{posterior} (whether \pkg{posterior} produced
  \code{stats}). When \pkg{posterior} is installed, \code{stats} is
  produced by its \code{summarise_draws} function, whose columns
  include \code{mean}, \code{median}, \code{sd}, \code{mad}, \code{q5},
  \code{q95}, \code{rhat}, \code{ess_bulk}, and \code{ess_tail};
  otherwise a built-in summary supplies \code{mean}, \code{sd},
  \code{q2.5}, \code{median}, and \code{q97.5}, with no convergence
  diagnostics. Its print method notes when \pkg{posterior} is missing,
  and separately when any R-hat exceeds 1.01; neither withholds the
  rest of the summary or errors.

  \code{as_draws_array} and \code{as_draws_df} return a \pkg{posterior}
  \code{draws_array}/\code{draws_df}: the fit's native (chain, sample[,
  variable]) storage, transposed to \pkg{posterior}'s (iteration,
  chain, variable) convention. The four own-class methods return the
  same convention over their own family-scoped \code{vars}.
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
