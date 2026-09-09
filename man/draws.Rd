\name{draws}
\alias{draws}
\alias{draws.bart}
\alias{draws.bartMultinomial}
\alias{draws.bartOrdinal}
\alias{draws.bartNegbin}
\alias{draws.bartHurdle}
\title{Chain-Dimensioned Draws From a BART Fit}
\description{
  Converts a fit's chain-dimensioned draws to a plain array with
  dimensions \code{(iteration, chain, variable)} and dimnames on the
  variable margin - the shape \code{posterior::as_draws_array} accepts
  unchanged from a caller who has that package installed
  (\code{as_draws_array(draws(fit))}). \code{\link{summary.bart}} is
  built on the same array.

  \code{\link{bart}}'s four own-class fits (\code{"bartMultinomial"},
  \code{"bartOrdinal"}, \code{"bartNegbin"}, \code{"bartHurdle"}; see
  \code{\link{bart}}) have their own \code{draws} methods, each
  exposing the scalar posterior parameters that family carries rather
  than its per-observation channels - never \code{yhat.train} itself.
}
\usage{
\method{draws}{bart}(x, vars = c("sigma", "k", "tau"), \dots)

\method{draws}{bartMultinomial}(x, vars = "meanProb", \dots)
\method{draws}{bartOrdinal}(x, vars = c("thresholds", "sigma", "k", "tau"), \dots)
\method{draws}{bartNegbin}(x, vars = c("dispersion", "sigma", "k", "tau"), \dots)
\method{draws}{bartHurdle}(x, vars = c("sigma", "k", "tau"), \dots)
}
\arguments{
  \item{x}{
    An object of class \code{bart}, as returned by \code{\link{bartBT}}
    or \code{\link{bart}}; for the four own-class methods, a \code{bart}
    fit of the matching class.
  }
  \item{vars}{
    Character vector of fields to gather. Requested fields absent from
    \code{x} (e.g. \code{k} when unmodeled, or \code{tau}, which no
    shipped family carries) are silently dropped. \code{sigma}, \code{k},
    and \code{tau} contribute one draws variable each; any other field
    (\code{varcount}, \code{varprobs}, \code{yhat.train}, \code{yhat.test})
    contributes one variable per column, named \code{"field[column]"}.
    \eqn{f(x)} draws (\code{yhat.*}) are reachable this way but are not
    summarized automatically by \code{\link{summary.bart}}, as they carry
    one variable per observation.

    On a heteroscedastic fit the \code{sigma} token resolves to
    \code{mean.s}: the mean over the training observations of that draw's
    variance surface \eqn{s(x)}, one value per draw. The \code{sigma} such
    a fit stores is the variance forest parameterization's fixed unit
    residual times the range of the response, a constant with no
    posterior content, so it is not reported as a parameter.
    \code{s.train} itself is reachable by name, contributing one variable
    per observation.

    For the four own-class methods, \code{vars} is scoped to that
    family's own vocabulary (see \code{\link{bart}}): a
    \code{"bartOrdinal"} fit's \code{"thresholds"} contributes
    \code{threshold[1]} (pinned at 0) through \code{threshold[K - 1]}; a
    \code{"bartNegbin"} fit's \code{"dispersion"} contributes the
    per-draw dispersion \eqn{r}; a \code{"bartMultinomial"} fit has a
    single channel, so its \code{vars} only ever means \code{"meanProb"}
    and it always reports \code{meanProb[<level>]}, its only scalar
    posterior parameter; a \code{"bartHurdle"} fit applies \code{vars} to
    both components and labels the result
    \code{occupancy.<field>}/\code{positive.<field>} (a dot, not a
    bracket, since \pkg{posterior} parses a bracket as an index).
  }
  \item{\dots}{Unused.}
}
\value{
  A numeric array with dimensions \code{(iteration, chain, variable)} and
  dimnames on the variable margin, matching the fit's native (chain,
  sample[, variable]) storage transposed to that convention. The four
  own-class methods return the same convention over their own
  family-scoped \code{vars}.
}
\seealso{
\code{\link{bart}}, \code{\link{bartBT}}, \code{\link{summary.bart}}
}
\examples{
\donttest{
fit <- bart(y ~ x, data.frame(y = rnorm(100), x = rnorm(100)), n.chains = 2L,
             n.samples = 20L, n.burn = 20L, n.trees = 5L, n.threads = 1L,
             verbose = FALSE)
d <- draws(fit)
dim(d)
}
}
\keyword{htest}
