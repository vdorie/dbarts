\name{survivalProbabilities}
\title{Survival Probability Draws from an AFT Fit}
\alias{survivalProbabilities}
\alias{survivalProbabilities.bart}
\alias{survivalProbabilities.rbart}
\description{
  Posterior draws of the survival probability \eqn{S(t \mid x)} from an
  accelerated failure time (AFT) log-normal fit produced by
  \code{\link{bart2}} (or \code{\link{bart}}) with \code{family = "aft"}.
}
\usage{
survivalProbabilities(object, ...)

\method{survivalProbabilities}{bart}(
  object,
  times,
  newdata = NULL,
  combineChains = TRUE,
  ...
)

\method{survivalProbabilities}{rbart}(object, ...)
}
\arguments{
  \item{object}{
    A fitted \code{bart} object from an AFT model (its \code{family} element
    equals \code{"aft"}).
  }
  \item{times}{
    A numeric vector of times (on the original, non-logarithmic scale) at
    which to evaluate the survival function.
  }
  \item{newdata}{
    Optional predictors at which to evaluate. When \code{NULL} the training
    observations are used; otherwise \code{object} must have been fit with
    \code{keepTrees = TRUE}.
  }
  \item{combineChains}{
    A logical determining whether the chain dimension is collapsed into the
    draw dimension, as elsewhere in the package.
  }
  \item{...}{
    Not used; for compatibility with the generic.
  }
}
\details{
  Under the log-normal AFT model \eqn{\log T = f(x) + \sigma \epsilon} with
  \eqn{\epsilon \sim N(0, 1)}, the survival function is
  \deqn{S(t \mid x) = 1 - \Phi\left(\frac{\log t - f(x)}{\sigma}\right),}
  where \eqn{f(x)} is the linear predictor \eqn{E[\log T \mid x]} - the
  log-time-scale quantity that \code{\link[=predict.bart]{predict}} and
  \code{\link[=extract.bart]{extract}} return for an \code{"aft"} fit - and
  \eqn{\sigma} the residual standard deviation. The probability is evaluated
  at every posterior draw of \eqn{f(x)} and \eqn{\sigma}, following the
  package's convention that draw-level functions return draws: take means or
  quantiles over the draw margin for point estimates and credible bands.

  There is no method for \code{rbart} fits: \code{\link{rbart_vi}} has no
  \code{family} argument, so a grouped (random-intercept) AFT model is not
  yet fittable from R, and a survival curve that ignored the fitted
  intercepts would be wrong. The \code{rbart} method fails with an error
  saying so.
}
\value{
  An array of survival probabilities. With \code{combineChains = TRUE}
  (or a single chain), the dimensions are draws \eqn{\times} times
  \eqn{\times} observations; with \code{combineChains = FALSE} on a
  multi-chain fit, a chain dimension precedes them (chains \eqn{\times}
  draws \eqn{\times} times \eqn{\times} observations). Observations are
  those of \code{newdata}, or the training data when \code{newdata} is
  \code{NULL}.
}
\seealso{
  \code{\link{bart2}}, \code{\link[=predict.bart]{predict}},
  \code{\link[=extract.bart]{extract}}
}
\examples{
set.seed(42)
n <- 100L
x <- matrix(runif(n * 2L), n, 2L)
eventTime <- exp(x[, 1L] + 0.5 * rnorm(n))
censorTime <- exp(x[, 1L] + 0.25 + 0.5 * rnorm(n))
status <- as.numeric(eventTime <= censorTime)
observedTime <- pmin(eventTime, censorTime)

fit <- bart2(
  x,
  cbind(observedTime, status),
  family = "aft",
  n.trees = 25L,
  n.burn = 50L,
  n.samples = 100L,
  n.chains = 1L,
  verbose = FALSE
)

surv <- survivalProbabilities(fit, times = c(0.5, 1, 2))
dim(surv) # draws x times x observations

## posterior-mean survival curve and 90\% band for the first observation
apply(surv[, , 1L], 2L, mean)
apply(surv[, , 1L], 2L, quantile, probs = c(0.05, 0.95))
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\keyword{nonparametric}
\keyword{survival}
\keyword{regression}
