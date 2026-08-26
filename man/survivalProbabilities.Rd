\name{survivalProbabilities}
\title{Survival Probability Draws from a Survival Fit}
\alias{survivalProbabilities}
\alias{survivalProbabilities.bart}
\alias{survivalProbabilities.rbart}
\alias{survivalProbabilities.bartMultinomial}
\alias{survivalProbabilities.bartOrdinal}
\alias{survivalProbabilities.bartNegbin}
\alias{survivalProbabilities.bartHurdle}
\description{
  Posterior draws of the survival probability \eqn{S(t \mid x)} from a survival
  fit: an accelerated failure time (AFT) log-normal fit produced by
  \code{\link{bart2}} (or \code{\link{bart}}) with \code{family = "aft"}, the
  grouped (random-intercept) AFT fit produced by \code{\link{rbart_vi}} with
  \code{family = "aft"}, or a discrete-time hazard fit produced by
  \code{\link{bart2}} with \code{family = "hazard"} (or
  \code{"hazard.logistic"}).

  \code{\link{bart2}}'s four own-class fits (\code{"bartMultinomial"},
  \code{"bartOrdinal"}, \code{"bartNegbin"}, \code{"bartHurdle"}) refuse
  by name instead: none of the four is a survival model, so none carries
  a hazard channel.
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

\method{survivalProbabilities}{rbart}(
  object,
  times,
  newdata = NULL,
  combineChains = TRUE,
  ...,
  group.by
)

\method{survivalProbabilities}{bartMultinomial}(object, ...)

\method{survivalProbabilities}{bartOrdinal}(object, ...)

\method{survivalProbabilities}{bartNegbin}(object, ...)

\method{survivalProbabilities}{bartHurdle}(object, ...)
}
\arguments{
  \item{object}{
    A fitted \code{bart} or \code{rbart} object from a survival model: an AFT
    model (its \code{family} element equals \code{"aft"}) or a discrete-time
    hazard model (it carries a \code{$periods} grid; its \code{family} records
    the binary link). A hazard fit must have been made with
    \code{keepTrees = TRUE}.
  }
  \item{times}{
    A numeric vector of times (on the original, non-logarithmic scale) at
    which to evaluate the survival function. For a discrete-time hazard fit the
    times are horizons on the period grid, and the default (when \code{times}
    is missing) is the training grid, \code{object$periods}.
  }
  \item{newdata}{
    Optional predictors at which to evaluate. When \code{NULL} the training
    observations are used; otherwise \code{object} must have been fit with
    \code{keepTrees = TRUE}.
  }
  \item{group.by}{
    For the \code{rbart} method with \code{newdata}, the grouping factor for
    the new observations, as in \code{\link[=predict.rbart]{predict}}. A group
    not seen in training draws its intercept from \eqn{N(0, \tau)}. Ignored
    (and unnecessary) when \code{newdata} is \code{NULL}. Supplied by name
    only - it follows \code{\dots} in the signature, so it is never matched
    positionally; a missing one is refused, naming itself. Refused by name on
    the \code{bart} method, which has no grouping.
  }
  \item{combineChains}{
    A logical determining whether the chain dimension is collapsed into the
    draw dimension, as elsewhere in the package.
  }
  \item{...}{
    Not a formal on either method: a name belonging to a sibling generic
    (\code{type}, \code{sample}, \code{ci.level}, \code{offset},
    \code{weights}, \code{n.threads}, \code{forest}, \code{contribution},
    \code{bases}) is refused by name rather than silently ignored, since this
    method returns the draws of \eqn{S(t \mid x)} at \code{times} alone.
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

  For a grouped fit from \code{\link{rbart_vi}} with \code{family = "aft"},
  the \code{rbart} method uses the same formula, with \eqn{f(x)} replaced by
  the linear predictor \eqn{E[\log T \mid x, g] = f(x) + \alpha_{g}} that
  includes the drawn random intercept for the observation's group. It is
  sourced from the expected-value (\code{"ev"}) channel, so dropping the
  intercepts (which would misplace every grouped curve) is not possible.

  For a discrete-time hazard fit, the survival function is the cumulative
  product
  \deqn{S(t \mid x) = \prod_{k \,:\, \mathrm{periods}[k] \le t} (1 - h(k \mid x)),}
  where \eqn{h(k \mid x) = g(f(x, k) + o)} is the per-period hazard through the
  fit's binary link \eqn{g}. Because the training design is ragged (each
  subject carries only its at-risk rows), the method ALWAYS re-expands its
  subjects onto the full grid and replays the trees - so it requires
  \code{keepTrees = TRUE} even for the training data (\code{newdata = NULL}).
  With \code{newdata}, each new subject is expanded to one row per period and
  its curve evaluated the same way.
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
  \code{\link{bart2}}, \code{\link{rbart_vi}},
  \code{\link[=predict.bart]{predict}}, \code{\link[=extract.bart]{extract}}
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
