\name{dbartsSpec}
\alias{dbartsSpec}
\title{Discrete Bayesian Additive Regression Trees Sampler Specification}
\description{
  Resolves the \code{control}/\code{model}/\code{data} triple and family token that describe a sampler, without constructing one. Intended for packages that embed dbarts - in particular \code{LinkingTo: dbarts} packages that create and hold a sampler through the C API in \file{dbarts.h} - and that supply their own design matrix.
}
\usage{
dbartsSpec(
    data, control = dbarts::dbartsControl(),
    tree.prior = cgm, node.prior = normal,
    resid.prior = chisq, resid.dist = gaussian,
    proposal.probs = c(
        birth_death = 0.5, swap = 0.1, change = 0.4, birth = 0.5),
    monotone = NULL, interactions = NULL, blocks = NULL,
    variance = NULL, n.trees.variance = 40L,
    power.variance = NULL, base.variance = NULL,
    sigma = NA_real_, seed = NA_integer_,
    family = c("auto", "gaussian", "probit", "logistic", "aft",
               "ordinal", "nbinom"),
    dispersion = NA_real_, survival = NULL,
    parentEnv = parent.frame())
}
\arguments{
  \item{data}{
    A \code{dbartsData} object, as built by \code{\link{dbartsData}}. Unlike \code{\link{dbarts}}, this function does no response ingestion of its own: the response must already be materialized.
  }
  \item{control, tree.prior, node.prior, resid.prior, resid.dist, proposal.probs, monotone, interactions, blocks, variance, n.trees.variance, power.variance, base.variance, sigma, seed, family, dispersion}{
    As in \code{\link{dbarts}}. The prior arguments are evaluated in dbarts's own prior vocabulary, so bare expressions such as \code{normal(k = chi(1.25, Inf))} resolve regardless of what the caller has attached.
  }
  \item{survival}{
    For \code{family = "aft"} only: a vector of per-observation status codes, 1 for an observed event and 0 for a right-censored observation, with the response holding log event/censoring times. \code{\link{dbarts}} reads this off a \code{Surv} or two-column response; a caller supplying log-times directly supplies it here. \code{NULL} for every other family.
  }
  \item{parentEnv}{
    The environment in which the prior arguments are evaluated; defaults to the caller's frame, which is almost always what is wanted. Pass explicitly when forwarding arguments from another function.
  }
}
\details{
  \code{\link{dbarts}} does two separable things: it ingests a formula or a pair of matrices into a \code{dbartsData}, and it then resolves that data into a full sampler specification - choosing the response family, coding the response, applying the family's weight policy, estimating a starting \code{sigma}, parsing the priors, and attaching the attributes that carry monotonicity, interaction and block constraints, a variance forest, residual degrees of freedom, and survival status. \code{\link{dbartsData}} exposes the first half. This function exposes the second, and returns the pieces rather than a sampler.

  Both entry points call the same resolution, so a family can never resolve two ways: a sampler built from this function's output draws the same samples as the equivalent \code{\link{dbarts}} call.

  \subsection{Use with the C API}{
    The \code{family} element of the returned list is the token to pass as the fourth argument of \code{dbarts_sampler_create}. Passing the empty string there instead asks the library to dispatch on the shape of the response, which is correct for \code{"gaussian"}, \code{"probit"}, \code{"ordinal"}, and \code{"nbinom"} - each is inferable from the response coding or from an attribute the library reads - but is \emph{wrong} for \code{"aft"} and \code{"logistic"}, which are indistinguishable by shape from \code{"gaussian"} and \code{"probit"} respectively. Pass the returned token and the question does not arise.
  }
  \subsection{Differences from dbarts()}{
    \code{sigma} defaults to leaving the data object's own value alone, rather than overwriting it: a caller that computed a starting estimate keeps it, and an unset (\code{NA}) value is still estimated during resolution. Cut points are taken from \code{control} only when the data does not already carry resolved per-column counts.

    Families that require ingestion this function does not perform are unavailable: \code{"hazard"} needs person-period expansion, and \code{"multinomial"} and \code{"hurdle.lognormal"} describe more than one sampler.
  }
}
\value{
  A list with components
  \item{control}{a \code{dbartsControl}, carrying any family and variance-forest attributes.}
  \item{model}{a \code{dbartsModel}, carrying the parsed priors and any monotone, interaction, block, and residual-degrees-of-freedom attributes.}
  \item{data}{a \code{dbartsData}, with the response coded and \code{sigma} resolved.}
  \item{family}{the resolved family token, never \code{"auto"}.}
}
\seealso{
  \code{\link{dbarts}}, \code{\link{dbartsData}}, \code{\link{dbartsControl}}, \code{\link{dbartsPriors}}
}
\examples{
f <- function(x) {
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
      10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

data <- dbartsData(x, y)
control <- dbartsControl(n.chains = 1L, n.samples = 20L, n.burn = 10L,
                         n.trees = 25L, updateState = FALSE)

## an ordinary specification, then the sampler it describes
spec <- dbartsSpec(data, control = control)
spec$family

sampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
samples <- sampler$run()

## additive-only fit, monotone in the fourth predictor
constrained <- dbartsSpec(data, control = control,
                          interactions = interactions(max.order = 1L),
                          monotone = c(rep(0, 3), 1, rep(0, 6)))
attr(constrained$model, "interaction.max.order")
}
