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
    variance = NULL,
    forests = NULL,
    sigma = NA_real_, seed = NA_integer_,
    family = c("auto", "gaussian", "probit", "logistic", "aft",
               "multinomial", "ordinal", "nbinom"),
    dispersion = NA_real_, survival = NULL,
    parentEnv = parent.frame())
}
\arguments{
  \item{data}{
    A \code{dbartsData} object, as built by \code{\link{dbartsData}}. Unlike \code{\link{dbarts}}, this function does no response ingestion of its own: the response must already be materialized. A data object carrying an \eqn{n \times K}{n x K} \code{counts} matrix resolves to \code{family = "multinomial"} even from \code{"auto"} - the slot IS the declaration - and is refused under any other family.
  }
  \item{control, tree.prior, node.prior, resid.prior, resid.dist, proposal.probs, monotone, interactions, blocks, variance, forests, sigma, seed, family, dispersion}{
    As in \code{\link{dbarts}}. The prior arguments are evaluated in dbarts's own prior vocabulary, so bare expressions such as \code{normal(k = chi(1.25, Inf))} resolve regardless of what the caller has attached. A \code{\link{forest}} \code{basis} given as a one-sided formula is evaluated in \code{parentEnv}, this surface performing no data ingestion of its own; the data object's rows are already whatever its own \code{subset} kept. A \code{basis} declared here always reaches the sampler: unlike \code{\link{dbarts}}, which refuses a basis declaration once \code{formula} is already a \code{dbartsData}, this function's first argument is always one, so the declaration is installed and REPLACES whatever bases the data object carried. \code{family}'s eight-token vocabulary is narrower than \code{\link{dbarts}}'s by design - multinomial's K-forest engine and \code{hurdle.lognormal}'s two-sampler composition build through \code{\link{bart2}}, never through this construct-only bypass; the wider set lives there.
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

  Both entry points call the same resolution on the \code{dbartsData} they end up with, so a sampler built from this function's output draws the same samples as the equivalent \code{\link{dbarts}} call given the SAME data object. The two can still resolve a given family differently, because \code{\link{dbarts}} performs its own response ingestion first - coding a factor response into a count matrix, for instance - while this function does none: a call this function refuses for the raw response \code{\link{dbarts}} can still accept, by ingesting it into the shape this function requires.

  \subsection{Use with the C API}{
    The \code{family} element of the returned list is the token to pass as the fourth argument of \code{dbarts_sampler_create}. Passing the empty string there instead asks the library to dispatch on the shape of the response, which is correct for \code{"gaussian"}, \code{"probit"}, \code{"ordinal"}, and \code{"nbinom"} - each is inferable from the response coding or from an attribute the library reads - but is \emph{wrong} for \code{"aft"} and \code{"logistic"}, which are indistinguishable by shape from \code{"gaussian"} and \code{"probit"} respectively. Pass the returned token and the question does not arise.
  }
  \subsection{Differences from dbarts()}{
    \code{sigma} defaults to leaving the data object's own value alone, rather than overwriting it: a caller that computed a starting estimate keeps it, and an unset (\code{NA}) value is still estimated during resolution. Cut points are taken from \code{control} only when the data does not already carry resolved per-column counts.

    Families that require ingestion this function does not perform are unavailable: \code{"hazard"} needs person-period expansion, and \code{"hurdle.lognormal"} describes more than one sampler. \code{"multinomial"} is in this function's own vocabulary, but only for a \code{data} object whose response is already an \eqn{n \times K}{n x K} \code{counts} matrix (see \code{data} above, and \code{\link{dbartsData}}'s \code{counts} argument); given a factor-coded or plain numeric response instead, it is refused with a message naming the counts route, even though \code{\link{dbarts}} accepts the same factor by coding it into that shape first.
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
