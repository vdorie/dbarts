\name{samplePriorPredictive}
\alias{samplePriorPredictive}
\title{Draw From The BART Prior}
\description{
  Draws repeatedly from a \code{\linkS4class{dbartsSampler}}'s prior - tree structure, end-node
  parameters, and optionally response noise - and evaluates the resulting forest at \code{x.test}.
  Useful for calibrating priors before ever fitting to data, e.g. inspecting the implied prior on a
  treatment effect \eqn{f(x_1) - f(x_0)}{f(x1) - f(x0)}. The sampler passed in is never modified.
}
\usage{
samplePriorPredictive(
  sampler, x.test = NULL, n.samples = 200L,
  type = c("ev", "ppd"), offset.test = NULL,
  n.threads = sampler$control@n.threads
)
}
\arguments{
  \item{sampler}{
    A \code{\linkS4class{dbartsSampler}}. Left unmodified: the function operates on a freshly
    constructed private sampler, built from this sampler's control, model, and data, with its own
    engine state.
  }
  \item{x.test}{
    A matrix of predictors at which to evaluate the drawn forests, in the form accepted by
    \code{\link[=dbartsSampler-class]{predict}}. \code{NULL} (the default) evaluates at the sampler's
    own training predictors.
  }
  \item{n.samples}{
    A positive integer giving the number of independent prior draws to return.
  }
  \item{type}{
    \code{"ev"} draws \eqn{E[y \mid x]}{E[y | x]} under the prior: the forest sum on the response
    scale for a gaussian sampler, or the link-transformed probability for a probit/logistic one -
    matching what \code{\link{extract}} and \code{predict} call \code{"ev"}. \code{"ppd"} adds the
    observation model on top: gaussian noise with sigma freshly drawn from its own prior for each
    sample, or a bernoulli draw from the \code{"ev"} probabilities for a binary sampler; on an aft
    (survival) sampler the draw is on the log-time scale (the model fits \eqn{\log T}{log T}).
    Unlike \code{\link{extract}}'s posterior \code{"ppd"}, prior predictive draws are weight-blind:
    they use unit weights rather than a weighted precision (\eqn{\sigma / \sqrt{w}}{sigma / sqrt(w)}
    for gaussian) or a \eqn{\mathrm{Binomial}(w, p)}{Binomial(w, p)} success count (weighted
    logistic); supply any weighted precision downstream.
  }
  \item{offset.test}{
    A numeric offset added to the forest sum, as in \code{\link[=dbartsSampler-class]{predict}};
    recycled to the number of rows of \code{x.test} if of length 1, or \code{NULL} (the default) for
    none.
  }
  \item{n.threads}{
    Forwarded to the underlying \code{\link[=dbartsSampler-class]{predict}} calls.
  }
}
\details{
  Each of the \code{n.samples} draws independently redraws every tree's structure
  (\code{sampleTreesFromPrior}) and end-node parameters (\code{sampleNodeParametersFromPrior}) from
  scratch, then evaluates the resulting forest at \code{x.test} - there is no MCMC dependence between
  rows of the result, unlike a posterior sample from \code{run}.

  For \code{type = "ppd"} on a gaussian sampler, sigma is drawn once per sample from the model's
  residual-variance prior: a scaled-inverse-chi-squared calibrated so that
  \eqn{P(\sigma < \code{sigest}) = \code{quantile}}{P(sigma < sigest) = quantile} when the prior is
  \code{chisq} (see \code{\link{dbartsPriors}}), or the sampler's fixed value when the prior is
  \code{fixed}.

  If \code{sampler} has more than one chain, prior draws are chain-free: each of the \code{n.samples}
  draws only keeps the first chain's stream (the other chains still advance their own independent
  engine RNGs, but their draws are discarded), since there is no mixing or convergence concept for
  draws that are already independent by construction.

  Seeding follows from the private sampler being constructed fresh on every call. When the
  \code{\link[=dbartsControl]{control}}'s \code{seed} is \code{NA} (the default), the engine RNG
  is seeded from R's random number stream at construction, so successive calls return independent
  draws, and a \code{\link{set.seed}} beforehand makes the whole draw - forests and observation
  layer alike - reproducible. When the control fixes \code{seed}, the engine stream is pinned:
  every call replays the identical forest sequence, with only the R-side \code{"ppd"} observation
  layer still governed by \code{set.seed}.

  To compare the prior at two or more covariate settings with paired forests - as for a treatment
  effect - stack the settings as rows of a single \code{x.test} and difference the corresponding
  columns of the result; see \sQuote{Examples}. Two separate calls draw independent sets of forests
  (unless \code{seed} is fixed, as above) and so cannot be differenced meaningfully.
}
\value{
  A matrix with \code{n.samples} rows and one column per row of \code{x.test} (or of the training
  predictors, when \code{x.test} is \code{NULL}) - the combined-draws convention \code{\link{extract}}
  uses for a single chain, with no chain dimension.
}
\seealso{
  \code{\linkS4class{dbartsSampler}} for \code{sampleTreesFromPrior}, \code{sampleNodeParametersFromPrior},
  and \code{predict}, the building blocks this function composes. \code{\link{dbarts}} for
  creating a sampler. \code{\link{extract}} for the same \code{"ev"}/\code{"ppd"} vocabulary applied to
  posterior draws.
}
\examples{
n <- 100
x <- matrix(runif(n * 2), n, 2)
y <- x[, 1] - x[, 2] + rnorm(n, 0, 0.2)

sampler <- dbarts(y ~ x, control = dbartsControl(n.chains = 1L, n.threads = 1L))

## prior-predictive draws of a treatment effect f(x1) - f(x0), holding the
## other covariate fixed at its midpoint: calibrates the prior before ever
## touching y. Stacking both settings into one x.test pairs the forests, so
## the two columns of the result come from the same drawn function.
xContrast <- rbind(c(0, 0.5), c(1, 0.5))
colnames(xContrast) <- colnames(sampler$data@x)

f <- samplePriorPredictive(sampler, x.test = xContrast, n.samples = 500L)
effect <- f[, 2] - f[, 1]
hist(effect, main = "prior on the treatment effect", xlab = "f(x1) - f(x0)")
}
