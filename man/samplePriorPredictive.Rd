\name{samplePriorPredictive}
\alias{samplePriorPredictive}
\title{Draw From The BART Prior}
\description{
  Draws repeatedly from the prior a \code{\linkS4class{dbartsSampler}} runs under - the
  hyperparameters that carry a prior, tree structure, leaf values, and optionally the observation
  model - and evaluates the resulting forest at \code{x.test}.
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
    \code{"ev"} draws \eqn{E[y \mid x]}{E[y | x]} under the prior, matching what
    \code{\link{extract}} and \code{predict} call \code{"ev"}: the forest sum on the response scale
    for a gaussian or Student-t sampler; the link-transformed probability for a probit or logistic
    one; the mean count \eqn{\exp(\eta)}{exp(eta)} for a negative-binomial one, the center and
    \code{offset.test} inside \eqn{\eta}{eta}; the \eqn{K}{K} category probabilities for an
    ordinal or multinomial one. \code{"ppd"} adds the observation model on top: gaussian noise
    with sigma freshly drawn from its own prior for each sample; for Student-t, sigma times
    \eqn{t_\nu}{t_nu} noise at the family's df, or at the sample's df drawn from its prior when the
    df is estimated; a bernoulli draw from the \code{"ev"} probabilities for a binary sampler; a
    negative-binomial count at the sample's shape, fixed or drawn; an integer category code
    \eqn{1, \ldots, K}{1, ..., K} for an ordinal sampler, and for a multinomial one fitted to one
    category per row; a count vector of one trial for a multinomial sampler fitted to count rows,
    as \code{predict} draws at a new row. On an aft (survival) sampler the draw is on the log-time
    scale (the model fits \eqn{\log T}{log T}), and on a heteroscedastic one the noise is drawn at
    a prior draw of that observation's own \eqn{s(x)}{s(x)} rather than at a scalar sigma.
    Unlike \code{\link{extract}}'s posterior \code{"ppd"}, prior predictive draws are weight-blind:
    they use unit weights rather than a weighted precision (\eqn{\sigma / \sqrt{w}}{sigma / sqrt(w)}
    for gaussian) or a \eqn{\mathrm{Binomial}(w, p)}{Binomial(w, p)} success count (weighted
    logistic); supply any weighted precision downstream.
  }
  \item{offset.test}{
    A numeric offset added to the forest sum, as in \code{\link[=dbartsSampler-class]{predict}};
    recycled to the number of rows of \code{x.test} if of length 1, or \code{NULL} (the default) for
    none. On a multinomial sampler it is \code{predict}'s matrix of category offsets, one row per
    row of \code{x.test} and one column per category; a sampler carrying a category offset requires
    it (an all-zero matrix for the offset-free surface).
  }
  \item{n.threads}{
    Forwarded to the underlying \code{\link[=dbartsSampler-class]{predict}} calls.
  }
}
\details{
  Each of the \code{n.samples} draws independently redraws the hyperparameters, every tree's
  structure (\code{sampleTreesFromPrior}) and leaf values (\code{sampleLeafParametersFromPrior})
  from scratch, then evaluates the resulting forest at \code{x.test} - there is no MCMC dependence
  between rows of the result, unlike a posterior sample from \code{run}.

  The hyperparameters drawn are those the sampler draws: \eqn{k}{k} from its \code{chi()}
  hyperprior (which \code{normal(sd = invchi())} is), a Student-t df from its prior on the
  family's grid when the df is estimated, a negative-binomial shape from its prior on its grid when
  the shape is estimated, and an ordinal sampler's cutpoints (the first at 0, each gap log-normal),
  one set per sample shared by every row. A fixed \eqn{k}{k}, df or shape is used as it stands.
  An improper leaf prior, \code{invchi()} with scale 0 or \code{chi()} with an infinite scale, has
  no prior draws and is refused.

  Under \code{split.probs = dart()} the prior the sampler runs is not the stated Dirichlet on the
  split probabilities alone: the sampler never holds a tree with an empty leaf, so its prior on the
  split probabilities and the trees together is the stated one restricted to forests in which no
  tree has an empty leaf. Each sample draws the concentration and the split probabilities from the
  stated prior and every tree from the tree prior given them, and redraws the whole, split
  probabilities included, while any tree has an empty leaf, with no limit on the attempts, so
  every sample is an exact draw of that prior. How often a whole draw is kept falls with the number
  of trees and rises with the number of rows: about one in a hundred at 60 rows, 8 columns and 200
  trees, where 1,000 samples take a few seconds. With few rows and many trees the draw can take
  long: at 8 columns, 1,000 samples take about a minute at 60 rows and 500 trees, and one and a
  half to two and a half hours at 60 rows and 1,000 trees or at 30 rows and 500 trees. Once the
  attempts made project a call longer than about a minute, a warning, printed when issued, gives the
  share of draws kept so far and the projected time. The draw can be interrupted as a run is; the
  sampler passed in is left as it was.

  For \code{type = "ppd"} on a gaussian sampler, sigma is drawn once per sample from the model's
  residual-variance prior: a scaled-inverse-chi-squared calibrated so that
  \eqn{P(\sigma < \code{sigest}) = \code{quant}}{P(sigma < sigest) = quant} when the prior is
  \code{chisq} (see \code{\link{dbartsPriors}}), or the sampler's fixed value when the prior is
  \code{fixed}.

  A heteroscedastic sampler (\code{\link{dbarts}}'s \code{variance}) has no single sigma to draw,
  so the noise is its own surface: each sample additionally redraws the variance forest from its
  prior (\code{sampleVarianceForestFromPrior}) and adds \eqn{s(x)\epsilon}{s(x) eps} at the rows
  being predicted, \eqn{s^2(x)}{s^2(x)} read there through
  \code{\link[=dbartsSampler-class]{getVariance}}. The mean draws, and every \code{type = "ev"}
  draw, are unaffected.

  Prior draws are chain-free: the private sampler runs one chain whatever \code{sampler}'s
  \code{n.chains}, since there is no mixing or convergence concept for draws that are already
  independent by construction.

  Seeding follows from the private sampler being constructed fresh on every call. When the
  \code{\link[=dbartsControl]{control}}'s \code{seed} is \code{NULL} (the default), the engine RNG
  is seeded from R's random number stream at construction, so successive calls return independent
  draws, and a \code{\link{set.seed}} beforehand makes the whole draw - forests and observation
  layer alike - reproducible. When the control fixes \code{seed}, the engine stream is pinned:
  every call replays the identical sequence of hyperparameters and forests, with only the R-side
  \code{"ppd"} observation layer (sigma, the noise, counts and categories) still governed by
  \code{set.seed}.

  To compare the prior at two or more covariate settings with paired forests - as for a treatment
  effect - stack the settings as rows of a single \code{x.test} and difference the corresponding
  columns of the result; see \sQuote{Examples}. Two separate calls draw independent sets of forests
  (unless \code{seed} is fixed, as above) and so cannot be differenced meaningfully.
}
\value{
  A matrix with \code{n.samples} rows and one column per row of \code{x.test} (or of the training
  predictors, when \code{x.test} is \code{NULL}) - the combined-draws convention \code{\link{extract}}
  uses for a single chain, with no chain dimension. For an ordinal or multinomial sampler
  \code{"ev"} is an \code{n.samples} by rows by \eqn{K}{K} array of category probabilities, as is
  \code{"ppd"} on a multinomial sampler fitted to count rows, holding the counts; category codes
  are integers.
}
\seealso{
  \code{\linkS4class{dbartsSampler}} for \code{sampleTreesFromPrior}, \code{sampleLeafParametersFromPrior},
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
