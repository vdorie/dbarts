#!/usr/bin/env Rscript

# Statistical-equivalence harness: verifies that two engine builds target the
# same posterior (docs/design/core-generalization.md). Data are fixed per
# scenario; only the MCMC seed varies, so across-seed spread is pure Monte
# Carlo variability. Record a baseline with build A, then compare from a tree
# with build B installed; per-summary Welch z-statistics should look standard
# normal if the posteriors agree.
#
# If both builds produce identical RNG streams the comparison reports an
# exact match instead, which doubles as a no-unintended-RNG-shift check for
# refactors of the current engine.
#
# Usage:
#   Rscript equivalence.R record [out.rds]
#   Rscript equivalence.R compare baseline.rds [--strict-coverage]
# Append 'quick' for a fast smoke test (not comparable to full runs).
# --strict-coverage fails compare (instead of warning) when the installed
# engine has scenarios the baseline predates.
# EQUIVALENCE_SCENARIOS (comma-separated names) restricts a run to a subset
# of scenarios, for targeted record/compare passes; unset runs everything.
# Compare treats scenarios missing on either side gracefully: baseline-only
# ones report "skipped", run-only ones report as uncovered (a warning unless
# --strict-coverage), and the exit code stays 0 when every compared scenario
# is clean.

suppressPackageStartupMessages(library(dbarts))

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args
args <- setdiff(args, "quick")
strictCoverage <- "--strict-coverage" %in% args
args <- setdiff(args, "--strict-coverage")
mode <- if (length(args) >= 1L) args[[1L]] else "record"

n.seeds <- if (quick) 3L else 20L
ndpost <- if (quick) 250L else 1000L
nskip <- if (quick) 100L else 500L
ntree <- if (quick) 50L else 200L
n.test <- 25L

friedman <- function(x) {
  10 *
    sin(pi * x[, 1L] * x[, 2L]) +
    20 * (x[, 3L] - 0.5)^2 +
    10 * x[, 4L] +
    5 * x[, 5L]
}

makeScenarios <- function() {
  result <- list()

  set.seed(5101L)
  x <- matrix(runif(500L * 10L), 500L)
  result$friedman <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE
  )

  set.seed(5102L)
  x <- matrix(runif(1000L * 10L), 1000L)
  result$probit <- list(
    x = x,
    y = rbinom(1000L, 1L, pnorm(scale(friedman(x)))),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE
  )

  set.seed(5103L)
  x <- matrix(runif(500L * 10L), 500L)
  weights <- sample(c(1, 2), 500L, replace = TRUE)
  result$weighted <- list(
    x = x,
    y = friedman(x) + rnorm(500L) / sqrt(weights),
    weights = weights,
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE
  )

  # exercises the non-uniform split-variable selection path (DART's seam)
  set.seed(5104L)
  x <- matrix(runif(500L * 10L), 500L)
  result$splitprobs <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    splitprobs = c(rep(0.15, 5L), rep(0.05, 5L))
  )

  # binary with a chi hyperprior on k; driven through the sampler
  # API (old engine) and the package's bartcore bridge (new engine), since
  # bart() cannot express an adaptive k. Scenarios pin ALL prior settings
  # explicitly so package default changes cannot shift the anchor.
  set.seed(5105L)
  x <- matrix(runif(600L * 10L), 600L)
  result$chik <- list(
    x = x,
    y = rbinom(600L, 1L, pnorm(scale(friedman(x)))),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    samplerApi = TRUE,
    samplerArgs = list(
      node.prior = dbarts:::normal(dbarts:::chi(1.5, Inf))
    )
  )

  # two chains pooled, driven through the sampler API of both engines; each
  # chain has its own generator seeded from R's stream
  set.seed(5107L)
  x <- matrix(runif(500L * 10L), 500L)
  result$chains <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    nChains = 2L
  )

  # whole-data replacement mid-chain: burn in on one draw of the process,
  # setData to a larger draw (cuts rebuild, splits remap, arrays resize),
  # then re-burn only briefly so the carried-over tree state matters
  set.seed(5108L)
  x <- matrix(runif(400L * 10L), 400L)
  x2 <- matrix(runif(500L * 10L), 500L)
  result$setdata <- list(
    x = x,
    y = friedman(x) + rnorm(400L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    setData = list(x = x2, y = friedman(x2) + rnorm(500L))
  )

  # weights swapped and a test offset installed mid-chain: exercises the
  # setWeights pointer swap (weighted node statistics and sigma draws) and
  # the test-offset addition to recorded test fits; created with a scalar
  # offset so the auto-synced test offset path runs too
  set.seed(5109L)
  x <- matrix(runif(500L * 10L), 500L)
  weights <- runif(500L, 0.5, 2)
  result$wtoffset <- list(
    x = x,
    y = friedman(x) + 0.4 + rnorm(500L) / sqrt(weights),
    weights = weights,
    offset = 0.4,
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    mutate = list(
      weights = runif(500L, 0.5, 2),
      offset.test = seq(-1, 1, length.out = n.test)
    )
  )

  # quantile cut points over a mix of continuous columns (thinned to numcut)
  # and discrete ones (columns 4-5 and 8-10 on an 11-level grid, inducing
  # 10 unique-midpoint cuts each)
  set.seed(5106L)
  x <- matrix(runif(500L * 10L), 500L)
  x.test <- matrix(runif(n.test * 10L), n.test)
  discrete <- c(4L, 5L, 8L, 9L, 10L)
  x[, discrete] <- round(x[, discrete], 1L)
  x.test[, discrete] <- round(x.test[, discrete], 1L)
  result$quants <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = x.test,
    binary = FALSE,
    usequants = TRUE
  )

  # --- 1.0-0 feature paths (all driven through the sampler API so node.prior,
  # tree.prior, family, factor data, missingness, and zero weights each reach
  # their own sampling code; the equivalence gate otherwise touches none) ---

  # factor predictors: categorical split rules and the mask machinery
  set.seed(5110L)
  xc <- data.frame(
    a = runif(500L),
    b = runif(500L),
    f = factor(sample(letters[1:4], 500L, replace = TRUE)),
    g = factor(sample(LETTERS[1:3], 500L, replace = TRUE))
  )
  fcat <- 10 *
    sin(pi * xc$a * xc$b) +
    5 * (as.integer(xc$f) - 2L) +
    3 * (as.integer(xc$g) - 1L)
  xc.test <- data.frame(
    a = runif(n.test),
    b = runif(n.test),
    f = factor(
      sample(letters[1:4], n.test, replace = TRUE),
      levels = levels(xc$f)
    ),
    g = factor(
      sample(LETTERS[1:3], n.test, replace = TRUE),
      levels = levels(xc$g)
    )
  )
  result$categorical <- list(
    x = xc,
    y = fcat + rnorm(500L),
    x.test = xc.test,
    binary = FALSE,
    samplerApi = TRUE
  )

  # missing predictors: MIA routing (reserved codes + the extra Bernoulli on
  # NA-bearing columns). y is built from the complete design, then NAs injected
  set.seed(5111L)
  xm.full <- matrix(runif(500L * 10L), 500L)
  ym <- friedman(xm.full) + rnorm(500L)
  xm <- xm.full
  xm[sample.int(500L * 10L, floor(500L * 10L * 0.08))] <- NA
  xm.test <- matrix(runif(n.test * 10L), n.test)
  xm.test[sample.int(n.test * 10L, floor(n.test * 10L * 0.08))] <- NA
  result$missing <- list(
    x = xm,
    y = ym,
    x.test = xm.test,
    binary = FALSE,
    samplerApi = TRUE,
    samplerArgs = list(missing = "incorporate")
  )

  # adaptive DART tree prior: the Dirichlet split-probability updates
  set.seed(5112L)
  x <- matrix(runif(500L * 10L), 500L)
  result$dart <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    samplerArgs = list(tree.prior = dbarts:::dart(2, 0.95))
  )

  # linear leaves with a chi hyperprior on k: conjugate ridge draws + sampled k
  set.seed(5113L)
  x <- matrix(runif(500L * 10L), 500L)
  result$linear <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    samplerArgs = list(
      node.prior = dbarts:::linear(c(1L, 4L), k = dbarts:::chi(1.25, Inf))
    )
  )

  # GP leaves with a chi hyperprior: marginal/Matheron draw + kernel cache + k
  set.seed(5114L)
  x <- matrix(runif(300L * 10L), 300L)
  result$gp <- list(
    x = x,
    y = friedman(x) + rnorm(300L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    samplerArgs = list(
      node.prior = dbarts:::gp(
        1L,
        k = dbarts:::chi(1.25, Inf),
        max.leaf.size = 100L
      )
    )
  )

  # logistic family: the Polya-Gamma latent path (distinct from probit)
  set.seed(5115L)
  x <- matrix(runif(600L * 10L), 600L)
  result$logistic <- list(
    x = x,
    y = rbinom(600L, 1L, plogis(scale(friedman(x)))),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    samplerApi = TRUE,
    samplerArgs = list(
      family = "logistic",
      node.prior = dbarts:::normal(dbarts:::chi(1.5, Inf))
    )
  )

  # logistic family with integer count weights: the weighted Polya-Gamma path,
  # where a count-w observation's omega is the sum of w PG(1, psi) draws
  set.seed(5121L)
  x <- matrix(runif(600L * 10L), 600L)
  result$wtlogistic <- list(
    x = x,
    y = rbinom(600L, 1L, plogis(scale(friedman(x)))),
    weights = sample(1:3, 600L, replace = TRUE),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    samplerApi = TRUE,
    samplerArgs = list(
      family = "logistic",
      node.prior = dbarts:::normal(dbarts:::chi(1.5, Inf))
    )
  )

  # exact-zero training weights: the constant-leaf no-likelihood path and the
  # collapseEmptyNodes weightTotal > 0 guard
  set.seed(5116L)
  x <- matrix(runif(500L * 10L), 500L)
  weights <- sample(c(0, 1, 2), 500L, replace = TRUE, prob = c(0.15, 0.5, 0.35))
  result$zeroweights <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    weights = weights,
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE
  )

  # active-row mask (docs/plans/latent-subset-mask.md) on a probit sampler:
  # the family that accepts no case weights at all still takes a
  # between-draws 0/1 row mask, exercised mid-chain through $setActiveRows
  # exactly as wtoffset exercises setWeights above - the counterpart to
  # zeroweights for a family zero weights cannot reach.
  set.seed(5132L)
  x <- matrix(runif(600L * 10L), 600L)
  result$maskprobit <- list(
    x = x,
    y = rbinom(600L, 1L, pnorm(scale(friedman(x)))),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    samplerApi = TRUE,
    samplerArgs = list(family = "probit"),
    mutate = list(activeRows = as.double(rbinom(600L, 1L, 0.75)))
  )

  # the same mask on an ordinal sampler: two family-specific sums
  # (computeScales, the cutpoint proposal; cutpointLogAcceptance, its
  # target) ride the mask beside the shared latent-skip machinery the probit
  # scenario above exercises, so this is the one guard on the cutpoint pass
  # moving under a mask installed mid-chain.
  set.seed(5133L)
  x <- matrix(runif(500L * 10L), 500L)
  latentMask <- as.vector(scale(friedman(x))) + rnorm(500L)
  result$maskordinal <- list(
    x = x,
    y = cut(
      latentMask,
      c(-Inf, -0.8, 0.2, 1.0, Inf),
      labels = c("a", "b", "c", "d"),
      ordered_result = TRUE
    ),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    samplerApi = TRUE,
    samplerArgs = list(family = "ordinal"),
    mutate = list(activeRows = as.double(rbinom(500L, 1L, 0.75)))
  )

  # sparse (dgCMatrix) predictors: the rank-bitmap column store and sparse
  # partition kernel. x.test stays dense (sparse test input is unsupported).
  # Skipped when Matrix is unavailable so the core gate still runs.
  if (requireNamespace("Matrix", quietly = TRUE)) {
    set.seed(5117L)
    dense <- matrix(runif(500L * 10L), 500L)
    mask <- matrix(runif(500L * 10L) < 0.15, 500L) # ~85% structural zeros
    dense[!mask] <- 0
    xs <- as(dense, "CsparseMatrix")
    result$sparse <- list(
      x = xs,
      y = friedman(dense) + rnorm(500L),
      x.test = matrix(runif(n.test * 10L), n.test),
      binary = FALSE,
      samplerApi = TRUE
    )
  }

  # --- gate-hardening-1.0 scenarios (poison-sweep blind spots; see
  # docs/plans/gate-blindspot-audit.md and gate-hardening-1.0.md) ---

  # GP leaves under NON-UNIT weights: the weighted GP nugget
  # residualVariance / w_i (poison 16 was blind to every R gate because the
  # gp scenario is unweighted). All weights strictly positive so the score
  # path's nugget - not the zero-weight fallback - is exercised.
  set.seed(5118L)
  x <- matrix(runif(300L * 10L), 300L)
  weights <- runif(300L, 0.5, 2)
  result$wtgp <- list(
    x = x,
    y = friedman(x) + rnorm(300L) / sqrt(weights),
    weights = weights,
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    samplerArgs = list(
      node.prior = dbarts:::gp(
        1L,
        k = dbarts:::chi(1.25, Inf),
        max.leaf.size = 100L
      )
    )
  )

  # grouped random intercepts (rbart_vi's in-core Gibbs blocks) under
  # NON-UNIT weights: the weighted group-precision accumulation (poison 12
  # was cpp-only; no grouped scenario existed and the tinytest "catch" was a
  # hang). The gamma tau prior, NOT the half-Cauchy default: an SBC
  # experiment found the half-Cauchy tail can stall the tau slice sampler at
  # extreme draws. Weights spread over [0.5, 3] so the per-group weight sum
  # differs materially from the member count.
  set.seed(5119L)
  x <- matrix(runif(400L * 10L), 400L)
  groups <- factor(sample(8L, 400L, replace = TRUE))
  groupEffects <- rnorm(8L, 0, 1.5)
  weights <- runif(400L, 0.5, 3)
  result$grouped <- list(
    x = x,
    y = friedman(x) +
      groupEffects[as.integer(groups)] +
      rnorm(400L) / sqrt(weights),
    weights = weights,
    group.by = groups,
    group.by.test = factor(
      sample(8L, n.test, replace = TRUE),
      levels = levels(groups)
    ),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    rbart = TRUE
  )

  # chi hyperprior with a df where the chi-k SHAPE term is statistically
  # separable (poison 8, the shape mislabel 0.5(M + 2 nu - 1) for
  # 0.5(M + nu), was invisible at nu = 1.5: a +0.25 shift against M ~ 500
  # leaves). With nu = 50 and 20 trees (M ~ 60 leaves at this n), the
  # mislabel's +0.5(nu - 1) = +24.5 on a posterior shape of ~55 lifts
  # k.mean from 4.14 (across-seed sd 0.16) to 5.2-6.3 - zero across-seed
  # overlap - and can run away entirely (the k feedback turns unstable), so
  # the compare's disjoint-range channel catches it even when a diverged
  # seed degenerates the Welch z.
  set.seed(5120L)
  x <- matrix(runif(400L * 10L), 400L)
  result$chik2 <- list(
    x = x,
    y = friedman(x) + rnorm(400L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    nTrees = 20L,
    samplerArgs = list(node.prior = dbarts:::normal(dbarts:::chi(50, Inf)))
  )

  # grouped AFT survival: GroupedResponse(AFTResponse), riAFTBART's model,
  # newly reachable from rbart_vi via family = "aft". Random intercepts and
  # right-censoring together, on the log-time scale (a modest signal keeps
  # exp(log T) finite). Gamma tau prior, in-core fast path, no weights (aft
  # refuses them). Added after equivalence-de67cbb.rds, so a compare against
  # that baseline reports it as skipped; the anchor re-records at landing.
  set.seed(5121L)
  x <- matrix(runif(400L * 10L), 400L)
  groups <- factor(sample(8L, 400L, replace = TRUE))
  groupEffects <- rnorm(8L, 0, 0.8)
  log.t <- 0.05 * friedman(x) + groupEffects[as.integer(groups)] + rnorm(400L)
  cens.t <- 0.05 * friedman(x) + quantile(rnorm(4000L), 0.7) + rnorm(400L)
  status <- as.numeric(log.t <= cens.t)
  result$grouped_aft <- list(
    x = x,
    y = exp(ifelse(status == 1, log.t, cens.t)),
    status = status,
    group.by = groups,
    group.by.test = factor(
      sample(8L, n.test, replace = TRUE),
      levels = levels(groups)
    ),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    rbart = TRUE,
    aft = TRUE
  )

  # Student-t (robust) residuals in the estimated-nu mode (TResponse,
  # docs/design/robust-errors.md), newly reachable via resid.dist = student().
  # Contaminated-normal data - a gaussian bulk with a 5% heavy-outlier tail -
  # so both mixture channels do real work: the per-observation lambda draws
  # downweight the contaminants and the capped-grid nu draw responds to the
  # tail. The default gaussian() error law adds no draws (the existing anchors
  # are untouched), so this scenario is new and a compare against an earlier
  # baseline reports it skipped; the anchor re-records at landing.
  set.seed(5122L)
  x <- matrix(runif(500L * 10L), 500L)
  noise <- rnorm(500L)
  contaminated <- sample(500L, 25L)
  noise[contaminated] <- noise[contaminated] +
    sample(c(-1, 1), 25L, replace = TRUE) * runif(25L, 8, 12)
  result$student <- list(
    x = x,
    y = friedman(x) + noise,
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    samplerArgs = list(resid.dist = dbarts:::student())
  )

  # ordinal (cumulative-probit) responses (OrdinalResponse,
  # docs/design/ordinal.md), newly reachable via family = "ordinal": an
  # ordered 4-level outcome cut from a latent continuous truth, so the
  # doubly-truncated latent draws, the marginal-MH cutpoint block, and the
  # category-probability reporting all do real work. Recorded channels: the
  # latent test fits (fhat.test), varcount, the K-1 cutpoint draws, and the
  # posterior-mean test category probabilities (fitViaOrdinal/fitSummaries).
  # binary = TRUE skips the sigma summaries (sigma is fixed at 1). New after
  # equivalence-31dc05a.rds, so a compare against that baseline reports it
  # skipped; the anchor re-records at landing.
  set.seed(5123L)
  x <- matrix(runif(400L * 10L), 400L)
  latent <- as.vector(scale(friedman(x))) + rnorm(400L)
  result$ordinal <- list(
    x = x,
    y = cut(
      latent,
      c(-Inf, -0.8, 0.2, 1.0, Inf),
      labels = c("a", "b", "c", "d"),
      ordered_result = TRUE
    ),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    ordinalFit = TRUE
  )

  # negative-binomial counts (NBResponse, docs/design/negative-binomial.md),
  # newly reachable via family = "nbinom": overdispersed counts drawn from an
  # NB with a nonlinear log-mean, estimated dispersion (the default). The
  # Polya-Gamma count augmentation, the grid r update, and the mean-count
  # reporting all do real work. Recorded channels: the latent test fits
  # (fhat.test = psi), varcount, the per-draw dispersion r, and the
  # posterior-mean test mean counts (fitViaNbinom/fitSummaries). The omega
  # augmentation drives the trees, so its stream is locked transitively through
  # these downstream channels (the ordinal precedent, which likewise locked its
  # z augmentation through the latent fit). binary = TRUE skips the sigma
  # summaries (sigma is fixed at 1). New after equivalence-227f46a.rds, so a
  # compare against that baseline reports it skipped; the anchor re-records at
  # landing.
  set.seed(6120L)
  x <- matrix(runif(300L * 10L), 300L)
  mu <- exp(0.6 * as.vector(scale(friedman(x))))
  result$nbinom <- list(
    x = x,
    y = as.double(rnbinom(300L, size = 4L, mu = mu)),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    nbinomFit = TRUE
  )

  # discrete-time hazard (person-period expansion + binary family,
  # docs/design/survival.md, "Discrete-time hazard"): the family adds NO engine
  # code - it remaps to probit before the sampler is built and runs an ordinary
  # binary fit on the expanded rows - so this scenario just drives the probit
  # path on person-period data, adding no draw to any existing family's stream.
  # A small integer grid (K periods) keeps the expansion modest. Recorded
  # channels: the per-draw survival probabilities on the held-out subjects at a
  # short horizon (fhat.test), varcount over the expanded design (its trailing
  # column is the ordinal period predictor), and the posterior-mean survival
  # surface over several horizons (surv.test). New after equivalence-1d9388d.rds,
  # so a compare against that baseline reports it skipped; the anchor re-records
  # at landing. binary = TRUE skips the sigma summaries.
  set.seed(6220L)
  x <- matrix(runif(300L * 10L), 300L)
  K.haz <- 6L
  baseline <- seq(-2.0, -0.6, length.out = K.haz)
  fHaz <- as.vector(scale(friedman(x)))
  eventPeriod <- rep(K.haz + 1L, 300L)
  for (i in seq_len(300L)) {
    for (k in seq_len(K.haz)) {
      if (runif(1L) < pnorm(baseline[k] + fHaz[i])) {
        eventPeriod[i] <- k
        break
      }
    }
  }
  censPeriod <- sample.int(K.haz, 300L, replace = TRUE)
  result$hazard <- list(
    x = x,
    time = as.double(pmin(eventPeriod, censPeriod, K.haz)),
    status = as.double(eventPeriod <= censPeriod & eventPeriod <= K.haz),
    x.test = matrix(runif(n.test * 10L), n.test),
    hazardTimes = seq_len(K.haz),
    binary = TRUE,
    hazardFit = TRUE
  )

  # hurdle.lognormal (semicontinuous two-part, docs/design/hurdle.md): family =
  # "hurdle.lognormal" composes a probit occupancy fit of 1{y > 0} over all n
  # with a gaussian fit of log(y) over the y > 0 subset - NO shared engine
  # code, so this scenario just drives the existing probit and gaussian paths
  # verbatim on a semicontinuous response, adding no draw to any existing
  # family's stream. A smaller n than its single-forest siblings, since each
  # seed drives two component fits. Recorded channels: the occupancy
  # probability pi(x.test) (occ.test), the positive-part log-scale linear
  # predictor f(x.test) (pos.test), and the combined natural-scale predict
  # E[y | x.test] (the standard yhat.test/fhat.test channel) - all three via
  # predict.bartHurdle's replay, which only touches the RNG for type = "ppd"
  # (unused here), so a given seed reproduces. New scenario (added after
  # equivalence-f494156.rds): a compare against that baseline reports it
  # skipped/uncovered; the anchor re-records at landing. binary = TRUE skips
  # the sigma summaries (the fit has no single top-level sigma).
  set.seed(6320L)
  x <- matrix(runif(200L * 10L), 200L)
  fHurdle <- as.vector(scale(friedman(x)))
  pi.true <- pnorm(0.6 * fHurdle - 0.3)
  mu.true <- 0.5 + 0.4 * fHurdle
  occupied <- rbinom(200L, 1L, pi.true) == 1L
  y <- numeric(200L)
  y[occupied] <- exp(rnorm(sum(occupied), mu.true[occupied], 0.5))
  result$hurdle <- list(
    x = x,
    y = y,
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    hurdleFit = TRUE
  )

  # --- predictor-mutation scenarios (docs/plans/multiforest-predictor-
  # mutation.md, "The harness this arc needs"). Before these, no anchor in the
  # trio drove any predictor mutation at all: this file's only mutation
  # surfaces were setData, setWeights and setTestOffset, so the transactional
  # two-phase path, the per-observation session and the forced refresh were
  # bitwise-uncovered. Each scenario burns in, mutates mid-chain, then runs a
  # SHORT second leg, so the carried-over tree state - exactly what a widened
  # revalidation moves - is what the recorded draws depend on. ---

  # transactional whole-matrix setPredictor: validate every tree against the
  # proposed codes, then rebuild the fits from the recovered parameters. The
  # replacement is a jitter of the design sized so the transaction is ACCEPTED
  # (5/5 probe seeds accept at this sd; 0.01 accepts 3/5), which is the arm the
  # rollback scenario below cannot reach.
  set.seed(5124L)
  x <- matrix(runif(500L * 10L), 500L)
  result$predswap <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    mutate = list(
      predictor = pmin(
        pmax(x + matrix(rnorm(500L * 10L, 0, 0.005), 500L), 0),
        1
      )
    )
  )

  # transactional single-column update: the same two-phase path entered through
  # updatePredictor, which touches one column's codes and leaves the rest of
  # the store alone (the subset the arc's pruning argument is stated over).
  # Column 3 is one friedman splits on, and the jitter is sized to be accepted.
  set.seed(5125L)
  x <- matrix(runif(500L * 10L), 500L)
  result$predcol <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    mutate = list(
      column = list(
        index = 3L,
        values = pmin(pmax(x[, 3L] + rnorm(500L, 0, 0.02), 0), 1)
      )
    )
  )

  # the per-observation update session: an independent draw for one column,
  # installed row by row under a scan permutation the ENGINE draws (from chain
  # 0's generator), with a row skipped whenever installing it would empty a
  # leaf. The only scenario that consumes that permutation, so it is the one
  # that pins the session refactor; the install mask rides into the second leg,
  # so both the permutation and the veto are locked by the recorded draws.
  set.seed(5126L)
  x <- matrix(runif(500L * 10L), 500L)
  result$predpartial <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    mutate = list(partial = list(index = 4L, values = runif(500L)))
  )

  # a transactional proposal built to be ROLLED BACK: a two-level replacement
  # column at updateCutPoints = FALSE collapses every observation onto two
  # values of the existing grid, which empties leaves in any tree splitting on
  # it, so the transaction is refused (0/5 probe seeds accept). LIMITATION,
  # stated because the scenario cannot state it itself: this gates that a
  # rejected transaction leaves the run BITWISE, and nothing more - it cannot
  # by itself distinguish an accept from a reject, since a build that wrongly
  # accepted would simply record different draws, which is what the accept
  # scenarios above are for.
  set.seed(5127L)
  x <- matrix(runif(500L * 10L), 500L)
  result$predreject <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    mutate = list(
      column = list(
        index = 5L,
        values = ifelse(seq_len(500L) %% 2L == 0L, 0.25, 0.75)
      )
    )
  )

  # the forced whole-matrix swap (the bartCause propensity pattern): no veto,
  # no rollback - every forest is re-routed against the new codes and any leaf
  # that empties is collapsed, so an arbitrary replacement design is legal.
  set.seed(5128L)
  x <- matrix(runif(500L * 10L), 500L)
  result$predforce <- list(
    x = x,
    y = friedman(x) + rnorm(500L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = FALSE,
    samplerApi = TRUE,
    mutate = list(forced = matrix(runif(500L * 10L), 500L))
  )

  # a HETEROSCEDASTIC sampler taking the forced whole-matrix swap. There was no
  # heteroscedastic scenario in this harness at all, so refreshVarianceForest -
  # the forced path's variance arm - had zero equivalence coverage. Recorded
  # channels include s2.test, the variance surface on the held-out rows, which
  # a change in the variance forest's routing moves directly rather than
  # transitively through the mean forest's weights. binary = TRUE skips the
  # sigma summaries: under a variance forest sigma is structurally pinned (the
  # bridge refuses to set it for exactly that reason), so its draws are a
  # constant and both summaries would be degenerate.
  set.seed(5129L)
  x <- matrix(runif(400L * 10L), 400L)
  result$hetforce <- list(
    x = x,
    y = friedman(x) + (0.5 + 2 * x[, 6L]) * rnorm(400L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    samplerApi = TRUE,
    samplerArgs = list(variance = varianceForest(n.trees = 40L)),
    mutate = list(forced = matrix(runif(400L * 10L), 400L))
  )

  # the two TRANSACTIONAL heteroscedastic paths, which a variance forest began
  # accepting when the two-phase revalidation and the per-observation session
  # reached it (docs/plans/multiforest-predictor-mutation.md, S3). New streams:
  # both were refused at the bridge before this tip, so they have no earlier
  # baseline and become the regression floor from here. Each records the
  # engine's verdict beside the draws - recordVerdict is opt-in, so the
  # scenarios above keep the summary vectors they recorded - so a build that
  # flipped an accept into a rollback, or moved one row's install decision,
  # fails on the verdict rather than only on the post-mutation draws. Their
  # channels include s2.test, the direct read on the variance forest's routing.
  # binary = TRUE for hetforce's reason: sigma is structurally pinned under a
  # variance forest.
  #
  # hetswap: the whole-matrix transaction, a jitter sized to be ACCEPTED, so
  # every variance tree re-routes against the new codes and its factors scatter
  # through the new partition.
  set.seed(5130L)
  x <- matrix(runif(400L * 10L), 400L)
  result$hetswap <- list(
    x = x,
    y = friedman(x) + (0.5 + 2 * x[, 6L]) * rnorm(400L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    samplerApi = TRUE,
    recordVerdict = TRUE,
    samplerArgs = list(variance = varianceForest(n.trees = 40L)),
    mutate = list(
      predictor = pmin(
        pmax(x + matrix(rnorm(400L * 10L, 0, 0.005), 400L), 0),
        1
      )
    )
  )

  # hetpartial: the per-observation session, whose cell guard now caches the
  # variance trees that split on the column beside the mean forest's. The
  # verdict channel is the whole install mask - the session answers per row -
  # and this is the only heteroscedastic scenario that consumes the engine's
  # scan permutation.
  set.seed(5131L)
  x <- matrix(runif(400L * 10L), 400L)
  result$hetpartial <- list(
    x = x,
    y = friedman(x) + (0.5 + 2 * x[, 6L]) * rnorm(400L),
    x.test = matrix(runif(n.test * 10L), n.test),
    binary = TRUE,
    samplerApi = TRUE,
    recordVerdict = TRUE,
    samplerArgs = list(variance = varianceForest(n.trees = 40L)),
    mutate = list(partial = list(index = 6L, values = runif(400L)))
  )

  # xbart's crossvalidation loop, the only scenario that drives it: each
  # replication draws a fold partition, every parameter cell is fit on the
  # training folds and scored on the held-out one, and the reported loss is
  # the fold average. The recorded channel is that loss array, so a build
  # that moved the draws shows up through the scoring rather than through a
  # fit's own summaries. Its budget is LITERAL, kept out of the guarded
  # ndpost/nskip/ntree, so settingsList() stays identical to the earlier
  # baselines and the neutrality compare against them still runs; the
  # k3counts/set_predictor precedent in the sister harnesses.
  set.seed(5134L)
  x <- matrix(runif(150L * 5L), 150L)
  result$xbart <- list(
    x = x,
    y = friedman(x) + rnorm(150L),
    binary = FALSE,
    xbartFit = TRUE
  )

  # bart2() gaussian on a formula (bart2-argument-consolidation.md section 7):
  # equivalence.R otherwise reaches bart2 only through the four
  # alternate-family fitters above, all called on the x/y matrix convenience
  # form with no data =, so the consolidated surface's own FORMULA parsing has
  # no equivalence eyes. Non-default args, both formula-interface-only and
  # untouched by every fitter above: weights and subset (350 of 400 rows
  # kept).
  set.seed(5135L)
  x <- matrix(runif(400L * 10L), 400L)
  d <- as.data.frame(x)
  d$y <- friedman(x) + rnorm(400L)
  result$bart2gauss <- list(
    data = d,
    formula = reformulate(colnames(d)[1:10], response = "y"),
    weights = runif(400L, 0.5, 2),
    subsetIdx = sort(sample.int(400L, 350L)),
    x.test = as.data.frame(matrix(runif(n.test * 10L), n.test)),
    binary = FALSE,
    bart2GaussFit = TRUE
  )

  # ditto, probit: offset (formula-interface-only, unreachable without
  # data =) and TWO pooled chains (n.chains = 2L, combineChains = TRUE) -
  # every fitter above pins n.chains = 1L, so bart2's own chain pooling has no
  # other anchor.
  set.seed(5136L)
  x <- matrix(runif(400L * 10L), 400L)
  d <- as.data.frame(x)
  d$y <- rbinom(400L, 1L, pnorm(scale(friedman(x))))
  x.test.bp <- matrix(runif(n.test * 10L), n.test)
  result$bart2probit <- list(
    data = d,
    formula = reformulate(colnames(d)[1:10], response = "y"),
    offset = as.vector(0.3 * scale(x[, 6L])),
    offsetTest = as.vector(0.3 * scale(x.test.bp[, 6L])),
    x.test = as.data.frame(x.test.bp),
    binary = TRUE,
    bart2ProbitFit = TRUE
  )

  # an amplitude-coupled two-forest bart2 fit through the S12 term route
  # (bart2-argument-consolidation.md section 5): the canonical
  # `zf:forest(x1 + x2)` sugar declares a second forest, modulated by a
  # 3-level factor's own indicator basis, alongside the main x1 + x2 forest -
  # the surface's only path to a K > 1 amplitude-coupled fit from bart2's
  # formula interface, and the B2 landing's own per-forest reporting channels
  # (forestFits, glue) have no equivalence anchor at all. Sizes are this
  # scenario's own literals, not the guarded ndpost/nskip/ntree.
  set.seed(5137L)
  n.tf <- 100L
  x1.tf <- runif(n.tf)
  x2.tf <- runif(n.tf)
  zf.tf <- factor(sample(c("a", "b", "c"), n.tf, replace = TRUE))
  d.tf <- data.frame(x1 = x1.tf, x2 = x2.tf, zf = zf.tf)
  d.tf$y <- 10 *
    sin(pi * x1.tf * x2.tf) +
    as.numeric(zf.tf == "b") * (1 + x2.tf) +
    rnorm(n.tf, sd = 0.3)
  result$bart2twoforest <- list(
    data = d.tf,
    formula = y ~ x1 + x2 + zf:forest(x1 + x2),
    binary = FALSE,
    twoforestFit = TRUE
  )

  # a dbartsMixedMatrix container (R/mixedMatrix.R, R/data.R;
  # inst/tinytest/test-data-mixed.R): a data frame holding dense numeric and
  # factor columns alongside a Matrix::sparseVector and a two-column
  # dgCMatrix, which dbartsData assembles into ONE mixed dense/sparse
  # container - the B2 landing's own predictor container, which no scenario
  # in this file drives (the "sparse" scenario above is all-sparse, never
  # mixed). Records both TRAIN (fhat.train, the container's column routing on
  # the rows the trees actually split on) and TEST (fhat.test, the
  # predict-time densification path) plus varcount, so the routing is
  # draw-guarded on both legs. Skipped when Matrix is unavailable, as
  # "sparse" already is.
  if (requireNamespace("Matrix", quietly = TRUE)) {
    set.seed(5138L)
    n.mm <- 200L
    x1.mm <- rnorm(n.mm)
    f.mm <- factor(sample(c("a", "b", "c"), n.mm, replace = TRUE))
    sv.mm <- Matrix::sparseVector(
      x = 0.5 + runif(30L),
      i = sort(sample.int(n.mm, 30L)),
      length = n.mm
    )
    sm.dense.mm <- matrix(0, n.mm, 2L)
    for (j in 1:2) {
      nz <- runif(n.mm) < 0.1
      sm.dense.mm[nz, j] <- 0.5 + runif(sum(nz))
    }
    sm.mm <- methods::as(sm.dense.mm, "CsparseMatrix")
    x.frame.mm <- data.frame(x1 = x1.mm, f = f.mm)
    x.frame.mm$sv <- sv.mm
    x.frame.mm$sm <- sm.mm
    y.mm <- 2 *
      as.double(sv.mm) -
      1.5 * x1.mm +
      0.5 * (f.mm == "b") +
      rowSums(sm.dense.mm) +
      rnorm(n.mm, sd = 0.3)
    x.test.mm <- data.frame(
      x1 = x1.mm[seq_len(n.test)],
      f = f.mm[seq_len(n.test)]
    )
    x.test.mm$sv <- as.double(sv.mm)[seq_len(n.test)]
    x.test.mm$sm <- sm.dense.mm[seq_len(n.test), , drop = FALSE]
    result$mixedmatrix <- list(
      x = x.frame.mm,
      y = y.mm,
      x.test = x.test.mm,
      binary = FALSE,
      samplerApi = TRUE,
      recordTrain = TRUE
    )
  }

  # bart2's multinomial path (docs/design/multinomial.md): the K-forest
  # softmax classifier, reachable only under an explicit family =
  # "multinomial" token (family = "auto" never reaches it). No scenario
  # above touches it, so this is the surface's first anchor in this file -
  # the engine-level bitwise gate lives in
  # benchmarks/R/multinomial-equivalence.R. Its own literal seed continues
  # the file's own literal-seed sequence, and every scenario re-seeds before
  # generating its own data, so its position perturbs nothing around it.
  set.seed(5139L)
  n.mn <- 200L
  x.mn <- matrix(runif(n.mn * 4L), n.mn, 4L)
  eta.mn <- cbind(
    2 * (x.mn[, 1L] - 0.5),
    x.mn[, 2L] - x.mn[, 3L],
    1.5 * (x.mn[, 4L] - 0.5)
  )
  probs.mn <- exp(eta.mn) / rowSums(exp(eta.mn))
  labels.mn <- vapply(
    seq_len(n.mn),
    function(i) sample.int(3L, 1L, prob = probs.mn[i, ]) - 1L,
    integer(1L)
  )
  y.mn <- factor(c("a", "b", "c")[labels.mn + 1L], levels = c("a", "b", "c"))
  result$bart2multinom <- list(
    x = x.mn,
    y = y.mn,
    x.test = matrix(runif(n.test * 4L), n.test, 4L),
    binary = FALSE,
    multinomialFit = TRUE
  )

  # --- predictor column-kind scenarios: the three predictor column shapes
  # that no scenario above carries, each guarding a different part of the
  # store's per-column typing (docs/design/data-store.md, "Cut grid"). Every
  # factor case above is built with plain factor(), all-dense and complete,
  # so an ordered factor's grid, a factor column's missing route, and a
  # CSC-backed factor's category count moved no recorded draw in this file.
  # Sizes are small and the forests are this group's own literal nTrees: the
  # corpus is a bitwise gate, not a fit. ---

  # an ORDERED-factor predictor, which this file had none of - the two
  # ordered_result = TRUE uses above are RESPONSES. R collapses is.ordered()
  # into the same ordinal type a numeric column takes, so the column reaches
  # the store as threshold splits over n.cuts uniform cut points spanning its
  # observed code range rather than one cut per level boundary. K = 150 levels
  # against the default n.cuts = 100 leaves only 101 reachable codes, so 49
  # adjacent level pairs share a code and no rule can separate them; K is also
  # above the quantile path's own thinning threshold (numCuts + 1), so a grid
  # rebuilt through that path is measured here too. Three interior levels are
  # deliberately absent from training and present in the test frame, so a grid
  # taken from observed values rather than declared ones is visible as well.
  # Two numeric nuisance columns only, so the trees reach the ordered column
  # often.
  set.seed(5140L)
  n.of <- 400L
  K.of <- 150L
  levels.of <- sprintf("L%03d", seq_len(K.of))
  absent.of <- c(37L, 88L, 131L)
  observed.of <- setdiff(seq_len(K.of), absent.of)
  codes.of <- sample(observed.of, n.of, replace = TRUE)
  x.of <- data.frame(
    a = runif(n.of),
    b = runif(n.of),
    o = factor(levels.of[codes.of], levels = levels.of, ordered = TRUE)
  )
  codes.test.of <- c(absent.of, sample(observed.of, n.test - 3L))
  x.test.of <- data.frame(
    a = runif(n.test),
    b = runif(n.test),
    o = factor(levels.of[codes.test.of], levels = levels.of, ordered = TRUE)
  )
  result$ordfactor <- list(
    x = x.of,
    y = 10 *
      sin(pi * x.of$a * x.of$b) +
      0.06 * (codes.of - K.of / 2) +
      rnorm(n.of),
    x.test = x.test.of,
    binary = FALSE,
    samplerApi = TRUE,
    nTrees = 50L
  )

  # an NA-bearing FACTOR predictor. The "missing" scenario above injects NAs
  # into a numeric matrix only, so MIA on a categorical column - the reserved
  # missing code beside the level codes, and the per-rule missing direction
  # drawn only on columns whose training values carried an NA
  # (docs/design/mia-missingness.md) - had no anchor: a build that flipped a
  # factor column's hasMissing, or dropped the extra Bernoulli there, moved
  # nothing. Missingness itself carries signal on both factor columns, so the
  # learned direction shows in the fits rather than only in the stream, and
  # the test frame carries NAs in both so the test-side routing is recorded
  # too. The numeric columns stay complete, so the two factor columns are the
  # only NA-bearing ones in the store.
  set.seed(5141L)
  n.nf <- 400L
  g.nf <- factor(sample(letters[1:4], n.nf, replace = TRUE))
  h.nf <- factor(sample(LETTERS[1:6], n.nf, replace = TRUE))
  gMissing.nf <- runif(n.nf) < 0.18
  hMissing.nf <- runif(n.nf) < 0.12
  g.nf[gMissing.nf] <- NA
  h.nf[hMissing.nf] <- NA
  x.nf <- data.frame(
    a = runif(n.nf),
    b = runif(n.nf),
    g = g.nf,
    h = h.nf
  )
  g.test.nf <- factor(
    sample(letters[1:4], n.test, replace = TRUE),
    levels = letters[1:4]
  )
  h.test.nf <- factor(
    sample(LETTERS[1:6], n.test, replace = TRUE),
    levels = LETTERS[1:6]
  )
  g.test.nf[c(3L, 11L, 19L)] <- NA
  h.test.nf[c(7L, 20L)] <- NA
  x.test.nf <- data.frame(
    a = runif(n.test),
    b = runif(n.test),
    g = g.test.nf,
    h = h.test.nf
  )
  result$nafactor <- list(
    x = x.nf,
    y = 10 *
      sin(pi * x.nf$a * x.nf$b) +
      2.5 * gMissing.nf -
      1.5 * hMissing.nf +
      ifelse(!is.na(g.nf) & g.nf == "b", 2, 0) +
      rnorm(n.nf),
    x.test = x.test.nf,
    binary = FALSE,
    samplerApi = TRUE,
    nTrees = 50L,
    samplerArgs = list(missing = "incorporate")
  )

  # a CSC-backed FACTOR column (a sparseFactor, docs/design/sparse-columns.md):
  # the one column shape whose category count comes from the container's own
  # declared level table rather than from a sweep of the values, since the
  # bridge's declared-count read skips a CSC-backed source. The "sparse"
  # scenario above is all-ordinal and "mixedmatrix" carries its factor DENSE,
  # so no scenario reached a categorical column over the sparse block at all:
  # its pooled mask tier, its code validation and its declared count were all
  # unguarded. Two of the eight levels are absent from training and present in
  # the test frame, which a count inferred from observed values would get
  # wrong - the TOP level among them, so an inferred count would come out 7
  # where the declared table says 8. Reaches the engine through the x/y
  # interface, the only entrance a sparseFactor has. Skipped when Matrix is unavailable, as "sparse" already
  # is.
  if (requireNamespace("Matrix", quietly = TRUE)) {
    set.seed(5142L)
    n.sf <- 300L
    K.sf <- 8L
    levels.sf <- paste0("s", seq_len(K.sf))
    observed.sf <- c(1L, 2L, 3L, 5L, 6L, 7L) # levels 4 and 8 unobserved
    codes.sf <- sample(
      observed.sf,
      n.sf,
      replace = TRUE,
      prob = c(0.7, rep(0.06, 5L))
    )
    x.sf <- data.frame(a = runif(n.sf), b = runif(n.sf))
    x.sf$f <- sparseFactor(
      factor(levels.sf[codes.sf], levels = levels.sf),
      reference = "s1"
    )
    x.test.sf <- data.frame(a = runif(n.test), b = runif(n.test))
    x.test.sf$f <- sparseFactor(
      factor(
        levels.sf[c(4L, 8L, sample(observed.sf, n.test - 2L, replace = TRUE))],
        levels = levels.sf
      ),
      reference = "s1"
    )
    result$sparsefactor <- list(
      x = x.sf,
      y = 10 *
        sin(pi * x.sf$a * x.sf$b) +
        3 * (codes.sf == 3L) -
        2 * (codes.sf == 6L) +
        rnorm(n.sf),
      x.test = x.test.sf,
      binary = FALSE,
      samplerApi = TRUE,
      nTrees = 50L,
      recordTrain = TRUE
    )
  }

  # --- typed-predictor-source scenarios: the store entrances that REBUILD or
  # RE-READ a factor column's storage, none of which any scenario above
  # records. The three column shapes just above enter the store once, at
  # creation; these five drive the test store on a replacement, a leaf model
  # reading a factor column's raw values, the per-observation write-through
  # into a retained dense block, and a fold view over a parent holding codes.
  # Sizes are small and the forests are this group's own literal nTrees: the
  # corpus is a bitwise gate, not a fit. ---

  # the test container REPLACED mid-chain. Every scenario above installs its
  # test set once, at creation: setTestPredictor and setTestPredictorAndOffset
  # were absent from the mutation vocabulary entirely, so the two bridge
  # funnels that rebuild the test store - and the store's own test build,
  # which re-reads the training cut grid and the training level tables against
  # fresh test values - moved no recorded draw. The design is the widest shape
  # the test side accepts: two numeric columns, a complete dense factor, an
  # NA-bearing dense factor and a CSC-backed sparseFactor, so one replacement
  # crosses the dense, the missing-code and the sparse arms of the test build
  # at once. THREE legs, each drawing its own samples: the creation container
  # (fhat.test.pre), the container setTestPredictor installs
  # (fhat.test.mid), and the container-plus-offset setTestPredictorAndOffset
  # installs (the standard fhat.test channel). A build that rebuilt one of the
  # three differently moves the leg reading it and leaves the others alone.
  # Skipped when Matrix is unavailable, as "sparsefactor" already is.
  if (requireNamespace("Matrix", quietly = TRUE)) {
    set.seed(5143L)
    n.ts <- 300L
    levels.ts <- letters[1:4]
    levelsNA.ts <- LETTERS[1:4]
    levels.sp.ts <- paste0("v", seq_len(6L))
    a.ts <- runif(n.ts)
    b.ts <- runif(n.ts)
    f.ts <- factor(sample(levels.ts, n.ts, replace = TRUE), levels = levels.ts)
    g.ts <- factor(
      sample(levelsNA.ts, n.ts, replace = TRUE),
      levels = levelsNA.ts
    )
    gMissing.ts <- runif(n.ts) < 0.15
    g.ts[gMissing.ts] <- NA
    codes.sp.ts <- sample.int(6L, n.ts, replace = TRUE)
    x.ts <- data.frame(a = a.ts, b = b.ts, f = f.ts, g = g.ts)
    x.ts$s <- sparseFactor(
      factor(levels.sp.ts[codes.sp.ts], levels = levels.sp.ts),
      reference = "v1"
    )
    # one container per leg, same shape, independent values
    testFrame.ts <- function() {
      frame <- data.frame(
        a = runif(n.test),
        b = runif(n.test),
        f = factor(
          sample(levels.ts, n.test, replace = TRUE),
          levels = levels.ts
        ),
        g = factor(
          sample(levelsNA.ts, n.test, replace = TRUE),
          levels = levelsNA.ts
        )
      )
      frame$g[c(4L, 13L)] <- NA
      frame$s <- sparseFactor(
        factor(
          levels.sp.ts[sample.int(6L, n.test, replace = TRUE)],
          levels = levels.sp.ts
        ),
        reference = "v1"
      )
      frame
    }
    result$testswap <- list(
      x = x.ts,
      y = 10 *
        sin(pi * a.ts * b.ts) +
        2 * (f.ts == "b") +
        1.5 * gMissing.ts +
        1.2 * (codes.sp.ts == 3L) +
        rnorm(n.ts),
      x.test = testFrame.ts(),
      binary = FALSE,
      samplerApi = TRUE,
      nTrees = 50L,
      samplerArgs = list(missing = "incorporate"),
      testSwap = list(
        predictor = testFrame.ts(),
        predictorAndOffset = testFrame.ts(),
        offset.test = seq(-0.5, 0.5, length.out = n.test)
      )
    )
  }

  # an ORDERED FACTOR designated as a LINEAR-LEAF covariate, on an all-dense
  # container. The linear and gp scenarios above designate plain numeric
  # columns, so a leaf covariate that is a FACTOR - admissible, since only a
  # categorical column is refused as a covariate - had no anchor at all: the
  # leaf model reads the column's raw values, standardizes them by the store's
  # own moments and solves the per-leaf ridge over them, on the training side
  # from the columns the store gathers at build and on the test side from the
  # test block. Two covariates, one numeric and one the ordered factor, so the
  # block solve mixes the two scales.
  set.seed(5144L)
  n.lf <- 300L
  K.lf <- 8L
  levels.lf <- sprintf("O%d", seq_len(K.lf))
  codes.lf <- sample.int(K.lf, n.lf, replace = TRUE)
  x.lf <- data.frame(
    a = runif(n.lf),
    b = runif(n.lf),
    o = factor(levels.lf[codes.lf], levels = levels.lf, ordered = TRUE)
  )
  x.test.lf <- data.frame(
    a = runif(n.test),
    b = runif(n.test),
    o = factor(
      levels.lf[sample.int(K.lf, n.test, replace = TRUE)],
      levels = levels.lf,
      ordered = TRUE
    )
  )
  y.lf <- 10 *
    sin(pi * x.lf$a * x.lf$b) +
    0.8 * (codes.lf - K.lf / 2) +
    rnorm(n.lf)
  result$leaffactor <- list(
    x = x.lf,
    y = y.lf,
    x.test = x.test.lf,
    binary = FALSE,
    samplerApi = TRUE,
    nTrees = 50L,
    samplerArgs = list(node.prior = dbarts:::linear(c(1L, 3L)))
  )

  # the SAME designation on a MIXED container - a two-column dgCMatrix beside
  # the dense numeric and ordered-factor columns. The container decides where
  # the leaf model's raw values come from: an all-dense container gathers the
  # designated columns at build, a mixed one gathers nothing and serves them
  # from the dense block it retains. Same covariates, same kinds, different
  # storage, and only this arm has a retained copy to lose. Skipped when
  # Matrix is unavailable, as "sparse" already is.
  if (requireNamespace("Matrix", quietly = TRUE)) {
    set.seed(5145L)
    sm.lf <- matrix(0, n.lf, 2L)
    for (j in 1:2) {
      nz <- runif(n.lf) < 0.1
      sm.lf[nz, j] <- 0.5 + runif(sum(nz))
    }
    x.mlf <- x.lf
    x.mlf$sm <- methods::as(sm.lf, "CsparseMatrix")
    x.test.mlf <- x.test.lf
    x.test.mlf$sm <- sm.lf[seq_len(n.test), , drop = FALSE]
    result$leaffactormixed <- list(
      x = x.mlf,
      y = y.lf + rowSums(sm.lf),
      x.test = x.test.mlf,
      binary = FALSE,
      samplerApi = TRUE,
      nTrees = 50L,
      samplerArgs = list(node.prior = dbarts:::linear(c(1L, 3L)))
    )
  }

  # a PER-OBSERVATION predictor update on a dense FACTOR column of a mixed
  # store. predpartial and hetpartial above drive the session on numeric
  # columns of a plain matrix, where the store retains no raw values at all,
  # so the write-through that installs an accepted row's value into the
  # store's own dense block - and the snapshot/restore pair a vetoed row rolls
  # back through - ran over no factor column anywhere in this file. A mixed
  # store is what makes that write-through live, and the bridge deliberately
  # keeps a mixed store's dense factor columns open to the session (the
  # latent-in-a-sparse-design case). The replacement collapses the column onto
  # two of its five levels, predreject's shape, which is what makes the VETO
  # fire: a uniform redraw installs every row, while this one leaves a
  # seed-dependent handful rolled back, so the snapshot/restore arm is
  # recorded too. The install mask rides the verdict channel, so a moved
  # row-level decision fails here rather than only downstream. Skipped when
  # Matrix is unavailable.
  if (requireNamespace("Matrix", quietly = TRUE)) {
    set.seed(5146L)
    n.fp <- 300L
    K.fp <- 5L
    a.fp <- runif(n.fp)
    b.fp <- runif(n.fp)
    codes.fp <- sample.int(K.fp, n.fp, replace = TRUE)
    sm.fp <- matrix(0, n.fp, 2L)
    for (j in 1:2) {
      nz <- runif(n.fp) < 0.1
      sm.fp[nz, j] <- 0.5 + runif(sum(nz))
    }
    x.fp <- data.frame(
      a = a.fp,
      b = b.fp,
      f = factor(letters[codes.fp], levels = letters[seq_len(K.fp)])
    )
    x.fp$sm <- methods::as(sm.fp, "CsparseMatrix")
    x.test.fp <- data.frame(
      a = runif(n.test),
      b = runif(n.test),
      f = factor(
        letters[sample.int(K.fp, n.test, replace = TRUE)],
        levels = letters[seq_len(K.fp)]
      )
    )
    x.test.fp$sm <- sm.fp[seq_len(n.test), , drop = FALSE]
    result$factorpartial <- list(
      x = x.fp,
      y = 10 *
        sin(pi * a.fp * b.fp) +
        2 * (codes.fp == 2L) +
        rowSums(sm.fp) +
        rnorm(n.fp),
      x.test = x.test.fp,
      binary = FALSE,
      samplerApi = TRUE,
      nTrees = 50L,
      recordVerdict = TRUE,
      mutate = list(
        partial = list(
          index = 3L,
          values = as.double(seq_len(n.fp) %% 2L)
        )
      )
    )
  }

  # xbart over a MIXED container carrying factors. The xbart scenario above is
  # a plain numeric matrix and is this file's only driver of the fold VIEW - a
  # store built from a parent store's rows rather than from a host's values -
  # so a view over a parent holding factor codes, a CSC-backed block and a
  # per-column source map had no anchor. Both factor kinds ride the parent, so
  # the view's per-column copy is exercised on each. Its budget is LITERAL,
  # the xbart scenario's own precedent, so settingsList() stays identical to
  # the earlier baselines. Skipped when Matrix is unavailable.
  if (requireNamespace("Matrix", quietly = TRUE)) {
    set.seed(5147L)
    n.xm <- 150L
    K.xm <- 4L
    levels.ord.xm <- sprintf("O%d", seq_len(5L))
    a.xm <- runif(n.xm)
    b.xm <- runif(n.xm)
    codes.xm <- sample.int(K.xm, n.xm, replace = TRUE)
    ord.xm <- sample.int(5L, n.xm, replace = TRUE)
    sm.xm <- matrix(0, n.xm, 2L)
    for (j in 1:2) {
      nz <- runif(n.xm) < 0.12
      sm.xm[nz, j] <- 0.5 + runif(sum(nz))
    }
    x.xm <- data.frame(
      a = a.xm,
      b = b.xm,
      f = factor(letters[codes.xm], levels = letters[seq_len(K.xm)]),
      o = factor(
        levels.ord.xm[ord.xm],
        levels = levels.ord.xm,
        ordered = TRUE
      )
    )
    x.xm$sm <- methods::as(sm.xm, "CsparseMatrix")
    result$xbartmixed <- list(
      x = x.xm,
      y = 10 *
        sin(pi * a.xm * b.xm) +
        2 * (codes.xm == 2L) +
        0.8 * ord.xm +
        rowSums(sm.xm) +
        rnorm(n.xm),
      binary = FALSE,
      xbartFit = TRUE
    )
  }

  result
}

# benign warnings the fitting paths emit for legitimate scenarios: test-data
# weights (no posterior-predictive use here), zero training weights
# (deliberate in the zero-weight scenario), and the rbart_vi formula path's
# positional test-column matching and test-weight handling (the grouped
# scenario's test data carries no weights by design). Muffle only these so
# real ones stay visible.
muffleBenignWarning <- function(w) {
  msg <- conditionMessage(w)
  if (
    grepl("'weights' are ignored for test data", msg) ||
      grepl("'weights' of 0 will be ignored", msg) ||
      grepl("columns of 'test' will be matched by position", msg) ||
      grepl("weights specified but not found in test data", msg)
  ) {
    invokeRestart("muffleWarning")
  }
}

# runs through the public dbartsSampler surface on the installed package's
# engine (the control's engine flag retired with the classic engine).
# scenario$samplerArgs (a named list) is spliced into the dbarts() call so a
# scenario can select node.prior/tree.prior/family/missing without new plumbing.
fitViaSamplerApi <- function(scenario) {
  n.chains <- if (!is.null(scenario$nChains)) scenario$nChains else 1L
  control <- dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = if (!is.null(scenario$nTrees)) scenario$nTrees else ntree,
    updateState = FALSE
  )
  dbartsArgs <- list(
    scenario$x,
    scenario$y,
    test = scenario$x.test,
    control = control
  )
  if (!is.null(scenario$weights)) {
    dbartsArgs$weights <- scenario$weights
  }
  if (!is.null(scenario$offset)) {
    dbartsArgs$offset <- scenario$offset
  }
  dbartsArgs <- c(dbartsArgs, scenario$samplerArgs)
  sampler <- withCallingHandlers(
    do.call(dbarts, dbartsArgs),
    warning = muffleBenignWarning
  )
  # [d, S] or [d, S, C] -> (S * C) x d, pooling chains
  poolChains <- function(a) {
    if (length(dim(a)) == 3L) t(matrix(a, nrow = dim(a)[1L])) else t(a)
  }
  verdict <- NULL
  testPre <- NULL
  testMid <- NULL
  r <- if (!is.null(scenario$setData)) {
    sampler$run(nskip, 0L)
    sampler$setData(dbartsData(
      scenario$setData$x,
      scenario$setData$y,
      test = scenario$x.test
    ))
    sampler$run(ceiling(nskip / 4), ndpost)
  } else if (!is.null(scenario$testSwap)) {
    # the test container is replaced mid-chain, once through each bridge
    # funnel that rebuilds the test store: setTestPredictor takes the
    # container alone, setTestPredictorAndOffset takes it with a test offset.
    # Every leg draws its own samples, so each container's test store is
    # recorded in its own channel rather than inferred from the last one.
    swap <- scenario$testSwap
    testPre <- poolChains(sampler$run(nskip, ndpost)$test)
    sampler$setTestPredictor(swap$predictor)
    testMid <- poolChains(sampler$run(ceiling(nskip / 4), ndpost)$test)
    sampler$setTestPredictorAndOffset(
      swap$predictorAndOffset,
      swap$offset.test
    )
    sampler$run(ceiling(nskip / 4), ndpost)
  } else if (!is.null(scenario$mutate)) {
    sampler$run(nskip, 0L)
    if (!is.null(scenario$mutate$weights)) {
      sampler$setWeights(scenario$mutate$weights)
    }
    if (!is.null(scenario$mutate$activeRows)) {
      sampler$setActiveRows(scenario$mutate$activeRows)
    }
    if (!is.null(scenario$mutate$offset.test)) {
      sampler$setTestOffset(scenario$mutate$offset.test)
    }
    # predictor keys, read with [[ ]] rather than $: $ partial-matches on
    # lists, so a scenario carrying only one of these would silently take
    # another's branch. The transactional flavors can be REFUSED by the engine
    # (a proposal that would empty a leaf rolls back), which is a legitimate
    # outcome the scenarios exercise deliberately. The return value - the
    # accept/rollback flag, or the session's per-row install mask - rides into
    # the summaries only for a scenario that asks for it (recordVerdict), so
    # the scenarios recorded before it existed keep their summary vectors.
    mutate <- scenario$mutate
    if (!is.null(mutate[["predictor"]])) {
      verdict <- sampler$setPredictor(
        mutate[["predictor"]],
        forceUpdate = FALSE
      )
    }
    if (!is.null(mutate[["column"]])) {
      verdict <- sampler$setPredictor(
        mutate[["column"]]$values,
        column = mutate[["column"]]$index,
        forceUpdate = FALSE,
        updateCutPoints = FALSE
      )
    }
    if (!is.null(mutate[["partial"]])) {
      verdict <- sampler$setPredictor(
        mutate[["partial"]]$values,
        column = mutate[["partial"]]$index,
        forceUpdate = "partial"
      )
    }
    if (!is.null(mutate[["forced"]])) {
      sampler$setPredictor(mutate[["forced"]], forceUpdate = TRUE)
    }
    sampler$run(ceiling(nskip / 4), ndpost)
  } else {
    sampler$run(nskip, ndpost)
  }
  list(
    yhat.test = poolChains(r$test),
    varcount = poolChains(r$varcount),
    sigma = as.vector(r$sigma),
    k = if (!is.null(r$k)) as.vector(r$k),
    # a variance forest additionally reports s^2(x) on the test rows; absent -
    # and so unsummarized - for every other sampler, which leaves the existing
    # scenarios' summary vectors untouched
    variance.test = if (!is.null(r$varianceTest)) poolChains(r$varianceTest),
    verdict = if (isTRUE(scenario$recordVerdict)) as.double(verdict),
    # the in-sample fit, for a scenario that asks for it (recordTrain): the
    # mixedmatrix scenario's dbartsMixedMatrix container routes dense and
    # sparse columns alike on the rows the trees actually split on, which the
    # densified test-row predict path below does not exercise. Absent - and
    # so unsummarized - for every scenario recorded before this channel
    # existed.
    yhat.train = if (isTRUE(scenario$recordTrain) && !is.null(r$train)) {
      poolChains(r$train)
    },
    # the test fits of the containers a testSwap scenario replaces, one field
    # per leg; NULL - and so unsummarized - for every other scenario
    yhat.test.pre = testPre,
    yhat.test.mid = testMid
  )
}

# runs bart2's ordinal (cumulative-probit) path (docs/design/ordinal.md):
# family = "ordinal" explicitly, so no auto announcement. yhat.test carries the
# LATENT eta draws on the test rows (the identified quantity on the unit latent
# scale, the standard fhat.test channel shape); the ordinal-only channels -
# per-draw cutpoints and the posterior-mean test category probabilities - ride
# their own fields, summarized by fitSummaries' guarded blocks.
fitViaOrdinal <- function(scenario) {
  fit <- bart2(
    scenario$x,
    scenario$y,
    test = scenario$x.test,
    family = "ordinal",
    n.samples = ndpost,
    n.burn = nskip,
    n.trees = ntree,
    n.chains = 1L,
    n.threads = 1L,
    combineChains = TRUE,
    verbose = FALSE
  )
  list(
    yhat.test = fit$latent.test,
    varcount = fit$varcount,
    cutpoints = fit$thresholds,
    probs.test = fit$yhat.test
  )
}

# runs bart2's negative-binomial count path (docs/design/negative-binomial.md):
# family = "nbinom" explicitly (a count response is never auto), estimated
# dispersion (the default). yhat.test carries the LATENT psi draws on the test
# rows (the identified log-odds quantity, the standard fhat.test channel shape);
# the nbinom-only channels - the per-draw dispersion r and the posterior-mean
# test mean counts - ride their own fields, summarized by fitSummaries' guarded
# blocks.
fitViaNbinom <- function(scenario) {
  fit <- bart2(
    scenario$x,
    scenario$y,
    test = scenario$x.test,
    family = "nbinom",
    n.samples = ndpost,
    n.burn = nskip,
    n.trees = ntree,
    n.chains = 1L,
    n.threads = 1L,
    combineChains = TRUE,
    verbose = FALSE
  )
  list(
    yhat.test = fit$latent.test,
    varcount = fit$varcount,
    dispersion = fit$dispersion,
    means.test = fit$yhat.test
  )
}

# runs bart2's discrete-time hazard path (docs/design/survival.md): family =
# "hazard" (probit link) on person-period-expanded (time, status) data,
# keepTrees so survivalProbabilities can replay the trees onto the held-out
# subjects. yhat.test carries the SURVIVAL probability draws at the first
# horizon (draws x n.test, the standard fhat.test channel shape); the
# hazard-only channel - the posterior-mean survival surface over several
# horizons - rides surv.test, summarized by fitSummaries' guarded block. Only
# type = "ev" prediction is used (draw-neutral), so a given seed reproduces.
fitViaHazard <- function(scenario) {
  fit <- bart2(
    scenario$x,
    cbind(scenario$time, scenario$status),
    family = "hazard",
    n.samples = ndpost,
    n.burn = nskip,
    n.trees = ntree,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    combineChains = TRUE,
    verbose = FALSE
  )
  sp <- survivalProbabilities(
    fit,
    scenario$hazardTimes,
    newdata = scenario$x.test
  )
  list(
    yhat.test = sp[, 1L, ],
    varcount = fit$varcount,
    surv.test = apply(sp, c(2L, 3L), mean)
  )
}

# runs bart2's hurdle.lognormal path (docs/design/hurdle.md): family =
# "hurdle.lognormal" composes an occupancy probit fit over all n with a
# gaussian fit of log(y) over the y > 0 subset - no shared engine code, so
# this drives the existing probit and gaussian paths verbatim on a
# semicontinuous response. keepTrees so predict.bartHurdle can replay both
# saved forests onto the held-out rows. Recorded channels: the occupancy
# probability pi(x.test) (occ.test), the positive-part log-scale linear
# predictor f(x.test) (pos.test), and the combined natural-scale predict
# E[y | x.test] (the standard yhat.test/fhat.test channel below); varcount
# sums both forests' per-draw split counts (both are built from the same
# matched call, so their p columns and draw counts agree). Only type = "ppd"
# touches the RNG (unused here), so a given seed reproduces.
fitViaHurdle <- function(scenario) {
  fit <- bart2(
    scenario$x,
    scenario$y,
    family = "hurdle.lognormal",
    n.samples = ndpost,
    n.burn = nskip,
    n.trees = ntree,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    combineChains = TRUE,
    verbose = FALSE
  )
  list(
    yhat.test = predict(fit, scenario$x.test, type = "ev"),
    varcount = fit$occupancy$varcount + fit$positive$varcount,
    occ.test = predict(fit, scenario$x.test, type = "prob"),
    pos.test = predict(fit, scenario$x.test, type = "log")
  )
}

# runs bart2's gaussian path through its own FORMULA interface
# (bart2-argument-consolidation.md section 7): the four family fitters above
# all call bart2 on the x/y matrix convenience form with no data =, so the
# consolidated surface's formula parsing has no equivalence eyes at all. Two
# non-default args, both formula-interface-only and untouched by every
# fitter above: weights and subset (350 of 400 rows kept).
fitViaBart2Gauss <- function(scenario) {
  # model.frame() evaluates 'weights =' and 'subset =' in the FORMULA's own
  # environment (its default 'env', per ?model.frame.default), which is
  # wherever makeScenarios() created it - not this frame - so scenario$...
  # would resolve there and fail; rebase it to this call site.
  formula <- scenario$formula
  environment(formula) <- environment()
  withCallingHandlers(
    bart2(
      formula,
      data = scenario$data,
      test = scenario$x.test,
      weights = scenario$weights,
      subset = scenario$subsetIdx,
      n.samples = ndpost,
      n.burn = nskip,
      n.trees = ntree,
      n.chains = 1L,
      n.threads = 1L,
      combineChains = TRUE,
      verbose = FALSE
    ),
    warning = muffleBenignWarning
  )
}

# ditto, probit: offset (formula-interface-only, unreachable without data =)
# and TWO pooled chains (n.chains = 2L, combineChains = TRUE) - every family
# fitter above pins n.chains = 1L, so bart2's own chain pooling has no other
# anchor.
fitViaBart2Probit <- function(scenario) {
  # see fitViaBart2Gauss: rebase the formula's environment to this call site
  # so 'offset =' resolves scenario here rather than in makeScenarios().
  formula <- scenario$formula
  environment(formula) <- environment()
  withCallingHandlers(
    bart2(
      formula,
      data = scenario$data,
      test = scenario$x.test,
      offset = scenario$offset,
      offset.test = scenario$offsetTest,
      family = "probit",
      n.samples = ndpost,
      n.burn = nskip,
      n.trees = ntree,
      n.chains = 2L,
      n.threads = 1L,
      combineChains = TRUE,
      verbose = FALSE
    ),
    warning = muffleBenignWarning
  )
}

# an amplitude-coupled two-forest bart2 fit through the S12 term route
# (bart2-argument-consolidation.md section 5): the canonical `zf:forest(x1 +
# x2)` sugar declares a second forest, modulated by a 3-level factor's own
# indicator basis (no reference dropped), alongside the main x1 + x2 forest -
# the only path in this file to a K > 1 amplitude-coupled fit off bart2's
# formula interface. Such a fit refuses test = by name (section 5.4), so
# there is no held-out leg; instead this records the in-sample per-forest
# channels the fit exposes - forestFits (each forest's own response-scale raw
# total), glue (the ragged per-observation multiplier) and the
# forest-widened varcount - none of which any other scenario reaches. Sizes
# are this scenario's own literals, not the guarded ndpost/nskip/ntree.
fitViaBart2TwoForest <- function(scenario) {
  fit <- bart2(
    scenario$formula,
    data = scenario$data,
    n.samples = ndpost,
    n.burn = nskip,
    n.trees = ntree,
    n.chains = 1L,
    n.threads = 1L,
    combineChains = TRUE,
    verbose = FALSE
  )
  forestMeans <- apply(fit$forestFits, c(2L, 3L), mean)
  glueMeans <- as.vector(colMeans(fit$glue))
  vprop <- function(f) {
    colMeans(fit$varcount[,, f] / rowSums(fit$varcount[,, f]))
  }
  c(
    setNames(
      as.vector(forestMeans),
      paste0("forestfit.", seq_along(forestMeans))
    ),
    setNames(glueMeans, paste0("glue.", seq_along(glueMeans))),
    setNames(vprop(1L), paste0("vprop.forest1.", seq_len(ncol(fit$varcount)))),
    setNames(vprop(2L), paste0("vprop.forest2.", seq_len(ncol(fit$varcount)))),
    sigma.mean = mean(fit$sigma),
    sigma.sd = sd(fit$sigma)
  )
}

# runs bart2's multinomial path (docs/design/multinomial.md): family =
# "multinomial" explicitly - the only route to a K-forest softmax classifier
# ("auto" never reaches it, section 1.2). Unlike the two-forest fit above,
# multinomial DOES accept test =, so this records the posterior-mean
# K-category test probabilities - the standard fhat.test channel's
# K-widened shape - and each category forest's own per-draw split-usage
# proportions; there is no sigma and no k, the softmax having neither. The
# shape does not match the common yhat.test/varcount pair fitSummaries
# accumulates below, so this fitter returns its own finished summary
# vector, the xbart/twoforest precedent.
fitViaBart2Multinomial <- function(scenario) {
  fit <- bart2(
    scenario$x,
    scenario$y,
    test = scenario$x.test,
    family = "multinomial",
    n.samples = ndpost,
    n.burn = nskip,
    n.trees = ntree,
    n.chains = 1L,
    n.threads = 1L,
    combineChains = TRUE,
    verbose = FALSE
  )
  probMeans <- apply(fit$yhat.test, c(2L, 3L), mean)
  K <- dim(fit$yhat.test)[3L]
  vpropForest <- function(k) {
    vc <- fit$varcount[,, k]
    colMeans(vc / rowSums(vc))
  }
  c(
    setNames(
      as.vector(probMeans),
      paste0("prob.test.", seq_along(probMeans))
    ),
    unlist(lapply(seq_len(K), function(k) {
      setNames(
        vpropForest(k),
        paste0("vprop.forest", k, ".", seq_len(ncol(fit$varcount)))
      )
    }))
  )
}

# runs xbart's k-fold crossvalidation over a (n.trees x k) grid. Sizes are
# this scenario's own literals rather than the guarded ndpost/nskip/ntree,
# and n.threads = 1L keeps the whole sweep in this process, under the
# harness's own fork-level parallelism. The result is a loss array, not a
# fit, so this fitter returns the finished summary vector itself.
fitViaXbart <- function(scenario) {
  loss <- xbart(
    scenario$x,
    scenario$y,
    n.samples = 25L,
    n.burn = c(20L, 5L),
    method = "k-fold",
    n.test = 4,
    n.reps = 2L,
    n.trees = c(10L, 25L),
    k = c(1, 3),
    n.threads = 1L,
    verbose = FALSE
  )
  setNames(as.vector(loss), paste0("loss.", seq_along(loss)))
}

# runs rbart_vi's in-core grouped path (built-in tau prior, no callback):
# single chain, no thinning, the harness's global budget. Draw matrices come
# back sample-major except varcount (predictor x sample), transposed here so
# fitSummaries sees the bart() orientation.
fitViaRbart <- function(scenario) {
  fit <- if (isTRUE(scenario$aft)) {
    # grouped AFT survival (riAFTBART's model): the response enters as a Surv
    # on the formula's left-hand side - inlined (as the gaussian branch's
    # response is) so no intermediate reads as unused - and aft refuses
    # weights, so this scenario carries none; status rides scenario$status and
    # the observed time scenario$y. Still the in-core fast path (gamma prior).
    withCallingHandlers(
      rbart_vi(
        structure(
          cbind(time = scenario$y, status = scenario$status),
          class = "Surv",
          type = "right"
        ) ~
          scenario$x,
        test = scenario$x.test,
        group.by = scenario$group.by,
        group.by.test = scenario$group.by.test,
        family = "aft",
        prior = gamma,
        n.samples = ndpost,
        n.burn = nskip,
        n.thin = 1L,
        n.chains = 1L,
        n.threads = 1L,
        n.trees = ntree,
        keepTrees = FALSE,
        keepTestFits = TRUE,
        verbose = FALSE
      ),
      warning = muffleBenignWarning
    )
  } else {
    withCallingHandlers(
      rbart_vi(
        scenario$y ~ scenario$x,
        test = scenario$x.test,
        weights = scenario$weights,
        group.by = scenario$group.by,
        group.by.test = scenario$group.by.test,
        prior = gamma,
        n.samples = ndpost,
        n.burn = nskip,
        n.thin = 1L,
        n.chains = 1L,
        n.threads = 1L,
        n.trees = ntree,
        keepTrees = FALSE,
        keepTestFits = TRUE,
        verbose = FALSE
      ),
      warning = muffleBenignWarning
    )
  }
  list(
    yhat.test = fit$yhat.test,
    varcount = t(fit$varcount),
    sigma = as.vector(fit$sigma),
    tau = as.vector(fit$tau),
    ranef = fit$ranef
  )
}

fitSummaries <- function(scenario, seed) {
  set.seed(seed)
  # Test-data weights are irrelevant here (no posterior-predictive use);
  # muffle only that warning so real ones stay visible.
  splitprobs <- scenario$splitprobs
  # xbart reports a loss array rather than a fit, and carries none of the
  # channels assembled below, so it returns its own summary vector whole
  if (!is.null(scenario$xbartFit)) {
    return(fitViaXbart(scenario))
  }
  # a two-forest fit refuses 'test =' (section 5.4), so it carries no
  # yhat.test/varcount pair in the common shape below; it returns its own
  # summary vector whole, the xbart precedent.
  if (!is.null(scenario$twoforestFit)) {
    return(fitViaBart2TwoForest(scenario))
  }
  # a multinomial fit's yhat.test carries a K margin the common
  # yhat.test/varcount pair below has no shape for; it returns its own
  # summary vector whole, the twoforest precedent just above.
  if (!is.null(scenario$multinomialFit)) {
    return(fitViaBart2Multinomial(scenario))
  }
  fit <- if (!is.null(scenario$rbart)) {
    fitViaRbart(scenario)
  } else if (!is.null(scenario$ordinalFit)) {
    fitViaOrdinal(scenario)
  } else if (!is.null(scenario$nbinomFit)) {
    fitViaNbinom(scenario)
  } else if (!is.null(scenario$hazardFit)) {
    fitViaHazard(scenario)
  } else if (!is.null(scenario$hurdleFit)) {
    fitViaHurdle(scenario)
  } else if (!is.null(scenario$bart2GaussFit)) {
    fitViaBart2Gauss(scenario)
  } else if (!is.null(scenario$bart2ProbitFit)) {
    fitViaBart2Probit(scenario)
  } else if (!is.null(scenario$samplerApi)) {
    fitViaSamplerApi(scenario)
  } else if (!is.null(splitprobs)) {
    withCallingHandlers(
      bart(
        scenario$x,
        scenario$y,
        x.test = scenario$x.test,
        weights = scenario$weights,
        splitprobs = splitprobs,
        usequants = isTRUE(scenario$usequants),
        ntree = ntree,
        ndpost = ndpost,
        nskip = nskip,
        nchain = 1L,
        nthread = 1L,
        verbose = FALSE
      ),
      warning = muffleBenignWarning
    )
  } else {
    withCallingHandlers(
      bart(
        scenario$x,
        scenario$y,
        x.test = scenario$x.test,
        weights = scenario$weights,
        usequants = isTRUE(scenario$usequants),
        ntree = ntree,
        ndpost = ndpost,
        nskip = nskip,
        nchain = 1L,
        nthread = 1L,
        verbose = FALSE
      ),
      warning = muffleBenignWarning
    )
  }
  # varcount width follows the fitted forest (one column per predictor, a factor
  # counting as one variable), which need not equal ncol(x) for data frames.
  vc <- fit$varcount
  result <- c(
    setNames(colMeans(fit$yhat.test), paste0("fhat.test.", seq_len(n.test))),
    setNames(colMeans(vc / rowSums(vc)), paste0("vprop.", seq_len(ncol(vc))))
  )
  if (!scenario$binary) {
    result <- c(result, sigma.mean = mean(fit$sigma), sigma.sd = sd(fit$sigma))
  }
  if (!is.null(fit$k)) {
    result <- c(result, k.mean = mean(fit$k), k.sd = sd(fit$k))
  }
  if (!is.null(fit$tau)) {
    result <- c(result, tau.mean = mean(fit$tau), tau.sd = sd(fit$tau))
  }
  if (!is.null(fit$ranef)) {
    result <- c(
      result,
      setNames(colMeans(fit$ranef), paste0("ranef.", seq_len(ncol(fit$ranef))))
    )
  }
  # ordinal-only channels (fitViaOrdinal); NULL - and so absent - for every
  # other fitter, leaving the existing scenarios' summary vectors untouched.
  # cutpoint.*.1 is the pinned gamma_1 = 0 (a constant across seeds: its Welch
  # z is NaN and ignored, but it anchors the identical-draws comparison).
  if (!is.null(fit[["cutpoints"]])) {
    cp <- as.matrix(fit[["cutpoints"]])
    result <- c(
      result,
      setNames(colMeans(cp), paste0("cutpoint.mean.", seq_len(ncol(cp)))),
      setNames(apply(cp, 2L, sd), paste0("cutpoint.sd.", seq_len(ncol(cp))))
    )
  }
  if (!is.null(fit[["probs.test"]])) {
    probMeans <- apply(fit[["probs.test"]], c(2L, 3L), mean)
    result <- c(
      result,
      setNames(
        as.vector(probMeans),
        paste0("prob.test.", seq_len(length(probMeans)))
      )
    )
  }
  # nbinom-only channels (fitViaNbinom); NULL - and so absent - for every other
  # fitter, leaving the existing scenarios' summary vectors untouched. The
  # per-draw dispersion r and the posterior-mean test mean counts.
  if (!is.null(fit[["dispersion"]])) {
    d <- as.vector(fit[["dispersion"]])
    result <- c(result, dispersion.mean = mean(d), dispersion.sd = sd(d))
  }
  if (!is.null(fit[["means.test"]])) {
    result <- c(
      result,
      setNames(
        colMeans(fit[["means.test"]]),
        paste0("mean.test.", seq_len(n.test))
      )
    )
  }
  # hazard-only channel (fitViaHazard); NULL - and so absent - for every other
  # fitter, leaving the existing scenarios' summary vectors untouched. The
  # posterior-mean survival surface (times x test subjects), flattened.
  if (!is.null(fit[["surv.test"]])) {
    st <- fit[["surv.test"]]
    result <- c(
      result,
      setNames(as.vector(st), paste0("surv.test.", seq_along(st)))
    )
  }
  # hurdle-only channels (fitViaHurdle); NULL - and so absent - for every
  # other fitter, leaving the existing scenarios' summary vectors untouched.
  # The occupancy probability pi(x.test) and the positive-part log-scale
  # linear predictor f(x.test); the combined predict rides the standard
  # yhat.test/fhat.test channel above.
  if (!is.null(fit[["occ.test"]])) {
    result <- c(
      result,
      setNames(
        colMeans(fit[["occ.test"]]),
        paste0("occ.test.", seq_len(n.test))
      )
    )
  }
  if (!is.null(fit[["pos.test"]])) {
    result <- c(
      result,
      setNames(
        colMeans(fit[["pos.test"]]),
        paste0("pos.test.", seq_len(n.test))
      )
    )
  }
  # heteroscedastic-only channel (fitViaSamplerApi under variance = TRUE);
  # NULL - and so absent - for every other sampler. The posterior-mean variance
  # surface on the test rows: the direct read on the variance forest's routing,
  # which the mean-side channels see only transitively (s^2 enters the mean
  # forest's likelihood as a per-observation weight).
  if (!is.null(fit[["variance.test"]])) {
    result <- c(
      result,
      setNames(
        colMeans(fit[["variance.test"]]),
        paste0("s2.test.", seq_len(n.test))
      )
    )
  }
  # the engine's verdict on a transactional mutation, for a scenario that asked
  # to record it (recordVerdict): the accept/rollback flag as one entry, or the
  # per-observation session's whole install mask, one entry per row. NULL - and
  # so absent - everywhere else, including the predictor-mutation scenarios
  # recorded before this channel existed.
  if (!is.null(fit[["verdict"]])) {
    v <- fit[["verdict"]]
    result <- c(result, setNames(v, paste0("verdict.", seq_along(v))))
  }
  # the in-sample fit, mixedmatrix only. Gated on recordTrain itself, not
  # merely on fit$yhat.train's presence: a raw bart2()/bart() fit object
  # (fitViaBart2Gauss/Probit, or the plain bart() branches) always carries
  # yhat.train under the default keepTrainingFits, so a presence-only guard
  # would have added this channel to scenarios that never asked for it.
  if (isTRUE(scenario$recordTrain) && !is.null(fit[["yhat.train"]])) {
    result <- c(
      result,
      setNames(
        colMeans(fit[["yhat.train"]]),
        paste0("fhat.train.", seq_len(ncol(fit[["yhat.train"]])))
      )
    )
  }
  # the pre-replacement test fits, testSwap only: the test store the sampler
  # was CREATED with, and the one setTestPredictor installs. The standard
  # fhat.test channel above carries the last leg, so all three containers are
  # recorded and a build that rebuilt one of them differently moves only that
  # leg. NULL - and so absent - everywhere else.
  if (!is.null(fit[["yhat.test.pre"]])) {
    result <- c(
      result,
      setNames(
        colMeans(fit[["yhat.test.pre"]]),
        paste0("fhat.test.pre.", seq_len(n.test))
      )
    )
  }
  if (!is.null(fit[["yhat.test.mid"]])) {
    result <- c(
      result,
      setNames(
        colMeans(fit[["yhat.test.mid"]]),
        paste0("fhat.test.mid.", seq_len(n.test))
      )
    )
  }
  result
}

# The scenario x seed units are independent and fully self-seeded
# (fitSummaries calls set.seed first), so fork-level parallelism reproduces
# the serial results bitwise. EQUIVALENCE_CORES overrides the width; Windows
# has no fork and runs serially.
numCores <- local({
  env <- Sys.getenv("EQUIVALENCE_CORES", "")
  if (nzchar(env)) {
    max(1L, as.integer(env))
  } else if (.Platform$OS.type == "windows") {
    1L
  } else {
    parallel::detectCores()
  }
})

runAll <- function(scenarios) {
  grid <- expand.grid(
    seed = seq_len(n.seeds),
    scenario = names(scenarios),
    stringsAsFactors = FALSE
  )
  fitOne <- function(i) {
    fitSummaries(scenarios[[grid$scenario[i]]], grid$seed[i])
  }
  rows <- if (numCores > 1L) {
    parallel::mclapply(
      seq_len(nrow(grid)),
      fitOne,
      mc.cores = numCores,
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_len(nrow(grid)), fitOne)
  }
  failed <- vapply(rows, inherits, TRUE, "try-error")
  if (any(failed)) {
    stop(
      "fit failed for ",
      paste(grid$scenario[failed], grid$seed[failed], collapse = ", "),
      ": ",
      conditionMessage(attr(rows[[which(failed)[1L]]], "condition"))
    )
  }
  lapply(setNames(names(scenarios), names(scenarios)), function(name) {
    do.call(rbind, rows[grid$scenario == name])
  })
}

scenarios <- makeScenarios()

scenarioFilter <- Sys.getenv("EQUIVALENCE_SCENARIOS", "")
if (nzchar(scenarioFilter)) {
  wanted <- trimws(strsplit(scenarioFilter, ",", fixed = TRUE)[[1L]])
  unknown <- setdiff(wanted, names(scenarios))
  if (length(unknown) > 0L) {
    stop("unknown scenario(s): ", paste(unknown, collapse = ", "))
  }
  scenarios <- scenarios[wanted]
}

if (mode == "record") {
  out.file <- if (length(args) >= 2L) args[[2L]] else "equivalence-baseline.rds"
  results <- runAll(scenarios)
  meta <- list(
    rev = system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
    date = format(Sys.Date()),
    quick = quick,
    n.seeds = n.seeds,
    ndpost = ndpost,
    nskip = nskip,
    ntree = ntree
  )
  saveRDS(list(meta = meta, results = results), out.file)
  cat(
    "wrote baseline for",
    length(results),
    "scenarios x",
    n.seeds,
    "seeds to",
    out.file,
    "\n"
  )
} else if (mode == "compare") {
  if (length(args) < 2L) {
    stop("usage: equivalence.R compare baseline.rds")
  }
  baseline <- readRDS(args[[2L]])
  settings <- baseline$meta[c("quick", "n.seeds", "ndpost", "nskip", "ntree")]
  if (
    !identical(
      settings,
      list(
        quick = quick,
        n.seeds = n.seeds,
        ndpost = ndpost,
        nskip = nskip,
        ntree = ntree
      )
    )
  ) {
    stop(
      "baseline was recorded with different settings: ",
      paste(names(settings), unlist(settings), sep = "=", collapse = ", ")
    )
  }

  results <- runAll(scenarios)
  anyFailure <- FALSE
  compared <- character(0L)
  for (name in names(baseline$results)) {
    a <- baseline$results[[name]]
    b <- results[[name]]
    if (is.null(b)) {
      # scenario absent this run (e.g. sparse when Matrix is not installed)
      cat(sprintf("%-10s skipped (not produced this run)\n", name))
      next
    }
    compared <- c(compared, name)
    if (identical(a, b)) {
      cat(sprintf("%-10s identical draws (same RNG stream)\n", name))
      next
    }
    z <- (colMeans(a) - colMeans(b)) /
      sqrt(apply(a, 2L, var) / nrow(a) + apply(b, 2L, var) / nrow(b))
    # The Welch z is DEGENERATE when a run diverges: one blown-up seed
    # inflates the variance so much that |z| caps near sqrt(n.seeds)
    # (observed while poisoning the chi-k shape: 19/20 seeds cleanly
    # shifted plus one at 1e153 passed the z gate). Disjoint seed ranges
    # are an exact, parameter-free complement: under exchangeability
    # P(disjoint) = 2 / choose(2n, n) per summary, negligible at the full
    # 20 seeds; skipped for small-n quick runs where it is uninformative.
    disjoint <- rep(FALSE, ncol(a))
    if (2 / choose(nrow(a) + nrow(b), nrow(a)) < 1e-8) {
      disjoint <- vapply(
        seq_len(ncol(a)),
        function(j) {
          isTRUE(max(a[, j]) < min(b[, j])) ||
            isTRUE(max(b[, j]) < min(a[, j]))
        },
        TRUE
      )
    }
    n.warn <- sum(abs(z) > 3, na.rm = TRUE)
    n.fail <- sum(abs(z) > 4, na.rm = TRUE)
    n.disjoint <- sum(disjoint)
    cat(sprintf(
      "%-10s %d summaries, max |z| = %.2f%s%s%s\n",
      name,
      length(z),
      max(abs(z), na.rm = TRUE),
      if (n.warn > 0L) sprintf(", %d with |z| > 3", n.warn) else "",
      if (n.fail > 0L) sprintf(", %d with |z| > 4 <- FAIL", n.fail) else "",
      if (n.disjoint > 0L) {
        sprintf(", %d with disjoint seed ranges <- FAIL", n.disjoint)
      } else {
        ""
      }
    ))
    if (n.fail > 0L) {
      anyFailure <- TRUE
      failed <- which(abs(z) > 4)
      cat(
        "  worst offenders:",
        paste0(
          names(z)[failed],
          " (z=",
          round(z[failed], 2L),
          ")",
          collapse = ", "
        ),
        "\n"
      )
    }
    if (n.disjoint > 0L) {
      anyFailure <- TRUE
      cat(
        "  disjoint ranges:",
        paste0(
          names(z)[disjoint],
          " (",
          signif(colMeans(a)[disjoint], 3L),
          " vs ",
          signif(colMeans(b)[disjoint], 3L),
          ")",
          collapse = ", "
        ),
        "\n"
      )
    }
  }
  # scenarios the installed engine produces but the baseline predates (added
  # after it was recorded); the loop above only ever walks baseline names, so
  # these would otherwise go uncompared with no trace.
  uncovered <- setdiff(names(results), names(baseline$results))
  cat(sprintf(
    "\ncoverage: %d compared / %d skipped (not in baseline)%s\n",
    length(compared),
    length(uncovered),
    if (length(uncovered) > 0L) {
      paste0(": ", paste(uncovered, collapse = ", "))
    } else {
      ""
    }
  ))
  if (length(uncovered) > 0L) {
    msg <- paste0(
      "baseline predates these scenarios, so they were not compared: ",
      paste(uncovered, collapse = ", ")
    )
    if (strictCoverage) {
      anyFailure <- TRUE
      cat(msg, "\n")
    } else {
      warning(msg, call. = FALSE)
    }
  }

  if (anyFailure) {
    quit(status = 1L)
  }
  cat("\nOK: posteriors statistically indistinguishable at |z| > 4\n")
} else {
  stop("unknown mode '", mode, "'; use record or compare")
}
