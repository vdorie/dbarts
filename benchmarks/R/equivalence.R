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
useNewEngine <- "engine=new" %in% args
args <- setdiff(args, "engine=new")
strictCoverage <- "--strict-coverage" %in% args
args <- setdiff(args, "--strict-coverage")
mode <- if (length(args) >= 1L) args[[1L]] else "record"

if (useNewEngine) {
  source(file.path(
    dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
    "bartcore-shim.R"
  ))
  loadBartcoreShim()
}

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
fitViaSamplerApi <- function(scenario, engineIsNew) {
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
  r <- if (!is.null(scenario$setData)) {
    sampler$run(nskip, 0L)
    sampler$setData(dbartsData(
      scenario$setData$x,
      scenario$setData$y,
      test = scenario$x.test
    ))
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
    verdict = if (isTRUE(scenario$recordVerdict)) as.double(verdict)
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
    cutpoints = fit$cutpoints,
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
  } else if (!is.null(scenario$samplerApi)) {
    fitViaSamplerApi(scenario, useNewEngine)
  } else if (useNewEngine) {
    bartcore_bart(
      scenario$x,
      scenario$y,
      x.test = scenario$x.test,
      weights = scenario$weights,
      splitprobs = splitprobs,
      usequants = isTRUE(scenario$usequants),
      ntree = ntree,
      ndpost = ndpost,
      nskip = nskip
    )
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
