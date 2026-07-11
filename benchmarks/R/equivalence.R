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
    if (!is.null(scenario$mutate$offset.test)) {
      sampler$setTestOffset(scenario$mutate$offset.test)
    }
    sampler$run(ceiling(nskip / 4), ndpost)
  } else {
    sampler$run(nskip, ndpost)
  }
  list(
    yhat.test = poolChains(r$test),
    varcount = poolChains(r$varcount),
    sigma = as.vector(r$sigma),
    k = if (!is.null(r$k)) as.vector(r$k)
  )
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
  fit <- if (!is.null(scenario$rbart)) {
    fitViaRbart(scenario)
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
