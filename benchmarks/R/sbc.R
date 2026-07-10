# Simulation-based calibration (Talts, Betancourt, Simpson, Vehtari, Gelman
# 2018) for the shipped dbarts sampler. Drives the INSTALLED package through
# its R API only (no engine hooks): for a model configuration it draws theta0
# from the sampler's OWN prior, simulates y | theta0 through the assumed
# likelihood, refits with L near-independent retained draws, and ranks each
# scalar functional's theta0 among its L posterior draws. Over R replications a
# calibrated sampler yields uniform ranks; the rank histogram, an ecdf-diff
# simultaneous band, and a chi-square test flag any non-uniformity.
#
# Run against the installed package (R CMD INSTALL . first):
#   Rscript benchmarks/R/sbc.R                 # baseline gaussian, R=200
#   Rscript benchmarks/R/sbc.R gaussian 200 200 30
#   Rscript benchmarks/R/sbc.R probit  200 200 30
# Positional args: config R L thin. Or source() the file to reuse the API:
#   source("benchmarks/R/sbc.R"); res <- runSbc(sbcConfig("gaussian"), R = 200)
#
# SELF-CONSISTENCY is the whole game: theta0 must come from the same prior the
# sampler assumes in its posterior. The forest/leaf draw uses the sampler's own
# sampleTreesFromPrior + sampleNodeParametersFromPrior; the sigma draw is the
# reported-scale scaled-inverse-chi-squared the engine calibrates (verified by
# a moment check, sbcCheckSigmaPrior). The data scale is fixed once at build
# (setResponse with updateScale = FALSE keeps it) so prior and posterior share
# it. A wrong prior draw makes SBC lie, so the harness self-checks before use.

suppressMessages(library(dbarts))

# --- sigma prior -----------------------------------------------------------

# The engine's residual-variance prior, on the REPORTED (original) scale, is
# scaled-inverse-chi-squared: sigma^2 ~ df * sigest^2 * rawScale / chisq(df),
# with rawScale = qchisq(1 - quant, df) / df (R_interface_bartcore.cpp: the
# internal scale = (sigest / range)^2 * rawScale, and range^2 cancels once the
# reported sigma = internal * range, so the reported-scale prior is
# range-independent). This calibration puts P(sigma < sigest) = quant exactly.
sbcSigmaDraw <- function(sigest, df, quant) {
  rawScale <- qchisq(1 - quant, df) / df
  function(nDraws = 1L) {
    sqrt(df * sigest^2 * rawScale / rchisq(nDraws, df))
  }
}

# Moment/calibration check the spec requires before trusting the sigma prior:
# a correct scaled-inv-chisq draw has P(sigma < sigest) = quant and a known
# median. Returns the empirical coverage and a pass flag.
sbcCheckSigmaPrior <- function(sigest, df, quant, nDraws = 2e5L) {
  draws <- sbcSigmaDraw(sigest, df, quant)(nDraws)
  coverage <- mean(draws < sigest)
  medianTheory <- sqrt(
    df * sigest^2 * (qchisq(1 - quant, df) / df) / qchisq(0.5, df)
  )
  list(
    coverage = coverage,
    coverageTarget = quant,
    medianEmpirical = median(draws),
    medianTheory = medianTheory,
    pass = abs(coverage - quant) < 0.005 &&
      abs(median(draws) / medianTheory - 1) < 0.02
  )
}

# --- likelihood ------------------------------------------------------------

# Simulate y0 | theta0 through the family's assumed likelihood. The latent
# f0Train is on the reported scale predict() returns; binary families threshold
# it through their link. Weighted gaussian scales the noise by 1/sqrt(weight)
# (the engine's per-row precision), zero-weight rows carry pure prior noise the
# fit ignores. Grouped models add the drawn per-row group intercept before the
# family draw.
sbcSimulate <- function(config, f0Train, sig0) {
  mu <- f0Train
  if (!is.null(config$groupIntercepts)) {
    mu <- mu + config$groupIntercepts[config$group]
  }
  switch(
    config$family,
    gaussian = {
      # per-row noise sd = sigma / sqrt(weight); zero-weight rows are dropped by
      # the likelihood, so their simulated value is a finite placeholder (w = 1)
      sd <- if (is.null(config$weights)) {
        rep_len(sig0, config$n)
      } else {
        sig0 / sqrt(ifelse(config$weights > 0, config$weights, 1))
      }
      mu + sd * rnorm(config$n)
    },
    probit = as.double(rbinom(config$n, 1L, pnorm(mu))),
    logistic = as.double(rbinom(config$n, 1L, plogis(mu)))
  )
}

# --- configuration ---------------------------------------------------------

# A configuration bundles everything a run needs: the fixed design, the model
# priors, and family-specific prior draw / likelihood / functional logic. New
# configurations (linear/gp leaf, grouped, BCF, DART, weighted) extend the same
# shape; this tier exercises gaussian, probit, logistic, DART, weighted.
sbcConfig <- function(
  family = c("gaussian", "probit", "logistic"),
  n = 150L,
  p = 3L,
  nTrees = 50L,
  k = 2,
  sigDf = 3,
  sigQuant = 0.9,
  nTest = 5L,
  nodePrior = NULL,
  dartAlpha = 1.0,
  weights = NULL,
  configSeed = 1L
) {
  family <- match.arg(family)
  set.seed(configSeed)
  x <- matrix(runif(n * p), n, p)
  colnames(x) <- paste0("x", seq_len(p))
  xTest <- matrix(runif(nTest * p), nTest, p)
  colnames(xTest) <- colnames(x)
  # A deterministic build response fixes the internal data scale (range,
  # centre) once; setResponse(updateScale = FALSE) then keeps it across
  # replications so the prior draw and the posterior share one scale. Binary
  # families need a 0/1 build vector; the probit latent scale is fixed by the
  # link, so no continuous range is involved.
  yBuild <- if (family == "gaussian") {
    seq(-2.5, 2.5, length.out = n)
  } else {
    as.double(rep_len(c(0L, 1L), n))
  }
  # Prior constructors are not bare-exported; they live in dbartsPriors (or in
  # the special evaluation env of dbarts()'s prior arguments). Build objects
  # ahead of time so the config carries a concrete prior.
  if (is.null(nodePrior)) {
    nodePrior <- dbartsPriors$normal(k)
  }
  # sigest anchors the sigma prior for gaussian; binary families fix sigma = 1.
  sigest <- if (family == "gaussian") 1.0 else 1.0
  list(
    family = family,
    n = n,
    p = p,
    nTrees = nTrees,
    k = k,
    sigDf = sigDf,
    sigQuant = sigQuant,
    nTest = nTest,
    x = x,
    xTest = xTest,
    yBuild = yBuild,
    sigest = sigest,
    nodePrior = nodePrior,
    dartAlpha = dartAlpha,
    weights = weights,
    hasSigma = family == "gaussian"
  )
}

# Add grouped random intercepts to a base config: a fixed grouping assigned
# independently of x, a fixed tau relative scale, and optionally a fraction of
# zero-weight rows (dropped by the likelihood -- a review-3 self-consistency
# target). Returns the extended config for runSbcGrouped.
sbcAddGrouping <- function(
  config,
  nGroups = 8L,
  relScale = 0.5,
  zeroWeightFrac = 0
) {
  set.seed(101L)
  config$nGroups <- nGroups
  config$group <- as.integer(sample.int(nGroups, config$n, replace = TRUE))
  config$relScale <- relScale
  if (zeroWeightFrac > 0) {
    w <- rep_len(1, config$n)
    w[sample.int(config$n, floor(zeroWeightFrac * config$n))] <- 0
    config$weights <- w
  }
  config
}

# Inject NA values into designated columns of the fixed design (missing =
# "incorporate", the default, handles them: rules carry a missing direction and
# missing leaf covariates enter at the standardized mean, model.hpp:173). The
# NA pattern is fixed across replications; watchRows records which rows carry
# an NA so per-row functionals can sit exactly on the imputation path.
sbcAddMissing <- function(config, columns, frac = 0.15) {
  set.seed(303L)
  naRows <- sort(sample.int(config$n, ceiling(frac * config$n)))
  for (j in columns) {
    config$x[naRows, j] <- NA_real_
  }
  config$watchRows <- naRows[seq_len(min(3L, length(naRows)))]
  config$naColumns <- columns
  config
}

# Add BCF glue to a base config: a fixed 0/1 treatment (assigned by a propensity
# in x1 so mu(x, pihat) has something to condition on) and the glue prior
# scales. The prognostic scalar prior is Cauchy(0, sd.control).
sbcAddBCF <- function(
  config,
  sdControl = 2,
  sdModerate = 1,
  bPriorVariance = 0.5
) {
  set.seed(202L)
  pi <- pnorm(0.8 * (config$x[, 1L] - 0.5))
  config$z <- as.double(rbinom(config$n, 1L, pi))
  config$sdControl <- sdControl
  config$sdModerate <- sdModerate
  config$bPriorVariance <- bPriorVariance
  config
}

# Build the reusable sampler for a configuration. One sampler serves all
# replications: the prior draw advances its internal RNG, setResponse swaps y
# in place, and the fixed build scale is never disturbed.
sbcMakeSampler <- function(config, L, thin, seed) {
  ctrl <- dbartsControl(
    n.trees = config$nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = L,
    n.thin = thin,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  family <- if (config$family == "gaussian") "gaussian" else config$family
  # The matrix (xy) interface DROPS NA rows even under missing = "incorporate"
  # (dbartsData warns "row(s) dropped"); only the formula interface keeps them
  # (na.action = na.pass). NA designs therefore build through a formula.
  if (anyNA(config$x)) {
    df <- as.data.frame(config$x)
    df$.sbc.y <- config$yBuild
    args <- list(
      formula = as.formula(paste(
        ".sbc.y ~",
        paste(colnames(config$x), collapse = " + ")
      )),
      data = df,
      test = as.data.frame(config$xTest),
      resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
      node.prior = config$nodePrior,
      sigma = config$sigest,
      control = ctrl,
      family = family,
      missing = "incorporate"
    )
  } else {
    args <- list(
      config$x,
      config$yBuild,
      test = config$xTest,
      resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
      node.prior = config$nodePrior,
      sigma = config$sigest,
      control = ctrl,
      family = family
    )
  }
  if (!is.null(config$weights)) {
    args$weights <- config$weights
  }
  # weights-on-test-data warns benignly (test predictions stay unweighted)
  suppressWarnings(do.call(dbarts, args))
}

# --- one replication -------------------------------------------------------

# The engine's INTERNAL-scale total fits at the current state (the stored leaf
# values the likelihood actually uses), and the affine internal -> reported
# map recovered from one recorded sweep (exact; cf. the BCF map). For GP
# (function-valued) leaves predict() re-krigs the stored values with jitter and
# differs from the recorded training fits by ~2e-3, so theta0's f0Train must
# read the stored fits instead; the TEST path is exactly shared with the
# in-run recorded test fits, so f* keeps predict().
sbcInternalFits <- function(sampler) {
  getFits <- getFromNamespace("C_dbarts_bartcore_getForestFits", "dbarts")
  .Call(getFits, sampler$getPointer(), 0L)[, 1]
}

sbcRecoverFitMap <- function(sampler, config) {
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  y0 <- sbcSimulate(config, as.numeric(sampler$predict(config$x)), 1.0)
  sampler$setResponse(y0)
  res <- sampler$run(0L, 1L)
  df <- data.frame(
    reported = res$train[, 1],
    internal = sbcInternalFits(sampler)
  )
  fit <- lm(reported ~ internal, data = df)
  list(
    shift = unname(coef(fit)[1L]),
    scale = unname(coef(fit)[2L]),
    maxResid = max(abs(residuals(fit)))
  )
}

# Draw theta0 from the prior, simulate y0, refit, and rank each functional's
# theta0 among its L posterior draws. The rank is #{posterior_l < theta0} in
# {0, ..., L}, uniform under calibration. The MCMC is re-initialised from a
# FRESH independent prior draw (not the truth) so a finite burn-in cannot leave
# the chain parked at theta0 and bias ranks toward the centre. fitMap non-NULL
# switches f0Train to the stored-fits path (GP leaves).
sbcReplication <- function(sampler, config, drawSigma, L, burn, fitMap = NULL) {
  # 1. theta0 from the prior (forest + leaves via the engine's own machinery)
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  f0Train <- if (is.null(fitMap)) {
    as.numeric(sampler$predict(config$x))
  } else {
    fitMap$shift + fitMap$scale * sbcInternalFits(sampler)
  }
  f0Test <- as.numeric(sampler$predict(config$xTest))
  sig0 <- if (config$hasSigma) drawSigma(1L) else 1.0
  avgF0 <- mean(f0Train)

  # 2. simulate y0 through the assumed likelihood
  y0 <- sbcSimulate(config, f0Train, sig0)

  # 3. overdispersed init: a second, independent prior draw, then refit
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  if (config$hasSigma) {
    sampler$setSigma(config$sigest)
  }
  sampler$setResponse(y0)
  res <- sampler$run(burn, L)

  # 4. rank each functional
  ranks <- c(
    avg.f = sum(colMeans(res$train) < avgF0)
  )
  for (j in seq_len(config$nTest)) {
    ranks[paste0("f.star", j)] <- sum(res$test[j, ] < f0Test[j])
  }
  if (config$hasSigma) {
    ranks["sigma"] <- sum(as.numeric(res$sigma) < sig0)
  }
  # watch rows: f at designated TRAINING rows, ranked through the recorded
  # training fits -- used to point functionals at NA-bearing rows so the
  # missing-covariate path is checked exactly where it acts
  for (w in seq_along(config$watchRows)) {
    i <- config$watchRows[w]
    ranks[paste0("f.row", i)] <- sum(res$train[i, ] < f0Train[i])
  }
  ranks
}

# Self-consistency of the prior draw's f0 with the fit's likelihood: the
# recorded training fits and predict(x) must agree at the SAME sampler state
# (both route NA by the rule's missing direction and impute missing leaf
# covariates at the standardized mean, so any disagreement means the SBC
# theta0 is not the f the likelihood uses). The TEST path matters equally:
# theta0's f(x*) comes from predict(xTest) while its posterior draws come from
# the in-run recorded test fits, so those two maps must agree at one state too
# (function-valued GP leaves ride a separate prediction path). One sweep, then
# compare both.
sbcCheckFitConsistency <- function(config, seed = 99L) {
  set.seed(seed)
  sampler <- sbcMakeSampler(config, 1L, 1L, seed)
  fitMap <- if (isTRUE(config$f0FromForestFits)) {
    sbcRecoverFitMap(sampler, config)
  } else {
    NULL
  }
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  f0 <- as.numeric(sampler$predict(config$x))
  y0 <- sbcSimulate(config, f0, 1.0)
  sampler$setResponse(y0)
  res <- sampler$run(0L, 1L)
  trainRef <- if (is.null(fitMap)) {
    as.numeric(sampler$predict(config$x))
  } else {
    fitMap$shift + fitMap$scale * sbcInternalFits(sampler)
  }
  predTest <- as.numeric(sampler$predict(config$xTest))
  maxDiff <- max(abs(res$train[, 1] - trainRef))
  maxDiffTest <- max(abs(res$test[, 1] - predTest))
  list(
    maxDiff = maxDiff,
    maxDiffTest = maxDiffTest,
    pass = maxDiff < 1e-8 && maxDiffTest < 1e-8
  )
}

# --- driver ----------------------------------------------------------------

# Run R replications for a configuration. Returns the R x numFunctional rank
# matrix, the chosen (L, thin, burn), and wall-clock timing. Progress prints
# every `report` replications (SBC is long wall-clock, not quiet-machine).
runSbc <- function(
  config,
  R = 200L,
  L = 200L,
  thin = 30L,
  burn = 50L,
  seed = 20260709L,
  report = 25L
) {
  set.seed(seed)
  sampler <- sbcMakeSampler(config, L, thin, seed)
  drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
  fitMap <- NULL
  if (isTRUE(config$f0FromForestFits)) {
    fitMap <- sbcRecoverFitMap(sampler, config)
    if (fitMap$maxResid > 1e-8) {
      stop("internal -> reported fit map is not exact: ", fitMap$maxResid)
    }
  }

  ranks <- NULL
  started <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    row <- sbcReplication(sampler, config, drawSigma, L, burn, fitMap)
    if (is.null(ranks)) {
      ranks <- matrix(
        NA_integer_,
        R,
        length(row),
        dimnames = list(NULL, names(row))
      )
    }
    ranks[r, ] <- row
    if (report > 0L && (r %% report == 0L || r == R)) {
      elapsed <- proc.time()[["elapsed"]] - started
      cat(sprintf(
        "  [%s] rep %d/%d  %.1fs elapsed  %.2fs/rep\n",
        config$family,
        r,
        R,
        elapsed,
        elapsed / r
      ))
    }
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    ranks = ranks,
    L = L,
    thin = thin,
    burn = burn,
    R = R,
    config = config,
    elapsed = elapsed,
    perRep = elapsed / R
  )
}

# --- DART variable-selection calibration -----------------------------------

# sampleTreesFromPrior grows trees under the CURRENT split probabilities
# (uniform at DART init), NOT a Dirichlet draw, so a self-consistent DART SBC
# uses two samplers: a generator that grows the forest under a fixed split
# vector s0 (a genuine Dirichlet-prior draw), and a DART fit whose posterior s
# should cover s0. The joint prior is s0 ~ Dirichlet(alpha/p), forest | s0 with
# splits ~ s0 -- exactly what the generator produces and the fit assumes.

# Draw split probs from the DART Dirichlet(alpha/p) prior, applying the same
# 1e-300 floor the engine's posterior update uses (model.hpp DartPrior::update)
# so the prior draw and the fit share the floor -- required for the sparsity
# probe to be self-consistent.
sbcDartDirichlet <- function(alpha, p) {
  raw <- rgamma(p, alpha / p, 1)
  g <- pmax(raw, 1e-300)
  s <- g / sum(g)
  attr(s, "nFloor") <- sum(raw <= 1e-300)
  s
}

# Moment check: Dirichlet(alpha/p) has E[s_j] = 1/p and
# Var[s_j] = (1/p)(1 - 1/p)/(alpha + 1). At small alpha the floor bites and the
# empirical variance drifts below theory -- reported, not failed, for the probe.
sbcCheckDirichlet <- function(alpha, p, nDraws = 2e4L) {
  draws <- matrix(NA_real_, nDraws, p)
  for (i in seq_len(nDraws)) {
    draws[i, ] <- sbcDartDirichlet(alpha, p)
  }
  meanEmp <- mean(colMeans(draws))
  varEmp <- mean(apply(draws, 2, var))
  varTheory <- (1 / p) * (1 - 1 / p) / (alpha + 1)
  list(
    meanEmp = meanEmp,
    meanTheory = 1 / p,
    varEmp = varEmp,
    varTheory = varTheory,
    pass = abs(meanEmp - 1 / p) < 0.01 &&
      abs(varEmp / varTheory - 1) < 0.15
  )
}

# Generator sampler: constant-leaf gaussian whose CGM tree prior splits under a
# fixed s0 (non-DART). One MCMC-free prior draw yields f0 with splits ~ s0.
sbcMakeDartGenerator <- function(config, s0) {
  ctrl <- dbartsControl(
    n.trees = config$nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 1L,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  dbarts(
    config$x,
    config$yBuild,
    test = config$xTest,
    resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
    node.prior = config$nodePrior,
    tree.prior = dbartsPriors$cgm(split.probs = s0),
    sigma = config$sigest,
    control = ctrl
  )
}

# DART fit sampler: fixed alpha (update.alpha = FALSE) isolates the s
# calibration; update.delay small so s equilibrates well before samples are
# kept. Rebuilt once and reused across replications via setResponse.
sbcMakeDartFit <- function(config, L, thin) {
  ctrl <- dbartsControl(
    n.trees = config$nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = L,
    n.thin = thin,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  dbarts(
    config$x,
    config$yBuild,
    test = config$xTest,
    resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
    node.prior = config$nodePrior,
    tree.prior = dbartsPriors$dart(
      alpha = config$dartAlpha,
      update.alpha = FALSE,
      update.delay = 100L
    ),
    sigma = config$sigest,
    control = ctrl
  )
}

# One DART replication: draw s0, generate f0 under s0, simulate, fit with DART,
# rank s0_j among posterior varprobs plus the usual sigma / avg f / f(x*).
runSbcDart <- function(
  config,
  R = 200L,
  L = 200L,
  thin = 30L,
  burn = 50L,
  seed = 20260709L,
  report = 25L
) {
  set.seed(seed)
  fit <- sbcMakeDartFit(config, L, thin)
  drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
  ranks <- NULL
  floorHits <- 0L # count of s0 components pinned at the 1e-300 floor
  started <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    s0 <- sbcDartDirichlet(config$dartAlpha, config$p)
    floorHits <- floorHits + attr(s0, "nFloor")
    gen <- sbcMakeDartGenerator(config, s0)
    gen$sampleTreesFromPrior()
    gen$sampleNodeParametersFromPrior()
    f0Train <- as.numeric(gen$predict(config$x))
    f0Test <- as.numeric(gen$predict(config$xTest))
    sig0 <- drawSigma(1L)
    avgF0 <- mean(f0Train)
    y0 <- f0Train + sig0 * rnorm(config$n)

    fit$sampleTreesFromPrior()
    fit$sampleNodeParametersFromPrior()
    fit$setSigma(config$sigest)
    fit$setResponse(y0)
    res <- fit$run(burn, L)

    row <- c(
      avg.f = sum(colMeans(res$train) < avgF0),
      sigma = sum(as.numeric(res$sigma) < sig0)
    )
    for (j in seq_len(config$nTest)) {
      row[paste0("f.star", j)] <- sum(res$test[j, ] < f0Test[j])
    }
    for (j in seq_len(config$p)) {
      row[paste0("s", j)] <- sum(res$varprobs[j, ] < s0[j])
    }
    if (is.null(ranks)) {
      ranks <- matrix(
        NA_integer_,
        R,
        length(row),
        dimnames = list(NULL, names(row))
      )
    }
    ranks[r, ] <- row
    if (report > 0L && (r %% report == 0L || r == R)) {
      elapsed <- proc.time()[["elapsed"]] - started
      cat(sprintf(
        "  [dart a=%.2f p=%d] rep %d/%d  %.1fs  %.2fs/rep\n",
        config$dartAlpha,
        config$p,
        r,
        R,
        elapsed,
        elapsed / r
      ))
    }
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    ranks = ranks,
    L = L,
    thin = thin,
    burn = burn,
    R = R,
    config = config,
    elapsed = elapsed,
    perRep = elapsed / R,
    floorFrac = floorHits / (R * config$p)
  )
}

# --- grouped random intercepts ---------------------------------------------

# Grouped-tau prior (rbart_vi's in-core path). The engine offers two (model.hpp
# logTauPrior): reported tau ~ half-Cauchy(0, 2.5*rel.scale) or Gamma(shape=2.5,
# scale=2.5*rel.scale); both are reported-scale and fitScale-independent (the
# internal fitScale cancels on report, like sigma), and group effects
# b_g ~ N(0, tau^2). SBC uses the GAMMA prior: the half-Cauchy's infinite-
# variance tail produces occasional astronomically large tau0 that inflate the
# response scale and stall the engine's tau slice sampler (stepping out by a
# fixed width over a range up to 1e6+), making brute-force SBC intractable. The
# gamma prior is light-tailed, engine-supported, and gives a clean well-posed
# calibration. A grouped sampler fixes its response at creation (the bridge
# refuses setResponse), so each replication rebuilds the fit. Groups are
# assigned independently of x so the smooth f and the categorical b_g stay
# orthogonal -- the f/b partition is then clean regardless of the f-prior scale.

sbcTauDraw <- function(relScale) {
  function(nDraws = 1L) rgamma(nDraws, shape = 2.5, scale = 2.5 * relScale)
}

# Moment check: Gamma(2.5, 2.5*rel.scale) has mean 6.25*rel.scale and
# variance 2.5*(2.5*rel.scale)^2.
sbcCheckTau <- function(relScale, nDraws = 2e5L) {
  d <- sbcTauDraw(relScale)(nDraws)
  meanTheory <- 2.5 * (2.5 * relScale)
  varTheory <- 2.5 * (2.5 * relScale)^2
  list(
    meanEmp = mean(d),
    meanTheory = meanTheory,
    varEmp = var(d),
    varTheory = varTheory,
    pass = abs(mean(d) / meanTheory - 1) < 0.02 &&
      abs(var(d) / varTheory - 1) < 0.05
  )
}

# Reusable f-only generator: a plain (ungrouped) sampler at the config's family
# and build scale. sampleTreesFromPrior + predict give a fresh prior f0.
sbcMakeGroupGenerator <- function(config) {
  ctrl <- dbartsControl(
    n.trees = config$nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 1L,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  family <- if (config$family == "gaussian") "gaussian" else config$family
  dbarts(
    config$x,
    config$yBuild,
    test = config$xTest,
    resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
    node.prior = config$nodePrior,
    sigma = config$sigest,
    control = ctrl,
    family = family
  )
}

# Grouped fit built fresh for one y0: the bartcore.groups attribute carries the
# fixed grouping and a fixed rel.scale (NOT recomputed from y0, for self-
# consistency). n.steps is the per-sweep tau slice count.
sbcMakeGroupedFit <- function(config, y0, L, thin) {
  ctrl <- dbartsControl(
    n.trees = config$nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = L,
    n.thin = thin,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  attr(ctrl, "bartcore.groups") <- list(
    indices = config$group,
    n.groups = config$nGroups,
    prior = "gamma",
    rel.scale = config$relScale,
    n.steps = 5L
  )
  family <- if (config$family == "gaussian") "gaussian" else config$family
  args <- list(
    config$x,
    y0,
    test = config$xTest,
    resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
    node.prior = config$nodePrior,
    sigma = config$sigest,
    control = ctrl,
    family = family
  )
  if (!is.null(config$weights)) {
    args$weights <- config$weights
  }
  # zero-weight rows and weights-on-test warn benignly for this construction
  suppressWarnings(do.call(dbarts, args))
}

runSbcGrouped <- function(
  config,
  R = 200L,
  L = 200L,
  thin = 30L,
  burn = 50L,
  seed = 20260709L,
  report = 25L
) {
  set.seed(seed)
  gen <- sbcMakeGroupGenerator(config)
  drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
  drawTau <- sbcTauDraw(config$relScale)
  ranks <- NULL
  started <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    gen$sampleTreesFromPrior()
    gen$sampleNodeParametersFromPrior()
    f0Train <- as.numeric(gen$predict(config$x))
    f0Test <- as.numeric(gen$predict(config$xTest))
    tau0 <- drawTau(1L)
    b0 <- rnorm(config$nGroups, 0, tau0)
    sig0 <- if (config$hasSigma) drawSigma(1L) else 1.0
    avgF0 <- mean(f0Train)

    cfgLocal <- config
    cfgLocal$groupIntercepts <- b0
    y0 <- sbcSimulate(cfgLocal, f0Train, sig0)

    fit <- sbcMakeGroupedFit(config, y0, L, thin)
    res <- fit$run(burn, L)

    row <- c(
      tau = sum(as.numeric(res$tau) < tau0),
      b1 = sum(res$ranef[1, ] < b0[1]),
      b2 = sum(res$ranef[2, ] < b0[2]),
      avg.f = sum(colMeans(res$train) < avgF0)
    )
    if (config$hasSigma) {
      row["sigma"] <- sum(as.numeric(res$sigma) < sig0)
    }
    for (j in seq_len(config$nTest)) {
      row[paste0("f.star", j)] <- sum(res$test[j, ] < f0Test[j])
    }
    if (is.null(ranks)) {
      ranks <- matrix(
        NA_integer_,
        R,
        length(row),
        dimnames = list(NULL, names(row))
      )
    }
    ranks[r, ] <- row
    if (report > 0L && (r %% report == 0L || r == R)) {
      elapsed <- proc.time()[["elapsed"]] - started
      cat(sprintf(
        "  [grouped-%s G=%d] rep %d/%d  %.1fs  %.2fs/rep\n",
        config$family,
        config$nGroups,
        r,
        R,
        elapsed,
        elapsed / r
      ))
    }
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    ranks = ranks,
    L = L,
    thin = thin,
    burn = burn,
    R = R,
    config = config,
    elapsed = elapsed,
    perRep = elapsed / R
  )
}

# --- BCF (Bayesian causal forest glue) -------------------------------------

# BCF (docs/design/bcf.md): y = a*mu(x,pihat) + b_{z}*tau(x) + eps. Prognostic
# scalar a ~ Cauchy(0, aPriorScale = sd.control) (a scale-mixture, chain.hpp
# drawGlue); treatment coefficients b0, b1 ~ N(0, bPriorVariance) so the effect
# is (b1 - b0)*tau(x). The a-glue prior precision was the one true gate
# survivor, so its calibration is the headline. BCF is an internal bartcore
# sampler (not the R5 surface); it has numGroups == 0 so setResponse works,
# letting one sampler serve all reps with a fixed scale (no rebuild mismatch).
# The glue and per-forest fits are only exposed as CURRENT state, so posterior
# draws are collected one sample at a time. Identification: a*mu and b_z*tau are
# each sign-invariant under (a, mu) -> (-a, -mu), so the CLEAN functionals are
# the identified functions a*mu(x*) and (b1-b0)*tau(x*) at fixed points plus
# sigma; the raw a and (b1-b0) are reported too but carry that sign caveat.

.bcfNew <- getFromNamespace("bartcoreBCFSampler", "dbarts")
.bcfRun <- getFromNamespace("bartcoreRun", "dbarts")
.bcfGlue <- getFromNamespace("bartcoreBCFGlue", "dbarts")
.bcfForest <- getFromNamespace("bartcoreForestFits", "dbarts")
.bcfSetResponse <- getFromNamespace("bartcoreSetResponse", "dbarts")
.bcfPriorTrees <- getFromNamespace(
  "C_dbarts_bartcore_sampleTreesFromPrior",
  "dbarts"
)
.bcfPriorNodes <- getFromNamespace(
  "C_dbarts_bartcore_sampleNodeParametersFromPrior",
  "dbarts"
)

sbcBCFGlueDraw <- function(aPriorScale, bPriorVariance) {
  function() {
    list(
      a = rcauchy(1L, 0, aPriorScale),
      b0 = rnorm(1L, 0, sqrt(bPriorVariance)),
      b1 = rnorm(1L, 0, sqrt(bPriorVariance))
    )
  }
}

# Moment check for the glue priors: a ~ Cauchy(0, s) (median 0, IQR 2s),
# b0/b1 ~ N(0, bVar) (sd sqrt(bVar)), b1 - b0 ~ N(0, 2 bVar).
sbcCheckBCFGlue <- function(aPriorScale, bPriorVariance, nDraws = 2e5L) {
  a <- rcauchy(nDraws, 0, aPriorScale)
  b0 <- rnorm(nDraws, 0, sqrt(bPriorVariance))
  b1 <- rnorm(nDraws, 0, sqrt(bPriorVariance))
  list(
    aIqrEmp = as.numeric(diff(quantile(a, c(0.25, 0.75)))),
    aIqrTheory = 2 * aPriorScale,
    bSdEmp = sd(b0),
    bSdTheory = sqrt(bPriorVariance),
    diffSdEmp = sd(b1 - b0),
    diffSdTheory = sqrt(2 * bPriorVariance),
    pass = abs(diff(quantile(a, c(0.25, 0.75))) / (2 * aPriorScale) - 1) <
      0.03 &&
      abs(sd(b0) / sqrt(bPriorVariance) - 1) < 0.02
  )
}

# Build a BCF sampler on a fixed design and recover the reported <- internal
# affine map (fitScale, fitShift) by regressing one run's reported combined
# fits on the internal a*mu + b_z*tau. Returns the sampler, z, and the map.
# fixedGlue = TRUE holds the glue at its initial values (a=1, b0=0, b1=1) via
# update.a = update.b = FALSE; the SBC generator then uses those constants,
# isolating the two-forest backfit from the glue draw for diagnosis.
sbcMakeBCF <- function(config, L, thin, fixedGlue = FALSE) {
  ctrl <- dbartsControl(
    n.trees = config$nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = L,
    n.thin = thin,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  base <- dbarts(
    config$x,
    config$yBuild,
    resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
    node.prior = config$nodePrior,
    sigma = config$sigest,
    control = ctrl
  )
  bcf <- .bcfNew(
    base,
    config$z,
    sd.control = config$sdControl,
    sd.moderate = config$sdModerate,
    b.prior.variance = config$bPriorVariance,
    update.a = !fixedGlue,
    update.b = !fixedGlue
  )
  # recover the affine reported/internal map from one warm sample
  res <- .bcfRun(bcf, 0L, 1L)
  glue <- .bcfGlue(bcf)
  mu <- .bcfForest(bcf, 0L)[, 1]
  tau <- .bcfForest(bcf, 1L)[, 1]
  bz <- ifelse(config$z != 0, glue[3L], glue[2L])
  map <- data.frame(
    reported = res$train[, 1],
    internal = glue[1L] * mu + bz * tau
  )
  fit <- lm(reported ~ internal, data = map)
  list(
    bcf = bcf,
    fitShift = unname(coef(fit)[1L]),
    fitScale = unname(coef(fit)[2L]),
    mapR2 = summary(fit)$r.squared
  )
}

# Collect one BCF posterior sample's glue + per-forest internal fits.
sbcBCFSample <- function(bcf) {
  res <- .bcfRun(bcf, 0L, 1L)
  glue <- .bcfGlue(bcf)
  list(
    a = glue[1L],
    b0 = glue[2L],
    b1 = glue[3L],
    mu = .bcfForest(bcf, 0L)[, 1],
    tau = .bcfForest(bcf, 1L)[, 1],
    sigma = as.numeric(res$sigma)[1L]
  )
}

runSbcBCF <- function(
  config,
  R = 200L,
  L = 200L,
  thin = 30L,
  burn = 50L,
  seed = 20260709L,
  report = 25L,
  fixedGlue = FALSE
) {
  set.seed(seed)
  built <- sbcMakeBCF(config, L, thin, fixedGlue = fixedGlue)
  bcf <- built$bcf
  drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
  drawGlue <- if (fixedGlue) {
    function() list(a = 1, b0 = 0, b1 = 1) # the engine's fixed initial glue
  } else {
    sbcBCFGlueDraw(config$sdControl, config$bPriorVariance)
  }
  idx <- seq_len(config$nTest) # first nTest training rows are the eval points
  ranks <- NULL
  started <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    # theta0: forests from the engine prior, glue + sigma drawn in R
    .Call(.bcfPriorTrees, bcf$ptr)
    .Call(.bcfPriorNodes, bcf$ptr)
    mu0 <- .bcfForest(bcf, 0L)[, 1]
    tau0 <- .bcfForest(bcf, 1L)[, 1]
    g0 <- drawGlue()
    sig0 <- drawSigma(1L)
    bz0 <- ifelse(config$z != 0, g0$b1, g0$b0)
    internal0 <- g0$a * mu0 + bz0 * tau0
    reported0 <- built$fitScale * internal0 + built$fitShift
    y0 <- reported0 + sig0 * rnorm(config$n)

    # identified theta0 functionals (internal-scale, matched to posterior)
    prog0 <- g0$a * mu0[idx]
    eff0 <- (g0$b1 - g0$b0) * tau0[idx]

    # overdispersed init (fresh prior forests), then fit and collect
    .Call(.bcfPriorTrees, bcf$ptr)
    .Call(.bcfPriorNodes, bcf$ptr)
    .bcfSetResponse(bcf, y0, FALSE)
    invisible(.bcfRun(bcf, burn, 0L))

    aDraws <- numeric(L)
    diffDraws <- numeric(L)
    sigmaDraws <- numeric(L)
    progDraws <- matrix(NA_real_, config$nTest, L)
    effDraws <- matrix(NA_real_, config$nTest, L)
    for (l in seq_len(L)) {
      s <- sbcBCFSample(bcf)
      aDraws[l] <- s$a
      diffDraws[l] <- s$b1 - s$b0
      sigmaDraws[l] <- s$sigma
      progDraws[, l] <- s$a * s$mu[idx]
      effDraws[, l] <- (s$b1 - s$b0) * s$tau[idx]
    }

    # a*mu and b_z*tau are invariant under a joint sign flip, so the raw a and
    # (b1-b0) posteriors are sign-symmetric (bimodal) and their SBC is ill-posed
    # -- reported to demonstrate that. The identified quantities are the
    # magnitudes |a|, |b1-b0| and the functions prog=a*mu, eff=(b1-b0)*tau.
    # Under fixedGlue the glue functionals are degenerate constants; drop them.
    row <- if (fixedGlue) {
      c(sigma = sum(sigmaDraws < sig0))
    } else {
      c(
        sigma = sum(sigmaDraws < sig0),
        a = sum(aDraws < g0$a),
        abs.a = sum(abs(aDraws) < abs(g0$a)),
        b1.minus.b0 = sum(diffDraws < (g0$b1 - g0$b0)),
        abs.diff = sum(abs(diffDraws) < abs(g0$b1 - g0$b0))
      )
    }
    for (j in seq_len(config$nTest)) {
      row[paste0("prog", j)] <- sum(progDraws[j, ] < prog0[j])
      row[paste0("eff", j)] <- sum(effDraws[j, ] < eff0[j])
    }
    if (is.null(ranks)) {
      ranks <- matrix(
        NA_integer_,
        R,
        length(row),
        dimnames = list(NULL, names(row))
      )
    }
    ranks[r, ] <- row
    if (report > 0L && (r %% report == 0L || r == R)) {
      elapsed <- proc.time()[["elapsed"]] - started
      cat(sprintf(
        "  [bcf] rep %d/%d  %.1fs  %.2fs/rep  (map R2=%.5f)\n",
        r,
        R,
        elapsed,
        elapsed / r,
        built$mapR2
      ))
    }
  }
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    ranks = ranks,
    L = L,
    thin = thin,
    burn = burn,
    R = R,
    config = config,
    elapsed = elapsed,
    perRep = elapsed / R,
    mapR2 = built$mapR2
  )
}

# --- diagnostics -----------------------------------------------------------

# Autocorrelation of a long unthinned chain, to justify the thinning choice.
# Fits one prior-drawn dataset, runs `nDraw` unthinned samples, and reports the
# ACF of sigma and a couple of f(x*) functionals plus the first lag under 0.1.
sbcThinningDiagnostic <- function(
  config,
  nDraw = 4000L,
  burn = 500L,
  seed = 3L
) {
  set.seed(seed)
  sampler <- sbcMakeSampler(config, nDraw, 1L, seed)
  drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  f0Train <- as.numeric(sampler$predict(config$x))
  sig0 <- if (config$hasSigma) drawSigma(1L) else 1.0
  y0 <- sbcSimulate(config, f0Train, sig0)
  sampler$sampleTreesFromPrior()
  sampler$sampleNodeParametersFromPrior()
  if (config$hasSigma) {
    sampler$setSigma(config$sigest)
  }
  sampler$setResponse(y0)
  res <- sampler$run(burn, nDraw)
  chains <- list(
    avg.f = colMeans(res$train),
    f.star1 = res$test[1, ],
    f.star3 = res$test[min(3L, config$nTest), ]
  )
  if (config$hasSigma) {
    chains$sigma <- as.numeric(res$sigma)
  }
  firstUnder <- function(v) {
    a <- acf(v, lag.max = 60L, plot = FALSE)$acf[,, 1]
    idx <- which(a < 0.1)[1]
    if (is.na(idx)) NA_integer_ else idx - 1L
  }
  lags <- c(1L, 2L, 4L, 8L, 15L, 30L, 45L)
  acfAt <- function(v) {
    a <- acf(v, lag.max = max(lags), plot = FALSE)$acf[,, 1]
    a[lags + 1L]
  }
  list(
    lags = lags,
    acf = lapply(chains, acfAt),
    firstUnder = vapply(chains, firstUnder, integer(1))
  )
}

# Uniformity verdict for one functional's ranks. The headline verdict is the
# ecdf-difference statistic against a simulation-based simultaneous 95% band
# (Talts fig. 1): already corrected for multiple looks across the rank grid, so
# it is the robust primary test. Chi-square goodness of fit on equal-width bins
# and a KS test against the discrete uniform (jitter handles rank ties) are
# reported as secondary signals -- at 20 bins over R=200 a lone chi-square
# p < 0.01 across many functionals is within multiple-comparison noise, so the
# verdict does not hinge on it.
rankUniformity <- function(
  ranks,
  L,
  nBins = 20L,
  nSim = 2000L,
  alpha = 0.05,
  seed = 1L
) {
  R <- length(ranks)
  # chi-square on nBins equal-width bins of {0, ..., L}
  edges <- seq(0, L + 1L, length.out = nBins + 1L)
  counts <- as.integer(table(cut(
    ranks,
    breaks = edges,
    include.lowest = TRUE,
    right = FALSE
  )))
  expected <- R / nBins
  chisqStat <- sum((counts - expected)^2 / expected)
  chisqP <- pchisq(chisqStat, df = nBins - 1L, lower.tail = FALSE)

  # KS against discrete uniform on {0, ..., L}: map rank -> (rank + U) / (L + 1)
  set.seed(seed)
  u <- (ranks + runif(R)) / (L + 1)
  ksP <- suppressWarnings(ks.test(u, "punif")$p.value)

  # ecdf-difference simultaneous band via simulation of the null
  grid <- 0:L
  ecdfDiff <- function(rk) {
    ec <- vapply(grid, function(g) mean(rk <= g), numeric(1))
    ec - (grid + 1) / (L + 1)
  }
  observed <- max(abs(ecdfDiff(ranks)))
  nullMax <- numeric(nSim)
  for (s in seq_len(nSim)) {
    nullMax[s] <- max(abs(ecdfDiff(sample.int(L + 1L, R, replace = TRUE) - 1L)))
  }
  band <- as.numeric(quantile(nullMax, 1 - alpha))
  list(
    counts = counts,
    nBins = nBins,
    expected = expected,
    chisqP = chisqP,
    ksP = ksP,
    ecdfDiff = observed,
    ecdfBand = band,
    pass = observed <= band,
    mean = mean(ranks),
    meanTarget = L / 2
  )
}

# A compact ASCII rank histogram with the +/- band around the uniform mean.
sbcAsciiHistogram <- function(ranks, L, nBins = 20L, width = 40L) {
  edges <- seq(0, L + 1L, length.out = nBins + 1L)
  counts <- as.integer(table(cut(
    ranks,
    breaks = edges,
    include.lowest = TRUE,
    right = FALSE
  )))
  expected <- length(ranks) / nBins
  scale <- width / max(counts, expected)
  lines <- character(nBins)
  for (b in seq_len(nBins)) {
    bar <- strrep("#", round(counts[b] * scale))
    lines[b] <- sprintf(
      "  %3d-%3d | %-*s %d",
      round(edges[b]),
      round(edges[b + 1L]) - 1L,
      width,
      bar,
      counts[b]
    )
  }
  paste(
    c(sprintf("  (expected %.1f per bin, uniform)", expected), lines),
    collapse = "\n"
  )
}

# Full report for a runSbc result: per-functional verdict table + histograms.
sbcReport <- function(fit, nBins = 20L) {
  cat(sprintf(
    "\nSBC report: family=%s n=%d p=%d nTrees=%d | R=%d L=%d thin=%d burn=%d\n",
    fit$config$family,
    fit$config$n,
    fit$config$p,
    fit$config$nTrees,
    fit$R,
    fit$L,
    fit$thin,
    fit$burn
  ))
  cat(sprintf(
    "wall-clock: %.1fs total, %.3fs/rep\n\n",
    fit$elapsed,
    fit$perRep
  ))
  funcs <- colnames(fit$ranks)
  cat(sprintf(
    "%-10s %8s %8s %9s %8s %6s\n",
    "functional",
    "chisqP",
    "ksP",
    "ecdfDiff",
    "band",
    "verdict"
  ))
  verdicts <- character(length(funcs))
  for (i in seq_along(funcs)) {
    u <- rankUniformity(fit$ranks[, funcs[i]], fit$L, nBins = nBins)
    verdicts[i] <- if (u$pass) "PASS" else "FLAG"
    cat(sprintf(
      "%-10s %8.3f %8.3f %9.4f %8.4f %6s\n",
      funcs[i],
      u$chisqP,
      u$ksP,
      u$ecdfDiff,
      u$ecdfBand,
      verdicts[i]
    ))
  }
  cat("\nRank histograms:\n")
  for (i in seq_along(funcs)) {
    cat(sprintf("\n[%s]\n", funcs[i]))
    cat(sbcAsciiHistogram(fit$ranks[, funcs[i]], fit$L, nBins = nBins), "\n")
  }
  invisible(verdicts)
}

# --- main ------------------------------------------------------------------

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  which <- if (length(args) >= 1L) args[1] else "gaussian"
  R <- if (length(args) >= 2L) as.integer(args[2]) else 200L
  L <- if (length(args) >= 3L) as.integer(args[3]) else 200L
  thin <- if (length(args) >= 4L) as.integer(args[4]) else 30L

  isDart <- which %in% c("dart", "dart-sparse")
  isWeighted <- which == "weighted"
  isGrouped <- which %in% c("grouped-gaussian", "grouped-probit")
  isBCF <- which %in% c("bcf", "bcf-weak")
  isLinear <- which %in%
    c("linear", "linear-na-leaf", "linear-na-split", "linear-weighted")
  isGP <- which %in% c("gp", "gp-na-leaf", "gp-weighted")

  config <- if (isDart) {
    sbcConfig(
      family = "gaussian",
      n = 200L,
      p = 10L,
      dartAlpha = if (which == "dart-sparse") 0.05 else 1.0
    )
  } else if (isWeighted) {
    cfg <- sbcConfig(family = "gaussian")
    set.seed(7L)
    cfg$weights <- rgamma(cfg$n, 2, 2) # known, positive, mean 1
    cfg
  } else if (isGrouped) {
    base <- if (which == "grouped-probit") "probit" else "gaussian"
    zw <- if (which == "grouped-gaussian") 0.2 else 0
    sbcAddGrouping(
      sbcConfig(family = base, n = 160L),
      nGroups = 8L,
      relScale = 0.2,
      zeroWeightFrac = zw
    )
  } else if (isBCF) {
    # prior-weak = small n so the a-glue prior term dominates the likelihood
    nBcf <- if (which == "bcf-weak") 40L else 200L
    sbcAddBCF(sbcConfig(family = "gaussian", n = nBcf))
  } else if (isLinear) {
    # columns 1:2 fit linearly inside leaves; column 3 is split-only
    cfg <- sbcConfig(
      family = "gaussian",
      nodePrior = dbartsPriors$linear(1:2, k = 2)
    )
    if (which == "linear-na-leaf") {
      cfg <- sbcAddMissing(cfg, columns = 1L) # NA in a designated leaf column
    } else if (which == "linear-na-split") {
      cfg <- sbcAddMissing(cfg, columns = 3L) # NA routed by splits only
    } else if (which == "linear-weighted") {
      set.seed(7L)
      cfg$weights <- rgamma(cfg$n, 2, 2)
    }
    cfg
  } else if (isGP) {
    # column 1 fits a GP inside leaves; max.leaf.size = 100 matches the
    # equivalence gp scenario's cap (n < 100 keeps every leaf a true GP leaf,
    # never the constant fallback). k is FIXED at 2: the equivalence
    # scenario's chi hyperprior samples k, but sampleNodeParametersFromPrior
    # draws at the CURRENT k with no API to install a hyperprior draw, so a
    # sampled-k SBC would be prior-mismatched (residual gap, recorded).
    # n/nTrees sized by measured cost: prior-drawn trees are shallow, so leaf
    # kernel solves scale with n^3 (31 ms/sweep at n=150 vs 1.8 at n=80).
    cfg <- sbcConfig(
      family = "gaussian",
      n = 80L,
      nTrees = 25L,
      nodePrior = dbartsPriors$gp(1L, k = 2, max.leaf.size = 100L)
    )
    cfg$f0FromForestFits <- TRUE
    if (which == "gp-na-leaf") {
      cfg <- sbcAddMissing(cfg, columns = 1L)
    } else if (which == "gp-weighted") {
      set.seed(7L)
      cfg$weights <- rgamma(cfg$n, 2, 2)
    }
    cfg
  } else {
    sbcConfig(family = which)
  }

  cat("== prior moment check ==\n")
  chk <- sbcCheckSigmaPrior(config$sigest, config$sigDf, config$sigQuant)
  cat(sprintf(
    "  sigma: P(sigma < sigest) = %.4f (target %.2f); median %.4f vs %.4f -> %s\n",
    chk$coverage,
    chk$coverageTarget,
    chk$medianEmpirical,
    chk$medianTheory,
    if (chk$pass) "PASS" else "FAIL"
  ))
  if (isDart) {
    d <- sbcCheckDirichlet(config$dartAlpha, config$p)
    cat(sprintf(
      "  dirichlet(a=%.2f,p=%d): mean %.4f vs %.4f; var %.2e vs %.2e -> %s\n",
      config$dartAlpha,
      config$p,
      d$meanEmp,
      d$meanTheory,
      d$varEmp,
      d$varTheory,
      if (d$pass) "PASS" else "NOTE(floor)"
    ))
  }
  if (isGrouped) {
    tc <- sbcCheckTau(config$relScale)
    cat(sprintf(
      "  tau: gamma mean %.4f vs %.4f; var %.4f vs %.4f -> %s\n",
      tc$meanEmp,
      tc$meanTheory,
      tc$varEmp,
      tc$varTheory,
      if (tc$pass) "PASS" else "FAIL"
    ))
  }
  if (isBCF) {
    gc <- sbcCheckBCFGlue(config$sdControl, config$bPriorVariance)
    cat(sprintf(
      "  glue: a Cauchy IQR %.4f vs %.4f; b sd %.4f vs %.4f -> %s\n",
      gc$aIqrEmp,
      gc$aIqrTheory,
      gc$bSdEmp,
      gc$bSdTheory,
      if (gc$pass) "PASS" else "FAIL"
    ))
  }
  if (isLinear || isGP) {
    fc <- sbcCheckFitConsistency(config)
    cat(sprintf(
      "  predict vs recorded fits: train %.2e, test %.2e -> %s\n",
      fc$maxDiff,
      fc$maxDiffTest,
      if (fc$pass) "PASS" else "FAIL"
    ))
  }

  cat(sprintf("\n== SBC run (%s R=%d L=%d thin=%d) ==\n", which, R, L, thin))
  fit <- if (isDart) {
    runSbcDart(config, R = R, L = L, thin = thin)
  } else if (isGrouped) {
    runSbcGrouped(config, R = R, L = L, thin = thin)
  } else if (isBCF) {
    runSbcBCF(config, R = R, L = L, thin = thin)
  } else {
    runSbc(config, R = R, L = L, thin = thin)
  }
  if (isDart) {
    cat(sprintf(
      "\nfloor incidence: %.3f of s0 components pinned at 1e-300\n",
      fit$floorFrac
    ))
  }
  sbcReport(fit)
}
