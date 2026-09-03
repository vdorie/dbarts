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
#   Rscript benchmarks/R/sbc.R ordinal 200 150 30   # family tiers, plan
#   Rscript benchmarks/R/sbc.R nbinom|t|multinom 200 150 30
#   Rscript benchmarks/R/sbc.R grouped-gaussian-swap 200 200 30  # the swap arm
#   Rscript benchmarks/R/sbc.R aft 200 200 30       # aft/survival, rebuilt
#   Rscript benchmarks/R/sbc.R gp-mixed 200 150 60 3000 # GP/constant mix
#   Rscript benchmarks/R/sbc.R discrete-selfcheck   # the discrete-rank gate
#   Rscript benchmarks/R/sbc.R burn-ordinal 20000 3 # the burn/cost ladder
# Positional args: config R L thin, plus an optional 5th, the burn in absolute
# sweeps, and an optional 6th, the driver seed. Or source() the file to reuse
# the API:
#   source("benchmarks/R/sbc.R"); res <- runSbc(sbcConfig("gaussian"), R = 200)
# SBC_FAIL_ON_FLAG=1 (env var, opt-in) makes the CLI exit status 1 if any
# functional's verdict is FLAG; unset, the CLI always exits 0 (unchanged
# default). Only affects Rscript use; source() usage is untouched.
#
# SELF-CONSISTENCY is the whole game: theta0 must come from the same prior the
# sampler assumes in its posterior. The forest/leaf draw uses the sampler's own
# sampleTreesFromPrior + sampleNodeParametersFromPrior; the sigma draw is the
# reported-scale scaled-inverse-chi-squared the engine calibrates (verified by
# a moment check, sbcCheckSigmaPrior). The data scale is fixed once at build
# (setResponse with updateScale = FALSE keeps it) so prior and posterior share
# it. A wrong prior draw makes SBC lie, so the harness self-checks before use.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

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

# --- discrete (grid) parameters --------------------------------------------

# Rank of theta0 among L posterior draws WHEN THE LAW HAS ATOMS. #{draws <
# theta0} is uniform only for an atomless law: an atom parks all its mass on one
# rank. Attach an iid Uniform(0, 1) tag to every draw AND to theta0 itself and
# rank the pairs lexicographically: theta0's tag is exchangeable with the tags
# of the tied draws, so the atom contributes a Uniform{0, ..., #ties} increment
# and the total rank is uniform on {0, ..., L} under calibration - exactly
# rankUniformity's null, unchanged.
#
# TWO kinds of atom need it, which is why the family driver applies it to every
# functional rather than only the declared-discrete ones (with no ties it is
# #{draws < theta0} exactly, and it consumes no rng, so an atomless functional
# is untouched). The obvious kind is a genuinely DISCRETE parameter - nbinom's
# grid dispersion r, the Student-t grid nu. The second is NUMERICAL: an ordinal
# top-category probability is mean_i (1 - Phi(gamma_K-1 - eta_i)), which
# UNDERFLOWS to exactly 0 whenever the prior draws the top cutpoint far out (a
# quarter of replications at K = 4, the empty-cell case ordinal.md section 9
# names), so theta0 and most of its posterior draws are all exactly 0. Without
# the tie-break those replications pile up at rank 0 and the functional flags -
# the same tie-degenerate artifact the DART 1e-300 floor probe recorded, not a
# calibration defect.
sbcDiscreteRank <- function(draws, theta0) {
  below <- sum(draws < theta0)
  ties <- sum(draws == theta0)
  if (ties == 0L) {
    return(below)
  }
  tag0 <- runif(1L)
  below + sum(runif(ties) < tag0)
}

# The engine's two DISCRETE grid priors, transcribed from src/bartcore/model.hpp
# (NBDispersionPrior, ResidualDfPrior): both normalize the same gamma(2, 0.1)
# kernel w_k propto grid_k * exp(-0.1 * grid_k) over a fixed capped grid, so a
# self-consistent prior draw must use the identical grid AND weights.
sbcNbGrid <- c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50)
sbcTGrid <- c(3, 4, 5, 6, 8, 10, 12, 15, 20)

sbcGridWeights <- function(grid) {
  w <- grid * exp(-0.1 * grid)
  w / sum(w)
}

sbcGridDraw <- function(grid) {
  w <- sbcGridWeights(grid)
  function(nDraws = 1L) sample(grid, nDraws, replace = TRUE, prob = w)
}

# Moment/calibration check for a grid prior: the empirical cell frequencies and
# the mean must match the normalized kernel.
sbcCheckGridPrior <- function(grid, nDraws = 2e5L) {
  w <- sbcGridWeights(grid)
  draws <- sbcGridDraw(grid)(nDraws)
  emp <- as.numeric(table(factor(draws, levels = grid))) / nDraws
  meanTheory <- sum(grid * w)
  list(
    maxCellDiff = max(abs(emp - w)),
    meanEmpirical = mean(draws),
    meanTheory = meanTheory,
    pass = max(abs(emp - w)) < 0.005 &&
      abs(mean(draws) / meanTheory - 1) < 0.02
  )
}

# Step-1 self-check for sbcDiscreteRank: a synthetic conjugate case whose
# posterior is available in CLOSED FORM, so the L "posterior draws" are exact
# and iid and any non-uniformity is the ranking rule's fault rather than a
# sampler's. The case mirrors the engine's own dispersion update - r0 from the
# nbinom grid prior, counts y_i ~ NB(r0, p) at a KNOWN p, posterior propto
# prior_k * prod_i dnbinom(y_i, r_k, p) over the same grid - so it also
# exercises the grid prior the nbinom arm draws from. n is small on purpose:
# the posterior must stay diffuse enough to tie often, which is the case the
# tie-breaker exists for.
sbcDiscreteSelfCheck <- function(
  R = 400L,
  L = 150L,
  n = 25L,
  prob = 0.5,
  seed = 20260804L
) {
  set.seed(seed)
  grid <- sbcNbGrid
  logPrior <- log(sbcGridWeights(grid))
  drawR <- sbcGridDraw(grid)
  ranks <- integer(R)
  tieFrac <- numeric(R)
  for (rep in seq_len(R)) {
    r0 <- drawR(1L)
    y <- rnbinom(n, size = r0, prob = prob)
    logPost <- logPrior +
      vapply(
        grid,
        function(rk) sum(dnbinom(y, size = rk, prob = prob, log = TRUE)),
        numeric(1)
      )
    post <- exp(logPost - max(logPost))
    draws <- sample(grid, L, replace = TRUE, prob = post)
    ranks[rep] <- sbcDiscreteRank(draws, r0)
    tieFrac[rep] <- mean(draws == r0)
  }
  uniformity <- rankUniformity(ranks, L)
  list(
    ranks = ranks,
    uniformity = uniformity,
    tieFrac = mean(tieFrac),
    pass = isTRUE(uniformity$pass)
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
# numCategories is K for the two categorical families (ordinal's ordered levels,
# multinomial's softmax categories) and is ignored elsewhere.
sbcConfig <- function(
  family = c(
    "gaussian",
    "probit",
    "logistic",
    "ordinal",
    "nbinom",
    "t",
    "multinomial"
  ),
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
  numCategories = 4L,
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
  # ordinal builds from an ordered factor over ALL K levels (its level set is
  # what fixes K, and a rebuilt fit must never re-derive a smaller K from a
  # replication whose simulated y happens to miss a category); nbinom builds
  # from a small count vector; t and multinomial build like their host family
  # (continuous gaussian).
  yBuild <- if (family %in% c("gaussian", "t", "multinomial")) {
    seq(-2.5, 2.5, length.out = n)
  } else if (family == "ordinal") {
    factor(
      rep_len(seq_len(numCategories), n),
      levels = seq_len(numCategories),
      ordered = TRUE
    )
  } else if (family == "nbinom") {
    as.double(rep_len(c(0L, 1L, 2L, 4L), n))
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
    K = as.integer(numCategories),
    hasSigma = family %in% c("gaussian", "t")
  )
}

# The dbarts() family token for a configuration: the Student-t and multinomial
# arms build a GAUSSIAN host (t adds resid.dist = student(); multinomial wraps
# the host in the K-forest softmax sampler), everything else names itself.
sbcSamplerFamily <- function(config) {
  switch(config$family, t = "gaussian", multinomial = "gaussian", config$family)
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

# The creation-time leaf-scale pin: what a plain sampler's own calibration
# resolves prior.scale to off a config's fixed build response. A rebuild that
# NAMES this as its own node.prior scale (dbartsPriors$normal(k, scale = )) is
# scored against the SAME response-unit leaf prior every replication, instead
# of one each rebuild's own y range would silently re-derive -- the transform
# is what moves per rebuild, never the leaf scale itself, so naming the
# transform's anchor is what holds the prior fixed (chain.hpp's
# resolvedNodeScale).
sbcAnchorScale <- function(config) {
  anchor <- dbarts(
    config$x,
    config$yBuild,
    control = dbartsControl(
      n.trees = config$nTrees,
      n.chains = 1L,
      n.threads = 1L,
      updateState = FALSE,
      verbose = FALSE
    )
  )
  anchor$getCalibration(1L)[1L, "prior.scale"]
}

# Build the reusable sampler for a configuration. One sampler serves all
# replications: the prior draw advances its internal RNG, setResponse swaps y
# in place, and the fixed build scale is never disturbed. `y` overrides the
# build response for the families whose fit is REBUILT per replication (ordinal
# and nbinom keep a slow-moving global - the cutpoints, the dispersion - across
# setResponse, which would break rank iid-ness); those families run at a fixed
# unit scale, so a rebuild re-anchors nothing.
sbcMakeSampler <- function(config, L, thin, seed, y = NULL) {
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
  family <- sbcSamplerFamily(config)
  if (is.null(y)) {
    y <- config$yBuild
  }
  # The matrix (xy) interface DROPS NA rows even under missing = "incorporate"
  # (dbartsData warns "row(s) dropped"); only the formula interface keeps them
  # (na.action = na.pass). NA designs therefore build through a formula.
  if (anyNA(config$x)) {
    df <- as.data.frame(config$x)
    df$.sbc.y <- y
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
      y,
      test = config$xTest,
      resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
      node.prior = config$nodePrior,
      sigma = config$sigest,
      control = ctrl,
      family = family
    )
  }
  # Student-t errors are a residual DISTRIBUTION on a gaussian response, not a
  # family; the constructor vocabulary is unexported, so reach it by namespace
  # exactly as the harness reaches the internal bartcore entry points.
  if (config$family == "t") {
    args$resid.dist <- getFromNamespace("dbartsResidDists", "dbarts")$student()
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
  # in absolute sweeps: the BCF sigma transient is tree-STRUCTURE mixing
  # under strong prognostic signal (settle ~72k sweeps at the Cauchy tail;
  # bcf-sigma-residual), so the default pins sweeps, not thinned units
  burn = as.integer(ceiling(72000 / thin)),
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
  # in absolute sweeps: the BCF sigma transient is tree-STRUCTURE mixing
  # under strong prognostic signal (settle ~72k sweeps at the Cauchy tail;
  # bcf-sigma-residual), so the default pins sweeps, not thinned units
  burn = as.integer(ceiling(72000 / thin)),
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
# calibration. The default arm REBUILDS the fit each replication; swap = TRUE
# instead keeps one fit and swaps the response into it, the arm that gates the
# bridge's grouped setResponse (see runSbcGrouped). Groups are
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

# swap = TRUE is the calibration gate on the bridge's grouped setResponse: one
# fit serves every replication, re-initialised through the STATE INSTALL rather
# than rebuilt. state0 is the pristine post-build state - two lines, because
# $storeState() returns invisible(NULL) and sbcMakeGroupedFit sets updateState =
# FALSE, so the field is filled by nothing else - and it is the only channel
# that returns b and tau to a y0-independent start (there is no $setTau). That
# start is CONSTANT rather than a fresh prior draw for tau: SBC needs
# independence from y0, not a prior draw, and the constancy is the one thing
# this arm does not randomize. The fixed build response also pins the internal
# scale for every replication, which is what setResponse(updateScale = FALSE)
# keeps -- the same self-consistency runSbc gets from its reused sampler. The
# state carries the chain rng too, so every replication restarts it at the same
# point; the draws still differ, since each conditions on its own y0, and each
# rank is marginally uniform, which is what rankUniformity tests.
runSbcGrouped <- function(
  config,
  R = 200L,
  L = 200L,
  thin = 30L,
  # in absolute sweeps: the BCF sigma transient is tree-STRUCTURE mixing
  # under strong prognostic signal (settle ~72k sweeps at the Cauchy tail;
  # bcf-sigma-residual), so the default pins sweeps, not thinned units
  burn = as.integer(ceiling(72000 / thin)),
  seed = 20260709L,
  report = 25L,
  swap = FALSE
) {
  set.seed(seed)
  gen <- sbcMakeGroupGenerator(config)
  drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
  drawTau <- sbcTauDraw(config$relScale)
  fit <- NULL
  state0 <- NULL
  if (swap) {
    fit <- sbcMakeGroupedFit(config, config$yBuild, L, thin)
    fit$storeState()
    state0 <- fit$state
  }
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

    if (swap) {
      # the same overdispersed init the rebuild arm gets from a fresh sampler:
      # a state install, then a second independent prior draw for the forest
      fit$setState(state0)
      fit$sampleTreesFromPrior()
      fit$sampleNodeParametersFromPrior()
      if (config$hasSigma) {
        fit$setSigma(config$sigest)
      }
      fit$setResponse(y0, updateScale = FALSE)
    } else {
      fit <- sbcMakeGroupedFit(config, y0, L, thin)
    }
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

# --- aft (accelerated failure time / survival) ------------------------------

# The engine's aft sampler (docs/design/survival.md) is a log-normal
# survival model: log T = f(x) + sigma*eps, uncensored rows enter as gaussian
# data on log T and right-censored ones contribute the upper normal tail past
# their log censoring time. Per-row censoring status is STRUCTURAL -- baked
# into the handle at construction (the bartcore.survival control attribute) --
# so unlike every other arm here the fit is REBUILT every replication rather
# than reused through setResponse; a drawn theta0 changes which rows censor.
# node.prior carries the ANCHOR scale (sbcAnchorScale, off the fixed build
# response) so every rebuild scores its draws against the prior theta0 came
# from rather than one its own y0 range would re-derive. The offset channel
# then pins the shift the same way: the build response is symmetric about 0,
# so the generator's own in-force prior.mean is 0, and zeroing each rebuild's
# prior.mean (setOffset(rep_len(-prior.mean, n), updateScale = FALSE), the
# documented recipe on $setCalibration) matches it. sigma is drawn
# conjugately exactly as gaussian's (chain.hpp), so it is this arm's log-time
# scale-parameter functional, ranked the same way avg.f/f.star are.

.aftMake <- getFromNamespace("bartcoreSampler", "dbarts")
.aftRun <- bartcoreRun
.aftPredict <- bartcorePredict
.aftCal <- bartcoreForestCalibration
.aftSetOffset <- bartcoreSetOffset
.aftSetTestOffset <- bartcoreSetTestOffset

# Fixed per-row right-censoring log-times: a design choice independent of any
# prior draw, pinned once (like sbcAddGrouping's grouping) so only the drawn
# log-time and the status it implies move across replications.
sbcAddCensoring <- function(config, shift = 1.0, sd = 1.2) {
  set.seed(505L)
  config$logC <- config$yBuild + shift + sd * rnorm(config$n)
  config
}

# The arm's fixed config: a continuous, symmetric-about-0 build response (so
# the shift pin's target is exactly 0), m = 50 trees, the anchor scale named
# as node.prior, and the fixed censoring fixture above.
sbcConfigAft <- function(n = 150L, nTrees = 50L) {
  config <- sbcConfig(family = "gaussian", n = n, nTrees = nTrees)
  config <- sbcAddCensoring(config)
  config$nodePrior <- dbartsPriors$normal(
    config$k,
    scale = sbcAnchorScale(config)
  )
  config$family <- "aft"
  config
}

# The f-only generator theta0 is drawn from: a plain sampler (response family
# is irrelevant -- only sampleTreesFromPrior/predict are used) at the
# config's pinned node.prior and fixed build response.
sbcMakeAftGenerator <- function(config) {
  dbarts(
    config$x,
    config$yBuild,
    test = config$xTest,
    resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
    node.prior = config$nodePrior,
    sigma = config$sigest,
    control = dbartsControl(
      n.trees = config$nTrees,
      n.chains = 1L,
      n.threads = 1L,
      n.samples = 1L,
      updateState = FALSE,
      verbose = FALSE,
      keepTrainingFits = TRUE
    )
  )
}

# A fresh aft handle for one y0/status0 pair (observed log-time, event
# indicator): pinned node.prior scale, then the shift pin above. Returns the
# handle with the applied offset, since bartcorePredict never reads a
# handle's own offset back (a caller must restate it) and the run's train and
# test channels are SEPARATE offset stores that both need the same value.
sbcMakeAftFit <- function(config, y, status, L, thin) {
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
  s <- dbarts(
    config$x,
    y,
    test = config$xTest,
    resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
    node.prior = config$nodePrior,
    sigma = config$sigest,
    control = ctrl
  )
  c2 <- s$control
  attr(c2, "bartcore.survival") <- as.numeric(status)
  s$control <- c2
  bc <- .aftMake(s, family = "aft")
  off <- -.aftCal(bc, 0L)[1L, "prior.mean"]
  .aftSetOffset(bc, rep_len(off, config$n), FALSE)
  .aftSetTestOffset(bc, rep_len(off, nrow(config$xTest)))
  list(bc = bc, off = off)
}

# predict() vs the recorded fit at one state -- the same wiring check every
# other arm runs before trusting its ranks (sbcCheckFitConsistency,
# sbcCheckLatentConsistency). y is the build response SHIFTED (not the build
# response itself, which is already centred at the pin's target) so the shift
# pin's offset is genuinely non-zero here, exercising the train/test offset
# pairing sbcMakeAftFit relies on rather than the degenerate zero-offset case.
sbcCheckAftFitConsistency <- function(config, seed = 99L) {
  set.seed(seed)
  built <- sbcMakeAftFit(
    config,
    config$yBuild + 3,
    rep_len(1, config$n),
    1L,
    1L
  )
  bc <- built$bc
  res <- .aftRun(bc, 0L, 1L)
  offTrain <- rep_len(built$off, config$n)
  offTest <- rep_len(built$off, nrow(config$xTest))
  maxDiff <- max(abs(res$train[, 1] - .aftPredict(bc, config$x, offTrain)))
  maxDiffTest <- max(abs(
    res$test[, 1] - .aftPredict(bc, config$xTest, offTest)
  ))
  list(
    maxDiff = maxDiff,
    maxDiffTest = maxDiffTest,
    pass = maxDiff < 1e-8 && maxDiffTest < 1e-8
  )
}

runSbcAft <- function(
  config,
  R = 200L,
  L = 200L,
  thin = 30L,
  burn = as.integer(ceiling(72000 / thin)),
  seed = 20260709L,
  report = 25L
) {
  set.seed(seed)
  gen <- sbcMakeAftGenerator(config)
  drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
  ranks <- NULL
  nCensored <- 0L
  started <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    gen$sampleTreesFromPrior()
    gen$sampleNodeParametersFromPrior()
    f0Train <- as.numeric(gen$predict(config$x))
    f0Test <- as.numeric(gen$predict(config$xTest))
    sig0 <- drawSigma(1L)
    avgF0 <- mean(f0Train)

    logT0 <- f0Train + sig0 * rnorm(config$n)
    status0 <- as.numeric(logT0 <= config$logC)
    y0 <- ifelse(status0 == 1, logT0, config$logC)
    nCensored <- nCensored + sum(status0 == 0)

    bc <- sbcMakeAftFit(config, y0, status0, L, thin)$bc
    res <- .aftRun(bc, burn, L)

    row <- c(
      avg.f = sum(colMeans(res$train) < avgF0),
      sigma = sum(as.numeric(res$sigma) < sig0)
    )
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
        "  [aft] rep %d/%d  %.1fs  %.2fs/rep  censoring %.3f\n",
        r,
        R,
        elapsed,
        elapsed / r,
        nCensored / (r * config$n)
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
    censoringRate = nCensored / (R * config$n)
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
.bcfRun <- bartcoreRun
.bcfGlue <- bartcoreForestAmplitudes
.bcfForest <- bartcoreForestFits
.bcfSetResponse <- bartcoreSetResponse
.bcfPriorTrees <- getFromNamespace(
  "C_dbarts_bartcore_sampleTreesFromPrior",
  "dbarts"
)
.bcfPriorNodes <- getFromNamespace(
  "C_dbarts_bartcore_sampleNodeParametersFromPrior",
  "dbarts"
)
.bcfStoreState <- bartcoreStoreState
.bcfSetState <- bartcoreSetState

# Install a drawn (a, b0, b1) as the sampler's LIVE glue. The tree prior is
# glue-dependent: each forest's prior trees are drawn conditioned on the
# no-empty-leaf set of that forest's own veto vector, w * b_z^2 for the
# treatment forest, so trees drawn before the glue is installed come from a
# different law than the theta0 that reports them. Round-tripped through the
# state, whose glue block is [K, q_1..q_K, amplitudes, K prior variances] - here
# K = 2 with widths (1, 2) - and which is re-installed exactly as stored in
# every other respect, the rng included.
sbcInstallBCFGlue <- function(bcf, glue) {
  state <- .bcfStoreState(bcf)
  state[[1L]][["glue"]][4:6] <- c(glue$a, glue$b0, glue$b1)
  .bcfSetState(bcf, state)
  invisible(NULL)
}

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
  # in absolute sweeps: the BCF sigma transient is tree-STRUCTURE mixing
  # under strong prognostic signal (settle ~72k sweeps at the Cauchy tail;
  # bcf-sigma-residual), so the default pins sweeps, not thinned units
  burn = as.integer(ceiling(72000 / thin)),
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
    # theta0: glue + sigma drawn in R, forests from the engine prior. The glue
    # is drawn and INSTALLED FIRST - the forest prior is conditioned on it (see
    # sbcInstallBCFGlue), so the order is what keeps theta0 a draw from one
    # joint prior rather than from two inconsistent ones.
    g0 <- drawGlue()
    sbcInstallBCFGlue(bcf, g0)
    .Call(.bcfPriorTrees, bcf$ptr)
    .Call(.bcfPriorNodes, bcf$ptr)
    mu0 <- .bcfForest(bcf, 0L)[, 1]
    tau0 <- .bcfForest(bcf, 1L)[, 1]
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

# --- family tiers: ordinal, nbinom, Student-t, multinomial -----------------

# The four remaining shipped families that admit a well-posed SBC
# (docs/plans/sbc-family-tiers.md). Each supplies the same four operations, so
# ONE driver ranks them and ONE diagnostic measures their burn ladder:
#
#   draw()        theta0 from the family's prior + the y it implies, as a named
#                 vector of scalar FUNCTIONALS plus the simulated response
#   fit(y)        the posterior sampler for that y, already re-initialised from
#                 an independent prior draw (never from theta0)
#   burnRun(f, b) b thinned units of burn-in, no samples kept
#   sample(f)     one retained draw's functionals, same names as theta
#
# sample() collects ONE draw at a time because the two grid parameters (nbinom's
# r, the Student-t nu) ride the STATE, not a run channel, so they are only
# readable between samples; the extra .Call per draw is microseconds against a
# thinned sweep block. Every functional is ranked by sbcDiscreteRank.
#
# Rebuild policy (plan step 3): ordinal, nbinom and multinomial REBUILD the fit
# per replication. Each keeps a slow-moving global across a response swap
# (OrdinalResponse::setResponse keeps gamma, NBResponse::setResponse keeps r,
# and a multinomial response refuses whole-data mutation outright), which would
# correlate consecutive replications and break rank iid-ness; all three run at a
# fixed unit scale, so a rebuild re-anchors nothing. Only the Student-t arm
# reuses one pinned sampler: setResponse cold-inits nu and lambda, and
# updateScale = FALSE keeps the build scale the prior draw shares.

# The category draw for a row-stochastic n x K probability matrix: category
# 1..K per row by inverse CDF.
sbcCategoricalDraw <- function(probs) {
  u <- runif(nrow(probs))
  1L + rowSums(t(apply(probs, 1L, cumsum)) < u)
}

sbcSoftmax <- function(f) {
  e <- exp(f - apply(f, 1L, max))
  e / rowSums(e)
}

# The cutpoint prior the ordinal engine assumes (docs/design/ordinal.md section
# 3, OrdinalResponse::logGapTarget): gamma_1 = 0 is pinned and the K-2 free
# interior cutpoints ride iid normal LOG-GAPS delta_j ~ N(0, 1.5^2), so
# gamma_{j+1} = gamma_j + exp(delta_j). The shipped constants are
# priorLogGapMean_ = 0 and priorLogGapSd_ = 1.5 (src/bartcore/model.hpp); a
# mismatch here would make ordinal SBC lie.
sbcOrdinalLogGapSd <- 1.5

sbcOrdinalGapDraw <- function(nDraws) {
  exp(rnorm(nDraws, 0, sbcOrdinalLogGapSd))
}

sbcOrdinalCutpointDraw <- function(K) {
  function() c(0, cumsum(sbcOrdinalGapDraw(K - 2L)))
}

# Moment check: the log-gaps are N(0, 1.5^2), so a gap has median exp(0) = 1.
sbcCheckOrdinalCutpointPrior <- function(nDraws = 2e5L) {
  gaps <- sbcOrdinalGapDraw(nDraws)
  logGaps <- log(gaps)
  list(
    sdEmpirical = sd(logGaps),
    sdTheory = sbcOrdinalLogGapSd,
    medianEmpirical = median(gaps),
    medianTheory = 1,
    pass = abs(sd(logGaps) / sbcOrdinalLogGapSd - 1) < 0.02 &&
      abs(median(gaps) - 1) < 0.03
  )
}

# The cumulative-probit category probabilities, P(y = k) = Phi(gamma_k - eta) -
# Phi(gamma_{k-1} - eta) with gamma_0 = -Inf and gamma_K = +Inf (the harness's
# own copy of the package's ordinalCategoryProbabilities).
sbcOrdinalProbs <- function(eta, gamma) {
  n <- length(eta)
  K <- length(gamma) + 1L
  bounds <- matrix(0, n, K + 1L)
  bounds[, K + 1L] <- 1
  for (j in seq_len(K - 1L)) {
    bounds[, j + 1L] <- pnorm(gamma[j] - eta)
  }
  bounds[, 2L:(K + 1L), drop = FALSE] - bounds[, 1L:K, drop = FALSE]
}

# A K-forest multinomial sampler over the configuration's design. The host
# gaussian sampler owns the data the wrapper borrows, so both are returned.
sbcMakeMultinomial <- function(config, labels, thin, seed) {
  host <- sbcMakeSampler(config, 1L, thin, seed)
  make <- getFromNamespace("bartcoreMultinomialSampler", "dbarts")
  list(host = host, mn = make(host, labels, config$K))
}

# The per-family operations the driver and the burn ladder share. `thin` is
# baked into the samplers the spec builds (the retained-draw spacing is a
# control setting, and the Student-t arm's pinned sampler is built once), so a
# spec is specific to one thinning.
sbcFamilySpec <- function(config, thin = 30L, seed = 20260709L) {
  storeState <- getFromNamespace("C_dbarts_bartcore_storeState", "dbarts")
  priorTrees <- getFromNamespace(
    "C_dbarts_bartcore_sampleTreesFromPrior",
    "dbarts"
  )
  priorNodes <- getFromNamespace(
    "C_dbarts_bartcore_sampleNodeParametersFromPrior",
    "dbarts"
  )
  bcRun <- bartcoreRun
  bcFits <- bartcoreForestFits
  K <- config$K

  if (config$family == "ordinal") {
    gen <- sbcMakeSampler(config, 1L, 1L, seed)
    drawGamma <- sbcOrdinalCutpointDraw(K)
    freeCuts <- seq_len(K - 2L) + 1L
    spec <- list(
      draw = function() {
        gen$sampleTreesFromPrior()
        gen$sampleNodeParametersFromPrior()
        eta0 <- as.numeric(gen$predict(config$x))
        eta0Test <- as.numeric(gen$predict(config$xTest))
        gamma0 <- drawGamma()
        p0 <- sbcOrdinalProbs(eta0, gamma0)
        y0 <- sbcCategoricalDraw(p0)
        theta <- c(
          setNames(gamma0[freeCuts], paste0("gamma", freeCuts)),
          avg.eta = mean(eta0),
          setNames(eta0Test, paste0("eta.star", seq_along(eta0Test))),
          setNames(colMeans(p0), paste0("p", seq_len(K)))
        )
        list(
          y = factor(y0, levels = seq_len(K), ordered = TRUE),
          theta = theta
        )
      },
      fit = function(y) {
        f <- sbcMakeSampler(config, 1L, thin, seed, y = y)
        f$sampleTreesFromPrior()
        f$sampleNodeParametersFromPrior()
        f
      },
      burnRun = function(f, burn) f$run(burn, 0L),
      sample = function(f) {
        res <- f$run(0L, 1L)
        gamma <- as.numeric(res$thresholds)
        eta <- res$train[, 1]
        p <- colMeans(sbcOrdinalProbs(eta, gamma))
        c(
          gamma[freeCuts],
          mean(eta),
          res$test[, 1],
          p
        )
      }
    )
  } else if (config$family == "nbinom") {
    gen <- sbcMakeSampler(config, 1L, 1L, seed)
    drawR <- sbcGridDraw(sbcNbGrid)
    spec <- list(
      draw = function() {
        gen$sampleTreesFromPrior()
        gen$sampleNodeParametersFromPrior()
        psi0 <- as.numeric(gen$predict(config$x))
        psi0Test <- as.numeric(gen$predict(config$xTest))
        r0 <- drawR(1L)
        # E[y | psi] = r exp(psi) under the engine's logit-p parameterization
        y0 <- rnbinom(config$n, size = r0, mu = r0 * exp(psi0))
        list(
          y = as.double(y0),
          theta = c(
            r = r0,
            avg.mu = mean(r0 * exp(psi0)),
            agg.psi = mean(psi0Test)
          )
        )
      },
      fit = function(y) {
        f <- sbcMakeSampler(config, 1L, thin, seed, y = y)
        f$sampleTreesFromPrior()
        f$sampleNodeParametersFromPrior()
        f
      },
      burnRun = function(f, burn) f$run(burn, 0L),
      sample = function(f) {
        res <- f$run(0L, 1L)
        r <- .Call(storeState, f$getPointer())[[1L]]$dispersion
        c(r, mean(r * exp(res$train[, 1])), mean(res$test[, 1]))
      }
    )
  } else if (config$family == "t") {
    sampler <- sbcMakeSampler(config, 1L, thin, seed)
    drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
    drawNu <- sbcGridDraw(sbcTGrid)
    spec <- list(
      draw = function() {
        sampler$sampleTreesFromPrior()
        sampler$sampleNodeParametersFromPrior()
        f0 <- as.numeric(sampler$predict(config$x))
        f0Test <- as.numeric(sampler$predict(config$xTest))
        sig0 <- drawSigma(1L)
        nu0 <- drawNu(1L)
        # r_i | lambda_i ~ N(0, sigma^2 / lambda_i), lambda_i ~ Gamma(nu/2,
        # nu/2) is exactly r_i = sigma * t_nu, so the mixture never needs
        # drawing: sigma is the CONDITIONAL scale the engine reports.
        list(
          y = f0 + sig0 * rt(config$n, nu0),
          theta = c(
            sigma = sig0,
            nu = nu0,
            avg.f = mean(f0),
            agg.f.star = mean(f0Test)
          )
        )
      },
      # ONE pinned sampler serves as generator and fit: rebuilding would
      # re-anchor the response scale the prior draw shares, and setResponse
      # (updateScale = FALSE) cold-inits nu and lambda, so the fresh prior draw
      # before it is a fully independent overdispersed start
      fit = function(y) {
        sampler$sampleTreesFromPrior()
        sampler$sampleNodeParametersFromPrior()
        sampler$setSigma(config$sigest)
        sampler$setResponse(y)
        sampler
      },
      burnRun = function(f, burn) f$run(burn, 0L),
      sample = function(f) {
        res <- f$run(0L, 1L)
        nu <- .Call(storeState, f$getPointer())[[1L]]$resid.df
        c(
          as.numeric(res$sigma)[1L],
          nu,
          mean(res$train[, 1]),
          mean(res$test[, 1])
        )
      }
    )
  } else if (config$family == "multinomial") {
    # eval points are the first nTest TRAINING rows: per-forest fits are exposed
    # for the training design only (bartcoreForestFits), and theta0's f_ik must
    # come from the SAME accessor the posterior draws do. The BCF arm's idx
    # convention exactly.
    idx <- seq_len(config$nTest)
    cells <- cbind(
      row = seq_len(min(3L, config$n)),
      cat = seq_len(min(3L, K))
    )
    cellNames <- paste0("f.", cells[, 1L], ".", cells[, 2L])
    buildLabels <- as.integer(rep_len(seq_len(K), config$n) - 1L)
    gen <- sbcMakeMultinomial(config, buildLabels, 1L, seed)
    forestFits <- function(handle) {
      vapply(
        seq_len(K),
        function(k) bcFits(handle$mn, k - 1L)[, 1],
        numeric(config$n)
      )
    }
    spec <- list(
      draw = function() {
        .Call(priorTrees, gen$mn$ptr)
        .Call(priorNodes, gen$mn$ptr)
        f0 <- forestFits(gen)
        p0 <- sbcSoftmax(f0)
        y0 <- sbcCategoricalDraw(p0)
        list(
          y = as.integer(y0 - 1L),
          theta = c(
            setNames(
              colMeans(p0[idx, , drop = FALSE]),
              paste0("p", seq_len(K))
            ),
            setNames(f0[cells], cellNames)
          )
        )
      },
      fit = function(y) {
        f <- sbcMakeMultinomial(config, y, thin, seed)
        .Call(priorTrees, f$mn$ptr)
        .Call(priorNodes, f$mn$ptr)
        f
      },
      burnRun = function(f, burn) bcRun(f$mn, burn, 0L),
      sample = function(f) {
        res <- bcRun(f$mn, 0L, 1L)
        probs <- array(res$train, c(config$n, K))
        c(colMeans(probs[idx, , drop = FALSE]), forestFits(f)[cells])
      }
    )
  } else {
    stop("no family spec for \"", config$family, "\"")
  }
  spec
}

# The configuration each family arm runs at, in one place so the burn ladder,
# the R=200 verdict run and the CI matrix cannot drift apart. Sizing notes:
# ordinal takes K = 4 because gamma_1 is pinned at 0 and only gamma_2..gamma_K-1
# are free, so K >= 4 is what makes the cutpoint block a real (multi-cutpoint)
# target; nbinom takes a TIGHTENED k = 8 (psi sd = node.scale/k = pi sqrt(3)/8
# ~ 0.68 rather than 2.7) because the Polya-Gamma draw loops sum(y_i + r) times
# per sweep and default-k psi draws are lognormal-tailed and unbudgetable - a
# tightened prior still validates NB; multinomial takes K = 3 forests.
sbcFamilyConfig <- function(family) {
  switch(
    family,
    ordinal = sbcConfig(family = "ordinal", numCategories = 4L, nTest = 3L),
    nbinom = sbcConfig(family = "nbinom", k = 8),
    t = sbcConfig(family = "t"),
    multinom = ,
    multinomial = sbcConfig(family = "multinomial", numCategories = 3L),
    stop("no family config for \"", family, "\"")
  )
}

# predict() vs the recorded latent channel at ONE state: theta0's latent (the
# ordinal eta, the nbinom psi, the Student-t f) is read with predict() while its
# posterior draws come from the run's train/test channels, so the two maps must
# agree exactly or the ranks compare different quantities.
sbcCheckLatentConsistency <- function(config, seed = 99L) {
  set.seed(seed)
  spec <- sbcFamilySpec(config, 1L, seed)
  drawn <- spec$draw()
  fit <- spec$fit(drawn$y)
  res <- fit$run(0L, 1L)
  maxDiff <- max(abs(res$train[, 1] - as.numeric(fit$predict(config$x))))
  maxDiffTest <- max(abs(res$test[, 1] - as.numeric(fit$predict(config$xTest))))
  list(
    maxDiff = maxDiff,
    maxDiffTest = maxDiffTest,
    pass = maxDiff < 1e-8 && maxDiffTest < 1e-8
  )
}

# The multinomial analogue: theta0's p_ik is the harness's softmax of the
# per-forest fits, while the posterior's p_ik rides the run's train channel, so
# those two maps must agree at one state (the GP/BCF fit-map precedent).
sbcCheckMultinomialProbs <- function(config, seed = 99L) {
  set.seed(seed)
  bcRun <- bartcoreRun
  bcFits <- bartcoreForestFits
  spec <- sbcFamilySpec(config, 1L, seed)
  drawn <- spec$draw()
  fit <- spec$fit(drawn$y)
  res <- bcRun(fit$mn, 0L, 1L)
  probs <- array(res$train, c(config$n, config$K))
  f <- vapply(
    seq_len(config$K),
    function(k) bcFits(fit$mn, k - 1L)[, 1],
    numeric(config$n)
  )
  maxDiff <- max(abs(sbcSoftmax(f) - probs))
  list(maxDiff = maxDiff, pass = maxDiff < 1e-10)
}

# The measured per-family burn floor, in absolute SWEEPS (plan step 2: 72000 was
# a BCF-specific number, so every arm re-measures). Read off sbcBurnLadder at
# 40000 sweeps x 3 datasets; the numbers recorded in
# docs/plans/sbc-family-tiers.md. The two categorical/count arms are set by a
# LIKELIHOOD RIDGE, not by a transient: ordinal's free cutpoints trade against
# the mean level (docs/design/ordinal.md section 9's f-vs-cutpoint-shift ridge -
# gamma2/gamma3 and the p2 that reads them stay autocorrelated past lag 200,
# while every eta functional clears 0.1 by lag ~16), and nbinom's r trades
# against the psi level because only mu = r exp(psi) is identified (r and
# agg.psi mirror each other block for block; avg.mu clears 0.1 at LAG 1). The
# Student-t settles in a couple of thousand sweeps with sigma/nu at lag ~40-60,
# and multinomial mixes fastest of all (every functional under lag 10).
sbcBurnSweeps <- c(
  ordinal = 36000,
  nbinom = 24000,
  t = 12000,
  multinomial = 6000
)

# Rank R replications of a family-spec configuration. The generic sibling of
# runSbc: same result shape, with every functional routed through the
# tie-breaking rank (these families' functionals carry atoms - see
# sbcDiscreteRank - while runSbc's gaussian/binary ones do not).
runSbcFamily <- function(
  config,
  R = 200L,
  L = 150L,
  thin = 30L,
  burnSweeps = sbcBurnSweeps[[config$family]],
  seed = 20260709L,
  report = 25L
) {
  burn <- as.integer(ceiling(burnSweeps / thin))
  set.seed(seed)
  spec <- sbcFamilySpec(config, thin, seed)
  ranks <- NULL
  started <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    drawn <- spec$draw()
    theta <- drawn$theta
    fit <- spec$fit(drawn$y)
    invisible(spec$burnRun(fit, burn))
    draws <- matrix(NA_real_, length(theta), L)
    for (l in seq_len(L)) {
      draws[, l] <- spec$sample(fit)
    }
    row <- integer(length(theta))
    names(row) <- names(theta)
    for (j in seq_along(theta)) {
      # every functional is ranked with the tie-break: it reduces to
      # #{draws < theta0} (and consumes no rng) unless the law has an atom, and
      # both a grid parameter and an underflowed tail probability do
      row[j] <- sbcDiscreteRank(draws[j, ], theta[[j]])
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
    burnSweeps = burnSweeps,
    R = R,
    config = config,
    elapsed = elapsed,
    perRep = elapsed / R
  )
}

# --- burn ladder (plan step 2) ---------------------------------------------

# Measure, per family, how long the chain takes to forget an overdispersed
# start and how fast it then mixes - the two numbers that set burn and thin.
# For each of nDataset prior-drawn datasets it runs nSweep UNTHINNED sweeps from
# an independent prior init and records every functional's trace, then reports
# (a) each block's mean as a z-score against the final block (the transient: the
# block where |z| stops exceeding ~1 is where the chain has settled), and (b)
# the first ACF lag under 0.1 on the trailing half (the thinning floor). It also
# times the sweeps, which is the per-sweep cost measurement the budget needs.
sbcBurnLadder <- function(
  config,
  nSweep = 20000L,
  nDataset = 3L,
  nBlock = 10L,
  seed = 20260804L
) {
  set.seed(seed)
  spec <- sbcFamilySpec(config, 1L, seed)
  blockSize <- nSweep %/% nBlock
  results <- vector("list", nDataset)
  totalElapsed <- 0
  for (d in seq_len(nDataset)) {
    drawn <- spec$draw()
    fit <- spec$fit(drawn$y)
    trace <- matrix(NA_real_, length(drawn$theta), nSweep)
    started <- proc.time()[["elapsed"]]
    for (s in seq_len(nSweep)) {
      trace[, s] <- spec$sample(fit)
    }
    totalElapsed <- totalElapsed + proc.time()[["elapsed"]] - started
    settled <- trace[, (nSweep %/% 2L + 1L):nSweep, drop = FALSE]
    z <- matrix(NA_real_, length(drawn$theta), nBlock)
    firstUnder <- integer(length(drawn$theta))
    for (j in seq_along(drawn$theta)) {
      scale <- sd(settled[j, ])
      if (!is.finite(scale) || scale <= 0) {
        scale <- 1
      }
      for (b in seq_len(nBlock)) {
        block <- trace[j, ((b - 1L) * blockSize + 1L):(b * blockSize)]
        z[j, b] <- (mean(block) - mean(settled[j, ])) /
          (scale / sqrt(blockSize))
      }
      a <- acf(settled[j, ], lag.max = 200L, plot = FALSE)$acf[,, 1]
      hit <- which(a < 0.1)[1L]
      firstUnder[j] <- if (is.na(hit)) NA_integer_ else hit - 1L
    }
    rownames(z) <- names(drawn$theta)
    names(firstUnder) <- names(drawn$theta)
    results[[d]] <- list(z = z, firstUnder = firstUnder, blockSize = blockSize)
  }
  list(
    family = config$family,
    nSweep = nSweep,
    nBlock = nBlock,
    datasets = results,
    elapsed = totalElapsed,
    perSweep = totalElapsed / (nDataset * nSweep)
  )
}

sbcReportBurnLadder <- function(ladder) {
  cat(sprintf(
    "\nburn ladder: family=%s  %d sweeps x %d datasets  %.1fs  %.1f us/sweep\n",
    ladder$family,
    ladder$nSweep,
    length(ladder$datasets),
    ladder$elapsed,
    1e6 * ladder$perSweep
  ))
  for (d in seq_along(ladder$datasets)) {
    res <- ladder$datasets[[d]]
    cat(sprintf(
      "\n dataset %d: block-mean z vs the final half (block = %d sweeps)\n",
      d,
      res$blockSize
    ))
    cat(sprintf(
      "  %-12s %s  acf<0.1\n",
      "functional",
      paste(sprintf("%6d", seq_len(ncol(res$z))), collapse = "")
    ))
    for (j in seq_len(nrow(res$z))) {
      cat(sprintf(
        "  %-12s %s  %6s\n",
        rownames(res$z)[j],
        paste(sprintf("%6.1f", res$z[j, ]), collapse = ""),
        format(res$firstUnder[j])
      ))
    }
  }
  invisible(NULL)
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

  # ecdf-difference simultaneous band via simulation of the null. The ecdf of
  # integer ranks is a cumulated tabulation, which is what makes a Bonferroni'd
  # alpha affordable: the band is a 1 - alpha quantile, so a small alpha needs
  # many more null draws to place stably (>= 20 in the tail), and the same RNG
  # calls in the same order keep every previously recorded band bit-identical.
  target <- seq_len(L + 1L) / (L + 1)
  ecdfDiff <- function(rk) {
    cumsum(tabulate(rk + 1L, L + 1L)) / length(rk) - target
  }
  observed <- max(abs(ecdfDiff(ranks)))
  nSim <- max(nSim, ceiling(20 / alpha))
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

# The CI matrix's admission level (plan step 4, replacing sbc-ci-gate's step 4):
# the ecdf band's alpha Bonferroni'd over the matrix's TOTAL functional count,
# so a full-matrix pass has probability ~0.95 on a fresh stream rather than each
# arm alarming independently at its own nominal 5%. M is
# gaussian 7 + ordinal 10 + nbinom 3 + t 4 + multinomial 6.
sbcMatrixConfigs <- c(
  "gaussian",
  "ordinal",
  "nbinom",
  "t",
  "multinom",
  "multinomial"
)
sbcMatrixFunctionals <- 7L + 10L + 3L + 4L + 6L
sbcMatrixAlpha <- 0.05 / sbcMatrixFunctionals

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
# alpha is the ecdf band's level; the CI matrix Bonferroni's it (see
# sbcMatrixAlpha) so that a whole matrix of arms passes with probability ~0.95
# rather than each arm alarming at its own nominal 5%.
sbcReport <- function(fit, nBins = 20L, alpha = 0.05) {
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
    "wall-clock: %.1fs total, %.3fs/rep; band alpha = %.5f\n\n",
    fit$elapsed,
    fit$perRep,
    alpha
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
    u <- rankUniformity(
      fit$ranks[, funcs[i]],
      fit$L,
      nBins = nBins,
      alpha = alpha
    )
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
  # optional 5th arg: the burn in absolute SWEEPS, which otherwise comes from
  # the measured sbcBurnSweeps (family tiers) or the driver default. It exists
  # so the chain-length diagnostic ladder (the A4e protocol: re-run a flagged
  # arm at several thin/burn points and see whether the bias SHRINKS into the
  # band or plateaus) is a recordable command rather than a scratch script, and
  # so an arm recorded at a shorter burn than today's default is reproducible
  # from the command line.
  burnSweeps <- if (length(args) >= 5L) as.numeric(args[5]) else NULL
  # optional 6th arg: the driver seed, so the A4e adjudication's replication
  # step can draw a FRESH stream at settings otherwise held fixed. Absent, the
  # driver's own pinned seed keeps every recorded run reproducible.
  runSeed <- if (length(args) >= 6L) as.integer(args[6]) else NULL

  # Step-1 self-check mode: the discrete rank against a closed-form conjugate
  # posterior. No engine involved, so it runs in seconds and gates the two grid
  # functionals (nbinom r, Student-t nu) before either arm is trusted.
  if (which == "discrete-selfcheck") {
    cat("== discrete-rank self-check (closed-form conjugate posterior) ==\n")
    for (grid in list(sbcNbGrid, sbcTGrid)) {
      g <- sbcCheckGridPrior(grid)
      cat(sprintf(
        "  grid prior (%d cells, max %g): max cell diff %.5f; mean %.4f vs %.4f -> %s\n",
        length(grid),
        max(grid),
        g$maxCellDiff,
        g$meanEmpirical,
        g$meanTheory,
        if (g$pass) "PASS" else "FAIL"
      ))
      if (!isTRUE(g$pass)) {
        stop("grid prior moment check failed")
      }
    }
    chk <- sbcDiscreteSelfCheck(R = if (length(args) >= 2L) R else 400L, L = L)
    u <- chk$uniformity
    cat(sprintf(
      "\n  ranks: mean %.1f (target %.1f); tied draws %.3f of L\n",
      u$mean,
      u$meanTarget,
      chk$tieFrac
    ))
    cat(sprintf(
      "  chisqP %.3f  ksP %.3f  ecdfDiff %.4f  band %.4f -> %s\n",
      u$chisqP,
      u$ksP,
      u$ecdfDiff,
      u$ecdfBand,
      if (chk$pass) "PASS" else "FLAG"
    ))
    cat("\n", sbcAsciiHistogram(chk$ranks, L), "\n", sep = "")
    if (nzchar(Sys.getenv("SBC_FAIL_ON_FLAG", "")) && !chk$pass) {
      quit(status = 1L, save = "no")
    }
    quit(status = 0L, save = "no")
  }

  # Step-2 burn-ladder mode: "burn-<family>" measures the transient and the
  # per-sweep cost instead of ranking. Positional args become nSweep, nDataset.
  if (startsWith(which, "burn-")) {
    family <- sub("^burn-", "", which)
    nSweep <- if (length(args) >= 2L) as.integer(args[2]) else 20000L
    nDataset <- if (length(args) >= 3L) as.integer(args[3]) else 3L
    sbcReportBurnLadder(sbcBurnLadder(
      sbcFamilyConfig(family),
      nSweep = nSweep,
      nDataset = nDataset
    ))
    quit(status = 0L, save = "no")
  }

  isDart <- which %in% c("dart", "dart-sparse")
  isWeighted <- which == "weighted"
  # "-swap" runs the same configuration through the response-swap arm; it is a
  # one-time local adjudication of the bridge's grouped setResponse and is
  # deliberately NOT in sbc.yaml's matrix
  isGrouped <- which %in%
    c(
      "grouped-gaussian",
      "grouped-probit",
      "grouped-gaussian-swap",
      "grouped-probit-swap"
    )
  isBCF <- which %in% c("bcf", "bcf-weak")
  isAft <- which == "aft"
  isLinear <- which %in%
    c("linear", "linear-na-leaf", "linear-na-split", "linear-weighted")
  isGP <- which %in% c("gp", "gp-na-leaf", "gp-weighted", "gp-mixed")
  isFamilyTier <- which %in%
    c("ordinal", "nbinom", "t", "multinom", "multinomial")

  config <- if (isFamilyTier) {
    sbcFamilyConfig(which)
  } else if (isDart) {
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
    base <- if (startsWith(which, "grouped-probit")) "probit" else "gaussian"
    zw <- if (startsWith(which, "grouped-gaussian")) 0.2 else 0
    cfg <- sbcAddGrouping(
      sbcConfig(family = base, n = 160L),
      nGroups = 8L,
      relScale = 0.2,
      zeroWeightFrac = zw
    )
    # pin: score every rebuild against the anchor's own leaf prior instead of
    # each replication's own y0 range (retires the sigma-functional FLAG)
    cfg$nodePrior <- dbartsPriors$normal(cfg$k, scale = sbcAnchorScale(cfg))
    cfg
  } else if (isBCF) {
    # prior-weak = small n so the a-glue prior term dominates the likelihood
    nBcf <- if (which == "bcf-weak") 40L else 200L
    sbcAddBCF(sbcConfig(family = "gaussian", n = nBcf))
  } else if (isAft) {
    sbcConfigAft()
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
    # "gp-mixed" caps at 30 instead, the MEASURED median leaf size of this
    # config's own prior draws (600 draws x 25 trees: median 29, 2.42 leaves
    # per tree), so ~49% of leaves exceed the cap and fall back to constant
    # fits and ~79% of trees carry both leaf kinds at once -- the only arm
    # that exercises the mixed path inside a single tree.
    cfg <- sbcConfig(
      family = "gaussian",
      n = 80L,
      nTrees = 25L,
      nodePrior = dbartsPriors$gp(
        1L,
        k = 2,
        max.leaf.size = if (which == "gp-mixed") 30L else 100L
      )
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
  selfCheckPass <- c(sigma = isTRUE(chk$pass))
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
    selfCheckPass["tau"] <- isTRUE(tc$pass)
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
    selfCheckPass["glue"] <- isTRUE(gc$pass)
  }
  if (isFamilyTier) {
    if (config$family == "ordinal") {
      oc <- sbcCheckOrdinalCutpointPrior()
      cat(sprintf(
        "  log-gap: sd %.4f vs %.4f; gap median %.4f vs %.4f -> %s\n",
        oc$sdEmpirical,
        oc$sdTheory,
        oc$medianEmpirical,
        oc$medianTheory,
        if (oc$pass) "PASS" else "FAIL"
      ))
      selfCheckPass["cutpoints"] <- isTRUE(oc$pass)
    }
    if (config$family %in% c("nbinom", "t")) {
      grid <- if (config$family == "nbinom") sbcNbGrid else sbcTGrid
      gp <- sbcCheckGridPrior(grid)
      cat(sprintf(
        "  grid prior: max cell diff %.5f; mean %.4f vs %.4f -> %s\n",
        gp$maxCellDiff,
        gp$meanEmpirical,
        gp$meanTheory,
        if (gp$pass) "PASS" else "FAIL"
      ))
      selfCheckPass["grid"] <- isTRUE(gp$pass)
    }
    if (config$family == "multinomial") {
      mc <- sbcCheckMultinomialProbs(config)
      cat(sprintf(
        "  softmax(forest fits) vs reported probabilities: %.2e -> %s\n",
        mc$maxDiff,
        if (mc$pass) "PASS" else "FAIL"
      ))
      selfCheckPass["softmax"] <- isTRUE(mc$pass)
    } else {
      lc <- sbcCheckLatentConsistency(config)
      cat(sprintf(
        "  predict vs recorded latent: train %.2e, test %.2e -> %s\n",
        lc$maxDiff,
        lc$maxDiffTest,
        if (lc$pass) "PASS" else "FAIL"
      ))
      selfCheckPass["latent"] <- isTRUE(lc$pass)
    }
  }
  if (isLinear || isGP) {
    fc <- sbcCheckFitConsistency(config)
    cat(sprintf(
      "  predict vs recorded fits: train %.2e, test %.2e -> %s\n",
      fc$maxDiff,
      fc$maxDiffTest,
      if (fc$pass) "PASS" else "FAIL"
    ))
    selfCheckPass["fit"] <- isTRUE(fc$pass)
  }
  if (isAft) {
    ac <- sbcCheckAftFitConsistency(config)
    cat(sprintf(
      "  predict vs recorded fits: train %.2e, test %.2e -> %s\n",
      ac$maxDiff,
      ac$maxDiffTest,
      if (ac$pass) "PASS" else "FAIL"
    ))
    selfCheckPass["fit"] <- isTRUE(ac$pass)
  }

  # Harness integrity: a failed self-check means the prior/fit reference is
  # miscalibrated, so the SBC result would be meaningless (or falsely clean).
  # Abort unconditionally - unlike the functional FLAG gate below, this is not
  # opt-in behind SBC_FAIL_ON_FLAG.
  if (any(!selfCheckPass)) {
    stop(
      "SBC harness self-check failed (",
      paste(names(selfCheckPass)[!selfCheckPass], collapse = ", "),
      "): reference is miscalibrated; SBC results are invalid."
    )
  }

  cat(sprintf("\n== SBC run (%s R=%d L=%d thin=%d) ==\n", which, R, L, thin))
  fit <- if (isFamilyTier) {
    if (is.null(burnSweeps)) {
      runSbcFamily(config, R = R, L = L, thin = thin)
    } else {
      runSbcFamily(config, R = R, L = L, thin = thin, burnSweeps = burnSweeps)
    }
  } else if (isDart) {
    runSbcDart(config, R = R, L = L, thin = thin)
  } else if (isGrouped) {
    runSbcGrouped(
      config,
      R = R,
      L = L,
      thin = thin,
      swap = endsWith(which, "-swap")
    )
  } else if (isBCF) {
    runSbcBCF(config, R = R, L = L, thin = thin)
  } else if (isAft) {
    runSbcAft(config, R = R, L = L, thin = thin)
  } else {
    plainArgs <- list(config, R = R, L = L, thin = thin)
    if (!is.null(burnSweeps)) {
      plainArgs$burn <- as.integer(ceiling(burnSweeps / thin))
    }
    if (!is.null(runSeed)) {
      plainArgs$seed <- runSeed
    }
    do.call(runSbc, plainArgs)
  }
  if (isDart) {
    cat(sprintf(
      "\nfloor incidence: %.3f of s0 components pinned at 1e-300\n",
      fit$floorFrac
    ))
  }
  if (isAft) {
    cat(sprintf("\ncensoring rate: %.3f\n", fit$censoringRate))
  }
  # matrix arms are admitted at the Bonferroni'd level; every other config keeps
  # the per-functional 5% band its recorded result was read at
  verdicts <- sbcReport(
    fit,
    alpha = if (which %in% sbcMatrixConfigs) sbcMatrixAlpha else 0.05
  )
  if (nzchar(Sys.getenv("SBC_FAIL_ON_FLAG", "")) && any(verdicts == "FLAG")) {
    quit(status = 1L, save = "no")
  }
}
