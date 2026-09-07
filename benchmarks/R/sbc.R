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
#   Rscript benchmarks/R/sbc.R aft 200 150 30 <burn> # aft/survival, reused
#   Rscript benchmarks/R/sbc.R gp-mixed 200 150 60 3000 # GP/constant mix
#   Rscript benchmarks/R/sbc.R bcf-probit 200 150 30 <burn> # latent BCF arms
#   Rscript benchmarks/R/sbc.R discrete-selfcheck   # the discrete-rank gate
#   Rscript benchmarks/R/sbc.R burn-ordinal 20000 3 # the burn/cost ladder
#   Rscript benchmarks/R/sbc.R burn-bcf-probit 40000 24 # its repriced ladder
#   Rscript benchmarks/R/sbc.R burn-aft 20000 3     # the aft arm's own
# Positional args: config R L thin, plus an optional 5th, the burn in absolute
# sweeps, and an optional 6th, the driver seed. Or source() the file to reuse
# the API:
#   source("benchmarks/R/sbc.R"); res <- runSbc(sbcConfig("gaussian"), R = 200)
# SBC_FAIL_ON_FLAG=1 (env var, opt-in) makes the CLI exit status 1 if any
# functional's verdict is FLAG; unset, the CLI always exits 0 (unchanged
# default). Only affects Rscript use; source() usage is untouched.
# SBC_EXPECTED_FLAGS (env var, opt-in) is a comma-separated list of functional
# names allowed to FLAG without failing SBC_FAIL_ON_FLAG -- e.g. nbinom's r
# and agg.psi trade off on an adjudicated identifiability ridge, not a defect
# (docs/plans/sbc-family-tiers.md Step 3; confirmed 2026-08-18). A listed
# functional that FLAGs prints "FLAG (expected)" and is excluded from the exit
# check; one that PASSES is unaffected; any FLAG outside the list still fails
# as before. Only the CLI reads this variable; source() usage is untouched.
# SBC_POISON (env var, opt-in) names one or more deliberate generator/sampler
# mismatches for the two latent BCF arms, comma-separated (sbcBCFPoisons):
# "link" simulates through the OTHER link, "glue-sd" draws the glue at
# gaussian's sd.control = 2 while the sampler runs at the family default 1, and
# "sigma" adds latent noise the fit cannot model. Each is a by-hand
# discrimination run whose functionals must FLAG, never a recorded verdict; an
# unknown name refuses the run. SBC_FIXED_GLUE (env var, opt-in) holds a BCF
# arm's glue at the engine's initial (1, 0, 1) - the control that isolates the
# two-forest backfit from the glue draw, and the one "glue-sd" carries.
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
# fit ignores.
sbcSimulate <- function(config, f0Train, sig0) {
  mu <- f0Train
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
# configurations (linear/gp leaf, BCF, DART, weighted) extend the same
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

# --- aft (accelerated failure time / survival) ------------------------------

# The engine's aft sampler (docs/design/survival.md) is a log-normal
# survival model: log T = f(x) + sigma*eps, uncensored rows enter as gaussian
# data on log T and right-censored ones contribute the upper normal tail past
# their log censoring time. ONE sampler serves every replication, the shape
# every other reused arm has: $setResponse(y0, status = status0) installs the
# observed log times and the censoring structure they imply in a single call,
# so a drawn theta0 censors whatever rows it likes without a rebuild.
#
# That retires the two pins the rebuilt arm carried, both of which existed
# only because a rebuild re-derived the response transform from range(y0): an
# anchor leaf scale named as node.prior held the leaf prior against it, and an
# offset zeroing each rebuild's prior.mean held the shift. One build fixes the
# transform once -- off the symmetric build response, response.shift 0 and
# response.scale 5, so prior.mean is already 0 and prior.scale already what
# the anchor named -- and updateScale = FALSE keeps it, so the prior draw and
# the posterior share one transform with nothing left to pin.
#
# sigma is drawn conjugately exactly as gaussian's (chain.hpp), so it is this
# arm's log-time scale functional, ranked the same way avg.f and f.star are.
# Two functionals are this family's alone. S(t0 | x*) = 1 - Phi((log t0 -
# f(x*)) / sigma) is the reported survival deliverable, at the t0 pinned below
# and the first test point. And logT0[i] at the LOWEST-INDEXED censored row is
# ranked against that row's own posterior latents, read with $getLatents one
# retained sample at a time: it is the only functional that ranks the
# truncated-normal imputation, every other one reading a channel the
# imputation moves only through the fit. The censored set is the replication's
# own draw and can be empty; such a replication contributes NO rank there (an
# NA the driver leaves alone) and that functional's own R is reported
# separately.

# Fixed per-row right-censoring log-times: a design choice independent of any
# prior draw, pinned once so only the drawn
# log-time and the status it implies move across replications.
sbcAddCensoring <- function(config, shift = 1.0, sd = 1.2) {
  set.seed(505L)
  config$logC <- config$yBuild + shift + sd * rnorm(config$n)
  config
}

# The arm's fixed config: a continuous, symmetric-about-0 build response,
# m = 50 trees, the censoring fixture above, and the survival deliverable's
# own time. t0 is the build response's median survival time exp(0) = 1, the
# transform's centre, so (log t0 - f(x*)) / sigma is centred at 0 under the
# prior draw and S(t0 | x*) spreads over (0, 1) rather than piling into a tail
# it would eventually underflow to an atom in.
sbcConfigAft <- function(n = 150L, nTrees = 50L) {
  config <- sbcConfig(family = "gaussian", n = n, nTrees = nTrees)
  config <- sbcAddCensoring(config)
  config$t0 <- 1
  config$family <- "aft"
  config
}

# The arm's one sampler, generator and fit both. An aft response is a
# (time, status) pair, so the build times are exp(yBuild) -- the engine fits
# the log times, which is what the transform is derived from -- and the build
# status is all events, the structure the first $setResponse replaces.
sbcMakeAftSampler <- function(config, thin) {
  ctrl <- dbartsControl(
    n.trees = config$nTrees,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 1L,
    n.thin = thin,
    updateState = FALSE,
    verbose = FALSE,
    keepTrainingFits = TRUE
  )
  dbarts(
    config$x,
    cbind(exp(config$yBuild), rep_len(1, config$n)),
    test = config$xTest,
    resid.prior = dbartsPriors$chisq(config$sigDf, config$sigQuant),
    node.prior = config$nodePrior,
    sigma = config$sigest,
    control = ctrl,
    family = "aft"
  )
}

# S(t | x) under the log-normal model, the arm's reported deliverable.
sbcAftSurvival <- function(t0, f, sigma) {
  1 - pnorm((log(t0) - f) / sigma)
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
#
# The two LATENT arms (bcf-probit, bcf-logistic) share the machinery below and
# differ in four places, one per thing a latent family makes new. Sigma is
# pinned at exactly 1, so its functional and its moment check go and
# sbcCheckBCFLatent replaces them; the reported/internal transform is the
# identity, so the affine map sbcMakeBCF regresses is a self-check rather than
# a conversion; y is Bernoulli at the link of the COMBINED index a mu + b_z tau,
# the location the latent refresh runs against, with no offset and no noise;
# and p_j, the link at each evaluation row's index, joins the ranked
# functionals - the reported deliverable, bounded, and rank-equivalent to the
# index itself because the link is increasing, where neither prog_j nor eff_j
# alone is the index. Thirteen at nTest = 3: 4 glue, 3 prog, 3 eff, 3 p. Their
# driver is runSbcFamily, off the sbcFamilySpec branch below, so one generator
# serves both the burn ladder and the R-replication run.

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

# A BCF config whose family is one of the two latent links, as against the
# gaussian arm's.
sbcBCFLatent <- function(config) {
  config$family %in% c("probit", "logistic")
}

# The link a latent arm simulates and reports through, and its opposite - the
# generator half of poison (i).
sbcBCFLink <- function(family) {
  switch(
    family,
    probit = pnorm,
    logistic = plogis,
    stop("no link for family \"", family, "\"")
  )
}

sbcBCFOtherLink <- function(family) {
  sbcBCFLink(switch(family, probit = "logistic", logistic = "probit"))
}

# The named discrimination poisons (docs/plans/bcf-latent-evidence.md Decision
# 4), each a deliberate mismatch between the generator and the sampler that
# must redden the arm. They are run once by hand and never recorded as a
# verdict, so the names are validated where they are read: a typo would
# otherwise score a clean arm and read as a pass.
sbcBCFPoisons <- c("link", "glue-sd", "sigma")

sbcBCFPoison <- function(poison) {
  if (is.null(poison)) {
    return(character(0))
  }
  poison <- trimws(poison[nzchar(trimws(poison))])
  unknown <- setdiff(poison, sbcBCFPoisons)
  if (length(unknown) > 0L) {
    stop(
      "unknown SBC poison(s): ",
      paste(unknown, collapse = ", "),
      "; known: ",
      paste(sbcBCFPoisons, collapse = ", ")
    )
  }
  poison
}

# The latent arm's ranked functionals, from one (glue, mu, tau) - theta0's or a
# posterior draw's, which is what makes the two comparable. Evaluation rows are
# the first nTest TRAINING rows, the arm's idx convention, so z is theirs too.
# Under a held glue the four glue functionals are degenerate constants and are
# dropped rather than ranked.
sbcBCFFunctionals <- function(
  config,
  glue,
  mu,
  tau,
  idx,
  link,
  fixedGlue = FALSE
) {
  diff <- glue$b1 - glue$b0
  out <- if (fixedGlue) {
    numeric(0)
  } else {
    c(
      a = glue$a,
      abs.a = abs(glue$a),
      b1.minus.b0 = diff,
      abs.diff = abs(diff)
    )
  }
  bz <- ifelse(config$z[idx] != 0, glue$b1, glue$b0)
  index <- glue$a * mu[idx] + bz * tau[idx]
  for (j in seq_along(idx)) {
    out[paste0("prog", j)] <- glue$a * mu[idx[j]]
    out[paste0("eff", j)] <- diff * tau[idx[j]]
    out[paste0("p", j)] <- link(index[j])
  }
  out
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
    # the family formal writes the link into the model copy the bridge reads;
    # NULL, the default, leaves the host gaussian sampler's own
    family = if (sbcBCFLatent(config)) config$family else NULL,
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

# What replaces the gaussian arm's sigma moment check, and the reason the
# regressed map is a self-check here. Three claims, read at the prior state and
# again after a response swap: sigma is pinned at EXACTLY 1, the
# reported/internal transform is the identity, and the recorded combined train
# fits are a mu + b_z tau - the location the latent refresh runs against, which
# is neither forest's own fits.
sbcCheckBCFLatent <- function(config, seed = 99L) {
  set.seed(seed)
  built <- sbcMakeBCF(config, 1L, 1L, fixedGlue = isTRUE(config$fixedGlue))
  bcf <- built$bcf
  combined <- function() {
    res <- .bcfRun(bcf, 0L, 1L)
    glue <- .bcfGlue(bcf)
    bz <- ifelse(config$z != 0, glue[3L], glue[2L])
    mu <- .bcfForest(bcf, 0L)[, 1]
    tau <- .bcfForest(bcf, 1L)[, 1]
    list(
      fits = res$train[, 1],
      index = glue[1L] * mu + bz * tau,
      sigma = as.numeric(res$sigma)[1L]
    )
  }
  prior <- combined()
  y0 <- as.double(rbinom(
    config$n,
    1L,
    sbcBCFLink(config$family)(prior$index)
  ))
  .bcfSetResponse(bcf, y0, FALSE)
  fitted <- combined()
  maxDiff <- max(
    abs(prior$fits - prior$index),
    abs(fitted$fits - fitted$index)
  )
  maxSigma <- max(abs(c(prior$sigma, fitted$sigma) - 1))
  list(
    fitScale = built$fitScale,
    fitShift = built$fitShift,
    mapR2 = built$mapR2,
    maxDiff = maxDiff,
    maxSigma = maxSigma,
    pass = maxDiff < 1e-12 &&
      maxSigma == 0 &&
      abs(built$fitScale - 1) < 1e-10 &&
      abs(built$fitShift) < 1e-10 &&
      abs(built$mapR2 - 1) < 1e-12
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

  if (!is.null(config$z)) {
    # A latent BCF arm - a config carrying a treatment vector is a BCF one.
    # One handle serves as generator and fit, as the gaussian arm's does: numGroups == 0 makes setResponse legal, so the build scale the
    # prior draw shares is never disturbed. Everything the arm needs beyond
    # sbcAddBCF rides the config - `fixedGlue` holds the glue at the engine's
    # initial (1, 0, 1), `poison` names a deliberate mismatch - so the burn
    # ladder and the R-replication driver share one generator.
    if (!sbcBCFLatent(config)) {
      stop("the gaussian BCF arm's driver is runSbcBCF, not the family one")
    }
    poison <- sbcBCFPoison(config$poison)
    fixedGlue <- isTRUE(config$fixedGlue)
    link <- sbcBCFLink(config$family)
    # poison (i): the generator's link, wrong on purpose
    simLink <- if ("link" %in% poison) {
      sbcBCFOtherLink(config$family)
    } else {
      link
    }
    # poison (ii): the generator's a-prior scale at gaussian's 2 while the
    # sampler runs at the family default. Inert under a held glue, which is the
    # control this poison carries and the other two do not
    glueScale <- if ("glue-sd" %in% poison) 2 else config$sdControl
    drawGlue <- if (fixedGlue) {
      function() list(a = 1, b0 = 0, b1 = 1) # the engine's fixed initial glue
    } else {
      sbcBCFGlueDraw(glueScale, config$bPriorVariance)
    }
    # poison (iii): latent noise the fit cannot model, at the gaussian arm's
    # own sigma prior
    drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
    built <- sbcMakeBCF(config, 1L, thin, fixedGlue = fixedGlue)
    bcf <- built$bcf
    idx <- seq_len(config$nTest)
    # `a` prescribes the prognostic scalar's MAGNITUDE (the ladder's strata);
    # its sign is unidentified, so the positive representative is drawn
    drawOne <- function(a = NULL) {
      g0 <- drawGlue()
      if (!is.null(a)) {
        g0$a <- a
      }
      # the glue is installed BEFORE the forests: each forest's prior trees are
      # drawn against its own veto vector, which the glue sets
      sbcInstallBCFGlue(bcf, g0)
      .Call(priorTrees, bcf$ptr)
      .Call(priorNodes, bcf$ptr)
      mu0 <- bcFits(bcf, 0L)[, 1]
      tau0 <- bcFits(bcf, 1L)[, 1]
      bz0 <- ifelse(config$z != 0, g0$b1, g0$b0)
      index0 <- g0$a * mu0 + bz0 * tau0
      if ("sigma" %in% poison) {
        index0 <- index0 + drawSigma(1L) * rnorm(config$n)
      }
      list(
        y = as.double(rbinom(config$n, 1L, simLink(index0))),
        theta = sbcBCFFunctionals(config, g0, mu0, tau0, idx, link, fixedGlue)
      )
    }
    spec <- list(
      draw = drawOne,
      drawAt = drawOne,
      fit = function(y) {
        .Call(priorTrees, bcf$ptr)
        .Call(priorNodes, bcf$ptr)
        .bcfSetResponse(bcf, y, FALSE)
        bcf
      },
      burnRun = function(f, burn) bcRun(f, burn, 0L),
      sample = function(f) {
        bcRun(f, 0L, 1L)
        glue <- .bcfGlue(f)
        sbcBCFFunctionals(
          config,
          list(a = glue[1L], b0 = glue[2L], b1 = glue[3L]),
          bcFits(f, 0L)[, 1],
          bcFits(f, 1L)[, 1],
          idx,
          link,
          fixedGlue
        )
      }
    )
  } else if (config$family == "ordinal") {
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
  } else if (config$family == "aft") {
    # ONE pinned sampler serves as generator and fit, the Student-t arm's
    # shape: $setResponse(y, status = ) replaces the censoring structure along
    # with the response, so nothing is rebuilt and updateScale = FALSE keeps
    # the build transform the prior draw shares. The fresh prior draw before
    # it is the overdispersed start.
    sampler <- sbcMakeAftSampler(config, thin)
    drawSigma <- sbcSigmaDraw(config$sigest, config$sigDf, config$sigQuant)
    # the row the censored-latent functional reads: closure state, since only
    # draw() knows the replication's status and only sample() reads latents.
    # NA_integer_ when nothing censored, which is the no-rank case
    censoredRow <- NA_integer_
    spec <- list(
      draw = function() {
        sampler$sampleTreesFromPrior()
        sampler$sampleNodeParametersFromPrior()
        f0 <- as.numeric(sampler$predict(config$x))
        f0Test <- as.numeric(sampler$predict(config$xTest))
        sig0 <- drawSigma(1L)
        # the latent log times, then what the pinned censoring fixture leaves
        # observable of them: an event row reports its own log time, a
        # censored one reports its censoring time and the fact it was passed
        logT0 <- f0 + sig0 * rnorm(config$n)
        status0 <- as.numeric(logT0 <= config$logC)
        y0 <- pmin(logT0, config$logC)
        censored <- which(status0 == 0)
        censoredRow <<- if (length(censored) > 0L) {
          censored[1L]
        } else {
          NA_integer_
        }
        list(
          y = cbind(y0, status0),
          theta = c(
            sigma = sig0,
            avg.f = mean(f0),
            setNames(f0Test, paste0("f.star", seq_along(f0Test))),
            S.star1 = sbcAftSurvival(config$t0, f0Test[1L], sig0),
            logT.cens = if (is.na(censoredRow)) {
              NA_real_
            } else {
              logT0[censoredRow]
            }
          )
        )
      },
      fit = function(y) {
        sampler$sampleTreesFromPrior()
        sampler$sampleNodeParametersFromPrior()
        sampler$setSigma(config$sigest)
        sampler$setResponse(y[, 1L], updateScale = FALSE, status = y[, 2L])
        sampler
      },
      burnRun = function(f, burn) f$run(burn, 0L),
      sample = function(f) {
        res <- f$run(0L, 1L)
        sigma <- as.numeric(res$sigma)[1L]
        c(
          sigma,
          mean(res$train[, 1]),
          res$test[, 1],
          sbcAftSurvival(config$t0, res$test[1L, 1L], sigma),
          # the imputed log survival time at the censored row, the draw the
          # generator's own logT0[i] is ranked among
          if (is.na(censoredRow)) NA_real_ else f$getLatents()[censoredRow]
        )
      }
    )
  } else {
    stop("no family spec for \"", config$family, "\"")
  }
  spec
}

# A latent BCF arm's configuration: the gaussian arm's n at nTest = 3, and
# sd.control at the FAMILY default 1 rather than gaussian's 2 (bcf.md's
# calibration section: under a latent family s is the link's own fixed error sd
# and sigma is pinned, so 2 would assert a median prognostic signal twice the
# noise). `arm` is what names it, since its family token is the link and the
# plain probit arm already owns that.
sbcBCFLatentConfig <- function(link) {
  config <- sbcAddBCF(
    sbcConfig(family = link, n = 200L, nTest = 3L),
    sdControl = 1
  )
  config$arm <- paste0("bcf-", link)
  # the ladder's PRESCRIBED |a| strata, run beside its prior-drawn datasets:
  # the settle time scales in |a| (sigma being pinned), and the half-Cauchy at
  # scale 1 reaches its own tail too rarely to measure the limit from draws
  # alone
  config$ladderStrata <- c(0.5, 2, 5, 10)
  config
}

# An arm's own name. Only where an arm does not name itself by its family does
# a config carry one: the two latent BCF arms share the plain probit arm's
# family token, so keying anything per arm on the family would collide.
sbcArmName <- function(config) {
  if (is.null(config$arm)) config$family else config$arm
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
    aft = sbcConfigAft(),
    "bcf-probit" = sbcBCFLatentConfig("probit"),
    "bcf-logistic" = sbcBCFLatentConfig("logistic"),
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

# The aft arm's own wiring check, on the channel only its censored-latent
# functional reads. $getLatents reports aft's imputed log survival time, so at
# an EVENT row it must be exactly the observed log time the response carries -
# data, not a draw - and at a CENSORED row it must sit strictly above that
# row's bound, which is the lower truncation the functional ranks theta0's
# logT0 against. A latent that was neither would make the rank compare two
# different quantities. Reports the drawn censored count too, the fixture's
# censoring rate at one replication.
sbcCheckAftLatents <- function(config, seed = 99L) {
  set.seed(seed)
  spec <- sbcFamilySpec(config, 1L, seed)
  drawn <- spec$draw()
  fit <- spec$fit(drawn$y)
  invisible(fit$run(0L, 1L))
  latents <- fit$getLatents()
  y <- drawn$y[, 1L]
  event <- drawn$y[, 2L] == 1
  eventDiff <- max(abs(latents[event] - y[event]))
  censoredGap <- if (any(!event)) min(latents[!event] - y[!event]) else NA_real_
  list(
    maxEventDiff = eventDiff,
    minCensoredGap = censoredGap,
    nCensored = sum(!event),
    pass = eventDiff == 0 && isTRUE(censoredGap > 0)
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
# The two latent BCF arms carry NO pre-registered burn. The gaussian BCF arm's
# 72000 is the (a, mu) amplitude ridge co-relaxing with tree-structure mixing,
# READ THROUGH sigma; pinning sigma removes the readout, not the ridge, and the
# misfit it absorbed lands in the index that these arms rank instead. Their
# ladder run fills these two in. The aft arm carries none either: a sweep is a
# gaussian one plus a truncated-normal draw per censored row, which argues for
# a gaussian-like cost and a t-like burn but measures neither, and the
# imputation is a second block the transient has to clear.
sbcBurnSweeps <- c(
  ordinal = 36000,
  nbinom = 24000,
  t = 12000,
  multinomial = 6000,
  aft = NA_real_,
  "bcf-probit" = NA_real_,
  "bcf-logistic" = NA_real_
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
  burnSweeps = sbcBurnSweeps[[sbcArmName(config)]],
  seed = 20260709L,
  report = 25L
) {
  if (!is.finite(burnSweeps)) {
    stop(
      "no measured burn for arm \"",
      sbcArmName(config),
      "\": run its ladder (sbc.R burn-",
      sbcArmName(config),
      ") and record the sweeps in sbcBurnSweeps, or pass the burn in sweeps ",
      "as the 5th positional argument"
    )
  }
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
      # a functional this replication does not define -- the aft arm's
      # censored latent when nothing censored -- contributes no rank and
      # stays NA; the report ranks the rest and says on how many replications
      if (is.na(theta[[j]])) {
        row[j] <- NA_integer_
        next
      }
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
        sbcArmName(config),
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
#
# `strata` appends datasets drawn at a PRESCRIBED value of the arm's own
# settle-driving parameter, after the prior-drawn ones - the latent BCF arms'
# |a|, whose worst stratum the prior reaches too rarely to measure from draws
# alone. Only a spec that offers drawAt can take them.
sbcBurnLadder <- function(
  config,
  nSweep = 20000L,
  nDataset = 3L,
  nBlock = 10L,
  seed = 20260804L,
  strata = config$ladderStrata
) {
  set.seed(seed)
  spec <- sbcFamilySpec(config, 1L, seed)
  if (length(strata) > 0L && is.null(spec$drawAt)) {
    stop("this family's spec cannot draw at a prescribed stratum")
  }
  blockSize <- nSweep %/% nBlock
  inputs <- c(rep(list(NULL), nDataset), as.list(strata))
  results <- vector("list", length(inputs))
  totalElapsed <- 0
  for (d in seq_along(inputs)) {
    drawn <- if (is.null(inputs[[d]])) {
      spec$draw()
    } else {
      spec$drawAt(inputs[[d]])
    }
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
    results[[d]] <- list(
      z = z,
      firstUnder = firstUnder,
      blockSize = blockSize,
      stratum = inputs[[d]]
    )
  }
  list(
    family = sbcArmName(config),
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
      "\n dataset %d%s: block-mean z vs the final half (block = %d sweeps)\n",
      d,
      if (is.null(res$stratum)) {
        ""
      } else {
        sprintf(" (prescribed |a| = %g)", res$stratum)
      },
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
  # a replication that does not define the functional (the aft arm's censored
  # latent, with nothing censored) carries NA and is not a rank: it leaves the
  # tabulation, the band's own R and the reported mean alike
  ranks <- ranks[!is.na(ranks)]
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
sbcReport <- function(
  fit,
  nBins = 20L,
  alpha = 0.05,
  expectedFlags = character(0)
) {
  cat(sprintf(
    "\nSBC report: family=%s n=%d p=%d nTrees=%d | R=%d L=%d thin=%d burn=%d\n",
    sbcArmName(fit$config),
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
    ranked <- fit$ranks[!is.na(fit$ranks[, funcs[i]]), funcs[i]]
    # a functional no replication defined has nothing to test; NONE is not a
    # FLAG, so the exit gate does not read it
    if (length(ranked) == 0L) {
      verdicts[i] <- "NONE"
      cat(sprintf(
        "%-10s %8s %8s %9s %8s %6s\n",
        funcs[i],
        "",
        "",
        "",
        "",
        "NONE"
      ))
      cat("  no replication defined this functional; nothing ranked\n")
      next
    }
    u <- rankUniformity(
      ranked,
      fit$L,
      nBins = nBins,
      alpha = alpha
    )
    verdicts[i] <- if (u$pass) {
      "PASS"
    } else if (funcs[i] %in% expectedFlags) {
      "FLAG (expected)" # pre-adjudicated; excluded from the SBC_FAIL_ON_FLAG gate
    } else {
      "FLAG"
    }
    cat(sprintf(
      "%-10s %8.3f %8.3f %9.4f %8.4f %6s\n",
      funcs[i],
      u$chisqP,
      u$ksP,
      u$ecdfDiff,
      u$ecdfBand,
      verdicts[i]
    ))
    if (length(ranked) < fit$R) {
      # the functional's own R, reported separately because it is not the
      # run's: the aft arm's censored latent is undefined whenever the
      # replication drew no censored row
      cat(sprintf(
        "  ranked on %d of %d replications; %d did not define it\n",
        length(ranked),
        fit$R,
        fit$R - length(ranked)
      ))
    }
    if (identical(verdicts[i], "FLAG (expected)")) {
      cat(
        "  pre-adjudicated flag, excluded from the exit check; the ",
        "adjudication is cited where SBC_EXPECTED_FLAGS is set\n",
        sep = ""
      )
    }
  }
  cat("\nRank histograms:\n")
  for (i in seq_along(funcs)) {
    ranked <- fit$ranks[!is.na(fit$ranks[, funcs[i]]), funcs[i]]
    if (length(ranked) == 0L) {
      next
    }
    cat(sprintf("\n[%s]\n", funcs[i]))
    cat(sbcAsciiHistogram(ranked, fit$L, nBins = nBins), "\n")
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
  isLatentBCF <- which %in% c("bcf-probit", "bcf-logistic")
  isBCF <- which %in% c("bcf", "bcf-weak", "bcf-probit", "bcf-logistic")
  isLinear <- which %in%
    c("linear", "linear-na-leaf", "linear-na-split", "linear-weighted")
  isGP <- which %in% c("gp", "gp-na-leaf", "gp-weighted", "gp-mixed")
  isFamilyTier <- which %in%
    c("ordinal", "nbinom", "t", "multinom", "multinomial", "aft")

  config <- if (isFamilyTier || isLatentBCF) {
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

  # The BCF arms' two by-hand controls, both opt-in and both off by default.
  # A poison must redden the arm it names, so an unknown name or an arm that
  # cannot carry one refuses the run rather than reporting a clean result.
  poison <- sbcBCFPoison(strsplit(Sys.getenv("SBC_POISON", ""), ",")[[1]])
  if (length(poison) > 0L && !isLatentBCF) {
    stop("SBC_POISON applies to the latent BCF arms (bcf-probit, bcf-logistic)")
  }
  if (isBCF) {
    config$poison <- poison
    config$fixedGlue <- nzchar(Sys.getenv("SBC_FIXED_GLUE", ""))
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
  if (isLatentBCF) {
    lb <- sbcCheckBCFLatent(config)
    cat(sprintf(
      "  transform: scale %.12f, shift %.2e, R2 %.12f; max |sigma - 1| %g\n",
      lb$fitScale,
      lb$fitShift,
      lb$mapR2,
      lb$maxSigma
    ))
    cat(sprintf(
      "  combined fits vs a mu + b_z tau: %.2e -> %s\n",
      lb$maxDiff,
      if (lb$pass) "PASS" else "FAIL"
    ))
    selfCheckPass["latent"] <- isTRUE(lb$pass)
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
    if (config$family == "aft") {
      al <- sbcCheckAftLatents(config)
      cat(sprintf(
        "  latents: max |event - y| %.2e; min censored gap %.4f over %d rows -> %s\n",
        al$maxEventDiff,
        al$minCensoredGap,
        al$nCensored,
        if (al$pass) "PASS" else "FAIL"
      ))
      selfCheckPass["latents"] <- isTRUE(al$pass)
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
  fit <- if (isFamilyTier || isLatentBCF) {
    if (is.null(burnSweeps)) {
      runSbcFamily(config, R = R, L = L, thin = thin)
    } else {
      runSbcFamily(config, R = R, L = L, thin = thin, burnSweeps = burnSweeps)
    }
  } else if (isDart) {
    runSbcDart(config, R = R, L = L, thin = thin)
  } else if (isBCF) {
    bcfArgs <- list(
      config,
      R = R,
      L = L,
      thin = thin,
      fixedGlue = isTRUE(config$fixedGlue)
    )
    if (!is.null(burnSweeps)) {
      bcfArgs$burn <- as.integer(ceiling(burnSweeps / thin))
    }
    do.call(runSbcBCF, bcfArgs)
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
  # matrix arms are admitted at the Bonferroni'd level; every other config keeps
  # the per-functional 5% band its recorded result was read at
  expectedFlags <- strsplit(Sys.getenv("SBC_EXPECTED_FLAGS", ""), ",")[[1]]
  expectedFlags <- trimws(expectedFlags[nzchar(trimws(expectedFlags))])
  verdicts <- sbcReport(
    fit,
    alpha = if (which %in% sbcMatrixConfigs) sbcMatrixAlpha else 0.05,
    expectedFlags = expectedFlags
  )
  if (nzchar(Sys.getenv("SBC_FAIL_ON_FLAG", "")) && any(verdicts == "FLAG")) {
    quit(status = 1L, save = "no")
  }
}
