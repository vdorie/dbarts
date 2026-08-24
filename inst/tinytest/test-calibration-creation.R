# Naming the leaf calibration at creation (docs/design/prior-defaults.md,
# "prior.scale"): a composed model states the prior sd of its forest total in
# RESPONSE units instead of inheriting it from the range of whatever vector the
# sampler happened to be constructed on. The oracles are two composed probit
# arms whose construction ranges are 16x apart - which today moves the posterior
# sd of f by about 5x - plus the refusals and the channels that must hold it.

set.seed(11)
n <- 100L
p <- 3L
x <- matrix(runif(n * p), n, p)
colnames(x) <- paste0("x", seq_len(p))
eta <- 1.5 * (x[, 1L] - x[, 2L])
yBinary <- rbinom(n, 1L, pnorm(eta))

# An Albert-Chib probit composed in R over a gaussian sampler: R draws the
# truncated-normal latents and the engine draws f | latents. The composed
# sampler's leaf prior is inherited from its CONSTRUCTION vector unless named.
# Returns the kept sweeps' fits, chain-major; keepFrom = sweeps keeps only the
# last, which is all the arm-vs-arm comparison needs.
composeProbit <- function(
  design,
  response,
  constructionRange,
  constructionShift,
  priorScale,
  sweeps,
  numChains = 2L,
  numTrees = 25L,
  keepFrom = sweeps,
  useRecipe = FALSE
) {
  rows <- length(response)
  yInit <- constructionShift +
    constructionRange * (seq(0, 1, length.out = rows) - 0.5)
  sampler <- dbarts(
    design,
    yInit,
    control = dbartsControl(
      n.chains = numChains,
      n.threads = 1L,
      n.trees = numTrees,
      n.samples = 1L,
      updateState = FALSE,
      seed = 29L,
      keepTrees = FALSE
    ),
    node.prior = if (is.na(priorScale)) {
      dbartsPriors$normal(k = 2)
    } else {
      dbartsPriors$normal(k = 2, scale = priorScale)
    },
    resid.prior = dbartsPriors$fixed(1),
    sigma = 1
  )
  if (useRecipe) {
    # the documented location lever: prior.mean is the response transform's
    # shift, and an offset of -prior.mean re-centers the modelled quantity.
    # ($getCalibration reports the same number once the mid-chain half lands.)
    priorMean <- constructionRange * 0.5 + min(yInit)
    sampler$setOffset(rep_len(-priorMean, rows), updateScale = FALSE)
  }
  latents <- ifelse(response == 1L, 0.5, -0.5)
  kept <- vector("list", sweeps - keepFrom + 1L)
  set.seed(5)
  for (sweep in seq_len(sweeps)) {
    sampler$setResponse(latents)
    fits <- drop(sampler$run(0L, 1L)$train)
    if (sweep >= keepFrom) {
      kept[[sweep - keepFrom + 1L]] <- fits
    }
    # the latents are driven off chain 1; every chain is compared, so a
    # conversion reaching only the first would still be caught
    mu <- if (is.null(dim(fits))) fits else fits[, 1L]
    u <- runif(rows)
    lower <- ifelse(response == 1L, pnorm(0, mu, 1), 0)
    upper <- ifelse(response == 1L, 1, pnorm(0, mu, 1))
    q <- pmin(pmax(lower + u * (upper - lower), 1e-12), 1 - 1e-12)
    latents <- qnorm(q, mu, 1)
  }
  if (length(kept) == 1L) kept[[1L]] else do.call(rbind, kept)
}

# --- arms built CENTERED (so the transform's shift matches by
# construction) at ranges 16x apart, naming the same prior.scale. The tolerance
# is pinned, not deferred: measured 1.0e-14 at the 120 sweeps shipped here and
# 4.3e-14 at 2000, growing like sqrt(sweeps), so 1e-12 holds with margin. ---
sweeps <- 120L
armA <- composeProbit(x, yBinary, 1.5, 0, 1.5, sweeps)
armB <- composeProbit(x, yBinary, 24, 0, 1.5, sweeps)
expect_equal(dim(armA), c(n, 2L))
expect_true(max(abs(armA - armB)) <= 1e-12)

# non-vacuity: unnamed, the same two arms are a different model entirely
bareA <- composeProbit(x, yBinary, 1.5, 0, NA_real_, sweeps)
bareB <- composeProbit(x, yBinary, 24, 0, NA_real_, sweeps)
expect_true(max(abs(bareA - bareB)) > 1)
expect_true(sd(bareB) / sd(bareA) > 3)

# Scope, stated here as it is in the help: both arms run at the same number
# of trees, so they see only fitScale-dependent errors and are BLIND to a
# calibration error that is a fixed function of the named value (one forgetting
# sqrt(m) passes this while reporting the wrong scale). That error is caught
# elsewhere, not left uncovered: test-calibration-prior-draws.R pins the
# ABSOLUTE prior sd at a known tree count, which a sqrt(m)-forgetting
# conversion fails outright, and test-calibration-midchain.R's set-then-get
# fidelity and static-m oracles pin the same error on the write path.

# --- the same arms NOT centered, each applying the documented offset
# recipe. This is the only test exercising the named scale and the location
# lever as one contract. ---
armA <- composeProbit(x, yBinary, 1.5, 3, 1.5, sweeps, useRecipe = TRUE)
armB <- composeProbit(x, yBinary, 24, -7, 1.5, sweeps, useRecipe = TRUE)
expect_true(max(abs(armA - armB)) <= 1e-12)

# without the recipe the shift is a function of the construction range and the
# arms part company by a wide margin
armA <- composeProbit(x, yBinary, 1.5, 3, 1.5, sweeps)
armB <- composeProbit(x, yBinary, 24, -7, 1.5, sweeps)
expect_true(max(abs(armA - armB)) > 1)

# --- a named composition targets the same posterior as the engine's own
# probit. The engine's probit anchor is node.scale 3.0 on a unit-scale latent,
# so prior.scale = 3.0 is the composition's statement of the same prior. ---
set.seed(21)
nO2 <- 200L
xO2 <- matrix(runif(nO2 * p), nO2, p)
colnames(xO2) <- paste0("x", seq_len(p))
etaO2 <- 2 * (xO2[, 1L] - xO2[, 2L])
yO2 <- rbinom(nO2, 1L, pnorm(etaO2))
o2Sweeps <- 400L
composeO2 <- function(constructionRange, priorScale) {
  composeProbit(
    xO2,
    yO2,
    constructionRange,
    0,
    priorScale,
    o2Sweeps,
    numChains = 1L,
    numTrees = 50L,
    keepFrom = o2Sweeps %/% 2L + 1L
  )
}
nativeSampler <- dbarts(
  xO2,
  yO2,
  control = dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 50L,
    n.samples = 1L,
    updateState = FALSE,
    seed = 29L,
    keepTrees = FALSE
  ),
  node.prior = dbartsPriors$normal(k = 2)
)
native <- t(nativeSampler$run(o2Sweeps %/% 2L, o2Sweeps %/% 2L)$train)
namedComposition <- composeO2(3, 3)

postSd <- function(draws) mean(apply(draws, 2L, sd))
rmse <- function(draws) sqrt(mean((colMeans(draws) - etaO2)^2))
expect_true(abs(postSd(namedComposition) / postSd(native) - 1) < 0.15)
expect_true(abs(rmse(namedComposition) / rmse(native) - 1) < 0.20)

# non-vacuity: the inherited calibration of a range-24 construction vector is
# not the probit prior, and both statistics say so loudly
inherited <- composeO2(24, NA_real_)
expect_true(postSd(inherited) / postSd(native) > 1.5)
expect_true(rmse(inherited) / rmse(native) > 1.5)

# --- the calibration refusal matrix, creation half. ---
set.seed(31)
nRef <- 120L
xRef <- matrix(runif(nRef * 4L), nRef, 4L)
colnames(xRef) <- paste0("x", seq_len(4L))
zRef <- rbinom(nRef, 1L, 0.5)
yRef <- 12 * (xRef[, 1L] - xRef[, 2L]) + 2 * zRef + rnorm(nRef)

refControl <- function(...) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 20L,
    n.samples = 10L,
    updateState = FALSE,
    seed = 17L,
    ...
  )
}

# (i) the two-forest model, refused R-side by name
expect_error(
  dbarts(
    xRef,
    yRef,
    forests = list(forest(), forest(basis = ~ factor(zRef))),
    control = refControl(),
    node.prior = normal(scale = 1.5)
  ),
  "a named 'prior.scale'"
)
# (ii) and again at the bridge, for a hand-built model the R layer never emits
resolved <- dbartsSpec(
  dbartsData(xRef, yRef),
  refControl(),
  forests = list(forest(), forest(basis = ~ factor(zRef)))
)
resolved$model@prior.scale <- 1.5
expect_error(
  new("dbartsSampler", resolved$control, resolved$model, resolved$data),
  "a named 'prior.scale'"
)
# (iii) the multinomial creation path, whose leaf scales come from the softmax
# calibration map: the low-level handle still refuses at the bridge (no
# R-level resolution runs on it), while bart2 now builds directly through
# dbarts()'s own resolveSamplerSpec, whose multinomial-specific check catches
# a named prior.scale earlier, with the same "a named 'prior.scale'" text
# (i)/(ii) above use, before any sampler is created at all
host <- dbarts(
  xRef,
  yRef,
  control = refControl(),
  node.prior = normal(scale = 1.5)
)
labels <- sample(0:2, nRef, replace = TRUE)
expect_error(
  dbarts:::bartcoreMultinomialSampler(host, labels, K = 3L),
  "softmax calibration map"
)
expect_error(
  bart2(
    xRef,
    factor(labels),
    family = "multinomial",
    prior.scale = 1.5,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.trees = 10L,
    verbose = FALSE
  ),
  "a named 'prior.scale'"
)

# a non-finite or non-positive value is an error, not a refusal
expect_error(dbartsPriors$normal(scale = -1), "'scale' must be positive")
expect_error(dbartsPriors$normal(sd = Inf), "'sd' must be positive")
expect_error(dbartsPriors$normal(sd = c(1, 2)), "single number")
expect_error(dbartsPriors$normal(sd = 1, scale = 1), "at most one")
expect_error(
  bart2(
    xRef,
    yRef,
    prior.scale = -1,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    verbose = FALSE
  ),
  "'prior.scale' must be positive"
)

# the sd spelling under a sampled k, with both remedies named. The binary
# default IS the sampled path, so this is the ordinary probit fit.
sdRefusal <- tryCatch(
  dbarts(
    xRef,
    zRef,
    control = refControl(),
    node.prior = normal(sd = 1.5)
  ),
  error = function(e) conditionMessage(e)
)
expect_true(is.character(sdRefusal))
expect_true(grepl("'scale'", sdRefusal, fixed = TRUE))
expect_true(grepl("fix 'k'", sdRefusal, fixed = TRUE))
# named 'scale' is honored under exactly the same hyperprior
expect_equal(
  dbarts(
    xRef,
    zRef,
    control = refControl(),
    node.prior = normal(scale = 1.5)
  )$model@prior.scale,
  1.5
)
# and the sd sugar resolves against a fixed k
expect_equal(
  dbarts(
    xRef,
    yRef,
    control = refControl(),
    node.prior = normal(k = 4, sd = 0.5)
  )$model@prior.scale,
  2.0
)

# a DART sampler is NOT refused: the named calibration and the Dirichlet split
# machinery are independent, and the fit runs clean
dartSampler <- dbarts(
  xRef,
  yRef,
  control = refControl(),
  tree.prior = dart(),
  node.prior = normal(scale = 1.5)
)
expect_equal(dartSampler$model@prior.scale, 1.5)
expect_true(all(is.finite(dartSampler$run(20L, 10L)$train)))

# monotone / linear / gp accept it too - the setter half is total over all four
# leaf models, and the creation half must be as well
expect_equal(
  dbarts(
    xRef,
    yRef,
    control = refControl(),
    monotone = c(x1 = 1),
    node.prior = normal(scale = 1.5)
  )$model@prior.scale,
  1.5
)
expect_equal(
  dbarts(
    xRef,
    yRef,
    control = refControl(),
    node.prior = linear("x1", scale = 1.5)
  )$model@prior.scale,
  1.5
)
expect_equal(
  dbarts(
    xRef,
    yRef,
    control = refControl(),
    node.prior = gp("x1", scale = 1.5)
  )$model@prior.scale,
  1.5
)

# --- the channels the named calibration must survive. ---

# the model slot records the INTENT and nothing else writes it
plain <- dbarts(xRef, yRef, control = refControl())
expect_true(is.na(plain$model@prior.scale))
named <- dbarts(
  xRef,
  yRef,
  control = refControl(),
  node.prior = normal(k = 2, scale = 1.5)
)
expect_equal(named$model@prior.scale, 1.5)
named$setResponse(yRef + 1, updateScale = TRUE)
expect_equal(named$model@prior.scale, 1.5)

# setModel's re-derivation: $setModel(sampler$model) is a documented no-op, and
# without the conversion it would revert the leaf scale to the family-keyed
# node.scale (an 8x move on this response). Bitwise, because a no-op is.
roundTripControl <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 10L,
  updateState = FALSE,
  seed = 17L,
  keepTrees = FALSE
)
withoutSetModel <- dbarts(
  xRef,
  yRef,
  control = roundTripControl,
  node.prior = normal(k = 2, scale = 1.5)
)
withSetModel <- dbarts(
  xRef,
  yRef,
  control = roundTripControl,
  node.prior = normal(k = 2, scale = 1.5)
)
withSetModel$setModel(withSetModel$model)
expect_identical(
  withoutSetModel$run(20L, 10L)$train,
  withSetModel$run(20L, 10L)$train
)
# non-vacuity: a model whose named intent has been stripped IS the revert, and
# the same round trip then moves the draws
reverted <- dbarts(
  xRef,
  yRef,
  control = roundTripControl,
  node.prior = normal(k = 2, scale = 1.5)
)
strippedModel <- reverted$model
strippedModel@prior.scale <- NA_real_
reverted$setModel(strippedModel)
expect_false(identical(
  withoutSetModel$run(20L, 10L)$train,
  reverted$run(20L, 10L)$train
))

# --- xbart rows: a named calibration is held across cells, in the branch
# that creates a sampler and in the branch that re-models one. xbart forces
# n.chains/n.threads/updateState itself (control = is no longer a formal),
# so nothing here needs to name them. ---
xbartArgs <- list(
  formula = xRef,
  data = yRef,
  n.samples = 20L,
  n.reps = 2L,
  n.burn = c(20L, 20L, 20L),
  n.trees = 25L,
  n.threads = 1L,
  seed = 5L,
  verbose = FALSE
)
# xbart reports one loss per replication per cell; average the replications
cellLoss <- function(result) {
  if (is.null(dim(result))) mean(result) else colMeans(result)
}
# the k grid: cell 1 creates the sampler, cells 2+ take the setModel branch
sweptLoss <- cellLoss(do.call(
  xbart,
  c(xbartArgs, list(k = c(1, 4), node.prior = quote(normal(scale = 1.5))))
))
# the same two cells run one at a time, so each takes the creation branch
singleLoss <- vapply(
  c(1, 4),
  function(kValue) {
    cellLoss(do.call(
      xbart,
      c(xbartArgs, list(k = kValue, node.prior = quote(normal(scale = 1.5))))
    ))
  },
  numeric(1L)
)
expect_equal(length(sweptLoss), 2L)
expect_true(all(is.finite(sweptLoss)))
expect_true(max(abs(sweptLoss / singleLoss - 1)) < 0.25)
# non-vacuity: the un-named grid is a materially different loss surface at the
# same cells, so the assertion above is not satisfied by any pair of numbers
inheritedLoss <- cellLoss(do.call(xbart, c(xbartArgs, list(k = c(1, 4)))))
expect_true(max(abs(inheritedLoss / sweptLoss - 1)) > 0.25)

# the k hyperprior arm: held across cells, and the sd spelling meets the
# sampled-k refusal here exactly as it does everywhere else
hyperLoss <- cellLoss(do.call(
  xbart,
  c(
    xbartArgs,
    list(
      k = dbartsPriors$chi(1.5, 2),
      node.prior = quote(normal(scale = 1.5))
    )
  )
))
expect_true(all(is.finite(hyperLoss)))
expect_error(
  do.call(
    xbart,
    c(
      xbartArgs,
      list(
        k = dbartsPriors$chi(1.5, 2),
        node.prior = quote(normal(sd = 1.5))
      )
    )
  ),
  "drawn every"
)
