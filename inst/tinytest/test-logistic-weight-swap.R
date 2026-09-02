# The logistic weight channel. A logistic sampler's case weights are the
# observation counts its Polya-Gamma latents are built from, so replacing them
# is a model change with a defined meaning: the engine redraws omega against
# the new counts on the spot, and an outer sampler can vary exposure between
# sweeps. Probit, ordinal, aft and nbinom decline by identification and keep
# their refusals.

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

set.seed(9021L, sample.kind = "Rejection")
n <- 80L
x <- matrix(runif(n * 2L), n)
f <- 2.5 * x[, 1L] - 1.2
y <- as.double(rbinom(n, 1L, plogis(f)))
w1 <- rep(1, n)
w2 <- as.double(sample(1:3, n, replace = TRUE))

logisticControl <- function(n.chains = 1L) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE,
    seed = 902L
  )
}
logisticSamplerWeightSwap <- function(weights, n.chains = 1L) {
  dbarts(
    x,
    y,
    weights = weights,
    family = "logistic",
    control = logisticControl(n.chains)
  )
}

# --- the build-vs-swap oracle -------------------------------------------
# Creation consumes no weight-dependent randomness (the cold start is
# deterministic) and leaves every chain generator at the same position, and
# setResponse is the already-shipped conduit that runs the SAME refresh against
# the same location, so a swapped sampler and one built with the new counts
# agree BITWISE rather than merely closely.
swapped <- logisticSamplerWeightSwap(w1)
built <- logisticSamplerWeightSwap(w2)
swapped$setWeights(w2)
built$setResponse(y)
expect_identical(swapped$getLatents(), built$getLatents())
expect_identical(swapped$run(0L, 2L)$train, built$run(0L, 2L)$train)

# not vacuous: a sampler left at the old counts is a different sampler
untouched <- logisticSamplerWeightSwap(w1)
expect_false(identical(swapped$getLatents(), untouched$getLatents()))

# the fan-out is per chain, each drawing from its own generator, which is what
# a broken one gets wrong on the second chain alone
swappedPair <- logisticSamplerWeightSwap(w1, 2L)
builtPair <- logisticSamplerWeightSwap(w2, 2L)
swappedPair$setWeights(w2)
builtPair$setResponse(y)
expect_identical(dim(swappedPair$getLatents()), c(n, 2L))
expect_identical(swappedPair$getLatents(), builtPair$getLatents())

# and off a grown forest, where the location the latents are drawn against is
# no longer the all-zero one creation leaves
ranA <- logisticSamplerWeightSwap(w1)
ranB <- logisticSamplerWeightSwap(w1)
invisible(ranA$run(5L, 1L))
invisible(ranB$run(5L, 1L))
grownBefore <- ranA$getLatents()
ranA$setWeights(w2)
ranB$setWeights(w2)
expect_identical(ranA$getLatents(), ranB$getLatents())
expect_false(identical(ranA$getLatents(), grownBefore))

# the shape, deterministically and with no run: at psi = 0 a count-w row's
# omega is PG(w, 0) with mean w/4, so the mean latent tracks the counts the
# swap installs. Maskless on purpose - see the mask cell below
unitCounts <- logisticSamplerWeightSwap(w1)
eightCounts <- logisticSamplerWeightSwap(w1)
unitCounts$setWeights(rep(1, n))
eightCounts$setWeights(rep(8, n))
expect_true(abs(mean(unitCounts$getLatents()) - 0.25) < 0.05)
expect_true(abs(mean(eightCounts$getLatents()) - 2) < 0.4)

# --- the outer Gibbs the channel exists for ------------------------------
# A one-shot swap self-heals: the next sweep's own refresh redraws omega from
# the new pointer, so a deferred refresh costs one sweep. ALTERNATING exposure
# never heals - every harvested sweep would move its trees under the other
# vector's precisions - so this is the arm with teeth. Outcome-tied counts
# separate the weighted curve from the unweighted one.
set.seed(410L, sample.kind = "Rejection")
nAlt <- 200L
xAlt <- matrix(runif(nAlt * 2L), nAlt)
yAlt <- as.double(rbinom(nAlt, 1L, plogis(3 * xAlt[, 1L] - 1.5)))
wFlat <- rep(1, nAlt)
wTied <- ifelse(yAlt == 1, 4, 1)
altControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 40L,
  updateState = FALSE,
  seed = 77L
)
altSampler <- function(weights) {
  dbarts(
    xAlt,
    yAlt,
    weights = weights,
    family = "logistic",
    control = altControl
  )
}
alternating <- altSampler(wTied)
invisible(alternating$run(100L, 1L))
harvest <- matrix(0, nAlt, 200L)
kept <- 0L
for (i in seq_len(400L)) {
  alternating$setWeights(if (i %% 2L == 1L) wFlat else wTied)
  draw <- alternating$run(0L, 1L)
  if (i %% 2L == 0L) {
    kept <- kept + 1L
    harvest[, kept] <- as.numeric(draw$train)
  }
}
expect_equal(kept, 200L)
altFit <- plogis(rowMeans(harvest))
tiedReference <- altSampler(wTied)
flatReference <- altSampler(wFlat)
invisible(tiedReference$run(100L, 1L))
invisible(flatReference$run(100L, 1L))
tiedFit <- plogis(rowMeans(tiedReference$run(0L, 200L)$train))
flatFit <- plogis(rowMeans(flatReference$run(0L, 200L)$train))
# the harvested sweeps sit on the weighted curve
expect_true(cor(altFit, tiedFit) > 0.9)
expect_true(mean(abs(altFit - tiedFit)) < 0.08)
# and the two curves are far enough apart for that to mean something
expect_true(mean(abs(tiedFit - flatFit)) > 0.12)
expect_true(mean(abs(altFit - flatFit)) > 0.12)

# --- the mask -------------------------------------------------------------
# A swap draws for ACTIVE rows only: consuming variates for an inactive row
# would desynchronize the stream against a sampler built on the retained rows.
# The inactive rows are not left carrying the OLD counts' omega either - they
# return to the deterministic cold start against the new count, 0.25 w.
masked <- logisticSamplerWeightSwap(w1)
mask <- rep(1, n)
mask[1:5] <- 0
masked$setActiveRows(mask)
invisible(masked$run(3L, 1L))
masked$setWeights(rep(4, n))
expect_equal(masked$getLatents()[1:5], rep(1, 5L))
expect_true(all(masked$getLatents()[-(1:5)] > 0))
expect_false(any(masked$getLatents()[-(1:5)] == 1))

# --- grouped random effects ------------------------------------------------
# rbart_vi's decorator delegates the swap to its base family and draws NOTHING
# of its own: b and tau are the sweep's blocks, not the swap's, exactly as they
# are under setResponse.
groupControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 15L,
  updateState = FALSE,
  seed = 41L
)
attr(groupControl, "bartcore.groups") <- list(
  indices = rep_len(1:3, n),
  n.groups = 3L,
  prior = "cauchy",
  rel.scale = 0.5,
  n.steps = 1L
)
grouped <- dbarts(
  x,
  y,
  weights = w1,
  family = "logistic",
  control = groupControl
)
invisible(grouped$run(3L, 1L))
grouped$storeState()
groupedBefore <- grouped$state[[1L]]
grouped$setWeights(w2)
grouped$storeState()
groupedAfter <- grouped$state[[1L]]
expect_identical(groupedAfter$ranef, groupedBefore$ranef)
expect_identical(groupedAfter$tau, groupedBefore$tau)
expect_false(identical(groupedAfter$latents, groupedBefore$latents))
expect_true(all(is.finite(grouped$run(0L, 1L)$train)))

# --- what a count may not be, on both conduits -----------------------------
countRefusal <- "must be positive integers"
pinned <- logisticSamplerWeightSwap(w2)
zeroWarnings <- captureWarnings(
  zeroData <- dbartsData(x, y, weights = replace(w2, 1L, 0))
)
# the data object only WARNS about a zero weight; the logistic sampler refuses
expect_equal(length(zeroWarnings), 1L)
expect_true(grepl("of 0 will be ignored", conditionMessage(zeroWarnings[[1L]])))
expect_inherits(zeroWarnings[[1L]], "dbartsIgnoredArgWarning")
expect_error(pinned$setWeights(replace(w2, 1L, 0)), countRefusal)
expect_error(pinned$setData(zeroData), countRefusal)
expect_error(pinned$setWeights(w2 + 0.5), countRefusal)
expect_error(pinned$setData(dbartsData(x, y, weights = w2 + 0.5)), countRefusal)
# a negative count is caught R-side first, on the rule every family shares
expect_error(pinned$setWeights(replace(w2, 1L, -1)), "must all be non-negative")
# every refusal leaves the mirrored vector where it was
expect_identical(pinned$data@weights, w2)

# replacement data given WITHOUT weights is unweighted, as at creation and as
# it is for gaussian - single-trial rows, unit counts. The latents are DRAWN
# against them rather than left at the cold start a data swap installs, so the
# same build-vs-swap oracle holds on this conduit too: a sampler created with
# counts and handed weightless data is bitwise one built unweighted.
reset <- logisticSamplerWeightSwap(w2)
unweighted <- logisticSamplerWeightSwap(w1)
reset$setData(dbartsData(x, y))
unweighted$setResponse(y)
expect_true(is.null(reset$data@weights))
expect_identical(reset$getLatents(), unweighted$getLatents())
expect_identical(reset$run(0L, 2L)$train, unweighted$run(0L, 2L)$train)
# not vacuous: the counts it was created with are gone, not carried through
expect_false(identical(
  reset$getLatents(),
  logisticSamplerWeightSwap(w2)$getLatents()
))
# and stated counts on the same conduit are drawn against, not cold-started
pinned$setData(dbartsData(x, y, weights = w2))
expect_false(isTRUE(all.equal(pinned$getLatents(), 0.25 * w2)))
expect_true(all(pinned$getLatents() > 0))

# --- the families that decline by identification ---------------------------
declined <- list(
  probit = list(
    y = y,
    text = "probit models do not support case weights"
  ),
  ordinal = list(
    y = as.double(1L + (seq_len(n) %% 3L)),
    text = "ordinal models do not support case weights"
  ),
  nbinom = list(
    y = as.double(seq_len(n) %% 5L),
    text = "nbinom \\(count\\) models do not support case weights"
  )
)
for (name in names(declined)) {
  cell <- declined[[name]]
  fit <- dbarts(x, cell$y, family = name, control = logisticControl())
  expect_error(fit$setWeights(w1), cell$text, info = name)
  expect_error(
    fit$setData(dbartsData(x, cell$y, weights = w1)),
    cell$text,
    info = name
  )
}
aftFit <- dbarts(
  x,
  cbind(exp(f + rnorm(n)), rep_len(c(1L, 0L), n)),
  family = "aft",
  control = logisticControl()
)
expect_error(
  aftFit$setWeights(w1),
  "aft \\(survival\\) models do not support case weights"
)
expect_error(aftFit$setData(dbartsData(x, y)), "fix the censoring structure")

rm(
  n,
  x,
  f,
  y,
  w1,
  w2,
  logisticControl,
  logisticSamplerWeightSwap,
  swapped,
  built,
  untouched,
  swappedPair,
  builtPair,
  ranA,
  ranB,
  grownBefore,
  unitCounts,
  eightCounts,
  nAlt,
  xAlt,
  yAlt,
  wFlat,
  wTied,
  altControl,
  altSampler,
  alternating,
  harvest,
  kept,
  i,
  draw,
  altFit,
  tiedReference,
  flatReference,
  tiedFit,
  flatFit,
  masked,
  mask,
  groupControl,
  grouped,
  groupedBefore,
  groupedAfter,
  countRefusal,
  pinned,
  zeroWarnings,
  zeroData,
  declined,
  name,
  cell,
  fit,
  aftFit
)
