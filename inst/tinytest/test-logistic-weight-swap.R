# The logistic weight channel. A logistic sampler's case weights are the
# observation counts its Polya-Gamma latents are built from, so replacing them
# is a model change with a defined meaning: setWeights stores the counts and
# rebuilds the working response from the omega in force, the next sweep drawing
# omega first against the new counts, and an outer sampler can vary exposure
# between sweeps. Probit, ordinal, aft and nbinom decline by identification and
# keep their refusals.

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
    updateState = FALSE,
    seed = 902L
  )
}
logisticSamplerWeightSwap <- function(weights, n.chains = 1L, offset = NULL) {
  dbarts(
    x,
    y,
    weights = weights,
    offset = offset,
    family = "logistic",
    control = logisticControl(n.chains),
    n.trees = 20L
  )
}

# --- the build-vs-swap oracle -------------------------------------------
# Creation consumes no weight-dependent randomness and leaves every chain
# generator at the same position, and the swap draws nothing, so a swapped
# sampler keeps its latents and its sweeps are BITWISE those of one built with
# the new counts, per chain.
for (chains in 1:2) {
  swapped <- logisticSamplerWeightSwap(w1, chains)
  built <- logisticSamplerWeightSwap(w2, chains)
  held <- swapped$getLatents()
  swapped$setWeights(w2)
  expect_identical(swapped$getLatents(), held)
  expect_identical(swapped$run(0L, 2L)$train, built$run(0L, 2L)$train)
}

# not vacuous: a sampler left at the old counts is a different sampler
untouched <- logisticSamplerWeightSwap(w1)
expect_false(identical(
  untouched$run(0L, 2L)$train,
  logisticSamplerWeightSwap(w2)$run(0L, 2L)$train
))

# the offset is the same conduit: one created with o1 and handed o2 is one
# created with o2
o1 <- 0.3 * x[, 2L]
o2 <- -0.5 + x[, 2L]
handed <- logisticSamplerWeightSwap(w2, offset = o1)
created <- logisticSamplerWeightSwap(w2, offset = o2)
handed$setOffset(o2)
expect_identical(handed$run(0L, 3L)$train, created$run(0L, 3L)$train)

# and off a grown forest, the swap leaves the latents where the sweep put them
grown <- logisticSamplerWeightSwap(w1)
invisible(grown$run(5L, 1L))
grownBefore <- grown$getLatents()
grown$setWeights(w2)
expect_identical(grown$getLatents(), grownBefore)

# the shape: at psi = 0 a count-w row's omega is PG(w, 0) with mean w/4, so the
# mean latent the first sweep draws tracks the counts the swap installed.
# Maskless on purpose - see the mask cell below
unitCounts <- logisticSamplerWeightSwap(w1)
eightCounts <- logisticSamplerWeightSwap(w1)
unitCounts$setWeights(rep(1, n))
eightCounts$setWeights(rep(8, n))
invisible(unitCounts$run(0L, 1L))
invisible(eightCounts$run(0L, 1L))
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
  updateState = FALSE,
  seed = 77L
)
altSampler <- function(weights) {
  dbarts(
    xAlt,
    yAlt,
    weights = weights,
    family = "logistic",
    control = altControl,
    n.trees = 40L
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
# A swap moves no latent, active or not; the next sweep draws for ACTIVE rows
# only, since consuming variates for an inactive row would desynchronize the
# stream against a sampler built on the retained rows.
masked <- logisticSamplerWeightSwap(w1)
mask <- rep(1, n)
mask[1:5] <- 0
masked$setActiveRows(mask)
invisible(masked$run(3L, 1L))
held <- masked$getLatents()
masked$setWeights(rep(4, n))
expect_identical(masked$getLatents(), held)
invisible(masked$run(0L, 1L))
expect_identical(masked$getLatents()[1:5], held[1:5])
expect_true(all(masked$getLatents()[-(1:5)] != held[-(1:5)]))

# --- what a count may not be, on both conduits -----------------------------
countRefusal <- "must be non-negative whole numbers"
pinned <- logisticSamplerWeightSwap(w2)
expect_error(pinned$setWeights(w2 + 0.5), countRefusal)
expect_error(pinned$setData(dbartsData(x, y, weights = w2 + 0.5)), countRefusal)
# a negative count is caught R-side first, on the rule every family shares
expect_error(pinned$setWeights(replace(w2, 1L, -1)), "must all be non-negative")
# every refusal leaves the mirrored vector where it was
expect_identical(pinned$data@weights, w2)
# a count of 0 is taken on both conduits, the row out of the likelihood; the
# data object warns of it, as on every family
zeroWarnings <- captureWarnings(
  zeroData <- dbartsData(x, y, weights = replace(w2, 1L, 0))
)
expect_equal(length(zeroWarnings), 1L)
expect_true(grepl(
  "of 0 are left out of the likelihood",
  conditionMessage(zeroWarnings[[1L]])
))
pinned$setWeights(replace(w2, 1L, 0))
expect_identical(pinned$data@weights, replace(w2, 1L, 0))
pinned$setData(zeroData)
expect_identical(pinned$data@weights, replace(w2, 1L, 0))
pinned$setWeights(w2)

# replacement data given WITHOUT weights is unweighted, as at creation and as
# it is for gaussian - single-trial rows, unit counts. setData starts the
# latents as creation does, so the build-vs-swap oracle holds on this conduit
# too: a sampler created with counts and handed weightless data is bitwise one
# built unweighted.
reset <- logisticSamplerWeightSwap(w2)
unweighted <- logisticSamplerWeightSwap(w1)
reset$setData(dbartsData(x, y))
expect_true(is.null(reset$data@weights))
expect_identical(reset$getLatents(), unweighted$getLatents())
expect_identical(reset$run(0L, 2L)$train, unweighted$run(0L, 2L)$train)
# not vacuous: the counts it was created with are gone, not carried through
expect_false(identical(
  reset$getLatents(),
  logisticSamplerWeightSwap(w2)$getLatents()
))
# and stated counts on the same conduit start where creation starts them
pinned$setData(dbartsData(x, y, weights = w2))
expect_identical(
  pinned$getLatents(),
  logisticSamplerWeightSwap(w2)$getLatents()
)

# --- the families that decline by identification ---------------------------
# probit and ordinal carry no precision, so a 0/1 vector is membership, kept
# on the data object, all ones as no weights; w2 is a count vector and stays
# refused on both conduits. nbinom reads 0/1 as membership too and takes w2 as
# counts on both.
masked <- list(
  probit = list(
    y = y,
    text = "probit models do not support case weights other than 0 and 1"
  ),
  ordinal = list(
    y = as.double(1L + (seq_len(n) %% 3L)),
    text = "ordinal models do not support case weights other than 0 and 1"
  )
)
for (name in names(masked)) {
  cell <- masked[[name]]
  fit <- dbarts(
    x,
    cell$y,
    family = name,
    control = logisticControl(),
    n.trees = 20L
  )
  expect_silent(fit$setWeights(w1))
  expect_true(is.null(fit$data@weights), info = name)
  expect_silent(fit$setData(dbartsData(x, cell$y, weights = w1)))
  expect_true(is.null(fit$data@weights), info = name)
  expect_error(fit$setWeights(w2), cell$text, info = name)
  expect_error(
    fit$setData(dbartsData(x, cell$y, weights = w2)),
    cell$text,
    info = name
  )
}
countFit <- dbarts(
  x,
  as.double(seq_len(n) %% 5L),
  family = "nbinom",
  control = logisticControl(),
  n.trees = 20L
)
expect_silent(countFit$setWeights(w1))
expect_true(is.null(countFit$data@weights))
expect_silent(countFit$setWeights(w2))
expect_identical(countFit$data@weights, w2)
expect_silent(countFit$setData(dbartsData(
  x,
  as.double(seq_len(n) %% 5L),
  weights = w2
)))
expect_identical(countFit$data@weights, w2)
expect_error(
  countFit$setWeights(w2 + 0.5),
  "nbinom weights are observation counts and must be non-negative whole numbers",
  fixed = TRUE
)
aftFit <- dbarts(
  x,
  cbind(exp(f + rnorm(n)), rep_len(c(1L, 0L), n)),
  family = "aft",
  control = logisticControl(),
  n.trees = 20L
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
  chains,
  swapped,
  built,
  held,
  untouched,
  countFit,
  o1,
  o2,
  handed,
  created,
  grown,
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
  countRefusal,
  pinned,
  reset,
  unweighted,
  zeroWarnings,
  zeroData,
  name,
  cell,
  fit,
  aftFit
)
