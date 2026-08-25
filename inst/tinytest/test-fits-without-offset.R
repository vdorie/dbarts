# $getFitsWithoutOffset(): the combined per-observation location on the
# response scale, without the installed offset.
#
# TWO INSTRUMENTS, because the storeSample refactor blinds one of them.
# (a) The IDENTITY cells below gate PLUMBING, not arithmetic: after the
# refactor both sides of getFitsWithoutOffset() == run()$train - offset
# evaluate the same Chain::fitsWithoutOffset, so a defect inside it moves both
# sides and the identity stays green. What they do gate is the bridge entry,
# the R-level method, the R/bartcore.R wrapper, the per-chain indexing, the return
# shape, the two refusals, and the OFFSET ADD - the one step the two sides do
# not share.
# (b) The ARITHMETIC gates are two, neither sharing the site: a tests/cpp cell
# in test_sampler.cpp pinning fitScale * combined + fitShift directly, and the
# two RECOMBINATION cells at the foot of this file, which rebuild the expected
# value through the independent read path ($getForestFits, a memcpy that never
# touches the accessor, plus $getForestAmplitudes and $getCalibration).

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

set.seed(3141L)
n <- 80L
p <- 3L
x <- matrix(runif(n * p), n, p)

samplerControl <- function(n.chains = 1L, ...) {
  dbartsControl(
    n.threads = 1L,
    n.trees = 20L,
    updateState = FALSE,
    seed = 61L,
    n.chains = n.chains,
    ...
  )
}

# Run one sweep, then read: storeSample is the last act of a sweep, so the
# accessor sees exactly the state the training channel recorded. The three
# assertions stay at TOP LEVEL - tinytest declines to register an expectation
# written as tinytest::expect_*, and an unqualified one inside a helper is a
# lintr finding - so the helper only gathers what they compare. Shape is
# asserted BEFORE the values: a NULL or a length mismatch would otherwise let
# the comparison pass vacuously.
fitIdentity <- function(sampler, n.obs = n) {
  r <- sampler$run(0L, 1L)
  nChains <- sampler$control@n.chains
  list(
    fits = sampler$getFitsWithoutOffset(),
    train = array(r$train, c(n.obs, nChains)),
    shape = c(n.obs, nChains)
  )
}

# --- gaussian: NON-UNIT response range and a NON-NULL offset, so the cell is
# not vacuous under either the internal-scale mutation or the added-offset one

yGauss <- 100 + 40 * x[, 1L] + 15 * x[, 2L] * x[, 3L] + rnorm(n, sd = 3)
offGauss <- rnorm(n, sd = 2)
expect_true(diff(range(yGauss)) > 10)
expect_true(any(offGauss != 0))
samplerGauss <- dbarts(
  x,
  yGauss,
  offset = offGauss,
  control = samplerControl(n.chains = 2L)
)
cell <- fitIdentity(samplerGauss)
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offGauss, cell$train)

# --- probit, logistic: the offset enters the latent scale

yBinary <- as.numeric(yGauss > median(yGauss))
offBinary <- rnorm(n, sd = 0.4)
cell <- fitIdentity(
  dbarts(x, yBinary, offset = offBinary, control = samplerControl(2L))
)
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offBinary, cell$train)

cell <- fitIdentity(dbarts(
  x,
  yBinary,
  offset = offBinary,
  family = "logistic",
  control = samplerControl()
))
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offBinary, cell$train)

# --- ordinal

yOrdinal <- cut(
  yGauss,
  quantile(yGauss, c(0, 1 / 3, 2 / 3, 1)),
  labels = FALSE,
  include.lowest = TRUE
)
cell <- fitIdentity(dbarts(
  x,
  yOrdinal,
  offset = offBinary,
  family = "ordinal",
  control = samplerControl()
))
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offBinary, cell$train)

# --- aft: a right-censored (time, status) response

timeAft <- exp(1 + x[, 1L] + rnorm(n, sd = 0.3))
statusAft <- rbinom(n, 1L, 0.8)
cell <- fitIdentity(dbarts(
  x,
  cbind(timeAft, statusAft),
  offset = offBinary,
  family = "aft",
  control = samplerControl()
))
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offBinary, cell$train)

# --- Student-t: a gaussian FAMILY whose residual distribution is t, which is
# also the sampler whose getLatents reports a PRECISION rather than a location

cell <- fitIdentity(dbarts(
  x,
  yGauss,
  offset = offGauss,
  resid.dist = student(df = 4),
  control = samplerControl()
))
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offGauss, cell$train)

# --- nbinom

yCount <- rnbinom(n, size = 4, mu = exp(0.5 + x[, 1L]))
cell <- fitIdentity(dbarts(
  x,
  yCount,
  offset = offBinary,
  family = "nbinom",
  control = samplerControl()
))
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offBinary, cell$train)

# --- Bayesian causal forest: the accessor's whole justification, since no
# other R route reaches the combined offset-free fit ($getForestFits gives one
# forest's internal-scale totals, $predict is refused, run()$train carries the
# offset)

z <- rbinom(n, 1L, 0.5)
yBcf <- 10 + 4 * x[, 1L] + z * (2 + 3 * x[, 2L]) + rnorm(n, sd = 0.5)
offBcf <- rnorm(n, sd = 0.5)
samplerBcf <- dbarts(
  x,
  yBcf,
  offset = offBcf,
  forests = list(forest(), forest(basis = ~ factor(z))),
  control = samplerControl()
)
cell <- fitIdentity(samplerBcf)
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offBcf, cell$train)

# --- grouped random intercepts (rbart_vi's internal control attribute)

groupedControl <- samplerControl()
attr(groupedControl, "bartcore.groups") <- list(
  indices = rep(seq_len(4L), length.out = n),
  n.groups = 4L,
  prior = "cauchy",
  rel.scale = sd(yGauss),
  n.steps = 1L
)
cell <- fitIdentity(
  dbarts(x, yGauss, offset = offGauss, control = groupedControl)
)
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offGauss, cell$train)

# --- heteroscedastic: the variance forest is a second reporting surface and
# must not leak into the location channel

cell <- fitIdentity(dbarts(
  x,
  yGauss,
  offset = offGauss,
  variance = varianceForest(n.trees = 10L),
  control = samplerControl()
))
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offGauss, cell$train)

# --- an active-row mask: an inactive row still receives a fitted value, so the
# accessor reports every row

samplerMasked <- dbarts(
  x,
  yGauss,
  offset = offGauss,
  control = samplerControl()
)
samplerMasked$setActiveRows(rep(c(1, 0), length.out = n))
cell <- fitIdentity(samplerMasked)
expect_true(!is.null(cell$fits))
expect_equal(dim(cell$fits), cell$shape)
expect_equal(cell$fits + offGauss, cell$train)

# --- the refusal -----------------------------------------------------------

# multinomial: exercised directly against the bridge-level refusal, through
# the low-level handle bartcoreMultinomialSampler builds (a bare $ptr/$x/$K
# environment) - the same C entry a real dbartsSampler's own
# $getFitsWithoutOffset() reaches, so this cell pins that backstop without
# needing a full bart2(family = "multinomial") fit.
makeMultinomial <- getFromNamespace("bartcoreMultinomialSampler", "dbarts")
fitsWithoutOffset <- bartcoreFitsWithoutOffset
multinomialHost <- dbarts(x, yGauss, control = samplerControl())
multinomial <- makeMultinomial(multinomialHost, rbinom(n, 2L, 0.5), 3L)
expect_error(
  fitsWithoutOffset(multinomial),
  pattern = "softmax probability channel per category"
)
# the caveat rides the message: predict reports the SAVED samples, not the
# current state, once the sampler keeps trees
expect_error(fitsWithoutOffset(multinomial), pattern = "keepTrees")

# --- the ARITHMETIC gates: recombination through the independent read path --

# single forest: response.scale * getForestFits(1) + response.shift. This
# moves under an internal-scale mutation, which every identity cell above
# survives.
calibrationGauss <- samplerGauss$getCalibration(1L)
forestFitsGauss <- samplerGauss$getForestFits(1L)
expect_equal(dim(forestFitsGauss), c(n, 2L))
expect_true(any(calibrationGauss[, "response.scale"] != 1))
expect_equal(
  sweep(
    sweep(forestFitsGauss, 2L, calibrationGauss[, "response.scale"], `*`),
    2L,
    calibrationGauss[, "response.shift"],
    `+`
  ),
  samplerGauss$getFitsWithoutOffset()
)

# BCF: response.scale * (a * mu + b_z * tau) + response.shift, the
# recombination refuseUndefinedTestFits's own message directs BCF consumers to.
# This is the ONLY cell here that moves when the accessor reports
# forest 0's totals rather than the combiner blend, so it must not be tidied
# away: the identity cells move together and the tests/cpp cell is
# single-forest, where forest 0's totals ARE the combined fit.
amplitudesBcf <- samplerBcf$getForestAmplitudes()
calibrationBcf <- samplerBcf$getCalibration(1L)
muBcf <- samplerBcf$getForestFits(1L)
tauBcf <- samplerBcf$getForestFits(2L)
expect_equal(dim(amplitudesBcf), c(3L, 1L))
expect_equal(dim(muBcf), c(n, 1L))
expect_true(any(tauBcf[, 1L] != 0))
expect_equal(
  calibrationBcf[1L, "response.scale"] *
    (amplitudesBcf[1L, 1L] *
      muBcf[, 1L] +
      ifelse(z == 1L, amplitudesBcf[3L, 1L], amplitudesBcf[2L, 1L]) *
        tauBcf[, 1L]) +
    calibrationBcf[1L, "response.shift"],
  as.vector(samplerBcf$getFitsWithoutOffset())
)
