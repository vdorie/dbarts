source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that works with non-standard models
x <- testData$x
y <- testData$y

k <- c(4, 8)
expect_silent(
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 3L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.test = 5,
    k = k,
    n.threads = 1L,
    resid.prior = chisq(2.5, 0.9)
  )
)

expect_silent(
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 3L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.test = 5,
    k = k,
    n.threads = 1L,
    resid.prior = fixed(2)
  )
)

# xbart's own resid.prior default is the literal chisq() object - it has no
# sigdf/sigquant shorthands for bart2's NULL-triggers-shorthand sentinel to
# build from - so an unsupplied resid.prior is byte-identical to naming the
# default explicitly
expect_identical(
  dbarts::xbart(
    x,
    y,
    n.reps = 2L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.threads = 1L,
    seed = 15L
  ),
  dbarts::xbart(
    x,
    y,
    n.reps = 2L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.threads = 1L,
    seed = 15L,
    resid.prior = chisq()
  )
)

n.trees <- c(5L, 10L)
expect_silent(
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 3L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.test = 5,
    n.trees = n.trees,
    n.threads = 1L
  )
)

# split.probs and dart reach the tree prior, which used to be hardcoded to a
# uniform cgm; both fix or sample the split-variable probabilities while power
# and base stay the swept grid
p <- ncol(x)
expect_silent(
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 2L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.test = 5,
    n.threads = 1L,
    split.probs = c(5, rep(1, p - 1L))
  )
)
expect_silent(
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 2L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.test = 5,
    n.threads = 1L,
    dart = TRUE
  )
)
expect_silent(
  dbarts::xbart(
    x,
    y,
    method = "k-fold",
    n.reps = 2L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.test = 5,
    n.threads = 1L,
    dart = dbarts::dbartsPriors$dart(a = 1)
  )
)
expect_error(
  dbarts::xbart(
    x,
    y,
    n.reps = 1L,
    n.threads = 1L,
    dart = TRUE,
    split.probs = c(5, rep(1, p - 1L))
  ),
  pattern = paste0(
    "'split.probs' cannot be combined with 'dart': a DART prior samples ",
    "its split probabilities"
  )
)
expect_error(
  dbarts::xbart(x, y, n.reps = 1L, n.threads = 1L, dart = "yes"),
  pattern = "'dart' must be TRUE, FALSE, or a prior created by dbartsPriors"
)
expect_error(
  dbarts::xbart(x, y, n.reps = 1L, n.threads = 1L, split.probs = c(1, 2)),
  pattern = "does not equal number of columns"
)

# The four flat knobs and tree.prior. Fits stay tiny; each assertion shows a
# knob demonstrably reaches the fit, a bad value is refused, or tree.prior's
# grid-override/collision rules hold.
quickXbart <- function(...) {
  dbarts::xbart(
    x,
    y,
    n.reps = 2L,
    n.samples = 6L,
    n.burn = c(10L, 5L, 1L),
    n.threads = 1L,
    ...
  )
}

# each knob moves the fit at a fixed seed: n.cuts changes the cut grid,
# useQuantiles changes rule placement, n.thin changes which sweeps are
# recorded, storage changes the residual's working precision
expect_false(identical(
  quickXbart(seed = 11L, n.cuts = 5L),
  quickXbart(seed = 11L, n.cuts = 90L)
))
expect_false(identical(
  quickXbart(seed = 12L, useQuantiles = FALSE),
  quickXbart(seed = 12L, useQuantiles = TRUE)
))
expect_false(identical(
  quickXbart(seed = 13L, n.thin = 1L),
  quickXbart(seed = 13L, n.thin = 3L)
))
# storage: xbart always creates its per-fold samplers over a shared data
# handle (bartcoreDataHandle/bartcoreSamplerFromHandle), a path the engine
# keeps fp64-only regardless of family (reduced-precision-storage.md sec 6);
# "single" reaching the fit is demonstrated by it reaching that refusal
# rather than being silently dropped, while the "double" default still runs
expect_silent(quickXbart(seed = 14L, storage = "double"))
expect_error(
  quickXbart(seed = 14L, storage = "single"),
  pattern = "storage = \"single\" is not supported"
)

# shape-check refusals (f2), mirroring dbartsControl's own validity messages
expect_error(
  quickXbart(n.cuts = 0L),
  pattern = "'n.cuts' must be a positive integer"
)
expect_error(
  quickXbart(n.cuts = NA_integer_),
  pattern = "'n.cuts' must be a positive integer"
)
expect_error(
  quickXbart(useQuantiles = NA),
  pattern = "'useQuantiles' must be TRUE/FALSE"
)
expect_error(
  quickXbart(n.thin = -1L),
  pattern = "'n.thin' must be a positive integer"
)
expect_error(quickXbart(storage = "half"), pattern = "should be one of")

# tree.prior grid-override (f4): power/base are xbart's grid axes, so a
# supplied cgm's own power/base are overridden every cell regardless of what
# it names - a call that sweeps the grid directly is byte-identical to the
# same sweep routed through an object naming DIFFERENT power/base
gridDirect <- quickXbart(seed = 21L, power = c(1.5, 3), base = c(0.6, 0.9))
gridViaObject <- quickXbart(
  seed = 21L,
  power = c(1.5, 3),
  base = c(0.6, 0.9),
  tree.prior = dbarts::dbartsPriors$cgm(power = 10, base = 0.1)
)
expect_identical(gridDirect, gridViaObject)

# dart-prior-via-tree.prior equals dart = <the same object>
dartObj <- dbarts::dbartsPriors$dart(a = 1)
expect_identical(
  quickXbart(seed = 22L, dart = dartObj),
  quickXbart(seed = 22L, tree.prior = dartObj)
)

# collision refusals: dart/split.probs would only duplicate what a supplied
# tree.prior already specifies; power/base/k stay legal alongside it, since
# they are grid axes, not duplicates (unlike bart2's tree.prior)
expect_error(
  quickXbart(tree.prior = dbarts::dbartsPriors$cgm(), dart = TRUE),
  pattern = paste0(
    "'tree.prior' cannot be combined with 'dart': supply the prior either ",
    "as an object or through its shorthand arguments, not both"
  )
)
expect_error(
  quickXbart(
    tree.prior = dbarts::dbartsPriors$cgm(),
    split.probs = c(5, rep(1, p - 1L))
  ),
  pattern = "'tree.prior' cannot be combined with 'split.probs': supply the"
)
expect_silent(quickXbart(
  tree.prior = dbarts::dbartsPriors$cgm(),
  power = c(1, 2),
  base = 0.9,
  k = 2
))

rm(quickXbart, gridDirect, gridViaObject, dartObj)

rm(p, n.trees, k, y, x)

rm(testData)

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that runs with binary data
x <- testData$X
z <- testData$Z

n.reps <- 3L
power <- c(1.5, 2)

xval <- dbarts::xbart(
  x,
  z,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = n.reps,
  power = power,
  n.threads = 2L
)

expect_inherits(xval, "matrix")
expect_equal(dim(xval), c(n.reps, length(power)))
expect_true(!anyNA(xval))
expect_equal(
  dimnames(xval),
  list(
    rep = NULL,
    power = as.character(power)
  )
)

rm(xval)

# the binary grid default is the same FIXED k the continuous arm defaults to
# (2), not bart2's chi hyperprior: a hyperprior is held rather than swept and
# is drawn every sweep, so it would collapse the k axis onto a single cell.
# A default binary run is therefore identical to the same run naming k = 2 -
# a hyperprior reaching a cell would draw and move the stream - and reports
# the k axis under drop = FALSE exactly as a continuous run does
binaryCV <- function(...) {
  dbarts::xbart(
    x,
    z,
    n.samples = 5L,
    n.burn = c(3L, 2L),
    method = "k-fold",
    n.test = 5,
    n.reps = 1L,
    n.threads = 1L,
    seed = 41L,
    drop = FALSE,
    ...
  )
}
binaryDefault <- binaryCV()
expect_identical(binaryDefault, binaryCV(k = 2))
expect_identical(dimnames(binaryDefault)[["k"]], "2")
expect_identical(dimnames(binaryCV(family = "gaussian"))[["k"]], "2")
rm(binaryCV, binaryDefault)

# family routes through to the folds: logistic and forced-gaussian fit 0/1
# responses, binary families reject continuous ones
xval.logistic <- dbarts::xbart(
  x,
  z,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = n.reps,
  n.threads = 1L,
  family = "logistic"
)
expect_equal(dim(xval.logistic), NULL)
expect_true(!anyNA(xval.logistic))

xval.gaussian <- dbarts::xbart(
  x,
  z,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = n.reps,
  n.threads = 1L,
  family = "gaussian"
)
expect_true(!anyNA(xval.gaussian))

expect_error(
  dbarts::xbart(
    x,
    rnorm(nrow(x)),
    n.reps = 1L,
    n.threads = 1L,
    family = "probit"
  ),
  pattern = "requires a response coded 0/1"
)

# the tree.prior grid-axis-override rule (f4) holds under a binary family
# too: a supplied cgm's own power/base are still overridden by the swept
# grid every cell, regardless of control@binary
gridDirect.bin <- dbarts::xbart(
  x,
  z,
  family = "probit",
  n.reps = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = c(3L, 2L, 1L),
  seed = 31L,
  power = c(1.5, 3),
  base = c(0.6, 0.9)
)
gridViaObject.bin <- dbarts::xbart(
  x,
  z,
  family = "probit",
  n.reps = 1L,
  n.threads = 1L,
  n.samples = 5L,
  n.burn = c(3L, 2L, 1L),
  seed = 31L,
  power = c(1.5, 3),
  base = c(0.6, 0.9),
  tree.prior = dbarts::dbartsPriors$cgm(power = 10, base = 0.1)
)
expect_identical(gridDirect.bin, gridViaObject.bin)
rm(gridDirect.bin, gridViaObject.bin)

rm(xval.logistic, xval.gaussian, power, n.reps, z, x)

rm(testData)
