# A draw the sampler has kept is the function that was drawn, and stays so:
# predict returns the same values before and after setData and before and
# after the response range is re-derived. A forest seeded from another fit
# starts at that fit's function even when the two standardize a leaf covariate
# differently, and the seeded sampler keeps its own standardization. Every
# identity is read from the engine - predictions at new rows and at the
# training rows, live fits, a freshly stored state - on samplers of several
# trees and two chains.

set.seed(20261007L)
n <- 120L
x <- cbind(x1 = runif(n), x2 = rnorm(n))
y <- 2 * x[, 2L] + sin(4 * x[, 1L]) + rnorm(n, 0, 0.2)
x.new <- cbind(x1 = runif(40L, -0.2, 1.2), x2 = rnorm(40L, 0, 2))
# the same rows read on another scale of the leaf covariate
x.stretched <- x
x.stretched[, 2L] <- 3 * x[, 2L] + 5

control <- function(
  n.chains = 2L,
  keepTrees = TRUE,
  seed = 7L,
  n.samples = 5L,
  n.burn = 30L
) {
  dbartsControl(
    n.chains = n.chains,
    n.threads = 1L,
    n.trees = 8L,
    n.samples = n.samples,
    n.burn = n.burn,
    keepTrees = keepTrees,
    updateState = FALSE,
    seed = seed
  )
}
# kept draws at rows the sampler never saw and at the ones it was fit to
keptDraws <- function(sampler) {
  drawn <- lapply(list(x.new, x), sampler$predict)
  if (is.list(drawn[[1L]])) {
    drawn <- unlist(drawn, recursive = FALSE)
  }
  drawn
}
worstGap <- function(a, b) {
  max(abs(unlist(a) - unlist(b)) / pmax(1, abs(unlist(a))))
}
liveFits <- function(sampler) {
  as.vector(sampler$getFitsWithoutOffset())
}
standardization <- function(sampler) {
  sampler$storeState()
  lapply(sampler$state, function(chain) {
    chain$forests[[1L]][c(
      "leaf.covariate.center",
      "leaf.covariate.scale",
      "leaf.lengthscales"
    )]
  })
}
liveTrees <- function(sampler) {
  sampler$storeState()
  lapply(sampler$state, function(chain) {
    chain$forests[[1L]][c("tree.vars", "tree.values", "tree.params")]
  })
}
cutGrid <- function(sampler) {
  sampler$storeState()
  attr(sampler$state, "cutPoints")
}

## --- 1. a leaf covariate without spread ------------------------------------
## Its scale is stored as NA beside a finite centre (it was written as the
## placeholder 1 the engine divides by), and such a sampler restores, copies
## and reloads.
x.flat <- cbind(x1 = x[, 1L], x2 = 1000)
x.flat.new <- cbind(x1 = runif(20L), x2 = 1000 + rnorm(20L))
flat <- dbarts(
  x.flat,
  y,
  control = control(),
  leaf.prior = linear(columns = 2L)
)
invisible(flat$run())
for (chain in standardization(flat)) {
  expect_identical(chain$leaf.covariate.center, 1000)
  expect_identical(chain$leaf.covariate.scale, NA_real_)
}
flat.kept <- flat$predict(x.flat.new)
expect_true(all(is.finite(flat.kept)))
# the covariate is divided by 1: moving it by one moves a draw by the sum of
# the slopes its rows reach, so two steps move it by twice one step
flat.step <- function(by) {
  shifted <- x.flat.new
  shifted[, 2L] <- shifted[, 2L] + by
  flat$predict(shifted) - flat.kept
}
expect_true(max(abs(flat.step(1))) > 1e-3)
expect_true(max(abs(flat.step(2) - 2 * flat.step(1))) < 1e-10)

flat.twin <- dbarts(
  x.flat,
  y,
  control = control(seed = 8L),
  leaf.prior = linear(columns = 2L)
)
flat$storeState()
expect_true(flat.twin$setState(flat$state))
expect_identical(flat.twin$predict(x.flat.new), flat.kept)
expect_identical(standardization(flat.twin), standardization(flat))

path <- tempfile(fileext = ".rds")
saveRDS(flat, path)
flat.reloaded <- readRDS(path)
unlink(path)
flat.copied <- flat$copy()
for (route in list(flat.copied, flat.reloaded)) {
  expect_identical(route$predict(x.flat.new), flat.kept)
  expect_identical(standardization(route), standardization(flat))
}
expect_equal(
  flat.copied$run(0L, 2L)$train,
  flat$run(0L, 2L)$train,
  tolerance = 1e-10
)
# the mark is the one value that is not a number a scale may hold
flat$storeState()
bad <- flat$state
bad[[2L]]$forests[[1L]]$leaf.covariate.scale <- 0
expect_error(flat$setState(bad), "not consistent with this sampler")
bad[[2L]]$forests[[1L]]$leaf.covariate.scale <- -1
expect_error(flat$setState(bad), "not consistent with this sampler")
bad[[2L]]$forests[[1L]]$leaf.covariate.scale <- Inf
expect_error(flat$setState(bad), "not consistent with this sampler")

# setData onto a covariate with spread. No observation informed a slope drawn
# against the constant, so the live ones are dropped and each leaf keeps the
# value it had there; the kept draws are only replayed, and stay the
# functions they were
x.varying <- cbind(x1 = x[, 1L], x2 = rnorm(n, 1000, 300))
slopes <- function(sampler) {
  unlist(lapply(liveTrees(sampler), `[[`, "tree.params"))
}
flat.kept <- flat$predict(x.flat.new)
flat.live <- liveFits(flat)
flat.donor <- flat$copy()
expect_true(max(abs(slopes(flat))) > 1e-3)
flat$setData(dbartsData(x.varying, y))
for (chain in standardization(flat)) {
  expect_true(chain$leaf.covariate.scale > 100)
}
expect_true(all(slopes(flat) == 0))
expect_true(max(abs(liveFits(flat) - flat.live)) < 1e-12)
expect_true(max(abs(flat$predict(x.flat.new) - flat.kept)) < 1e-10)
# and a warm start from such a donor into a sampler whose covariate varies
varied <- dbarts(
  x.varying,
  y,
  control = control(keepTrees = FALSE, seed = 9L),
  leaf.prior = linear(columns = 2L)
)
varied$installTrees(flat.donor, samples = c(5L, 10L))
expect_true(all(slopes(varied) == 0))
expect_true(max(abs(liveFits(varied) - flat.live)) < 1e-12)

## --- 2. kept draws across setData -------------------------------------------
## The same raw points, replayed before and after a setData whose leaf
## covariate is 3 x + 5: a linear leaf's kept draws came back off by about 8.
for (n.chains in 1:2) {
  linear.fit <- dbarts(
    x,
    y,
    control = control(n.chains = n.chains),
    leaf.prior = linear(columns = 2L)
  )
  invisible(linear.fit$run())
  before <- keptDraws(linear.fit)
  centers <- standardization(linear.fit)
  linear.fit$setData(dbartsData(x.stretched, y))
  expect_true(
    abs(
      standardization(linear.fit)[[n.chains]]$leaf.covariate.center -
        centers[[n.chains]]$leaf.covariate.center
    ) >
      4,
    info = n.chains
  )
  expect_true(diff(range(before[[1L]])) > 1, info = n.chains)
  expect_true(worstGap(before, keptDraws(linear.fit)) < 1e-12, info = n.chains)
  expect_true(all(is.finite(linear.fit$run(0L, 2L)$train)), info = n.chains)

  # a constant leaf beside it holds no coefficient to restate
  constant.fit <- dbarts(x, y, control = control(n.chains = n.chains))
  invisible(constant.fit$run())
  before <- keptDraws(constant.fit)
  constant.fit$setData(dbartsData(x.stretched, y))
  expect_identical(keptDraws(constant.fit), before, info = n.chains)
}

# rows appended inside every column's range leave the cut grid alone, so the
# live fit on the rows that stayed is what it was (it moved by about 1)
x.more <- rbind(
  x,
  cbind(
    x1 = runif(40L, min(x[, 1L]), max(x[, 1L])),
    x2 = runif(40L, 1.2, max(x[, 2L]))
  )
)
y.more <- c(y, runif(40L, min(y), max(y)))
appended <- dbarts(
  x,
  y,
  control = control(keepTrees = FALSE),
  leaf.prior = linear(columns = 2L)
)
invisible(appended$run())
grid <- cutGrid(appended)
before <- matrix(liveFits(appended), n)
centers <- standardization(appended)
appended$setData(dbartsData(x.more, y.more))
expect_identical(cutGrid(appended), grid)
expect_true(
  abs(
    standardization(appended)[[2L]]$leaf.covariate.center -
      centers[[2L]]$leaf.covariate.center
  ) >
    0.2
)
after <- matrix(liveFits(appended), n + 40L)[seq_len(n), ]
expect_true(max(abs(after - before)) < 1e-12)

## --- 3. kept draws across a re-anchor ---------------------------------------
## setResponse and setOffset with updateScale = TRUE, and setData with a
## response on another range, returned every kept draw in the new units:
## 3 f + 10 for a response of 3 y + 10.
y.moved <- 3 * y + 10
reanchored <- list(
  constant = function() dbarts(x, y, control = control()),
  student = function() dbarts(x, y, family = "student", control = control()),
  monotone = function() {
    dbarts(x, y, monotone = c(x2 = 1), control = control())
  },
  linear = function() {
    dbarts(x, y, leaf.prior = linear(columns = 2L), control = control())
  },
  variance = function() dbarts(x, y, variance = ~x1, control = control())
)
for (kind in names(reanchored)) {
  calls <- list(
    setResponse = function(sampler) {
      sampler$setResponse(y.moved, updateScale = TRUE)
    },
    setOffset = function(sampler) {
      sampler$setOffset(-5 + x[, 1L], updateScale = TRUE)
    },
    setData = function(sampler) sampler$setData(dbartsData(x, y.moved))
  )
  for (call in names(calls)) {
    info <- paste(kind, call)
    sampler <- reanchored[[kind]]()
    invisible(sampler$run())
    before <- keptDraws(sampler)
    live <- liveFits(sampler)
    calls[[call]](sampler)
    expect_true(worstGap(before, keptDraws(sampler)) < 1e-12, info = info)
    if (call != "setOffset") {
      # the live chain sits at the same internal values, as it always did
      expect_true(max(abs(liveFits(sampler) - (3 * live + 10))) < 1e-10, info)
    }
    expect_true(all(is.finite(sampler$run(0L, 2L)$train)), info = info)
  }
  # holding the range rewrites nothing
  sampler <- reanchored[[kind]]()
  invisible(sampler$run())
  before <- keptDraws(sampler)
  sampler$setResponse(y.moved, updateScale = FALSE)
  expect_identical(keptDraws(sampler), before, info = kind)
}

# a variance forest returns both channels, each held to the identity above
variance.fit <- reanchored$variance()
invisible(variance.fit$run())
expect_identical(names(variance.fit$predict(x.new)), c("mean", "variance"))
expect_identical(length(keptDraws(variance.fit)), 4L)

# aft: the range is that of the log times
time <- exp(y / 3)
status <- rbinom(n, 1L, 0.8)
aft.fit <- dbarts(x, cbind(time, status), family = "aft", control = control())
invisible(aft.fit$run())
before <- keptDraws(aft.fit)
aft.fit$setResponse(aft.fit$data@y + 0.5, updateScale = TRUE)
expect_true(worstGap(before, keptDraws(aft.fit)) < 1e-12)

# nbinom: the shift is the log of the mean count
counts <- rpois(n, exp(0.5 * sin(4 * x[, 1L]) + 0.5 * x[, 2L] + 1))
count.fit <- dbarts(x, counts, family = "nbinom", control = control())
invisible(count.fit$run())
before <- keptDraws(count.fit)
count.fit$setResponse(3L * counts, updateScale = TRUE)
expect_true(worstGap(before, keptDraws(count.fit)) < 1e-12)

# a store the runs have part filled: the draws a later run adds are the new
# range's, and the earlier ones are still what they were
partial <- dbarts(x, y, control = control(n.samples = 6L))
invisible(partial$run(30L, 2L))
before <- keptDraws(partial)
partial$setResponse(y.moved, updateScale = TRUE)
invisible(partial$run(0L, 2L))
after <- partial$predict(x.new)
expect_identical(dim(after), c(40L, 4L, 2L))
expect_true(max(abs(after[, 1:2, ] - before[[1L]])) < 1e-12)
expect_true(abs(mean(after[, 3:4, ]) - mean(3 * before[[1L]] + 10)) < 2)

## --- 4. a warm start between standardizations ------------------------------
## A recipient made on the stretched covariate and then given the donor's
## rows keeps the standardization it was made with. Its seeded fit was off the
## donor's by about 6, and a donor on another grid replaced its standardization.
warm.control <- function(...) control(keepTrees = FALSE, n.burn = 60L, ...)
donor <- dbarts(
  x,
  y,
  control = control(n.burn = 60L, n.samples = 4L),
  leaf.prior = linear(columns = 2L)
)
invisible(donor$run())
donor.grid <- cutGrid(donor)
donor.fits <- matrix(liveFits(donor), n)
offStandard <- function(seed, grid = donor.grid) {
  recipient <- dbarts(
    x.stretched,
    y,
    control = warm.control(seed = seed),
    leaf.prior = linear(columns = 2L)
  )
  invisible(recipient$setPredictor(x, forceUpdate = TRUE))
  recipient$setCutPoints(grid)
  recipient
}

# the donor's grid, another standardization; the donor's live trees are its
# newest kept draw, so the last of the pool seeds both chains
recipient <- offStandard(21L)
own <- standardization(recipient)
expect_identical(cutGrid(recipient), donor.grid)
expect_true(abs(own[[2L]]$leaf.covariate.center - 5) < 1)
donor.kept <- donor$predict(x)
recipient$installTrees(donor, samples = c(4L, 8L))
expect_true(max(abs(matrix(liveFits(recipient), n) - donor.fits)) < 1e-12)
expect_identical(standardization(recipient), own)
expect_true(all(is.finite(recipient$run(0L, 2L)$train)))

# a kept draw as the seed: each chain starts at the draw it was given
recipient <- offStandard(22L)
recipient$installTrees(donor, samples = c(6L, 1L))
seeded <- matrix(liveFits(recipient), n)
expect_true(max(abs(seeded[, 1L] - donor.kept[, 2L, 2L])) < 1e-12)
expect_true(max(abs(seeded[, 2L] - donor.kept[, 1L, 1L])) < 1e-12)
expect_identical(standardization(recipient), own)

# a grid that holds every donor split point and more
finer <- lapply(donor.grid, function(points) {
  sort(c(points, points[-1L] - diff(points) / 2))
})
recipient <- offStandard(23L, finer)
own <- standardization(recipient)
expect_false(identical(cutGrid(recipient), donor.grid))
recipient$installTrees(donor, samples = c(4L, 8L))
expect_true(max(abs(matrix(liveFits(recipient), n) - donor.fits)) < 1e-12)
expect_identical(standardization(recipient), own)

# equal standardizations: the donor's coefficients arrive as they are
twin <- dbarts(
  x,
  y,
  control = warm.control(seed = 24L),
  leaf.prior = linear(columns = 2L)
)
expect_identical(standardization(twin), standardization(donor))
twin$installTrees(donor, samples = c(4L, 8L))
expect_identical(liveTrees(twin), liveTrees(donor))

# a donor record a hand has broken is refused, the recipient as it was
recipient <- offStandard(25L)
before <- liveTrees(recipient)
donor$storeState()
broken <- donor$state
broken[[2L]]$forests[[1L]]$leaf.covariate.scale <- 0
expect_error(
  recipient$installTrees(broken),
  "malformed parameters in warm-start donor"
)
expect_identical(liveTrees(recipient), before)

# through bart: one donor forest written under two standardizations, the
# second restated here by the rule (slope s' / s, and the intercept gaining
# slope (m' - m) / s), is one function, so it seeds the same draws
restated <- function(state, shift, stretch) {
  for (chain in seq_along(state)) {
    forest <- state[[chain]]$forests[[1L]]
    values <- readBin(forest$tree.values, "double", length(forest$tree.vars))
    leaves <- forest$tree.vars < 0L
    values[leaves] <- values[leaves] +
      forest$tree.params * shift / forest$leaf.covariate.scale
    forest$tree.values <- writeBin(values, raw())
    forest$tree.params <- forest$tree.params * stretch
    forest$leaf.covariate.center <- forest$leaf.covariate.center + shift
    forest$leaf.covariate.scale <- forest$leaf.covariate.scale * stretch
    state[[chain]]$forests[[1L]] <- forest
  }
  state
}
donor.fit <- bart(
  x,
  y,
  leaf.prior = linear(columns = 2L),
  n.trees = 8L,
  n.samples = 4L,
  n.burn = 60L,
  n.chains = 1L,
  n.threads = 1L,
  keepSampler = TRUE,
  verbose = FALSE,
  seed = 1L
)
donor.fit$fit$storeState()
recorded <- donor.fit$fit$state
moved <- restated(recorded, shift = 0.7, stretch = 2.5)
warmed <- function(state) {
  bart(
    x.more,
    y.more,
    leaf.prior = linear(columns = 2L),
    n.trees = 8L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 2L,
    n.threads = 1L,
    warm.start = state,
    verbose = FALSE,
    seed = 2L
  )$yhat.train
}
from.recorded <- warmed(recorded)
expect_true(diff(range(from.recorded)) > 1)
expect_true(max(abs(warmed(moved) - from.recorded)) < 1e-8)
# and before any sweep, through installTrees
recipient <- offStandard(26L)
recipient$installTrees(recorded)
before <- liveFits(recipient)
recipient <- offStandard(26L)
recipient$installTrees(moved)
expect_true(max(abs(liveFits(recipient) - before)) < 1e-12)

## --- 5. gp leaves at a warm start -------------------------------------------
## Per-row fits hold no kernel: copied on the donor's grid, zero on another,
## and the recipient's centre, scale and lengthscale its own on both.
gp.donor <- dbarts(
  x,
  y,
  control = warm.control(),
  leaf.prior = gp(columns = 2L)
)
invisible(gp.donor$run())
gpOffStandard <- function(seed, grid) {
  recipient <- dbarts(
    x.stretched,
    y,
    control = warm.control(seed = seed),
    leaf.prior = gp(columns = 2L)
  )
  invisible(recipient$setPredictor(x, forceUpdate = TRUE))
  recipient$setCutPoints(grid)
  recipient
}
recipient <- gpOffStandard(31L, cutGrid(gp.donor))
own <- standardization(recipient)
expect_false(identical(own, standardization(gp.donor)))
recipient$installTrees(gp.donor)
expect_true(max(abs(liveFits(recipient) - liveFits(gp.donor))) < 1e-12)
expect_identical(standardization(recipient), own)

recipient <- gpOffStandard(32L, finer)
own <- standardization(recipient)
recipient$installTrees(gp.donor)
expect_identical(standardization(recipient), own)
expect_true(diff(range(liveFits(recipient))) == 0)

## --- 6. live coefficients across a covariate without spread ----------------
## A constant leaf covariate is centred at its one value and divided by 1, so
## a leaf reads zero for it on every training row and the intercept alone is
## the fit: between two constants a live leaf does not move, off a constant
## its intercept stays, onto one it takes the function's value there. Whether
## a covariate is such a column is read from the values and the centre and
## scale in force, so one given values, moved or made constant after creation
## converts as what it then is. Each case reads the fitted function at the
## training rows and at appended or fresh ones before and after, to 1e-10
## relative to max(1, |value|), and requires the next draws to stay in the
## range of a copy taken before the call, widened by twice its width either
## way.
tolerance <- 1e-10
newest <- c(5L, 10L) # each chain's last kept draw, which is its live forest
withLeafColumn <- function(rows, value) cbind(x1 = rows[, 1L], x2 = value)
madeOn <- function(rows, seed = 7L, response = y, ...) {
  dbarts(
    rows,
    response,
    control = control(n.burn = 60L, seed = seed, ...),
    leaf.prior = linear(columns = 2L)
  )
}
liveFunction <- function(sampler, rows) sampler$predict(rows)[, 5L, ]
liveMatrix <- function(sampler) matrix(liveFits(sampler), ncol = 2L)
nextDrawsWithin <- function(sampler, bounds) {
  drawn <- range(sampler$run(0L, 3L)$train)
  all(is.finite(drawn)) &&
    drawn[1L] >= bounds[1L] - 2 * diff(bounds) &&
    drawn[2L] <= bounds[2L] + 2 * diff(bounds)
}
drawsLikeTwin <- function(sampler, twin) {
  nextDrawsWithin(sampler, range(twin$run(0L, 3L)$train))
}
old.rows <- seq_len(n)
added <- x.more[-old.rows, ]

# neither side has spread: one constant replaced by another, at a value whose
# mean is exact and at one whose mean rounds. The live fit ran to thousands.
for (constant in c(1000, 0.1)) {
  x.constant <- withLeafColumn(x, constant)
  x.fresh <- withLeafColumn(x.new, constant + x.new[, 2L])
  constant.fit <- madeOn(x.constant)
  invisible(constant.fit$run())
  for (chain in standardization(constant.fit)) {
    expect_identical(chain$leaf.covariate.center, constant, info = constant)
    expect_identical(chain$leaf.covariate.scale, NA_real_, info = constant)
  }
  expect_true(max(abs(slopes(constant.fit))) > 1e-3, info = constant)
  twin <- constant.fit$copy()
  donor <- constant.fit$copy()
  live <- liveMatrix(constant.fit)
  trees <- liveTrees(constant.fit)
  kept <- constant.fit$predict(x.fresh)
  constant.fit$setData(dbartsData(withLeafColumn(x.more, 2 * constant), y.more))
  for (chain in standardization(constant.fit)) {
    expect_identical(chain$leaf.covariate.center, 2 * constant, info = constant)
    expect_identical(chain$leaf.covariate.scale, NA_real_, info = constant)
  }
  expect_true(
    worstGap(live, liveMatrix(constant.fit)[old.rows, ]) < tolerance,
    info = constant
  )
  expect_identical(liveTrees(constant.fit), trees, info = constant)
  expect_true(
    worstGap(kept, constant.fit$predict(x.fresh)) < tolerance,
    info = constant
  )
  expect_true(drawsLikeTwin(constant.fit, twin), info = constant)

  # the same pair through a warm start, and through bart, which runs a sweep
  # before anything can be read and is held to the donor's range
  recipient <- madeOn(withLeafColumn(x, 2 * constant), 41L, keepTrees = FALSE)
  recipient$installTrees(donor, samples = newest)
  expect_true(
    worstGap(live, liveMatrix(recipient)) < tolerance,
    info = constant
  )
  donor$storeState()
  seeded <- bart(
    withLeafColumn(x, 2 * constant),
    y,
    leaf.prior = linear(columns = 2L),
    n.trees = 8L,
    n.samples = 3L,
    n.burn = 0L,
    n.chains = 2L,
    n.threads = 1L,
    warm.start = donor$state,
    verbose = FALSE,
    seed = 2L
  )$yhat.train
  expect_true(drawsLikeTwin(recipient, donor), info = constant)
  theirs <- range(live)
  expect_true(
    min(seeded) >= theirs[1L] - 2 * diff(theirs) &&
      max(seeded) <= theirs[2L] + 2 * diff(theirs),
    info = constant
  )
}

# a placeholder of zeros given values after creation: its slopes are informed
# from then on and convert by the formula. They were all set to zero, and the
# live fit on the old rows moved by about 6.
placeholder <- function(seed = 7L, ...) {
  sampler <- madeOn(withLeafColumn(x, 0), seed, ...)
  expect_true(sampler$setPredictor(x[, 2L], 2L))
  invisible(sampler$run())
  sampler
}
given <- placeholder()
for (chain in standardization(given)) {
  expect_identical(chain$leaf.covariate.center, 0)
  expect_identical(chain$leaf.covariate.scale, 1)
}
expect_true(max(abs(slopes(given))) > 0.1)
twin <- given$copy()
donor <- given$copy()
copied <- given$copy()
live <- liveMatrix(given)
before <- liveFunction(given, x.more)
expect_true(worstGap(live, before[old.rows, ]) < tolerance)
kept <- keptDraws(given)
given$setData(dbartsData(x.more, y.more))
expect_true(standardization(given)[[2L]]$leaf.covariate.scale != 1)
expect_true(max(abs(slopes(given))) > 0.1)
# the old rows and the appended ones
expect_true(worstGap(before, liveMatrix(given)) < tolerance)
expect_true(worstGap(kept, keptDraws(given)) < tolerance)
expect_true(drawsLikeTwin(given, twin))
# a copy holds what the sampler does, so it converts the same way
copied$setData(dbartsData(x.more, y.more))
expect_true(worstGap(before, liveMatrix(copied)) < tolerance)
# as a donor into a sampler made on the values
recipient <- madeOn(x, 42L, keepTrees = FALSE)
recipient$installTrees(donor, samples = newest)
expect_true(worstGap(live, liveMatrix(recipient)) < tolerance)
expect_true(max(abs(slopes(recipient))) > 0.1)
expect_true(drawsLikeTwin(recipient, donor))
# and through bart: the donor's record and the same forest restated by the
# rule under a centre and scale of its own are one function
given <- placeholder(43L, keepTrees = FALSE)
given$storeState()
from.recorded <- warmed(given$state)
expect_true(diff(range(from.recorded)) > 1)
expect_true(
  max(abs(warmed(restated(given$state, 0.7, 2.5)) - from.recorded)) < 1e-8
)

# made on one constant and moved to another: every row reads the same value
# off the centre, so the function is kept where the next data spreads the
# covariate about the new constant. The slope's term was dropped: off by 7.
moved <- madeOn(withLeafColumn(x, 1000))
expect_true(moved$setPredictor(rep(2000, n), 2L))
invisible(moved$run())
for (chain in standardization(moved)) {
  expect_identical(chain$leaf.covariate.center, 1000)
  expect_identical(chain$leaf.covariate.scale, 1)
}
x.spread <- withLeafColumn(x.more, 2000 + x.more[, 2L])
before <- liveFunction(moved, x.spread)
moved$setData(dbartsData(x.spread, y.more))
expect_true(worstGap(before, liveMatrix(moved)) < tolerance)

# values made constant after creation keep the centre and scale they were
# given until a replacement re-derives them; onto the constant the slopes are
# dropped and each leaf takes its function's value there
held <- madeOn(x)
own <- standardization(held)
expect_true(held$setPredictor(rep(0.25, n), 2L))
invisible(held$run())
expect_identical(standardization(held), own)
expect_false(own[[2L]]$leaf.covariate.scale %in% c(1, NA))
x.held <- withLeafColumn(x.more, 0.25)
before <- liveFunction(held, x.held)
held$setData(dbartsData(x.held, y.more))
for (chain in standardization(held)) {
  expect_identical(chain$leaf.covariate.center, 0.25)
  expect_identical(chain$leaf.covariate.scale, NA_real_)
}
expect_true(all(slopes(held) == 0))
expect_true(worstGap(before, liveMatrix(held)) < tolerance)
# held at its centre exactly, the covariate reads zero under a real scale,
# which is not the placeholder: the scale stays in the state, so a copy
# predicts off the centre what the sampler does, and the coefficients convert
# by the formula
centred <- madeOn(x, 44L)
expect_true(centred$setPredictor(rep(own[[2L]]$leaf.covariate.center, n), 2L))
invisible(centred$run())
expect_identical(standardization(centred), own)
expect_identical(centred$copy()$predict(x.new), centred$predict(x.new))
before <- liveFunction(centred, x.more)
centred$setData(dbartsData(x.more, y.more))
expect_true(worstGap(before, liveMatrix(centred)) < tolerance)

# rows appended that give a constant covariate spread, at a constant whose
# mean rounds: the slopes are dropped and each leaf keeps the value it had at
# the constant, on the old rows and on the appended ones. The live fit ran
# to 1e16. The appended rows bring the covariate a copy never saw, so the
# next draws are held to the response's range instead of a copy's.
grown <- madeOn(withLeafColumn(x, 0.1))
invisible(grown$run())
live <- liveMatrix(grown)
before <- liveFunction(grown, withLeafColumn(x.more, 0.1))
expect_true(worstGap(live, before[old.rows, ]) < tolerance)
x.grown <- rbind(
  withLeafColumn(x, 0.1),
  withLeafColumn(added, 0.1 + added[, 2L])
)
x.fresh <- withLeafColumn(x.new, 0.1 + x.new[, 2L])
kept <- grown$predict(x.fresh)
grown$setData(dbartsData(x.grown, y.more))
expect_true(standardization(grown)[[2L]]$leaf.covariate.scale > 0.1)
expect_true(all(slopes(grown) == 0))
expect_true(worstGap(before, liveMatrix(grown)) < tolerance)
expect_true(worstGap(kept, grown$predict(x.fresh)) < tolerance)
expect_true(nextDrawsWithin(grown, range(y.more)))
