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
