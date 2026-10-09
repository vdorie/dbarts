# A copy installs the stored state as a reload does: a state whose saved-tree
# store has another capacity than the control now names (n.samples edited
# after the state was stored) is taken with its own capacity, not refused. A
# sampler with no stored state (updateState = FALSE) has its current state
# read for the copy, so every copy continues its source.

set.seed(0)
n <- 80L
x <- matrix(rnorm(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n)

make <- function(seed = 3L, updateState = TRUE, keepTrees = TRUE, ...) {
  dbarts::dbarts(
    x,
    y,
    control = dbarts::dbartsControl(
      keepTrees = keepTrees,
      updateState = updateState,
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 10L,
      n.samples = 20L,
      n.burn = 5L,
      verbose = FALSE,
      seed = seed
    ),
    ...
  )
}

# the control slot edited in place, the state stored before it
sampler <- make()
sampler$run(2L, 5L)
sampler$control@n.samples <- 8L
dupe <- sampler$copy()
expect_true(inherits(dupe, "dbartsSampler"))
expect_identical(dupe$state, sampler$state)
expect_equal(dupe$run(0L, 3L)$train, sampler$run(0L, 3L)$train)

# the supported route, setControl, after the store
sampler <- make()
sampler$run(2L, 5L)
newControl <- sampler$control
newControl@n.samples <- 8L
sampler$setControl(newControl)
dupe <- sampler$copy()
expect_identical(dupe$state, sampler$state)
expect_equal(dupe$run(0L, 3L)$train, sampler$run(0L, 3L)$train)

# the copy caches the very state it installed, also when the install had to
# convert it (here, into the units of a response re-anchored after the store)
sampler <- make()
sampler$run(1L, 3L)
sampler$setResponse(y * 3 + 1, updateScale = TRUE, updateState = FALSE)
dupe <- sampler$copy()
expect_identical(dupe$state, sampler$state)

# a named leaf sd stays what it was named after the response is replaced,
# also in the copy
sdOf <- function(s) {
  .Call(dbarts:::C_dbarts_bartcore_getLeafPrior, s$pointer, 0L)[[
    1L,
    "prior.scale"
  ]]
}
sampler <- make(leaf.prior = dbarts:::normal(sd = 0.9))
sampler$setResponse(y * 4 + 2)
sampler$run(1L, 3L)
dupe <- sampler$copy()
expect_equal(sdOf(dupe), sdOf(sampler))

# the per-forest weights and the active-row mask ride neither the state nor
# the data: the copy takes them from the source's mirrors
w <- runif(n, 0.2, 2)
z <- rbinom(n, 1L, 0.5)
build2 <- function(updateState) {
  dbarts::dbarts(
    x,
    y + z,
    forests = list(
      dbarts:::forest(),
      dbarts:::forest(basis = ~ factor(z))
    ),
    control = dbarts::dbartsControl(
      updateState = updateState,
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 10L,
      n.samples = 6L,
      verbose = FALSE,
      seed = 5L
    )
  )
}
for (updateState in c(TRUE, FALSE)) {
  sampler <- build2(updateState)
  sampler$setForestWeights(2L, w)
  sampler$run(1L, 2L)
  dupe <- sampler$copy()
  expect_equal(dupe$run(0L, 3L)$train, sampler$run(0L, 3L)$train)

  sampler <- make(updateState = updateState)
  sampler$setActiveRows(rep(c(TRUE, FALSE), length.out = n))
  sampler$run(1L, 2L)
  dupe <- sampler$copy()
  expect_equal(dupe$run(0L, 3L)$train, sampler$run(0L, 3L)$train)
}

# updateState = FALSE: seeded and unseeded, untouched and run, the copy
# continues its source, and the source continues as if never copied
for (seed in list(3L, NULL)) {
  for (burn in c(0L, 4L)) {
    sampler <- make(seed = seed, updateState = FALSE)
    if (burn > 0L) {
      sampler$run(burn, 3L)
    }
    expect_null(sampler$state)
    dupe <- sampler$copy()
    expect_null(sampler$state)
    expect_null(dupe$state)
    d <- dupe$run(2L, 5L)
    s <- sampler$run(2L, 5L)
    expect_equal(d$train, s$train)
    expect_equal(d$sigma, s$sigma)
  }
}

twinA <- make(updateState = FALSE)
twinB <- make(updateState = FALSE)
twinA$run(3L, 3L)
twinB$run(3L, 3L)
twinB$copy()
expect_identical(twinB$run(2L, 5L)$train, twinA$run(2L, 5L)$train)

# an unchanged updateState = TRUE copy still matches
sampler <- make()
sampler$run(3L, 3L)
dupe <- sampler$copy()
expect_equal(dupe$run(2L, 5L)$train, sampler$run(2L, 5L)$train)
