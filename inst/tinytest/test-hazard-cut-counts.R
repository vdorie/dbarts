# A hazard fit's appended period column takes n.cuts only by its own name, and
# is otherwise at the default 100, at every door that reaches a hazard fit.

set.seed(81)
n <- 80L
x <- cbind(x1 = runif(n), x2 = runif(n), x3 = runif(n))
time <- ceiling(rexp(n) * 6)
status <- rbinom(n, 1L, 0.7)
Y <- cbind(time, status)
frame <- data.frame(x, time = time, status = status)
form <- survival::Surv(time, status) ~ x1 + x2 + x3
run <- function(...) {
  dbartsControl(
    n.samples = 4L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    updateState = FALSE,
    ...
  )
}
# a count written on a control is read with a warning, which is not what this
# file tests
written <- function(...) {
  suppressWarnings(dbartsControl(
    n.trees = 5L,
    n.samples = 4L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE,
    updateState = FALSE,
    ...
  ))
}

doors <- list(
  dbarts = function(v) {
    dbarts(x, Y, family = "hazard", n.cuts = v, control = run(), n.trees = 5L)
  },
  control = function(v) {
    dbarts(x, Y, family = "hazard", control = written(n.cuts = v))
  },
  formula = function(v) {
    dbarts(
      form,
      frame,
      family = "hazard",
      n.cuts = v,
      control = run(),
      n.trees = 5L
    )
  },
  bart = function(v) {
    bart(
      x,
      Y,
      family = "hazard",
      n.cuts = v,
      n.trees = 5L,
      n.samples = 4L,
      n.burn = 0L,
      n.chains = 1L,
      n.threads = 1L,
      verbose = FALSE,
      keepSampler = TRUE
    )$fit
  }
)

for (door in names(doors)) {
  read <- function(v) doors[[door]](v)$data@n.cuts
  # a refusal of a count that should be taken fails the expectation
  got <- function(v) tryCatch(read(v), error = conditionMessage)
  expect_identical(got(c(5L, 7L, 9L)), c(5L, 7L, 9L, 100L), info = door)
  expect_identical(got(5L), c(5L, 5L, 5L, 100L), info = door)
  expect_identical(
    got(c(x1 = 5L, x2 = 7L, x3 = 9L)),
    c(5L, 7L, 9L, 100L),
    info = door
  )
  expect_identical(
    got(c(x1 = 5L, .default = 9L)),
    c(5L, 9L, 9L, 100L),
    info = door
  )
  expect_identical(
    got(c(period = 5L, .default = 5L)),
    c(5L, 5L, 5L, 5L),
    info = door
  )
  expect_identical(
    got(c(x1 = 5L, x2 = 7L, x3 = 9L, period = 11L)),
    c(5L, 7L, 9L, 11L),
    info = door
  )
  # four unnamed values stay by position
  expect_identical(got(c(5L, 7L, 9L, 11L)), c(5L, 7L, 9L, 11L), info = door)
  if (door != "bart") {
    expect_identical(
      doors[[door]](NULL)$data@n.cuts,
      rep(100L, 4L),
      info = door
    )
  }
  # the refusals name the caller's predictors, never the period
  expect_error(
    read(c(x2 = 7L)),
    "missing values for columns 'x1', 'x3'",
    info = door
  )
  expect_error(
    read(c(5L, 7L)),
    "'n.cuts' has 2 values but the model has 3",
    info = door
  )
  expect_error(
    read(1:5),
    "'n.cuts' has 5 values but the model has 3",
    info = door
  )
}

# a count that agrees with the sampler's is the sampler's
sampler <- dbarts(
  x,
  Y,
  family = "hazard",
  n.cuts = 5L,
  control = run(),
  n.trees = 5L
)
expect_silent(sampler$setControl(written(n.cuts = 5L)))
expect_silent(sampler$setControl(written(n.cuts = c(5L, 5L, 5L, 100L))))
expect_error(
  sampler$setControl(written(n.cuts = c(period = 5L, .default = 5L))),
  "changing 'n.cuts' is not available on an existing sampler"
)
expect_identical(sampler$copy()$data@n.cuts, c(5L, 5L, 5L, 100L))
reloaded <- unserialize(serialize(sampler, NULL))
expect_identical(reloaded$data@n.cuts, c(5L, 5L, 5L, 100L))

# a data object's count is for the columns it was built on; the appended
# period takes 100
hazardData <- dbartsData(form, frame, n.cuts = c(5L, 7L, 9L))
expect_identical(
  dbarts(
    hazardData,
    family = "hazard",
    control = run(),
    n.trees = 5L
  )$data@n.cuts,
  c(5L, 7L, 9L, 100L)
)

# the oracle: a hazard fit at a stated count draws as the binary fit on the
# expanded design at the counts that gives, and at the period named as well to
# one at the single count
drawsOf <- function(fit, counts, seed = 5L) {
  d <- fit$data
  link <- fit$model@family
  twin <- dbarts(
    d@x,
    d@y,
    family = link,
    n.cuts = counts,
    control = run(),
    seed = seed,
    n.trees = 5L
  )
  list(
    hazard = fit$run(0L, 4L)$train,
    binary = twin$run(0L, 4L)$train
  )
}
fiveAll <- dbarts(
  x,
  Y,
  family = "hazard",
  n.cuts = 5L,
  control = run(),
  seed = 5L,
  n.trees = 5L
)
oracle <- drawsOf(fiveAll, c(5L, 5L, 5L, 100L))
expect_identical(oracle$hazard, oracle$binary)
named <- dbarts(
  x,
  Y,
  family = "hazard",
  n.cuts = c(period = 5L, .default = 5L),
  control = run(),
  seed = 5L,
  n.trees = 5L
)
oracle <- drawsOf(named, 5L)
expect_identical(oracle$hazard, oracle$binary)
# and the two differ, the period's grid being coarser in the second
expect_false(identical(
  drawsOf(fiveAll, 5L)$binary,
  drawsOf(fiveAll, c(5L, 5L, 5L, 100L))$binary
))

# split.probs keeps its refusal, which names the period
expect_error(
  dbarts(
    x,
    Y,
    family = "hazard",
    split.probs = c(x1 = 0.5, x2 = 0.5, x3 = 0.5),
    control = run(),
    n.trees = 5L
  ),
  "period"
)

# a count stated to dbartsData() on a hazard sampler's person-period matrix is
# read as creation reads it, at $setData and at the direct door: the period is
# named or else 100
hazardFit <- dbarts(
  x,
  Y,
  family = "hazard",
  n.cuts = 5L,
  control = run(),
  n.trees = 5L
)
personPeriod <- hazardFit$data@x
atDoors <- function(value) {
  given <- dbartsData(personPeriod, hazardFit$data@y, n.cuts = value)
  swapped <- hazardFit$copy()
  swapped$setData(given)
  by.hand <- new("dbartsSampler", run(), hazardFit$model, given)
  list(setData = swapped$data@n.cuts, hand = by.hand$data@n.cuts)
}
for (case in list(
  list(7L, c(7L, 7L, 7L, 100L)),
  list(c(x1 = 3L, .default = 7L), c(3L, 7L, 7L, 100L)),
  list(c(period = 3L, .default = 7L), c(7L, 7L, 7L, 3L)),
  list(c(7L, 8L, 9L, 10L), c(7L, 8L, 9L, 10L)),
  list(c(7L, 8L, 9L), c(7L, 8L, 9L, 100L))
)) {
  got <- atDoors(case[[1L]])
  expect_identical(got$setData, case[[2L]], info = deparse(case[[1L]]))
  expect_identical(got$hand, case[[2L]], info = deparse(case[[1L]]))
}
# the count taken from the sampler's own data object keeps its values
back <- hazardFit$copy()
back$setData(hazardFit$data)
expect_identical(back$data@n.cuts, c(5L, 5L, 5L, 100L))

# the oracle: draws as a sampler created with the same count
seven <- dbarts(
  x,
  Y,
  family = "hazard",
  n.cuts = 7L,
  control = run(),
  seed = 5L,
  n.trees = 5L
)
expected <- seven$run(0L, 4L)$train
fiveFit <- dbarts(
  x,
  Y,
  family = "hazard",
  n.cuts = 5L,
  control = run(),
  seed = 5L,
  n.trees = 5L
)
fiveFit$setData(dbartsData(personPeriod, hazardFit$data@y, n.cuts = 7L))
expect_identical(fiveFit$run(0L, 4L)$train, expected)
by.hand <- new(
  "dbartsSampler",
  run(seed = 5L),
  seven$model,
  dbartsData(personPeriod, hazardFit$data@y, n.cuts = 7L)
)
expect_identical(by.hand$run(0L, 4L)$train, expected)

# without hazard the same data object is read over all its columns, and a value
# per predictor only is refused
plainFit <- dbarts(
  personPeriod,
  hazardFit$data@y,
  control = run(),
  n.trees = 5L
)
plainFit$setData(dbartsData(personPeriod, hazardFit$data@y, n.cuts = 7L))
expect_identical(plainFit$data@n.cuts, rep(7L, 4L))
expect_error(
  plainFit$setData(
    dbartsData(personPeriod, hazardFit$data@y, n.cuts = c(7L, 8L, 9L))
  ),
  "predictor column"
)

# a count stated on the control beside the same one given flat is no mismatch
notUsed <- function(e) {
  w <- character()
  withCallingHandlers(
    e,
    dbartsControlNotUsedWarning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    }
  )
  w
}
expect_identical(
  notUsed(dbarts(
    x,
    Y,
    family = "hazard",
    n.cuts = 7L,
    control = written(n.cuts = 7L)
  )),
  character()
)

# the direct door beside a control that states a count
expect_identical(
  notUsed(new(
    "dbartsSampler",
    written(n.cuts = 5L),
    hazardFit$model,
    hazardFit$data
  )),
  character()
)
expect_identical(
  notUsed(new(
    "dbartsSampler",
    written(n.cuts = c(5L, 5L, 5L)),
    hazardFit$model,
    hazardFit$data
  )),
  character()
)
expect_true(grepl(
  "(9, 100)",
  notUsed(new(
    "dbartsSampler",
    written(n.cuts = 9L),
    hazardFit$model,
    hazardFit$data
  )),
  fixed = TRUE
))

# a single value on a hazard fit with no predictors leaves the period alone
expect_identical(
  dbarts(
    matrix(numeric(), n, 0L),
    Y,
    family = "hazard",
    n.cuts = 5L,
    control = run(),
    n.trees = 5L
  )$data@n.cuts,
  100L
)

# the controls above write counts and the rule on purpose; give the keys of the
# warnings that read them back to the next file
onceState <- dbarts:::onceWarnState
for (key in paste0(
  "tombstone.control.",
  c("n.trees", "n.cuts", "useQuantiles")
)) {
  onceState[[key]] <- NULL
}
