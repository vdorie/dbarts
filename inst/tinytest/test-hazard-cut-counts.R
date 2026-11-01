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
# the tree count written on a control is read with a warning, which is not
# what this file tests
run <- function(...) {
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
written <- run

doors <- list(
  dbarts = function(v) {
    dbarts(x, Y, family = "hazard", n.cuts = v, control = run())
  },
  control = function(v) {
    dbarts(x, Y, family = "hazard", control = written(n.cuts = v))
  },
  formula = function(v) {
    dbarts(form, frame, family = "hazard", n.cuts = v, control = run())
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
sampler <- dbarts(x, Y, family = "hazard", n.cuts = 5L, control = run())
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
  dbarts(hazardData, family = "hazard", control = run())$data@n.cuts,
  c(5L, 7L, 9L, 100L)
)
expanded <- dbarts(form, frame, family = "hazard", control = run())$data
expanded@n.cuts <- NA_integer_
carried <- expanded
carried@n.cuts <- c(5L, 7L, 9L)
expect_identical(
  dbarts:::applyCutSettings(carried, list(stated = character()), TRUE)@n.cuts,
  c(5L, 7L, 9L, 100L)
)
expect_identical(
  dbarts:::applyCutSettings(carried, list(stated = character()), FALSE)@n.cuts,
  c(5L, 7L, 9L, 5L)
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
    seed = seed
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
  seed = 5L
)
oracle <- drawsOf(fiveAll, c(5L, 5L, 5L, 100L))
expect_identical(oracle$hazard, oracle$binary)
named <- dbarts(
  x,
  Y,
  family = "hazard",
  n.cuts = c(period = 5L, .default = 5L),
  control = run(),
  seed = 5L
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
    control = run()
  ),
  "period"
)
