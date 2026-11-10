# The collapsed rescaling move under logistic, from R: that it sits ahead of
# every recorded channel, runs per chain on any number of threads, runs finite
# under an offset, counts, a mask and a mask with every row out, reaches the
# forests that share one prior total and a binomial response, and that a
# sampler's state carries what it needs. Each expectation holds whether or not
# a fit takes the move, so none pins which do. The move's conditional, its mapping, its schedule, its
# declines and its place in the sweep are tests/cpp gates; its posterior is
# the probit-k-scale-exact gate; which fits take it is the equivalence
# baseline's record.

set.seed(101)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("x1", "x2", "x3")))
eta <- 1.5 * sin(3 * x[, 1L]) + x[, 2L] - 0.7
yBinary <- as.numeric(eta + rnorm(n) > 0)
xTest <- matrix(runif(30L), 10L, 3L, dimnames = list(NULL, colnames(x)))
counts <- rep(c(1, 2, 3), length.out = n)
mask <- rep(c(1, 1, 1, 1, 0), length.out = n)
d <- data.frame(yBinary, x)

moveControl <- function(n.chains = 1L, n.threads = 1L) {
  dbarts::dbartsControl(
    n.chains = n.chains,
    n.threads = n.threads,
    n.burn = 0L,
    n.samples = 20L,
    seed = 11L
  )
}
# the default leaf prior draws k under chi(1.5, 2), the fits the move reaches;
# shared puts the trees in two forests that share one prior total
moveSampler <- function(shared = FALSE, control = moveControl(), ...) {
  if (shared) {
    dbarts::dbarts(
      yBinary ~ forest(x1) + forest(x2 + x3),
      d,
      family = "logistic",
      n.trees = 10L,
      control = control,
      ...
    )
  } else {
    dbarts::dbarts(
      x,
      yBinary,
      family = "logistic",
      n.trees = 10L,
      control = control,
      ...
    )
  }
}
# a run that recorded k, with it and the fits finite
finiteWithK <- function(run) {
  !is.null(run$k) && all(is.finite(run$k)) && all(is.finite(run$train))
}

for (shared in c(FALSE, TRUE)) {
  # ---- per chain, on any number of threads ----

  twoChains <- function(n.threads) {
    moveSampler(
      shared,
      moveControl(n.chains = 2L, n.threads = n.threads)
    )$run(10L, 40L)
  }
  oneThread <- twoChains(1L)
  expect_true(finiteWithK(oneThread), info = shared)
  expect_identical(oneThread, twoChains(2L), info = shared)
  rm(oneThread)

  # ---- an offset, a mask installed mid-run and a mask with every row out ----

  masked <- moveSampler(shared, offset = 0.3 * x[, 3L])
  invisible(masked$run(20L, 1L))
  masked$setActiveRows(mask)
  maskedRun <- masked$run(20L, 40L)
  expect_true(finiteWithK(maskedRun), info = shared)
  masked$setActiveRows(rep(0, n))
  maskedRun <- masked$run(0L, 40L)
  expect_true(finiteWithK(maskedRun), info = shared)
  rm(masked, maskedRun)

  # ---- the state round-trip ----

  # a sampler restored from the store, a reload of the same store and a copy
  # continue identically, k and the fits both
  stored <- moveSampler(shared)
  invisible(stored$run(10L, 1L))
  stored$storeState()
  serialized <- tempfile(fileext = ".rds")
  saveRDS(stored, serialized)
  restored <- moveSampler(shared)
  restored$setState(stored$state)
  reloaded <- readRDS(serialized)
  unlink(serialized)
  copied <- stored$copy()
  continued <- restored$run(0L, 40L)
  expect_identical(continued, reloaded$run(0L, 40L), info = shared)
  expect_identical(continued, copied$run(0L, 40L), info = shared)
  rm(stored, restored, reloaded, copied, continued, serialized)
}

# ---- ahead of every recorded channel ----

# the saved trees replayed reproduce each kept draw's training and test fits,
# which holds only if the leaves are rescaled before any of them is written
treesFit <- dbarts::bart(
  x,
  yBinary,
  test = xTest,
  family = "logistic",
  n.trees = 10L,
  n.samples = 40L,
  n.burn = 30L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  seed = 11L,
  verbose = FALSE
)
expect_true(
  max(abs(predict(treesFit, x, type = "bart") - treesFit$yhat.train)) < 1e-10
)
expect_true(
  max(abs(predict(treesFit, xTest, type = "bart") - treesFit$yhat.test)) < 1e-10
)
rm(treesFit)

# ---- count weights, an infinite prior scale, binomial, hazard.logistic ----

counted <- moveSampler(weights = counts, offset = 0.3 * x[, 3L])
countedRun <- counted$run(20L, 40L)
expect_true(finiteWithK(countedRun))
rm(counted, countedRun)

flat <- moveSampler(leaf.prior = normal(k = chi(1.5, Inf)))
flatRun <- flat$run(20L, 40L)
expect_true(all(is.finite(flatRun$train)))
rm(flat, flatRun)

# a binomial response, successes out of 1 to 3 trials under count weights, on
# one forest and on the sharers
trials <- rep(1:3, length.out = n)
successes <- rbinom(n, trials, plogis(eta))
dBinomial <- data.frame(successes, failures = trials - successes, x)
for (shared in c(FALSE, TRUE)) {
  binomialSampler <- dbarts::dbarts(
    if (shared) {
      cbind(successes, failures) ~ forest(x1) + forest(x2 + x3)
    } else {
      cbind(successes, failures) ~ x1 + x2 + x3
    },
    dBinomial,
    family = binomial,
    weights = rep(1:2, length.out = n),
    n.trees = 10L,
    control = moveControl()
  )
  expect_true(finiteWithK(binomialSampler$run(20L, 40L)), info = shared)
}
rm(binomialSampler, trials, successes, dBinomial)

# three periods, an event in each and subjects censored along the way
periods <- rep(1:3, length.out = n)
events <- as.numeric(runif(n) < 0.6)
hazardFit <- dbarts::bart(
  x,
  cbind(periods, events),
  family = "hazard.logistic",
  n.trees = 10L,
  n.samples = 40L,
  n.burn = 20L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 11L,
  verbose = FALSE
)
expect_true(
  !is.null(hazardFit$k) &&
    all(is.finite(hazardFit$k)) &&
    all(is.finite(hazardFit$yhat.train))
)
rm(hazardFit)
rm(n, x, eta, yBinary, xTest, counts, mask, d, periods, events)
rm(moveControl, moveSampler, finiteWithK, twoChains, shared)
