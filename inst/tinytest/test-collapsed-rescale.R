# The collapsed rescaling move, from R, under both binary links: that it sits
# ahead of every recorded channel, runs per chain on any number of threads,
# runs finite under an offset, counts, a mask and a mask with every row out,
# and that a sampler's state carries what it needs. Each expectation holds
# whether or not a link's fits take the move, so none pins which do. The
# move's conditional, its mapping, its schedule and its declines are tests/cpp
# gates; its posterior is the probit-k-scale-exact gate; which fits take it
# is the equivalence baseline's record.

set.seed(101)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("x1", "x2", "x3")))
eta <- 1.5 * sin(3 * x[, 1L]) + x[, 2L] - 0.7
yBinary <- as.numeric(eta + rnorm(n) > 0)
xTest <- matrix(runif(30L), 10L, 3L, dimnames = list(NULL, colnames(x)))
counts <- rep(c(1, 2, 3), length.out = n)
mask <- rep(c(1, 1, 1, 1, 0), length.out = n)

moveControl <- function(n.chains = 1L, n.threads = 1L) {
  dbarts::dbartsControl(
    n.trees = 10L,
    n.chains = n.chains,
    n.threads = n.threads,
    n.burn = 0L,
    n.samples = 20L,
    seed = 11L
  )
}
# the default leaf prior draws k under chi(1.5, 2), the fits the move reaches
moveSampler <- function(family, control = moveControl(), ...) {
  dbarts::dbarts(x, yBinary, family = family, control = control, ...)
}
# a run that recorded k, with it and the fits finite
finiteWithK <- function(run) {
  !is.null(run$k) && all(is.finite(run$k)) && all(is.finite(run$train))
}

for (family in c("probit", "logistic")) {
  # ---- per chain, on any number of threads ----

  twoChains <- function(n.threads) {
    moveSampler(
      family,
      moveControl(n.chains = 2L, n.threads = n.threads)
    )$run(10L, 40L)
  }
  oneThread <- twoChains(1L)
  expect_true(finiteWithK(oneThread), info = family)
  expect_identical(oneThread, twoChains(2L), info = family)
  rm(oneThread)

  # ---- ahead of every recorded channel ----

  # the saved trees replayed reproduce each kept draw's training and test
  # fits, which holds only if the leaves are rescaled before any of them is
  # written
  treesFit <- dbarts::bart(
    x,
    yBinary,
    test = xTest,
    family = family,
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
    max(abs(predict(treesFit, x, type = "bart") - treesFit$yhat.train)) < 1e-10,
    info = family
  )
  expect_true(
    max(abs(predict(treesFit, xTest, type = "bart") - treesFit$yhat.test)) <
      1e-10,
    info = family
  )
  rm(treesFit)

  # ---- an offset, a mask installed mid-run and a mask with every row out ----

  masked <- moveSampler(family, offset = 0.3 * x[, 3L])
  invisible(masked$run(20L, 1L))
  masked$setActiveRows(mask)
  maskedRun <- masked$run(20L, 40L)
  expect_true(finiteWithK(maskedRun), info = family)
  masked$setActiveRows(rep(0, n))
  maskedRun <- masked$run(0L, 40L)
  expect_true(finiteWithK(maskedRun), info = family)
  rm(masked, maskedRun)

  # ---- the state round-trip ----

  # a sampler restored from the store, a reload of the same store and a copy
  # continue identically, k and the fits both
  stored <- moveSampler(family)
  invisible(stored$run(10L, 1L))
  stored$storeState()
  serialized <- tempfile(fileext = ".rds")
  saveRDS(stored, serialized)
  restored <- moveSampler(family)
  restored$setState(stored$state)
  reloaded <- readRDS(serialized)
  unlink(serialized)
  copied <- stored$copy()
  continued <- restored$run(0L, 40L)
  expect_identical(continued, reloaded$run(0L, 40L), info = family)
  expect_identical(continued, copied$run(0L, 40L), info = family)
  rm(stored, restored, reloaded, copied, continued, serialized)
}

# ---- logistic: count weights, an infinite prior scale, hazard.logistic ----

counted <- moveSampler("logistic", weights = counts, offset = 0.3 * x[, 3L])
countedRun <- counted$run(20L, 40L)
expect_true(finiteWithK(countedRun))
rm(counted, countedRun)

flat <- moveSampler("logistic", leaf.prior = normal(k = chi(1.5, Inf)))
flatRun <- flat$run(20L, 40L)
expect_true(all(is.finite(flatRun$train)))
rm(flat, flatRun)

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

# ---- probit: 0/1 weights are the mask ----

weighted <- moveSampler("probit", weights = mask)$run(20L, 40L)
byMask <- moveSampler("probit")
byMask$setActiveRows(mask)
expect_identical(weighted, byMask$run(20L, 40L))
rm(weighted, byMask)
