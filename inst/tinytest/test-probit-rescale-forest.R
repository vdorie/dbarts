# The probit rescaling step: its control switch, where it acts and where it is
# inert, that it sits ahead of every recorded channel, and that a sampler's
# state carries what it needs. The step's conditional, its mapping and its
# declines are tests/cpp gates; its posterior is the probit-k-scale-exact gate.

source(
  system.file("common", "stateContinuation.R", package = "dbarts"),
  local = TRUE
)

set.seed(99)
n <- 120L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("x1", "x2", "x3")))
eta <- 1.5 * sin(3 * x[, 1L]) + x[, 2L] - 0.7
yBinary <- as.numeric(eta + rnorm(n) > 0)
xTest <- matrix(runif(30L), 10L, 3L, dimnames = list(NULL, colnames(x)))

rescaleControl <- function(rescale, n.chains = 1L, n.threads = 1L) {
  dbarts::dbartsControl(
    probitRescaleForest = rescale,
    n.trees = 10L,
    n.chains = n.chains,
    n.threads = n.threads,
    n.burn = 0L,
    n.samples = 20L,
    seed = 7L
  )
}

# ---- the switch ----

expect_true(dbarts::dbartsControl()@probitRescaleForest)
expect_false(
  dbarts::dbartsControl(probitRescaleForest = FALSE)@probitRescaleForest
)
for (bad in list(NA, "not-a-logical")) {
  expect_error(
    dbarts::dbartsControl(probitRescaleForest = bad),
    "'probitRescaleForest' must be TRUE/FALSE"
  )
}
expect_error(
  dbarts::dbartsControl(probitRescaleForest = c(TRUE, FALSE)),
  "'probitRescaleForest' must be of length 1"
)

# fixed when the sampler is created, restating it accepted
sampler <- dbarts::dbarts(x, yBinary, control = rescaleControl(TRUE))
changed <- sampler$control
changed@probitRescaleForest <- FALSE
expect_error(
  sampler$setControl(changed),
  "changing 'probitRescaleForest' is not available"
)
expect_null(sampler$setControl(sampler$control))

# a control rebuilt by a front door keeps it
viaBart <- dbarts::bart(
  x,
  yBinary,
  control = dbarts::dbartsControl(probitRescaleForest = FALSE),
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepSampler = TRUE
)
expect_false(viaBart$fit$control@probitRescaleForest)
rm(viaBart, changed)

# ---- where it acts ----

# Two samplers that differ only in the switch, both installing one state,
# then run on: the draws part where the step is taken and agree to the
# bit where it is not. TRUE against FALSE from one state isolates the step
# from every other way two fits could differ.
runPairFrom <- function(onSampler, numSamples = 20L) {
  offControl <- onSampler$control
  offControl@probitRescaleForest <- FALSE
  offSampler <- methods::new(
    "dbartsSampler",
    offControl,
    onSampler$model,
    onSampler$data
  )
  invisible(onSampler$run(10L, 1L))
  onSampler$storeState()
  # an install re-derives the fits from the leaves, so both sides take one
  onSampler$setState(onSampler$state)
  offSampler$setState(onSampler$state)
  list(
    on = onSampler$run(0L, numSamples),
    off = offSampler$run(0L, numSamples)
  )
}
stepTaken <- function(pair) !identical(pair$on$train, pair$off$train)

# a drawn-k probit fit through dbarts, bart and bartBT given a drawn k
expect_true(stepTaken(runPairFrom(
  dbarts::dbarts(x, yBinary, control = rescaleControl(TRUE))
)))
expect_true(stepTaken(runPairFrom(
  dbarts::bart(
    x,
    yBinary,
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 7L,
    verbose = FALSE,
    samplerOnly = TRUE
  )
)))
bartBTFit <- dbarts::bartBT(
  x,
  yBinary,
  k = chi(1.5, 2),
  ntree = 10L,
  ndpost = 5L,
  nskip = 5L,
  verbose = FALSE,
  keepsampler = TRUE
)
expect_true(bartBTFit$fit$control@probitRescaleForest)
expect_true(stepTaken(runPairFrom(bartBTFit$fit)))
rm(bartBTFit)

# the hazard fit and the hurdle's zero part
survTime <- sample.int(5L, n, replace = TRUE)
survStatus <- rbinom(n, 1L, 0.7)
hazardFit <- dbarts::bart(
  x,
  cbind(survTime, survStatus),
  family = "hazard",
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 7L,
  verbose = FALSE,
  samplerOnly = TRUE
)
expect_true(stepTaken(runPairFrom(hazardFit)))
yHurdle <- ifelse(yBinary > 0, exp(rnorm(n)), 0)
hurdleFit <- dbarts::bart(
  x,
  yHurdle,
  family = "hurdle.lognormal",
  n.trees = 10L,
  n.samples = 20L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 7L,
  verbose = FALSE,
  control = dbarts::dbartsControl(probitRescaleForest = TRUE)
)
hurdleOff <- dbarts::bart(
  x,
  yHurdle,
  family = "hurdle.lognormal",
  n.trees = 10L,
  n.samples = 20L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 7L,
  verbose = FALSE,
  control = dbarts::dbartsControl(probitRescaleForest = FALSE)
)
expect_false(identical(
  dbarts::extract(hurdleFit, type = "ev"),
  dbarts::extract(hurdleOff, type = "ev")
))
rm(hazardFit, hurdleFit, hurdleOff, survTime, survStatus, yHurdle)

# ---- where it is inert ----

yGaussian <- eta + rnorm(n)
yCount <- rpois(n, exp(0.5 * eta + 1))
z <- rbinom(n, 1L, 0.5)
inertSamplers <- list(
  gaussian = dbarts::dbarts(x, yGaussian, control = rescaleControl(TRUE)),
  logistic = dbarts::dbarts(
    x,
    yBinary,
    family = "logistic",
    control = rescaleControl(TRUE)
  ),
  nbinom = dbarts::dbarts(
    x,
    yCount,
    family = "nbinom",
    control = rescaleControl(TRUE)
  ),
  fixedK = dbarts::dbarts(
    x,
    yBinary,
    leaf.prior = normal(k = 2),
    control = rescaleControl(TRUE)
  ),
  linearLeaf = dbarts::dbarts(
    x,
    yBinary,
    leaf.prior = linear("x2", k = chi(1.5, 2)),
    control = rescaleControl(TRUE)
  ),
  multinomial = dbarts::dbarts(
    x,
    factor(cut(eta + rnorm(n), 3L, labels = c("a", "b", "c"))),
    family = "multinomial",
    control = rescaleControl(TRUE)
  ),
  bcf = dbarts::dbarts(
    x,
    yBinary,
    family = "probit",
    forests = list(
      dbarts::dbartsForests$forest(),
      dbarts::dbartsForests$forest(basis = ~ factor(z), n.trees = 5L)
    ),
    control = rescaleControl(TRUE)
  )
)
for (name in names(inertSamplers)) {
  expect_false(stepTaken(runPairFrom(inertSamplers[[name]])), info = name)
}
rm(inertSamplers, name, yGaussian, yCount, z)

# ---- per chain, and on any number of threads ----

twoChains <- function(rescale, n.threads) {
  dbarts::dbarts(
    x,
    yBinary,
    control = rescaleControl(rescale, n.chains = 2L, n.threads = n.threads)
  )$run(10L, 20L)$train
}
onOne <- twoChains(TRUE, 1L)
offOne <- twoChains(FALSE, 1L)
for (chain in 1:2) {
  expect_false(
    identical(onOne[,, chain], offOne[,, chain]),
    info = paste("chain", chain)
  )
}
expect_identical(onOne, twoChains(TRUE, 2L))
rm(onOne, offOne, chain)

# ---- ahead of every recorded channel ----

# the saved trees replayed reproduce each kept draw's training and test fits,
# which holds only if the step rescales the leaves before any of them is
# written
treesFit <- dbarts::bart(
  x,
  yBinary,
  test = xTest,
  n.trees = 10L,
  n.samples = 30L,
  n.burn = 30L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  seed = 7L,
  verbose = FALSE
)
expect_true(
  max(abs(predict(treesFit, x, type = "bart") - treesFit$yhat.train)) < 1e-10
)
expect_true(
  max(abs(predict(treesFit, xTest, type = "bart") - treesFit$yhat.test)) < 1e-10
)
rm(treesFit)

# ---- an offset, a mask installed mid-run and an infinite prior scale ----

masked <- dbarts::dbarts(
  x,
  yBinary,
  offset = 0.3 * x[, 3L],
  leaf.prior = normal(k = chi(1.5, Inf)),
  control = rescaleControl(TRUE)
)
invisible(masked$run(20L, 1L))
masked$setActiveRows(rep(c(1, 1, 1, 1, 0), length.out = n))
maskedRun <- masked$run(20L, 20L)
expect_true(all(is.finite(maskedRun$train)) && all(is.finite(maskedRun$k)))
rm(masked, maskedRun)

# ---- the setting and the state round-trip ----

stateSampler <- dbarts::dbarts(x, yBinary, control = rescaleControl(FALSE))
invisible(stateSampler$run(10L, 1L))
stateSampler$storeState()
copied <- stateSampler$copy()
expect_false(copied$control@probitRescaleForest)
serialized <- tempfile(fileext = ".rds")
saveRDS(stateSampler, serialized)
reloaded <- readRDS(serialized)
unlink(serialized)
expect_false(reloaded$control@probitRescaleForest)
reloaded$storeState()
statesAgree(reloaded$state, stateSampler$state)

# a sampler restored from the store and a reload of the same store continue
# identically
onSampler <- dbarts::dbarts(x, yBinary, control = rescaleControl(TRUE))
invisible(onSampler$run(10L, 1L))
onSampler$storeState()
serialized <- tempfile(fileext = ".rds")
saveRDS(onSampler, serialized)
restored <- dbarts::dbarts(x, yBinary, control = rescaleControl(TRUE))
restored$setState(onSampler$state)
reloadedOn <- readRDS(serialized)
unlink(serialized)
expect_true(reloadedOn$control@probitRescaleForest)
expect_identical(restored$run(0L, 10L)$train, reloadedOn$run(0L, 10L)$train)
rm(stateSampler, copied, reloaded, onSampler, restored, reloadedOn, serialized)
