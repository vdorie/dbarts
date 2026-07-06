source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that rng cooperates with native generator
n.trees <- 5L

oldSeed <- if (exists(".Random.seed")) {
  .Random.seed
} else {
  runif(1L)
  .Random.seed
}
control <- dbarts::dbartsControl(
  n.trees = n.trees,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  updateState = FALSE
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
invisible(sampler$run(0L, 5L))
expect_true(any(.Random.seed != oldSeed))

rm(sampler, control, oldSeed, n.trees)


# test that rng with fixed seed matches native, one chain, one thread
set.seed(1234L, kind = "Mersenne-Twister", normal.kind = "Inversion")
control <- dbarts::dbartsControl(
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  updateState = FALSE
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
builtInResults <- sampler$run(0L, 5L)

control <- dbarts::dbartsControl(
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  updateState = FALSE,
  rngSeed = 1234L
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
seedOnlyResults <- sampler$run(0L, 5L)

expect_equal(builtInResults$train, seedOnlyResults$train)

# the dbarts() seed argument mirrors dbartsControl(rngSeed = ) and, when
# given, overrides a seed already set in the control
control <- dbarts::dbartsControl(
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  updateState = FALSE
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control, seed = 1234L)
seedArgResults <- sampler$run(0L, 5L)
expect_equal(seedArgResults$train, seedOnlyResults$train)

control@rngSeed <- 999L
sampler <- dbarts::dbarts(y ~ x, testData, control = control, seed = 1234L)
overrideResults <- sampler$run(0L, 5L)
expect_equal(overrideResults$train, seedOnlyResults$train)

rm(
  sampler,
  control,
  seedOnlyResults,
  builtInResults,
  seedArgResults,
  overrideResults
)


# test that bart with fixed seed is reproducible
fit1 <- dbarts::bart(
  y ~ x,
  testData,
  ntree = 5,
  nskip = 0,
  ndpost = 3,
  seed = 12345L,
  verbose = FALSE,
  nchain = 2L,
  nthread = 1L
)
fit2 <- dbarts::bart(
  y ~ x,
  testData,
  ntree = 5,
  nskip = 0,
  ndpost = 3,
  seed = 12345L,
  verbose = FALSE,
  nchain = 2L,
  nthread = 1L
)

expect_equal(fit1$yhat.train, fit2$yhat.train)

fit3 <- dbarts::bart(
  y ~ x,
  testData,
  ntree = 5,
  nskip = 0,
  ndpost = 3,
  seed = 12345L,
  verbose = FALSE,
  nchain = 2L,
  nthread = 2L
)
fit4 <- dbarts::bart(
  y ~ x,
  testData,
  ntree = 5,
  nskip = 0,
  ndpost = 3,
  seed = 12345L,
  verbose = FALSE,
  nchain = 2L,
  nthread = 2L
)

expect_equal(fit3$yhat.train, fit4$yhat.train)

# chain generators are seeded identically regardless of the thread count,
# so seeded results are thread-count independent
expect_equal(fit1$yhat.train, fit3$yhat.train)

rm(fit4, fit3, fit2, fit1)

rm(testData)
