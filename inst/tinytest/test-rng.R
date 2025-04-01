source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that rng cooperates with native generator
n.trees  <- 5L

oldSeed <- if (exists(".Random.seed")) .Random.seed else { runif(1L); .Random.seed }
control <- dbarts::dbartsControl(
  n.trees = n.trees, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE,
  rngKind = "default", rngNormalKind = "default"
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
invisible(sampler$run(0L, 5L))
expect_true(any(.Random.seed != oldSeed))

oldSeed <- .Random.seed
control <- dbarts::dbartsControl(
  n.trees = n.trees, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE,
  rngKind = "Mersenne-Twister", rngNormalKind = "Inversion"
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
invisible(sampler$run(0L, 5L))
expect_equal(.Random.seed, oldSeed)

## run once for 10 iterations
control <- dbarts::dbartsControl(
  n.trees = n.trees, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE,
  rngKind = "Mersenne-Twister", rngNormalKind = "Inversion"
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
sampler$storeState()
origSeed <- .Call(dbarts:::C_dbarts_deepCopy, sampler$state[[1L]]@rng.state)
invisible(sampler$run(0L, 10L, updateState = TRUE))
oldSeed <- sampler$state[[1L]]@rng.state

## run twice for 5 iterations, writing out state and restoring from that
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
sampler$storeState()
tempState <- sampler$state
tempState[[1L]]@rng.state <- .Call(dbarts:::C_dbarts_deepCopy, origSeed)
sampler$setState(tempState)

invisible(sampler$run(0L, 5L, updateState = TRUE))
tempState <- sampler$state

sampler <- dbarts::dbarts(y ~ x, testData, control = control)
sampler$setState(tempState)
invisible(sampler$run(0L, 5L, updateState = TRUE))

expect_equal(sampler$state[[1L]]@rng.state, oldSeed)

rm(sampler, tempState, oldSeed, origSeed, control, n.trees)


# test that rng with fixed seed matches native, one chain, one thread
set.seed(1234L, kind = "Mersenne-Twister", normal.kind = "Inversion")
control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
builtInResults <- sampler$run(0L, 5L)

control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE,
  rngSeed = 1234L
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
seedOnlyResults <- sampler$run(0L, 5L)

control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE,
  rngKind = "Mersenne-Twister", rngNormalKind = "Inversion",
  rngSeed = 1234L
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
allRngResults <- sampler$run(0L, 5L)

expect_equal(builtInResults$train, seedOnlyResults$train)
expect_equal(builtInResults$train, allRngResults$train)

rm(allRngResults, sampler, control, seedOnlyResults, builtInResults)


# test that rng with fixed seed matches native, two chains, one threads
set.seed(1234L, kind = "Mersenne-Twister", normal.kind = "Inversion")
control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 2L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
builtInResults <- sampler$run(0L, 5L)

control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE,
  rngSeed = 1234L
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
sequentialResults_1 <- sampler$run(0L, 5L)

control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
sequentialResults_2 <- sampler$run(0L, 5L)

expect_equal(
  as.vector(builtInResults$train),
  c(sequentialResults_1$train, sequentialResults_2$train)
)

rm(sequentialResults_2, sampler, control, sequentialResults_1, builtInResults)


# test that rng with fixed seed matches native, two chains, one threads
control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 2L, n.threads = 2L,
  verbose = FALSE, updateState = FALSE,
  rngSeed = 1234L
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
seedOnlyResults <- sampler$run(0L, 5L)


set.seed(1234L, kind = "Mersenne-Twister", normal.kind = "Inversion")
seed_1 <- runif(1L)
seed_2 <- runif(1L)

to_unsigned_int <- function(unif) {
  unif <- unif * (.Machine$integer.max + 0.5)
  msb <- floor(unif)
  one <- unif - msb
  bitwShiftL(msb, 1L) + 1 * (unif - msb >= 0.5)
}

control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE,
  rngKind = "Mersenne-Twister", rngNormalKind = "Inversion",
  rngSeed = to_unsigned_int(seed_1)
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
sequentialResults_1 <- sampler$run(0L, 5L)
control <- dbarts::dbartsControl(
  n.trees = 5L, n.chains = 1L, n.threads = 1L,
  verbose = FALSE, updateState = FALSE,
  rngKind = "Mersenne-Twister", rngNormalKind = "Inversion",
  rngSeed = to_unsigned_int(seed_2)
)
sampler <- dbarts::dbarts(y ~ x, testData, control = control)
sequentialResults_2 <- sampler$run(0L, 5L)

expect_equal(
  as.vector(seedOnlyResults$train),
  c(sequentialResults_1$train, sequentialResults_2$train)
)

rm(sequentialResults_2, sampler, control, sequentialResults_1)
rm(to_unsigned_int, seed_2, seed_1, seedOnlyResults)


# test that bart with fixed seed is reproducible
fit1 <- dbarts::bart(
  y ~ x, testData, ntree = 5, nskip = 0, ndpost = 3, seed = 12345L,
  verbose = FALSE,
  nchain = 2L, nthread = 1L
)
fit2 <- dbarts::bart(
  y ~ x, testData, ntree = 5, nskip = 0, ndpost = 3, seed = 12345L,
  verbose = FALSE,
  nchain = 2L, nthread = 1L
)

expect_equal(fit1$yhat.train, fit2$yhat.train)

fit3 <- dbarts::bart(
  y ~ x, testData, ntree = 5, nskip = 0, ndpost = 3, seed = 12345L,
  verbose = FALSE,
  nchain = 2L, nthread = 2L
)
fit4 <- dbarts::bart(
  y ~ x, testData, ntree = 5, nskip = 0, ndpost = 3, seed = 12345L,
  verbose = FALSE,
  nchain = 2L, nthread = 2L
)

expect_equal(fit3$yhat.train, fit4$yhat.train)

expect_true(any(fit1$yhat.train != fit3$yhat.train))

rm(fit4, fit3, fit2, fit1)

rm(testData)

