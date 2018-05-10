context("dbartsControl")

test_that("logical arguments are identified and 'bounded'", {
  expect_error(dbartsControl(verbose = NA))
  expect_error(dbartsControl(verbose = "not-a-logical"))

  expect_error(dbartsControl(keepTrainingFits = NA))
  expect_error(dbartsControl(keepTrainingFits = "not-a-logical"))

  expect_error(dbartsControl(useQuantiles = NA))
  expect_error(dbartsControl(useQuantiles = "not-a-logical"))

  expect_error(dbartsControl(updateState = NA))
  expect_error(dbartsControl(updateState = "not-a-logical"))
})

test_that("integer arguments are identified and bounded", {
  expect_error(dbartsControl(n.samples = "not-an-integer"))
  expect_error(dbartsControl(n.samples = -1L))

  expect_error(dbartsControl(n.burn = "not-an-integer"))
  expect_error(dbartsControl(n.burn = NA_integer_))
  expect_error(dbartsControl(n.burn = -1L))

  expect_error(dbartsControl(n.trees = "not-an-integer"))
  expect_error(dbartsControl(n.trees = NA_integer_))
  expect_error(dbartsControl(n.trees = 0L))
  
  expect_error(dbartsControl(n.chains = "not-an-integer"))
  expect_error(dbartsControl(n.chains = NA_integer_))
  expect_error(dbartsControl(n.chains = 0L))

  expect_error(dbartsControl(n.threads = "not-an-integer"))
  expect_error(dbartsControl(n.threads = NA_integer_))
  expect_error(dbartsControl(n.threads = 0L))

  expect_error(dbartsControl(n.thin = "not-an-integer"))
  expect_error(dbartsControl(n.thin = NA_integer_))
  expect_error(dbartsControl(n.thin = -1L))

  expect_error(dbartsControl(printEvery = "not-an-integer"))
  expect_error(dbartsControl(printEvery = NA_integer_))
  expect_error(dbartsControl(printEvery = -1L))

  expect_error(dbartsControl(printCutoffs = "not-an-integer"))
  expect_error(dbartsControl(printCutoffs = NA_integer_))
  expect_error(dbartsControl(printCutoffs = -1L))

  expect_error(dbartsControl(n.cuts = "not-an-integer"))
  expect_error(dbartsControl(n.cuts = NA_integer_))
  expect_error(dbartsControl(n.cuts = -1L))
  
  expect_error(dbartsControl(rngKind = "not-an-rng"))
  expect_error(dbartsControl(rngNormalKind = "not-an-rng"))
})

source(system.file("common", "friedmanData.R", package = "dbarts"))

test_that("control argument works", {
  n.samples <- 500L
  n.trees <- 50L
  n.cuts <- 50L
  n.chains <- 2L
  n.threads <- 2L
  control <- dbartsControl(n.samples = n.samples, verbose = FALSE, n.trees = n.trees, n.cuts = n.cuts,
                           n.chains = n.chains, n.threads = n.threads)
  sampler <- dbarts(y ~ x, testData, control = control)

  expect_equal(sampler$control@n.trees, n.trees)
  expect_equal(sampler$control@verbose, FALSE)
  expect_equal(sampler$control@n.samples, n.samples)
  expect_true(all(sampler$data@n.cuts == n.cuts))
  expect_equal(sampler$control@n.chains, n.chains)
  expect_equal(sampler$control@n.threads, n.threads)
})

test_that("keepevery behaves as it did in BayesTree", {
  n.burn <- 0L
  n.sims <- 200L
  keepevery <- 10L
  
  bartFit <- bart(testData$x, testData$y, ndpost = n.sims, nskip = n.burn, ntree = 25L, verbose = FALSE, keepevery = keepevery)
  
  expect_equal(nrow(bartFit$yhat.train), n.sims %/% keepevery)
})

test_that("call is propagated", {
  n.samples <- 5L
  n.burn    <- 1L
  n.trees   <- 5L
  
  control <- dbartsControl(n.samples = n.samples, verbose = FALSE, n.trees = n.trees, n.chains = 1L, n.threads = 1L)
  sampler <- dbarts(y ~ x, testData, control = control)
  expect_equal(sampler$control@call[[1L]], quote(dbarts))
  
  control@call <- call("NULL")
  sampler <- dbarts(y ~ x, testData, control = control)
  expect_equal(sampler$control@call, control@call)
  
  bartFit <- bart(y ~ x, testData, verbose = FALSE, ndpost = n.samples, nskip = n.burn, ntree = n.trees)
  expect_equal(bartFit$call[[1L]], quote(bart))
  
  bartFit <- bart(y ~ x, testData, verbose = FALSE, ndpost = n.samples, nskip = n.burn, ntree = n.trees, keepcall = FALSE)
  expect_equal(bartFit$call, call("NULL"))
})

test_that("rng cooperates with native generator", {
  n.trees  <- 5L
  
  oldSeed <- if (exists(".Random.seed")) .Random.seed else { runif(1L); .Random.seed }
  control <- dbartsControl(verbose = FALSE, n.trees = n.trees, n.chains = 1L, n.threads = 1L,
                           updateState = FALSE,
                           rngKind = "default", rngNormalKind = "default")
  sampler <- dbarts(y ~ x, testData, control = control)
  invisible(sampler$run(0L, 5L))
  expect_true(any(.Random.seed != oldSeed))
  
  oldSeed <- .Random.seed
  control <- dbartsControl(verbose = FALSE, n.trees = n.trees, n.chains = 1L, n.threads = 1L,
                           updateState = FALSE,
                           rngKind = "Mersenne-Twister", rngNormalKind = "Inversion")
  sampler <- dbarts(y ~ x, testData, control = control)
  invisible(sampler$run(0L, 5L))
  expect_equal(.Random.seed, oldSeed)
  
  ## run once for 10 iterations
  control <- dbartsControl(verbose = FALSE, n.trees = n.trees, n.chains = 1L, n.threads = 1L,
                           updateState = FALSE,
                           rngKind = "Mersenne-Twister", rngNormalKind = "Inversion", )
  sampler <- dbarts(y ~ x, testData, control = control)
  sampler$storeState()
  origSeed <- .Call(dbarts:::C_dbarts_deepCopy, sampler$state[[1L]]@rng.state)
  invisible(sampler$run(0L, 10L, TRUE))
  oldSeed <- sampler$state[[1L]]@rng.state
  
  ## run twice for 5 iterations, writing out state and restoring from that
  sampler <- dbarts(y ~ x, testData, control = control)
  sampler$storeState()
  tempState <- sampler$state
  tempState[[1L]]@rng.state <- .Call(dbarts:::C_dbarts_deepCopy, origSeed)
  sampler$setState(tempState)
  
  invisible(sampler$run(0L, 5L, TRUE))
  tempState <- sampler$state
  
  sampler <- dbarts(y ~ x, testData, control = control)
  sampler$setState(tempState)
  invisible(sampler$run(0L, 5L, TRUE))
  
  expect_equal(sampler$state[[1L]]@rng.state, oldSeed)
})
