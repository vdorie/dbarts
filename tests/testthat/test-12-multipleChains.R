context("bart with multiple chains")

source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

test_that("multiple chains single threads runs correctly", {
  set.seed(99)
  oldSeed <- .Random.seed
  fit <- bart(testData$x, testData$y, ndpost = 120L, nskip = 100L, ntree = 50L, nthread = 1L, nchain = 2L, combinechains = FALSE, verbose = FALSE)
  expect_true(any(oldSeed != .Random.seed))
})

test_that("multiple chains multiple threads runs correctly", {
  set.seed(99)
  oldSeed <- .Random.seed
  fit <- bart(testData$x, testData$y, ndpost = 120L, nskip = 100L, ntree = 50L, nthread = 2L, nchain = 2L, combinechains = FALSE, verbose = FALSE)
  
  expect_identical(dim(fit$yhat.train), c(2L, 120L, nrow(testData$x)))
  ## check that the random seeds didn't end up the same/the chains got different results
  expect_true(mean(abs(fit$yhat.train[1L,,] - fit$yhat.train[2L,,])) > 1.0e-6)
  expect_equal(oldSeed, .Random.seed)
})

