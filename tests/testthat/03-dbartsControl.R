source(system.file("common", "friedmanData.R", package = "dbarts"))

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
  expect_error(dbartsControl(n.samples = 0L))

  expect_error(dbartsControl(n.burn = "not-an-integer"))
  expect_error(dbartsControl(n.burn = NA_integer_))
  expect_error(dbartsControl(n.burn = -1L))

  expect_error(dbartsControl(n.trees = "not-an-integer"))
  expect_error(dbartsControl(n.trees = NA_integer_))
  expect_error(dbartsControl(n.trees = 0L))

  expect_error(dbartsControl(n.threads = "not-an-integer"))
  expect_error(dbartsControl(n.threads = NA_integer_))
  expect_error(dbartsControl(n.threads = 0L))

  expect_error(dbartsControl(n.thin = "not-an-integer"))
  expect_error(dbartsControl(n.thin = NA_integer_))
  expect_error(dbartsControl(n.thin = 0L))

  expect_error(dbartsControl(printEvery = "not-an-integer"))
  expect_error(dbartsControl(printEvery = NA_integer_))
  expect_error(dbartsControl(printEvery = -1L))

  expect_error(dbartsControl(printCutoffs = "not-an-integer"))
  expect_error(dbartsControl(printCutoffs = NA_integer_))
  expect_error(dbartsControl(printCutoffs = -1L))

  expect_error(dbartsControl(n.cuts = "not-an-integer"))
  expect_error(dbartsControl(n.cuts = NA_integer_))
  expect_error(dbartsControl(n.cuts = -1L))
})
