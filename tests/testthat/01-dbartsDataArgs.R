source(system.file("common", "friedmanData.R", package = "dbarts"))

context("dbarts data arguments")

test_that("formula specification raises errors", {
  expect_error(dbarts("not-a-formula", testData))
  expect_error(dbarts(y ~ x))
  expect_error(dbarts(NULL, testData))
  expect_error(dbarts(y ~ 0, testData))
})

test_that("extra arguments for formula specification raises errors", {
  testData <- as.data.frame(testData)
  testData$weights <- runif(nrow(testData))
  testData$offset  <- rnorm(nrow(testData))

  modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10

  expect_error(dbarts(modelFormula, testData, subset = "not-a-number"))

  expect_error(dbarts(modelFormula, testData, weights = "not-a-number"))
  expect_error(dbarts(modelFormula, testData, weights = rep("not-a-number", nrow(testData))))
  expect_error(dbarts(modelFormula, testData, weights = offset))
  
  expect_error(dbarts(modelFormula, testData, offset = offset[1]))
})

test_that("compatibility specification raises errors", {
  expect_error(dbarts(testData$x, "not-a-number"))
  expect_error(dbarts(testData$x, testData$y[1]))
})

test_that("formula specification creates valid objects", {
  testData <- as.data.frame(testData)
  testData$weights <- runif(nrow(testData))
  testData$offset  <- rnorm(nrow(testData))

  modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10
  
  expect_is(dbarts(modelFormula, testData), "dbartsSampler")
  expect_is(dbarts(modelFormula, testData, weights = weights), "dbartsSampler")
  expect_is(dbarts(modelFormula, testData, offset = offset), "dbartsSampler")
  expect_is(dbarts(modelFormula, testData, weights = weights, offset = offset), "dbartsSampler")
  expect_is(dbarts(modelFormula, testData, subset = 1:10, weights = weights, offset = offset), "dbartsSampler")
})

test_that("compatibility specification creates valid objects", {
  testData$weights <- runif(length(testData$y))
  testData$offset  <- rnorm(length(testData$y))
  attach(testData)
  
  expect_is(dbarts(x, y), "dbartsSampler")
  expect_is(dbarts(x, y, weights = weights), "dbartsSampler")
  expect_is(dbarts(x, y, offset = offset), "dbartsSampler")
  expect_is(dbarts(x, y, weights = weights, offset = offset), "dbartsSampler")
  expect_is(dbarts(x, y, subset = 1:10, weights = weights, offset = offset), "dbartsSampler")

  detach(testData)
})

test_that("test argument raises errors", {
  expect_error(dbarts(y ~ x, testData, testData$x[11:20, 1:9]))
  expect_error(dbarts(y ~ x, testData, "not-a-matrix"))
  expect_error(dbarts(y ~ x, testData, outOfScope))

  test <- testData$x[11:20,]
  colnames(test) <- paste0("x", c(1:9, 11))
  expect_warning(dbarts(y ~ x, testData, test))
})

test_that("test argument creates valid objects", {
  ## test when is embedded in passed data
  testData$test <- testData$x[11:20,]
  expect_is(dbarts(y ~ x, testData, test), "dbartsSampler")
  expect_is(dbarts(y ~ x, testData, testData$test), "dbartsSampler")

  ## test when is in environment of formula
  test <- testData$test
  testData$test <- NULL
  expect_is(dbarts(y ~ x, testData, test), "dbartsSampler")
  expect_is(dbarts(y ~ x, testData, testData$x[11:20,]), "dbartsSampler")
})

