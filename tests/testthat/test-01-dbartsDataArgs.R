context("dbarts data arguments")

source(system.file("common", "friedmanData.R", package = "dbarts"))

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

## rm(testData)

source(system.file("common", "probitData.R", package = "dbarts"))

test_that("test offset fills in control logicals depending on specification", {
  sampler <- dbarts(Z ~ X, testData, testData$X)
  
  expect_that(sampler$data@offset,      is_null())
  expect_that(sampler$data@offset.test, is_null())
  expect_that(sampler$data@testUsesRegularOffset, equals(NA))

  
  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2)
  
  expect_that(sampler$data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(sampler$data@offset.test[1:5], equals(rep(0.2, 5)))
  expect_that(sampler$data@testUsesRegularOffset, equals(TRUE))


  otherOffset <- 0.2 + 0.1
  sampler <- dbarts(Z ~ X, testData, testData$X, offset = otherOffset)

  expect_that(sampler$data@offset[1:5],      equals(rep(0.3, 5)))
  expect_that(sampler$data@offset.test[1:5], equals(rep(0.3, 5)))
  expect_that(sampler$data@testUsesRegularOffset, equals(TRUE))
  

  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2, offset.test = NULL)

  expect_that(sampler$data@offset[1:5], equals(rep(0.2, 5)))
  expect_that(sampler$data@offset.test, is_null())
  expect_that(sampler$data@testUsesRegularOffset, equals(FALSE))
  

  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2, offset.test = 0.1)

  expect_that(sampler$data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(sampler$data@offset.test[1:5], equals(rep(0.1, 5)))
  expect_that(sampler$data@testUsesRegularOffset, equals(FALSE))

  
  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2, offset.test = offset + 0.1)
  
  expect_that(sampler$data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(sampler$data@offset.test[1:5], equals(rep(0.3, 5)))
  expect_that(sampler$data@testUsesRegularOffset, equals(FALSE))


  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2, offset.test = offset)
  
  expect_that(sampler$data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(sampler$data@offset.test[1:5], equals(rep(0.2, 5)))
  expect_that(sampler$data@testUsesRegularOffset, equals(TRUE))


  otherOffset <- runif(nrow(testData$X))
  sampler <- dbarts(Z ~ X, testData, testData$X, offset = otherOffset)
  
  expect_that(sampler$data@offset[1:5],      equals(otherOffset[1:5]))
  expect_that(sampler$data@offset.test[1:5], equals(otherOffset[1:5]))
  expect_that(sampler$data@testUsesRegularOffset, equals(TRUE))


  expect_error(dbarts(Z ~ X, testData, testData$X[-1,], offset = otherOffset))

  
  sampler <- dbarts(Z ~ X, testData, testData$X, offset = otherOffset, offset.test = NULL)
  
  expect_that(sampler$data@offset[1:5], equals(otherOffset[1:5]))
  expect_that(sampler$data@offset.test, is_null())
  expect_that(sampler$data@testUsesRegularOffset, equals(FALSE))


  sampler <- dbarts(Z ~ X, testData, testData$X, offset = otherOffset, offset.test = 0.2)

  expect_that(sampler$data@offset[1:5],      equals(otherOffset[1:5]))
  expect_that(sampler$data@offset.test[1:5], equals(rep(0.2, 5)))
  expect_that(sampler$data@testUsesRegularOffset, equals(FALSE))


  sampler <- dbarts(Z ~ X, testData, testData$X, offset = otherOffset, offset.test = offset + 0.1)

  expect_that(sampler$data@offset[1:5],      equals(otherOffset[1:5]))
  expect_that(sampler$data@offset.test[1:5], equals(otherOffset[1:5] + 0.1))
  expect_that(sampler$data@testUsesRegularOffset, equals(FALSE))


  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2, offset.test = otherOffset)
  
  expect_that(sampler$data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(sampler$data@offset.test[1:5], equals(otherOffset[1:5]))
  expect_that(sampler$data@testUsesRegularOffset, equals(FALSE))

  
  sampler <- dbarts(Z ~ X, testData, testData$X, offset = NULL, offset.test = otherOffset)
  
  expect_that(sampler$data@offset,           is_null())
  expect_that(sampler$data@offset.test[1:5], equals(otherOffset[1:5]))
  expect_that(sampler$data@testUsesRegularOffset, equals(FALSE))
})
