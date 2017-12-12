context("dbarts data arguments")

source(system.file("common", "friedmanData.R", package = "dbarts"))

test_that("formula specification raises errors", {
  expect_error(dbartsData("not-a-formula", testData))
  expect_error(dbartsData(y ~ x))
  expect_error(dbartsData(NULL, testData))
  expect_error(dbartsData(y ~ 0, testData))
})

test_that("extra arguments for formula specification raises errors", {
  testData <- as.data.frame(testData)
  testData$weights <- runif(nrow(testData))
  testData$offset  <- rnorm(nrow(testData))

  modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10

  expect_error(dbartsData(modelFormula, testData, subset = "not-a-number"))

  expect_error(dbartsData(modelFormula, testData, weights = "not-a-number"))
  expect_error(dbartsData(modelFormula, testData, weights = rep("not-a-number", nrow(testData))))
  expect_error(dbartsData(modelFormula, testData, weights = offset))
})

test_that("compatibility specification raises errors", {
  expect_error(dbartsData(testData$x, "not-a-number"))
  expect_error(dbartsData(testData$x, testData$y[1]))
})

test_that("formula specification creates valid objects", {
  testData <- as.data.frame(testData)
  testData$weights <- runif(nrow(testData))
  testData$offset  <- rnorm(nrow(testData))

  modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10
  
  expect_is(dbartsData(modelFormula, testData), "dbartsData")
  expect_is(dbartsData(modelFormula, testData, weights = weights), "dbartsData")
  expect_is(dbartsData(modelFormula, testData, offset = offset), "dbartsData")
  expect_is(dbartsData(modelFormula, testData, weights = weights, offset = offset), "dbartsData")
  expect_is(dbartsData(modelFormula, testData, subset = 1:10, weights = weights, offset = offset), "dbartsData")
})

test_that("compatibility specification creates valid objects", {
  testData$weights <- runif(length(testData$y))
  testData$offset  <- rnorm(length(testData$y))
  attach(testData)
  
  expect_is(dbartsData(x, y), "dbartsData")
  expect_is(dbartsData(x, y, weights = weights), "dbartsData")
  expect_is(dbartsData(x, y, offset = offset), "dbartsData")
  expect_is(dbartsData(x, y, weights = weights, offset = offset), "dbartsData")
  expect_is(dbartsData(x, y, subset = 1:10, weights = weights, offset = offset), "dbartsData")
  
  detach(testData)
})

test_that("compatibility specification works with dimnames", {
  x <- testData$x
  y <- testData$y
  
  colnames(x) <- paste0("x.", seq_len(ncol(x)))
  expect_is(dbartsData(x, y), "dbartsData")
  
  x <- x[,1L,drop=FALSE]
  expect_is(dbartsData(x, y), "dbartsData")
})

test_that("compatibility specification works with duplicated dimnames", {
  x <- testData$x
  y <- testData$y
  
  colnames(x) <- c(paste0("x.", seq_len(ncol(x) - 2L)), "", "")
  x.test <- x
  x.test[,ncol(x.test)] <- x.test[,ncol(x.test)] + 1
   
  data <- dbartsData(x, y, x.test)
  expect_equal(data@x.test, x.test)
})

test_that("test argument raises errors", {
  expect_error(dbartsData(y ~ x, testData, testData$x[11:20, 1:9]))
  expect_error(dbartsData(y ~ x, testData, "not-a-matrix"))
  expect_error(dbartsData(y ~ x, testData, outOfScope))

  test <- testData$x[11:20,]
  colnames(test) <- paste0("x", c(1:9, 11))
  expect_warning(dbartsData(y ~ x, testData, test))
})

test_that("test argument creates valid objects", {
  ## test when is embedded in passed data
  testData$test <- testData$x[11:20,]
  expect_is(dbartsData(y ~ x, testData, test), "dbartsData")
  expect_is(dbartsData(y ~ x, testData, testData$test), "dbartsData")

  ## test when is in environment of formula
  test <- testData$test
  testData$test <- NULL
  expect_is(dbartsData(y ~ x, testData, test), "dbartsData")
  expect_is(dbartsData(y ~ x, testData, testData$x[11:20,]), "dbartsData")
})

## rm(testData)

source(system.file("common", "probitData.R", package = "dbarts"))

test_that("test offset fills in control logicals depending on specification", {
  data <- dbartsData(Z ~ X, testData, testData$X)
  
  expect_that(data@offset,      is_null())
  expect_that(data@offset.test, is_null())
  expect_that(data@testUsesRegularOffset, equals(NA))

  
  data <- dbartsData(Z ~ X, testData, testData$X, offset = 0.2)
  
  expect_that(data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(data@offset.test[1:5], equals(rep(0.2, 5)))
  expect_that(data@testUsesRegularOffset, equals(TRUE))


  otherOffset <- 0.2 + 0.1
  data <- dbartsData(Z ~ X, testData, testData$X, offset = otherOffset)

  expect_that(data@offset[1:5],      equals(rep(0.3, 5)))
  expect_that(data@offset.test[1:5], equals(rep(0.3, 5)))
  expect_that(data@testUsesRegularOffset, equals(TRUE))
  

  data <- dbartsData(Z ~ X, testData, testData$X, offset = 0.2, offset.test = NULL)

  expect_that(data@offset[1:5], equals(rep(0.2, 5)))
  expect_that(data@offset.test, is_null())
  expect_that(data@testUsesRegularOffset, equals(FALSE))
  

  data <- dbartsData(Z ~ X, testData, testData$X, offset = 0.2, offset.test = 0.1)

  expect_that(data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(data@offset.test[1:5], equals(rep(0.1, 5)))
  expect_that(data@testUsesRegularOffset, equals(FALSE))

  
  data <- dbartsData(Z ~ X, testData, testData$X, offset = 0.2, offset.test = offset + 0.1)
  
  expect_that(data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(data@offset.test[1:5], equals(rep(0.3, 5)))
  expect_that(data@testUsesRegularOffset, equals(FALSE))


  data <- dbartsData(Z ~ X, testData, testData$X, offset = 0.2, offset.test = offset)
  
  expect_that(data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(data@offset.test[1:5], equals(rep(0.2, 5)))
  expect_that(data@testUsesRegularOffset, equals(TRUE))


  otherOffset <- runif(nrow(testData$X))
  data <- dbartsData(Z ~ X, testData, testData$X, offset = otherOffset)
  
  expect_that(data@offset[1:5],      equals(otherOffset[1:5]))
  expect_that(data@offset.test[1:5], equals(otherOffset[1:5]))
  expect_that(data@testUsesRegularOffset, equals(TRUE))


  expect_error(dbartsData(Z ~ X, testData, testData$X[-1,], offset = otherOffset))

  
  data <- dbartsData(Z ~ X, testData, testData$X, offset = otherOffset, offset.test = NULL)
  
  expect_that(data@offset[1:5], equals(otherOffset[1:5]))
  expect_that(data@offset.test, is_null())
  expect_that(data@testUsesRegularOffset, equals(FALSE))


  data <- dbartsData(Z ~ X, testData, testData$X, offset = otherOffset, offset.test = 0.2)

  expect_that(data@offset[1:5],      equals(otherOffset[1:5]))
  expect_that(data@offset.test[1:5], equals(rep(0.2, 5)))
  expect_that(data@testUsesRegularOffset, equals(FALSE))


  data <- dbartsData(Z ~ X, testData, testData$X, offset = otherOffset, offset.test = offset + 0.1)

  expect_that(data@offset[1:5],      equals(otherOffset[1:5]))
  expect_that(data@offset.test[1:5], equals(otherOffset[1:5] + 0.1))
  expect_that(data@testUsesRegularOffset, equals(FALSE))


  data <- dbartsData(Z ~ X, testData, testData$X, offset = 0.2, offset.test = otherOffset)
  
  expect_that(data@offset[1:5],      equals(rep(0.2, 5)))
  expect_that(data@offset.test[1:5], equals(otherOffset[1:5]))
  expect_that(data@testUsesRegularOffset, equals(FALSE))

  
  data <- dbartsData(Z ~ X, testData, testData$X, offset = NULL, offset.test = otherOffset)
  
  expect_that(data@offset,           is_null())
  expect_that(data@offset.test[1:5], equals(otherOffset[1:5]))
  expect_that(data@testUsesRegularOffset, equals(FALSE))
})

source(system.file("common", "almostLinearBinaryData.R", package = "dbarts"))

test_that("bart creates viable sampler with formula, data specification", {
  data <- data.frame(y = testData$y, x = testData$x)
  modelFormula <- y ~ x.1 + x.2 + x.3
  
  expect_is(bart(modelFormula, data, nskip = 0L, ndpost = 1L, verbose = FALSE), "bart")
  expect_is(bart(modelFormula, data[1:100,], data[101:200,], nskip = 0L, ndpost = 1L, verbose = FALSE), "bart")
})

