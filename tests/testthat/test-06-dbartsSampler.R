context("dbarts sampler as a discrete object")

source(system.file("common", "hillData.R", package = "dbarts"))

test_that("dbarts sampler settors raise errors", {
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test  <- data.frame(x = testData$x, z = 1 - testData$z)

  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L)
  sampler <- dbarts(y ~ x + z, train, test, control = control)

  expect_error(sampler$setControl("not-a-control"))
  
  expect_error(sampler$setResponse(numeric(0)))
  expect_error(sampler$setOffset(numeric(0)))
  expect_error(sampler$setPredictor(numeric(0), 1))
  expect_error(sampler$setPredictor(testData$z, 3))
  expect_error(sampler$setTestPredictor(numeric(0), 1))
  expect_error(sampler$setTestPredictor(numeric(0)))
  expect_error(sampler$setTestPredictor(testData$z, 3))
})

test_that("dbarts sampler updates correctly", {
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test <- data.frame(x = testData$x, z = 1 - testData$z)
  
  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L)
  sampler <- dbarts(y ~ x + z, train, test, control = control)

  n <- testData$n
  z <- testData$z

  sampler$setOffset(numeric(n))
  expect_equal(sampler$data@offset, numeric(n))
  sampler$setOffset(NULL)
  expect_identical(sampler$data@offset, NULL)
  
  sampler$setPredictor(x = z, column = 2)
  expect_equal(as.numeric(sampler$data@x[,2]), z)
  
  sampler$setTestPredictor(x = 1 - z, column = 2)
  expect_equal(as.numeric(sampler$data@x.test[,2]), 1 - z)

  sampler$setTestPredictors(NULL)
  expect_identical(sampler$data@x.test, NULL)

  sampler$setTestPredictors(test)
  expect_equal(sampler$data@x.test, as.matrix(test))
})

test_that("dbarts sampler shallow/deep copies", {
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test <- data.frame(x = testData$x, z = 1 - testData$z)
  
  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L)
  sampler <- dbarts(y ~ x + z, train, test, control = control)

  shallowCopy <- sampler$copy(TRUE)

  n <- testData$n
  
  sampler$setPredictor(numeric(n), 2)
  expect_equal(sampler$data@x, shallowCopy$data@x)
  
  rm(shallowCopy)
  gc(verbose = FALSE)

  deepCopy <- sampler$copy(FALSE)

  sampler$setPredictor(train$z, 2)
  expect_false(all(sampler$data@x[,2] == deepCopy$data@x[,2]))
})

source(system.file("common", "probitData.R", package = "dbarts"))

test_that("dbarts sampler correctly updates R test offsets only when applicable", {
  n <- nrow(testData$X)
    
  sampler <- dbarts(Z ~ X, testData, testData$X)
  
  sampler$setOffset(0.2)
  expect_equal(sampler$data@offset, rep_len(0.2, n))
  expect_null(sampler$data@offset.test)

  sampler$setOffset(NULL)
  expect_null(sampler$data@offset)
  expect_null(sampler$data@offset.test)

  sampler$setOffset(runif(n))
  expect_null(sampler$data@offset.test)
  
  
  sampler <- dbarts(Z ~ X, testData, offset = 0.2)

  expect_null(sampler$data@offset.test)

  sampler$setOffset(-0.1)
  expect_equal(sampler$data@offset, rep_len(-0.1, n))
  expect_null(sampler$data@offset.test)
  

  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2)

  sampler$setOffset(0.1)
  expect_equal(sampler$data@offset, rep_len(0.1, n))
  expect_equal(sampler$data@offset.test, rep_len(0.1, n))
  
  sampler$setTestOffset(0.2)
  expect_equal(sampler$data@offset.test, rep_len(0.2, n))

  sampler$setOffset(-0.1)
  expect_equal(sampler$data@offset, rep_len(-0.1, n))
  expect_equal(sampler$data@offset.test, rep_len(0.2, n))

  
  sampler <- dbarts(Z ~ X, testData, testData$X[-1,], offset = 0.2)

  expect_equal(sampler$data@offset.test, rep_len(0.2, n - 1))
  
  sampler$setOffset(0.1)
  expect_equal(sampler$data@offset, rep_len(0.1, n))
  expect_equal(sampler$data@offset.test, rep_len(0.1, n - 1))

  sampler$setOffset(rep_len(-0.1, n))
  expect_equal(sampler$data@offset, rep_len(-0.1, n))
  expect_null(sampler$data@offset.test)

  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2, offset.test = -0.1)
  sampler$setOffset(0.3)

  expect_equal(sampler$data@offset, rep_len(0.3, n))
  expect_equal(sampler$data@offset.test, rep_len(-0.1, n))
})

source(system.file("common", "hillData.R", package = "dbarts"))

test_that("dbarts sampler runs", {
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test <- data.frame(x = testData$x, z = 1 - testData$z)
  
  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L)
  sampler <- dbarts(y ~ x + z, train, test, control = control)

  n <- testData$n
  y <- testData$y
  z <- testData$z
  p <- testData$p

  set.seed(0)
  for (i in 1:5) {
    samples <- sampler$run()
    
    mu0 <- ifelse(z == 0, samples$train[,1], samples$test[,1])
    mu1 <- ifelse(z == 1, samples$train[,1], samples$test[,1])
    
    p0 <- dnorm(y, mu0, samples$sigma[1]) * (1 - p)
    p1 <- dnorm(y, mu1, samples$sigma[1]) * p
    p.z <- p1 / (p0 + p1)

    z <- rbinom(n, 1, p.z)

    n1 <- sum(z); n0 <- n - n1
    p <- rbeta(1, 1 + n0, 1 + n1)

    sampler$setPredictor(x = z, column = 2)
    sampler$setTestPredictor(x = 1 - z, column = 2)
  }

  expect_equal(z[1:20], c(0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0))
  expect_equal(p, 0.459823235825831)
  expect_equal(samples$train[1:5], c(94.4144689198833, 89.4520172313702, 93.2821083074909, 96.499552417188, 97.0400773954802))
  expect_equal(samples$test[1:5], c(89.5592665825945, 94.4797604402427, 94.4227067142213, 92.7638379016566, 91.7194156602015))
})

source(system.file("common", "probitData.R", package = "dbarts"))

test_that("dbarts sampler updates offsets in C++", {
  control <- dbartsControl(updateState = FALSE, n.burn = 0L, n.samples = 1L, verbose = FALSE)
  sampler <- dbarts(Z ~ X, testData, testData$X[1:200,], 1:200,
                    control = control)
  
  sampler$setOffset(0.5)
  set.seed(0)
  samples <- sampler$run(25, 1)
  expect_equal(as.double(samples$train - samples$test), rep_len(0.5, 200))

  sampler$setTestOffset(0.5)
  samples <- sampler$run(0, 1)
  expect_equal(samples$train, samples$test)
})
