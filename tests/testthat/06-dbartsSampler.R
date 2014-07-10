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

  sampler$setTestPredictor(NULL)
  expect_identical(sampler$data@x.test, NULL)

  sampler$setTestPredictor(test)
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

  sampler$setPredictor(z, 2)
  expect_false(all(sampler$data@x[,2] == deepCopy$data@x[,2]))
})
  
  

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

  expect_equal(z[1:20], c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(p, 0.852304686322099)
  expect_equal(samples$train[1:5], c(92.5282515859516, 92.5282515859516, 93.4420183143099, 90.9300810722605, 92.8149376208215))
  expect_equal(samples$test[1:5], c(96.4889413170199, 96.4889413170199, 97.3120703909065, 94.8907708033289, 96.7676289468002))
})
