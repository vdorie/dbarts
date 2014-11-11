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

  n <- length(testData$y)
  expect_error(sampler$setPredictor(matrix(numeric(n * 3), n)))
  expect_error(sampler$setPredictor(matrix(numeric((n - 1) * 2), n - 1)))
  expect_error(sampler$setPredictor(matrix(numeric(n * 2), n), 1))
})

test_that("dbarts sampler updates predictors correctly", {
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
  
  invisible(sampler$setPredictor(numeric(n), 2))
  expect_equal(as.numeric(sampler$data@x[,2]), numeric(n))
  
  invisible(sampler$setPredictor(x = z, column = 2))
  expect_equal(as.numeric(sampler$data@x[,2]), z)
  
  sampler$setTestPredictor(x = 1 - z, column = 2)
  expect_equal(as.numeric(sampler$data@x.test[,2]), 1 - z)

  sampler$setTestPredictor(NULL)
  expect_identical(sampler$data@x.test, NULL)

  sampler$setTestPredictor(test)
  expect_equal(sampler$data@x.test, as.matrix(test))

  new.x <- rnorm(n)
  new.z <- as.double(rbinom(n, 1, 0.5))
  new.data <- cbind(new.x, new.z)
  invisible(sampler$setPredictor(new.data))

  
  expect_equal(as.numeric(sampler$data@x), as.numeric(new.data))
})

test_that("dbarts sampler shallow/deep copies", {
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test <- data.frame(x = testData$x, z = 1 - testData$z)
  
  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L)
  sampler <- dbarts(y ~ x + z, train, test, control = control)

  shallowCopy <- sampler$copy(TRUE)

  n <- testData$n
  
  invisible(sampler$setPredictor(numeric(n), 2))
  expect_equal(sampler$data@x, shallowCopy$data@x)
  
  rm(shallowCopy)
  gc(verbose = FALSE)

  deepCopy <- sampler$copy(FALSE)

  invisible(sampler$setPredictor(1 - train$z, 2))
  expect_false(all(sampler$data@x[,2] == deepCopy$data@x[,2]))
  
  invisible(sampler$setPredictor(deepCopy$data@x[,2], 2))
  expect_equal(sampler$data@x, deepCopy$data@x)
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

    new.z <- rbinom(n, 1, p.z)
    while (sampler$setPredictor(x = new.z, column = 2) == FALSE) {
      new.z <- rbinom(n, 1, p.z)
    }
    z <- new.z
    sampler$setTestPredictor(x = 1 - z, column = 2)

    n1 <- sum(z); n0 <- n - n1
    p <- rbeta(1, 1 + n0, 1 + n1)
  }

  expect_equal(z[1:20], c(0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0))
  expect_equal(p, 0.470788237121872)
  expect_equal(samples$train[1:5], c(89.6662721124284, 89.6662721124284, 94.4644684264674, 92.3429201084874, 91.1812768131237))
  expect_equal(samples$test[1:5], c(94.6496627941143, 94.6496627941143, 93.0602674319221, 91.9088151870825, 95.6941612719439))
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
