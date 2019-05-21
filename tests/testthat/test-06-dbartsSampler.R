context("dbarts sampler as a discrete object")

source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

test_that("dbarts sampler settors raise errors", {
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test  <- data.frame(x = testData$x, z = 1 - testData$z)

  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L,
                           n.chains = 1L, n.threads = 1L)
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
                           n.burn = 0L, n.samples = 1L, n.thin = 5L,
                           n.chains = 1L, n.threads = 1L)
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

test_that("dbarts sampler with matrix specification doesn't change variables in parent frame", {
  x.train <- dbarts::makeModelMatrixFromDataFrame(data.frame(x = testData$x, z = testData$z))
  y.train <- testData$y
  
  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L,
                           n.chains = 1L, n.threads = 1L)
  sampler <- dbarts(x.train, y.train, control = control)
  
  n <- testData$n
  z <- testData$z
  
  invisible(sampler$setPredictor(numeric(n), 2))
  expect_equal(as.numeric(sampler$data@x[,2]), numeric(n))
  expect_equal(as.numeric(x.train[,2]), z)
})

test_that("dbarts sampler setData yields valid model", {
  n <- 105L

  x <- runif(n, -10, 10)
  x.test <- seq.int(min(x), max(x), length.out = 101L)

  beta.1 <- -0.75
  beta.2 <-  0.5

  n.cuts <- 10L
  cutoffs <- min(x) + seq_len(n.cuts) * (max(x) - min(x)) / (n.cuts + 1)

  y <- ifelse(x <= cutoffs[1L] | x > cutoffs[n.cuts], beta.1, beta.2) + rnorm(n, 0, 0.15)
  
  control <- dbartsControl(n.trees = 1L, n.cuts = n.cuts, n.chains = 1L, n.threads = 1L, updateState = FALSE)
  sampler <- dbarts(y ~ x, control = control)

  samples1 <- sampler$run(500L, 1000L)

  x.new <- x + diff(cutoffs[1L:2L])

  sampler$setData(dbartsData(y ~ x.new))
  samples2 <- sampler$run(500L, 1000L)
  
  expect_equal(sd(apply(samples1$train, 1, mean) - y),
               sd(apply(samples2$train, 1, mean) - y),
               tol = 1.0e-2)
})

test_that("dbarts sampler shallow/deep copies", {
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test <- data.frame(x = testData$x, z = 1 - testData$z)
  
  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L,
                           n.chains = 1L, n.threads = 1L)
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

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

test_that("dbarts sampler correctly updates R test offsets only when applicable", {
  n <- nrow(testData$X)
  control <- dbartsControl(n.chains = 1L, n.threads = 2L, updateState = FALSE)
  
  sampler <- dbarts(Z ~ X, testData, testData$X, control = control)
  
  sampler$setOffset(0.2)
  expect_equal(sampler$data@offset, rep_len(0.2, n))
  expect_null(sampler$data@offset.test)

  sampler$setOffset(NULL)
  expect_null(sampler$data@offset)
  expect_null(sampler$data@offset.test)

  sampler$setOffset(runif(n))
  expect_null(sampler$data@offset.test)
  
  
  sampler <- dbarts(Z ~ X, testData, offset = 0.2, control = control)

  expect_null(sampler$data@offset.test)

  sampler$setOffset(-0.1)
  expect_equal(sampler$data@offset, rep_len(-0.1, n))
  expect_null(sampler$data@offset.test)
  

  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2, control = control)

  sampler$setOffset(0.1)
  expect_equal(sampler$data@offset, rep_len(0.1, n))
  expect_equal(sampler$data@offset.test, rep_len(0.1, n))
  
  sampler$setTestOffset(0.2)
  expect_equal(sampler$data@offset.test, rep_len(0.2, n))

  sampler$setOffset(-0.1)
  expect_equal(sampler$data@offset, rep_len(-0.1, n))
  expect_equal(sampler$data@offset.test, rep_len(0.2, n))

  
  sampler <- dbarts(Z ~ X, testData, testData$X[-1,], offset = 0.2, control = control)

  expect_equal(sampler$data@offset.test, rep_len(0.2, n - 1))
  
  sampler$setOffset(0.1)
  expect_equal(sampler$data@offset, rep_len(0.1, n))
  expect_equal(sampler$data@offset.test, rep_len(0.1, n - 1))

  sampler$setOffset(rep_len(-0.1, n))
  expect_equal(sampler$data@offset, rep_len(-0.1, n))
  expect_null(sampler$data@offset.test)

  sampler <- dbarts(Z ~ X, testData, testData$X, offset = 0.2, offset.test = -0.1, control = control)
  sampler$setOffset(0.3)

  expect_equal(sampler$data@offset, rep_len(0.3, n))
  expect_equal(sampler$data@offset.test, rep_len(-0.1, n))
})

source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

fitTinyModel <- function(testData, n.samples) {  
  train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
  test <- data.frame(x = testData$x, z = 1 - testData$z)
  
  control <- dbartsControl(updateState = FALSE, verbose = FALSE,
                           n.burn = 0L, n.samples = 1L, n.thin = 5L,
                           n.chains = 1L, n.threads = 1L, n.trees = 5)
  sampler <- dbarts(y ~ x + z, train, test, control = control)
  
  n <- testData$n
  y <- testData$y
  z <- testData$z
  p <- testData$p
  
  set.seed(0L)
  for (i in seq_len(n.samples)) {
    samples <- sampler$run()
    
    mu0 <- ifelse(z == 0L, samples$train[,1L], samples$test[,1L])
    mu1 <- ifelse(z == 1L, samples$train[,1L], samples$test[,1L])
    
    p0 <- dnorm(y, mu0, samples$sigma[1L]) * (1 - p)
    p1 <- dnorm(y, mu1, samples$sigma[1L]) * p
    p.z <- p1 / (p0 + p1)

    new.z <- rbinom(n, 1L, p.z)
    sampler$setPredictor(x = new.z, column = 2L, forceUpdate = TRUE)
    z <- new.z
    sampler$setTestPredictor(x = 1L - z, column = 2L)

    n1 <- sum(z); n0 <- n - n1
    p <- rbeta(1L, 1L + n0, 1L + n1)
  }
  
  dbarts:::namedList(sampler, z, p, samples)
}

test_that("dbarts sampler runs", {
  massign[sampler, z, p, samples] <- fitTinyModel(testData, 5L)
  
  expect_equal(z[1:20], c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0))
  expect_equal(p, 0.939263201710153)
  expect_equal(samples$train[1:5], c(92.2954347050941, 92.2954347050941, 93.2401525161323, 92.2954347050941, 93.2708029591151))
  expect_equal(samples$test[1:5], c(93.6109481580171, 93.6109481580171, 94.5556659690554, 93.6109481580171, 94.5863164120381))
  
  tree <- sampler$getTrees(5L)
  expect_equal(tree$n, c(120, 61, 59, 54, 5))
  expect_equal(tree$var, c(1, -1, 1, -1, -1))
  expect_equal(tree$value, c(34.8679947872694, -0.0436339769593335, 57.3007502897521, -0.0220356180367407, 0.0966367066202154))
})

test_that("sampler sets cut points correctly", {
  sampler <- fitTinyModel(testData, 5L)$sampler

  tree.orig <- sampler$getTrees(1L)
  
  sampler$storeState()
  sampler$setCutPoints(attr(sampler$state, "cutPoints")[[1L]][1L:87L], 1L)
  tree.new <- sampler$getTrees(1L)
  
  expect_equal(nrow(tree.orig), 5L)
  expect_equal(nrow(tree.new), 3L)
  expect_equal(tree.orig$value[4L], tree.new$value[3L])
  
  
  sampler <- fitTinyModel(testData, 5L)$sampler
  
  textCon <- textConnection("output.orig", open = "w")
  sink(textCon)
  sampler$printTrees(1L)
  sink()
  close(textCon)
  
  sampler$storeState()
  sampler$setCutPoints(attr(sampler$state, "cutPoints")[[1L]][1L:88L], 1L)
  
  textCon <- textConnection("output.new", open = "w")
  sink(textCon)
  sampler$printTrees(1L)
  sink()
  close(textCon)
  
  expect_equal(length(output.orig), length(output.new))
  expect_true(grepl("Avail: 01", output.new[5L]))
})

test_that("sampler set data collapses nodes correctly", {
  sampler <- fitTinyModel(testData, 15L)$sampler
  
  tree.orig <- sampler$getTrees(2L)
  
  x.new <- sampler$data@x
  
  valuesToSwap <- which(x.new[,1L] <= tree.orig$value[1L] & x.new[,2L] <= tree.orig$value[2L])
  targetValues <- which(x.new[,1L] > tree.orig$value[1L] & x.new[,2L] > tree.orig$value[2L])[seq_along(valuesToSwap)]
  
  temp <- x.new[valuesToSwap,2L]
  x.new[valuesToSwap,2L] <- x.new[targetValues,2L]
  x.new[targetValues,2L] <- temp
  
  data.new <- sampler$data
  data.new@x <- x.new
  sampler$setData(data.new)
  
  tree.new <- sampler$getTrees(2L)
  
  expect_equal(tree.new$value[2L], tree.orig$value[4L])
})

rm(fitTinyModel)

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

test_that("dbarts sampler updates offsets in C++", {
  control <- dbartsControl(updateState = FALSE, n.burn = 0L, n.samples = 1L, verbose = FALSE,
                           n.chains = 1L, n.threads = 1L)
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

