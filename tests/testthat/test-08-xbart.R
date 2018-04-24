context("xbart")

source(system.file("common", "friedmanData.R", package = "dbarts"))

test_that("random subsample runs correctly with valid inputs", {
  x <- testData$x
  y <- testData$y
  
  n.reps  <- 8L
  n.trees <- c(5L, 7L)
  k       <- c(1, 2, 4)
  power   <- c(1.5, 2)
  base    <- c(0.75, 0.8, 0.95)

  xval <- xbart(x, y, n.samples = 15L, n.burn = c(10L, 3L, 1L), method = "random subsample",
                n.reps = n.reps,
                n.trees = n.trees,
                k = k,
                power = power,
                base = base,
                n.threads = 2L)
  
  expect_is(xval, "array")
  expect_equal(dim(xval), c(n.reps, length(n.trees), length(k), length(power), length(base)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)))
})

test_that("k-fold runs correctly with valid inputs", {
  x <- testData$x
  y <- testData$y
  
  n.reps  <- 8L
  n.trees <- c(5L, 7L)
  k       <- c(1, 2, 4)
  power   <- c(1.5, 2)
  base    <- c(0.75, 0.8, 0.95)

  xval <- xbart(x, y, n.samples = 15L, n.burn = c(10L, 3L, 1L), method = "k-fold", n.test = 6,
                n.reps = n.reps,
                n.trees = n.trees,
                k = k,
                power = power,
                base = base,
                n.threads = 2L)
  
  expect_is(xval, "array")
  expect_equal(dim(xval), c(n.reps, length(n.trees), length(k), length(power), length(base)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)))
})

test_that("k-fold and random subsample are roughly similar", {
  x <- testData$x
  y <- testData$y
  
  k <- c(0.5, 8)
  
  set.seed(0)
  xval.kf <- xbart(x, y, method = "k-fold",
                   n.reps = 20L, n.samples = 100L, n.burn = c(100L, 50L, 10L),
                   k = k,
                   n.threads = 1L)

  xval.rs <- xbart(x, y, method = "random subsample", 
                   n.reps = 100L, n.samples = 100L, n.burn = c(100L, 50L, 10L),
                   k = k, 
                   n.threads = 1L)
  
  expect_true(sign(diff(apply(xval.rs, 2, mean))) == sign(diff(apply(xval.kf, 2, mean))))
})

test_that("works with custom loss", {
  x <- testData$x
  y <- testData$y
  
  n.reps  <- 8L
  n.trees <- c(5L, 7L)
  k       <- c(1, 2, 4)
  power   <- c(1.5, 2)
  base    <- c(0.75, 0.8, 0.95)
  
  mad <- function(y.train, y.train.hat) 
    mean(abs(y.train - apply(y.train.hat, 1L, mean)))

  xval <- xbart(x, y, n.samples = 15L, n.burn = c(10L, 3L, 1L), method = "k-fold", n.test = 6,
                n.reps = n.reps,
                n.trees = n.trees,
                k = k,
                power = power,
                base = base, loss = mad,
                n.threads = 1L)
  
  expect_is(xval, "array")
  expect_equal(dim(xval), c(n.reps, length(n.trees), length(k), length(power), length(base)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)))
  
  
  xval <- xbart(x, y, n.samples = 15L, n.burn = c(10L, 3L, 1L), method = "k-fold", n.test = 6,
                n.reps = n.reps,
                n.trees = n.trees,
                k = k,
                power = power,
                base = base, loss = mad,
                n.threads = 2L)
  
  expect_is(xval, "array")
  expect_equal(dim(xval), c(n.reps, length(n.trees), length(k), length(power), length(base)))
  expect_true(!anyNA(xval))
  expect_equal(dimnames(xval), list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)))
})

test_that("fails with invalid inputs", {
  x <- testData$x
  y <- testData$y
  
  expect_error(xbart(y ~ x, method = "not-a-method"))
  expect_error(xbart(y ~ x, method = NULL))
  expect_error(xbart(y ~ x, method = NA_character_))
  
  expect_error(xbart(y ~ x, n.samples = 0L))
  expect_error(xbart(y ~ x, n.samples = "not-a-integer"))
  expect_error(xbart(y ~ x, n.samples = NULL))
  expect_error(xbart(y ~ x, n.samples = NA_integer_))
  
  expect_error(xbart(y ~ x, method = "k-fold", n.test = 1))
  expect_error(xbart(y ~ x, method = "k-fold", n.test = length(testData$y) + 1))
  expect_error(xbart(y ~ x, method = "random subsample", n.test = 0))
  expect_error(xbart(y ~ x, method = "random subsample", n.test = length(testData$y) + 1))
  expect_error(xbart(y ~ x, n.test = "not-a-numeric"))
  expect_error(xbart(y ~ x, n.test = NULL))
  expect_error(xbart(y ~ x, n.test = NA_real_))
  
  expect_error(xbart(y ~ x, n.reps = 0L))
  expect_error(xbart(y ~ x, n.reps = "not-a-integer"))
  expect_error(xbart(y ~ x, n.reps = NULL))
  expect_error(xbart(y ~ x, n.reps = NA_integer_))
  
  expect_error(xbart(y ~ x, n.burn = c(200L, -1L, 50L)))
  expect_error(xbart(y ~ x, n.burn = "not-a-integer"))
  expect_error(xbart(y ~ x, n.burn = NULL))
  expect_error(xbart(y ~ x, n.burn = NA_integer_))
  
  expect_error(xbart(y ~ x, loss = "unknown-loss"))
  expect_error(xbart(y ~ x, loss = 2))
  expect_error(xbart(y ~ x, loss = function(x) x))
  expect_error(xbart(y ~ x, loss = list(function(x, y) NULL, "not-an-environment")))
  
  expect_error(xbart(y ~ x, n.threads = 0L))
  expect_error(xbart(y ~ x, n.threads = "not-a-integer"))
  expect_error(xbart(y ~ x, n.threads = NULL))
  
  expect_error(xbart(y ~ x, n.trees = 0L))
  expect_error(xbart(y ~ x, n.trees = "not-a-integer"))
  expect_error(xbart(y ~ x, n.trees = NULL))
  expect_error(xbart(y ~ x, n.trees = NA_integer_))
  
  expect_error(xbart(y ~ x, k = c(-0.5, 1)))
  expect_error(xbart(y ~ x, k = "not-a-numeric"))
  expect_error(xbart(y ~ x, k = NULL))
  expect_error(xbart(y ~ x, k = NA_real_))
  
  expect_error(xbart(y ~ x, power = c(0, 0.5)))
  expect_error(xbart(y ~ x, power = "not-a-numeric"))
  expect_error(xbart(y ~ x, power = NULL))
  expect_error(xbart(y ~ x, power = NA_real_))
  
  expect_error(xbart(y ~ x, base = c(0.5, 1)))
  expect_error(xbart(y ~ x, base = "not-a-numeric"))
  expect_error(xbart(y ~ x, base = NULL))
  expect_error(xbart(y ~ x, base = NA_real_))
  
  expect_error(xbart(y ~ x, verbose = "not-a-logical"))
  expect_error(xbart(y ~ x, verbose = NULL))
  expect_error(xbart(y ~ x, verbose = NA))
  
  expect_error(xbart(y ~ x, resid.prior = NULL))
  
  expect_error(xbart(y ~ x, sigma = -1))
  expect_error(xbart(y ~ x, sigma = "not-a-numeric"))
})
