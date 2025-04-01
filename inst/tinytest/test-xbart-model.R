source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that works with non-standard models
x <- testData$x
y <- testData$y

k <- c(4, 8)
expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    k = k, n.threads = 1L, resid.prior = chisq(2.5, 0.9)
  )
)

expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    k = k, n.threads = 1L, resid.prior = fixed(2)
  )
)

n.trees <- c(5L, 10L)
expect_silent(
  dbarts::xbart(
    x, y, method = "k-fold",
    n.reps = 3L, n.samples = 6L, n.burn = c(10L, 5L, 1L), n.test = 5,
    n.trees = n.trees, n.threads = 1L
  )
)


rm(n.trees, k, y, x)

rm(testData)

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that runs with binary data and k hyperprior
x <- testData$X
z <- testData$Z

n.reps  <- 3L
power   <- c(1.5, 2)

xval <- dbarts::xbart(
  x, z, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
  n.reps = n.reps,
  power = power,
  n.threads = 2L
)

expect_inherits(xval, "matrix")
expect_equal(dim(xval), c(n.reps, length(power)))
expect_true(!anyNA(xval))
expect_equal(
  dimnames(xval),
  list(
    rep     = NULL,
    power   = as.character(power)
  )
)

rm(xval, power, n.reps, z, x)

rm(testData)

