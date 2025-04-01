source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that random subsample runs correctly with valid inputs
x <- testData$x
y <- testData$y

n.reps  <- 3L
n.trees <- c(5L, 7L)
k       <- c(1, 2, 4)
power   <- c(1.5, 2)
base    <- c(0.75, 0.8, 0.95)

xval <- dbarts::xbart(
  x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "random subsample",
  n.reps = n.reps,
  n.trees = n.trees,
  k = k,
  power = power,
  base = base,
  n.threads = 2L
)

expect_inherits(xval, "array")
expect_equal(
  dim(xval),
  c(n.reps, length(n.trees), length(k), length(power), length(base))
)
expect_true(!anyNA(xval))
expect_equal(
  dimnames(xval),
  list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)
  )
)

rm(xval, base, power, k, n.trees, n.reps, y, x)


# test that k-fold runs correctly with valid inputs
x <- testData$x
y <- testData$y

n.reps  <- 3L
n.trees <- c(5L, 7L)
k       <- c(1, 2, 4)
power   <- c(1.5, 2)
base    <- c(0.75, 0.8, 0.95)

xval <- dbarts::xbart(
  x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
  n.reps = n.reps,
  n.trees = n.trees,
  k = k,
  power = power,
  base = base,
  n.threads = 2L
)

expect_inherits(xval, "array")
expect_equal(
  dim(xval),
  c(n.reps, length(n.trees), length(k), length(power), length(base))
)
expect_true(!anyNA(xval))
expect_equal(
  dimnames(xval),
  list(
    rep     = NULL,
    n.trees = as.character(n.trees),
    k       = as.character(k),
    power   = as.character(power),
    base    = as.character(base)
  )
)

rm(xval, base, power, k, n.trees, n.reps, y, x)

# test that k-fold runs correctly with one input
x <- testData$x
y <- testData$y

xval <- dbarts::xbart(
  x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L),
  n.reps = 3, n.test = 5,
  k = 2,
  n.threads = 2L
)

expect_equal(length(xval), 3L)

rm(xval, y, x)


# test that k-fold subdivides data correctly when data do not divide evenly by k
x <- testData$x[1L:24L,]
y <- testData$y[1L:24L]

k <- c(2, 4)

xval <- dbarts::xbart(
  x, y, n.samples = 6L, n.burn = c(5L, 3L, 1L), method = "k-fold", n.test = 5,
  n.reps = 3L,
  k = k,
  n.threads = 1L
)

expect_inherits(xval, "array")

rm(testData)

