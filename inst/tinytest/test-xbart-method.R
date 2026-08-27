source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "checkXvalShape.R", package = "dbarts"),
  local = TRUE
)

# test that random subsample runs correctly with valid inputs
x <- testData$x
y <- testData$y

n.reps <- 3L
n.trees <- c(5L, 7L)
k <- c(1, 2, 4)
power <- c(1.5, 2)
base <- c(0.75, 0.8, 0.95)

xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "random subsample",
  n.reps = n.reps,
  n.trees = n.trees,
  k = k,
  power = power,
  base = base,
  n.threads = 2L
)

checkXvalShape(xval, n.reps, n.trees, k, power, base)

rm(xval, base, power, k, n.trees, n.reps, y, x)


# test that k-fold runs correctly with valid inputs
x <- testData$x
y <- testData$y

n.reps <- 3L
n.trees <- c(5L, 7L)
k <- c(1, 2, 4)
power <- c(1.5, 2)
base <- c(0.75, 0.8, 0.95)

xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = n.reps,
  n.trees = n.trees,
  k = k,
  power = power,
  base = base,
  n.threads = 2L
)

checkXvalShape(xval, n.reps, n.trees, k, power, base)

rm(xval, base, power, k, n.trees, n.reps, y, x)

# test that k-fold runs correctly with one input
x <- testData$x
y <- testData$y

xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  n.reps = 3,
  n.test = 5,
  k = 2,
  n.threads = 2L
)

expect_equal(length(xval), 3L)

rm(xval, y, x)


# test that k-fold subdivides data correctly when data do not divide evenly by k
x <- testData$x[1L:24L, ]
y <- testData$y[1L:24L]

k <- c(2, 4)

xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = 3L,
  k = k,
  n.threads = 1L
)

expect_inherits(xval, "array")

# the fold sizes themselves: 24 rows over 5 folds is 5, 5, 5, 4, ... - a loss
# that reports its own held-out row count averages to 24 / 5 in every cell,
# which no fold plan that drops the remainder rows can reach
foldSizeXbartMethod <- function(y.test, y.test.hat, weights) length(y.test)
xvalFolds <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L, 1L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  k = k,
  loss = foldSizeXbartMethod,
  n.threads = 1L
)
expect_equal(as.vector(xvalFolds), rep_len(length(y) / 5, length(xvalFolds)))
rm(xvalFolds, foldSizeXbartMethod)

rm(testData)


# control = is removed: the driver runs fine without it, and passing it now
# hits R's own unused-argument error, since xbart has no dots
xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  k = c(1, 4),
  n.threads = 1L
)
expect_equal(dim(xval), c(2L, 2L))
expect_true(all(is.finite(xval)))
expect_error(
  dbarts::xbart(
    x,
    y,
    n.reps = 1L,
    n.threads = 1L,
    control = dbarts::dbartsControl()
  ),
  pattern = "unused argument"
)

rm(xval)

# k grid: the sweep order does not depend on the order k is listed in. A
# fixed seed and n.threads = 1 make the whole sweep bitwise deterministic,
# so aligning the k axis by name must give identical values either way,
# and each result's k axis must stay in the order it was given
seed <- 37L
xvalAscending <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  k = c(2, 8),
  n.trees = 5L,
  n.threads = 1L,
  seed = seed
)
xvalDescending <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  k = c(8, 2),
  n.trees = 5L,
  n.threads = 1L,
  seed = seed
)
expect_identical(dimnames(xvalAscending)$k, c("2", "8"))
expect_identical(dimnames(xvalDescending)$k, c("8", "2"))
expect_identical(xvalAscending, xvalDescending[, dimnames(xvalAscending)$k])
rm(xvalAscending, xvalDescending, seed)

# k as a hyperprior is held and drawn every sweep, never sorted, and
# contributes no k axis to the result
xvalHyperprior <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  k = dbarts:::chi(1.5, 2),
  n.trees = c(5L, 7L),
  n.threads = 1L
)
expect_true(!("k" %in% names(dimnames(xvalHyperprior))))
expect_equal(dim(xvalHyperprior), c(2L, 2L))
expect_true(all(is.finite(xvalHyperprior)))
rm(xvalHyperprior)

rm(k, y, x)
