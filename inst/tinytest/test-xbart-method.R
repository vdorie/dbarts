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
  n.burn = c(5L, 3L),
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
  n.burn = c(5L, 3L),
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
  n.burn = c(5L, 3L),
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
  n.burn = c(5L, 3L),
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
  n.burn = c(5L, 3L),
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


# control = is removed: the driver runs fine without it, and the name is
# refused by a tombstone naming the flat arguments its settings became
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
  pattern = "'control' has left 'xbart'"
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

# a bare k hyperprior is one modelled cell: k is drawn every sweep inside it,
# so there is nothing to sort, and the one-cell axis drops under drop = TRUE
xvalHyperprior <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  k = chi(1.5, 2),
  n.trees = c(5L, 7L),
  n.threads = 1L
)
expect_true(!("k" %in% names(dimnames(xvalHyperprior))))
expect_equal(dim(xvalHyperprior), c(2L, 2L))
expect_true(all(is.finite(xvalHyperprior)))
rm(xvalHyperprior)

# the k grid is a LIST when it mixes fixed values with modelled ones. Every
# cell is fit: the modelled cell is labelled by the call that rebuilds it and
# differs from both fixed cells, and the two fixed cells reproduce the same
# call made with the numeric grid alone, since a fixed cell's own sweep
# position is unchanged by a modelled one sorting after it.
xvalMixed <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  k = list(2, 8, chi(1.5, 2)),
  n.trees = 5L,
  n.threads = 1L,
  seed = 61L
)
expect_identical(dimnames(xvalMixed)$k, c("2", "8", "chi(1.5, 2)"))
expect_true(all(is.finite(xvalMixed)))
expect_true(all(xvalMixed[, "chi(1.5, 2)"] != xvalMixed[, "2"]))
expect_true(all(xvalMixed[, "chi(1.5, 2)"] != xvalMixed[, "8"]))
xvalFixedOnly <- dbarts::xbart(
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
  seed = 61L
)
expect_equal(xvalMixed[, c("2", "8")], xvalFixedOnly)
rm(xvalMixed, xvalFixedOnly)

# the chi() the calls above hold is the prior vocabulary's own constructor,
# resolved inside the 'k' argument; the name is not exported and is not
# reachable from here
expect_error(chi(1.5, 2), pattern = "could not find function")

rm(k, y, x)
