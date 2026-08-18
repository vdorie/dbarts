# Oracles for xbart's loss plumbing: every value asserted here is derived
# OUTSIDE xbart. Two arms run at the same seed - one under the loss being
# checked, one under a capturing loss function that records the (y.test,
# testSamples, weights) triple xbart hands its loss and returns a constant -
# and the built-in formula is then transcribed by hand from the captured
# triples and averaged over folds here. A loss draws nothing, so the two arms
# see the same splits and the same draws; no assertion below takes its
# expected value from xbart's own arithmetic.

set.seed(2718L)
n <- 12L
x <- matrix(runif(n * 2L), n, 2L)
y <- 3 * x[, 1L] + rnorm(n, 0, 0.5)
y.bin <- as.numeric(x[, 1L] + rnorm(n, 0, 0.5) > 0.5)
w <- rep(c(1, 5), length.out = n)

# 12 rows in 3 folds hold out 4 apiece; 2 replications over a single grid
# cell, so a run reports one fold-averaged loss per replication.
xval <- function(y, n.test = 3, ...) {
  dbarts::xbart(
    x,
    y,
    n.samples = 5L,
    n.burn = c(3L, 1L),
    method = "k-fold",
    n.test = n.test,
    n.reps = 2L,
    n.trees = 3L,
    n.threads = 1L,
    seed = 7L,
    ...
  )
}

captured <- new.env(parent = emptyenv())
captureLoss <- function(y.test, testSamples, weights) {
  captured$calls <- c(
    captured$calls,
    list(list(y = y.test, samples = testSamples, weights = weights))
  )
  0
}
foldsUnder <- function(y, ...) {
  captured$calls <- list()
  captured$reported <- xval(y, ..., loss = captureLoss)
  captured$calls
}
# replication-major, folds within: the first 3 calls are replication 1's
foldAverage <- function(values) c(mean(values[1L:3L]), mean(values[4L:6L]))

# (a) continuous, built-in rmse. The captured shapes pin the plumbing the
# formula reads: one row per held-out observation, one column per sample,
# no weights, and a replication's folds partitioning the data.
folds <- foldsUnder(y)
expect_equal(captured$reported, c(0, 0))
expect_equal(length(folds), 6L)
expect_equal(
  vapply(folds, function(f) dim(f$samples), integer(2L)),
  matrix(c(4L, 5L), 2L, 6L)
)
expect_true(all(vapply(folds, function(f) is.null(f$weights), NA)))
expect_equal(sort(unlist(lapply(folds[1L:3L], `[[`, "y"))), sort(y))
expect_equal(
  xval(y, loss = "rmse"),
  foldAverage(vapply(
    folds,
    function(f) sqrt(mean((f$y - rowMeans(f$samples))^2)),
    0
  ))
)

# (b) the weighted rmse branch, and the weights are not decoration: a
# non-unit vector moves the reported loss, a unit one leaves it alone
foldsWeighted <- foldsUnder(y, weights = w)
expect_equal(
  sort(unlist(lapply(foldsWeighted[1L:3L], `[[`, "weights"))),
  sort(w)
)
expect_equal(
  xval(y, weights = w, loss = "rmse"),
  foldAverage(vapply(
    foldsWeighted,
    function(f) {
      sqrt(sum(f$weights * (f$y - rowMeans(f$samples))^2) / sum(f$weights))
    },
    0
  ))
)
expect_false(isTRUE(all.equal(
  xval(y, weights = w, loss = "rmse"),
  xval(y, loss = "rmse")
)))
# a unit weight vector is the same fit, but not bitwise: the weighted
# sufficient statistics and the weighted loss re-associate the same sums, so
# the two arms part in the last ulp (measured 4e-16 on the draws, 7e-16
# relative on the reported loss) rather than agreeing exactly. The fold
# partition itself is untouched, which is exact.
expect_equal(
  xval(y, weights = rep(1, n), loss = "rmse"),
  xval(y, loss = "rmse"),
  tolerance = 1e-10
)
expect_identical(
  lapply(foldsUnder(y, weights = rep(1, n)), `[[`, "y"),
  lapply(folds, `[[`, "y")
)

# (c) the binary losses, both transcribed from one capture arm: the samples
# are latent, so the oracle applies the probit link itself
foldsBinary <- foldsUnder(y.bin)
probs <- lapply(foldsBinary, function(f) rowMeans(pnorm(f$samples)))
expect_equal(
  xval(y.bin, loss = "log"),
  foldAverage(vapply(
    seq_along(foldsBinary),
    function(i) {
      -mean(ifelse(
        foldsBinary[[i]]$y > 0,
        log(probs[[i]]),
        log1p(-probs[[i]])
      ))
    },
    0
  ))
)
expect_equal(
  xval(y.bin, loss = "mcr"),
  foldAverage(vapply(
    seq_along(foldsBinary),
    function(i) mean((probs[[i]] > 0.5) != foldsBinary[[i]]$y),
    0
  ))
)
expect_identical(xval(y.bin), xval(y.bin, loss = "log"))

# (d) the fold averaging on its own, by hand rather than by capture: a loss
# that reports its fold's size averages the fold sizes, and 12 rows split
# 4, 4, 4 in 3 folds and 3, 3, 2, 2, 2 in 5
foldSize <- function(y.test, testSamples, weights) length(y.test)
expect_equal(xval(y, loss = foldSize), c(4, 4))
expect_equal(xval(y, n.test = 5, loss = foldSize), c(2.4, 2.4))
