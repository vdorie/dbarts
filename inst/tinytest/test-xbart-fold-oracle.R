# Oracles for the parts of xbart upstream of its loss call: which rows a
# fold holds out, whether a fold's fitted values line up with the y it was
# scored against, whether a fold's training set ever saw its own held-out
# rows, and whether a labelled parameter cell in the result array is the
# cell that was actually fit. The loss arithmetic itself (test-xbart-oracle.R)
# is not retested here.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# the fold split is drawn via R's sample(); pin the sampling kind so the
# hand reconstruction below tracks xbart's own draws regardless of what
# earlier test files left behind
oldSampleKind <- RNGkind()[3L]
suppressWarnings(RNGkind(sample.kind = "Rejection"))

## fold assembly. Every seed the run uses is drawn from the call's seed in
## one pass - one split seed per replication, then one seed per sampler each
## (replication, fold) unit creates (one per distinct tree count; n.trees is
## a single value below, so one per unit) - so the permutation is
## reconstructible outside xbart from (seed, n, n.reps, fold count) alone, at
## any thread count. y is set to the row index, so a capturing loss's y.test
## IS that fold's row numbers directly.
n <- 24L
seed <- 7L
set.seed(4441L)
x <- matrix(runif(n * 2L), n, 2L)
y <- as.numeric(seq_len(n))

n.test <- 4L
n.reps <- 1L
foldSizes <- rep.int(n %/% n.test, n.test) +
  rep.int(c(1L, 0L), c(n %% n.test, n.test - n %% n.test))
set.seed(seed)
seeds <- sample.int(.Machine$integer.max, n.reps + n.reps * n.test)
set.seed(seeds[1L])
permutation <- sample.int(n)
expectedFolds <- vector("list", n.test)
foldOffset <- 0L
for (fold in seq_len(n.test)) {
  expectedFolds[[fold]] <-
    sort(permutation[foldOffset + seq_len(foldSizes[fold])])
  foldOffset <- foldOffset + foldSizes[fold]
}

captured <- new.env(parent = emptyenv())
captureRows <- function(y.test, testSamples, weights) {
  captured$calls <- c(captured$calls, list(as.integer(y.test)))
  0
}
invisible(dbarts::xbart(
  x,
  y,
  n.samples = 5L,
  n.burn = c(3L, 1L),
  method = "k-fold",
  n.test = n.test,
  n.reps = n.reps,
  n.trees = 3L,
  n.threads = 1L,
  seed = seed,
  loss = captureRows
))
actualFolds <- captured$calls

# the folds arrive in fold order even though each is an independent unit of
# work now, and each holds exactly the rows the seed's own permutation gives
# it
expect_equal(actualFolds, expectedFolds)

pairs <- combn(seq_len(n.test), 2L, simplify = FALSE)
overlaps <- vapply(
  pairs,
  function(p) length(intersect(actualFolds[[p[1L]]], actualFolds[[p[2L]]])),
  0L
)
names(overlaps) <- vapply(pairs, paste, "", collapse = "-")
expect_equal(overlaps, `names<-`(integer(length(pairs)), names(overlaps)))
expect_equal(sort(unlist(actualFolds)), seq_len(n))
# the capturing loss only ever receives held-out rows, so the training set
# a fold actually fit on is not observable through it directly; disjointness
# and coverage above pin it exactly, since seq_len(n)[-testRows] over folds
# that partition 1:n is necessarily the union of every other fold. Whether
# a fold's fit ever leaked from its own held-out rows is checked below.

rm(
  n,
  seed,
  x,
  y,
  n.test,
  n.reps,
  foldSizes,
  seeds,
  permutation,
  expectedFolds,
  foldOffset,
  fold,
  captured,
  captureRows,
  actualFolds,
  pairs,
  overlaps
)

## row alignment. A fold's posterior mean on its held-out rows should
## track the y it was scored against; a misordered gather of the test
## channel decorrelates the two while leaving shapes untouched. The
## permutation null below reshuffles the observed pairing to get that
## decorrelated baseline directly, rather than assuming a fixed number.
x <- testData$x
y <- testData$y

captured <- new.env(parent = emptyenv())
captureFit <- function(y.test, testSamples, weights) {
  captured$calls <- c(
    captured$calls,
    list(list(y = y.test, fit = rowMeans(testSamples)))
  )
  0
}
invisible(dbarts::xbart(
  x,
  y,
  n.samples = 40L,
  n.burn = c(20L, 10L),
  method = "k-fold",
  n.test = 4L,
  n.reps = 1L,
  n.trees = 50L,
  n.threads = 1L,
  seed = 31L,
  loss = captureFit
))
folds <- captured$calls

observed <- vapply(folds, function(f) cor(f$y, f$fit), 0)
pooledY <- unlist(lapply(folds, `[[`, "y"))
pooledFit <- unlist(lapply(folds, `[[`, "fit"))
set.seed(1L)
null <- replicate(200L, cor(pooledY, sample(pooledFit)))
threshold <- max(abs(null)) + 0.15

for (i in seq_along(folds)) {
  expect_true(
    observed[i] > threshold,
    info = paste0(
      "fold ",
      i,
      ": observed ",
      round(observed[i], 3L),
      " vs null threshold ",
      round(threshold, 3L)
    )
  )
}

rm(captured, captureFit, folds, observed, pooledY, pooledFit, null, threshold)

## leakage. y is independent of x, so a fold that never saw its own
## held-out rows can do no better than predict near the training mean;
## its rmse must sit at or above sd(y). A fit that (wrongly) trained on
## the rows it is scored against fits noise instead, and its rmse falls
## well short of sd(y) - which is exactly what an in-sample fit at the
## same hyperparameters looks like. The folds are deliberately two rows
## wide: a single leaked row then halves the fold's mean squared error,
## which clears the run-to-run spread of the ratio, whereas in a five-row
## fold one leaked row can only move the rmse by sqrt(4/5) and the leaked
## and clean cases overlap. Five replications average over the split.
set.seed(11L)
n <- 60L
x <- matrix(runif(n * 4L), n, 4L)
y <- rnorm(n)
sdY <- sd(y)

heldOut <- dbarts::xbart(
  x,
  y,
  n.samples = 30L,
  n.burn = c(20L, 10L),
  method = "k-fold",
  n.test = 30L,
  n.reps = 5L,
  n.trees = 50L,
  k = 0.5,
  n.threads = 1L,
  seed = 1011L,
  loss = "rmse"
)
fit <- dbarts::bart(
  y ~ x,
  n.trees = 50L,
  k = 0.5,
  n.samples = 30L,
  n.burn = 20L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 2011L,
  verbose = FALSE,
  keepTrainingFits = TRUE
)
heldOutRmse <- mean(heldOut)
inSampleRmse <- sqrt(mean((y - fit$yhat.train.mean)^2))

expect_true(
  heldOutRmse >= sdY,
  info = paste(
    "held-out rmse",
    round(heldOutRmse, 3L),
    "vs sd(y)",
    round(sdY, 3L)
  )
)
expect_true(
  inSampleRmse <= 0.65 * sdY,
  info = paste(
    "in-sample rmse",
    round(inSampleRmse, 3L),
    "vs 0.65*sd(y)",
    round(0.65 * sdY, 3L)
  )
)
expect_true(
  heldOutRmse >= 1.75 * inSampleRmse,
  info = paste(
    "held-out rmse",
    round(heldOutRmse, 3L),
    "in-sample rmse",
    round(inSampleRmse, 3L)
  )
)

rm(n, x, y, sdY, heldOut, fit, heldOutRmse, inSampleRmse)

## axis placement. k = 20 shrinks every tree toward the training mean
## regardless of n.trees, so that column's rmse must sit near sd(y.test) at
## both tree counts; a small k with 50 trees can use x and lands well
## below it. n.trees and k are both length-2 grids here, so a transposed
## pair of equal-length axes would place values in the wrong cell without
## triggering an array-shape error.
x <- testData$x
y <- testData$y
sdY <- sd(y)

xval <- dbarts::xbart(
  x,
  y,
  n.samples = 40L,
  n.burn = c(20L, 10L),
  method = "k-fold",
  n.test = 5L,
  n.reps = 2L,
  n.trees = c(1L, 50L),
  k = c(20, 2),
  n.threads = 1L,
  seed = 17L
)
k20 <- colMeans(xval[,, "20"])
kSmall50 <- mean(xval[, "50", "2"])

for (treeCount in names(k20)) {
  expect_true(
    k20[treeCount] >= 0.85 * sdY && k20[treeCount] <= 1.05 * sdY,
    info = paste0(
      "n.trees=",
      treeCount,
      ": k=20 rmse ",
      round(k20[treeCount], 3L),
      " vs band [",
      round(0.85 * sdY, 3L),
      ", ",
      round(1.05 * sdY, 3L),
      "]"
    )
  )
}
expect_true(
  kSmall50 < 0.65 * sdY,
  info = paste(
    "k=2, n.trees=50 rmse",
    round(kSmall50, 3L),
    "vs 0.65*sd(y)",
    round(0.65 * sdY, 3L)
  )
)

rm(x, y, sdY, xval, k20, kSmall50, treeCount)

## hand-rebuilt cell. The oracles above show a seeded xbart call is
## thread-count invariant; they do not show its VALUES are right. Rebuild one
## cell entirely from outside: reconstruct the fold split and the unit's
## sampler seed through the documented derivation - one split seed per
## replication, then one seed per sampler a unit creates - fit that fold with
## an ordinary dbarts() sampler at the seed and the cell's hyperparameters,
## score it by hand, and check it against xbart's own reported loss.
##
## sigest is pinned explicitly on both sides rather than left to the engine's
## own linear-model fallback, so the residual prior calibrates identically
## without reproducing that fallback here. x's global min and max sit at rows
## 1 and 2, and the seed is one for which the held-out fold never draws them,
## so the training-only cut grid dbarts() builds from x[trainRows, ] matches
## xbart's shared, full-data grid exactly (both default to
## useQuantiles = FALSE, a uniform grid over each column's range) - the one
## condition under which a sampler over the training rows alone reproduces a
## fold view over the full data bit for bit.
n <- 20L
numTest <- 4L
seed <- 2L
n.trees <- 5L
n.samples <- 8L
n.burn <- c(6L, 3L)
sigest <- 1.0

set.seed(9182L)
x <- matrix(runif(n), n, 1L)
x[1L] <- 0
x[2L] <- 1
y <- 3 * x[, 1L] + rnorm(n)

cellLoss <- dbarts::xbart(
  x,
  y,
  method = "random subsample",
  n.reps = 1L,
  n.test = numTest,
  n.trees = n.trees,
  n.threads = 1L,
  seed = seed,
  n.samples = n.samples,
  n.burn = n.burn,
  sigest = sigest
)

# useUnitSeed = FALSE reproduces an off-by-one in the seed index: the split
# seed handed to the sampler instead of the unit's own
rebuildCell <- function(useUnitSeed = TRUE) {
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, 2L)
  splitSeed <- seeds[1L]
  unitSeed <- seeds[2L]
  set.seed(splitSeed)
  testRows <- sort(sample.int(n, numTest))
  trainRows <- setdiff(seq_len(n), testRows)
  sampler <- dbarts::dbarts(
    x[trainRows, , drop = FALSE],
    y[trainRows],
    test = x[testRows, , drop = FALSE],
    control = dbarts::dbartsControl(
      n.chains = 1L,
      n.threads = 1L,
      n.trees = n.trees,
      n.samples = n.samples,
      n.cuts = 100L,
      useQuantiles = FALSE,
      keepTrees = FALSE,
      keepTrainingFits = FALSE,
      updateState = FALSE,
      verbose = FALSE,
      seed = if (useUnitSeed) unitSeed else splitSeed
    ),
    sigest = sigest
  )
  samples <- sampler$run(n.burn[1L], n.samples)
  sqrt(mean((y[testRows] - rowMeans(samples$test))^2))
}

expect_equal(as.vector(cellLoss), rebuildCell())
expect_false(isTRUE(all.equal(
  as.vector(cellLoss),
  rebuildCell(useUnitSeed = FALSE)
)))

rm(
  n,
  numTest,
  seed,
  n.trees,
  n.samples,
  n.burn,
  sigest,
  x,
  y,
  cellLoss,
  rebuildCell
)

suppressWarnings(RNGkind(sample.kind = oldSampleKind))
rm(oldSampleKind, testData)


# The default sigest is estimated per fold from the fold's training rows: the
# fold that holds a row out fits exactly the same whatever that row's
# response, while a fixed sigest stays one value for every fold.
foldDraws <- function(y, ...) {
  record <- new.env()
  record$draws <- list()
  recordLoss <- function(y.test, testSamples, weights) {
    record$draws[[length(record$draws) + 1L]] <- list(
      y = y.test,
      draws = testSamples
    )
    0
  }
  set.seed(5)
  x <- matrix(runif(60L), 30L, 2L)
  xbart(
    x,
    y,
    n.samples = 5L,
    n.burn = c(5L, 3L),
    method = "k-fold",
    n.test = 3L,
    n.reps = 1L,
    n.trees = 5L,
    k = 2,
    loss = recordLoss,
    seed = 2L,
    verbose = FALSE,
    n.threads = 1L,
    ...
  )
  record$draws
}
set.seed(6)
foldY <- rnorm(30L)
changedY <- foldY
changedY[1L] <- foldY[1L] + 50
original <- foldDraws(foldY)
changed <- foldDraws(changedY)
holdsRowOne <- which(vapply(
  changed,
  function(fold) any(fold$y == changedY[1L]),
  FALSE
))
expect_equal(length(holdsRowOne), 1L)
expect_identical(changed[[holdsRowOne]]$draws, original[[holdsRowOne]]$draws)
# a fold that trains on the changed row does move
trainsOnRowOne <- setdiff(seq_along(changed), holdsRowOne)[1L]
expect_false(identical(
  changed[[trainsOnRowOne]]$draws,
  original[[trainsOnRowOne]]$draws
))
rm(foldDraws, foldY, changedY, original, changed, holdsRowOne, trainsOnRowOne)
