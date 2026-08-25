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

## fold assembly. At n.threads = 1, the only draws before replication
## 1's row permutation are the chunk seed itself, so the permutation is
## reconstructible outside xbart from (seed, n) alone. y is set to the row
## index, so a capturing loss's y.test IS that fold's row numbers directly.
n <- 24L
seed <- 7L
set.seed(4441L)
x <- matrix(runif(n * 2L), n, 2L)
y <- as.numeric(seq_len(n))

n.test <- 4L
foldSizes <- rep.int(n %/% n.test, n.test) +
  rep.int(c(1L, 0L), c(n %% n.test, n.test - n %% n.test))
set.seed(seed)
chunkSeed <- sample.int(.Machine$integer.max, 1L)
set.seed(chunkSeed)
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
  n.reps = 1L,
  n.trees = 3L,
  n.threads = 1L,
  seed = seed,
  loss = captureRows
))
actualFolds <- captured$calls

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
  foldSizes,
  chunkSeed,
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
fit <- dbarts::bart2(
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

suppressWarnings(RNGkind(sample.kind = oldSampleKind))
rm(oldSampleKind, testData)
