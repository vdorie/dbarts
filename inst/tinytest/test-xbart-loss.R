source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)
source(
  system.file("common", "checkXvalShape.R", package = "dbarts"),
  local = TRUE
)

# test that works with custom loss
x <- testData$x
y <- testData$y

n.reps <- 3L
n.trees <- c(5L, 7L)
k <- c(1, 2, 4)
power <- c(1.5, 2)
base <- c(0.75, 0.8, 0.95)

mad <- function(y.train, y.train.hat, weights) {
  mean(abs(y.train - apply(y.train.hat, 1L, mean)))
}

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
  loss = mad,
  n.threads = 1L
)

checkXvalShape(xval, n.reps, n.trees, k, power, base)


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
  loss = mad,
  n.threads = 2L
)

checkXvalShape(xval, n.reps, n.trees, k, power, base)

# the reported cell IS the loss function's value, averaged over the folds
# rather than summed: a loss that always returns 3 makes every cell 3
constantLoss <- function(y.test, y.test.hat, weights) 3.0
xval <- dbarts::xbart(
  x,
  y,
  n.samples = 6L,
  n.burn = c(5L, 3L),
  method = "k-fold",
  n.test = 5,
  n.reps = 2L,
  n.trees = n.trees,
  loss = constantLoss,
  n.threads = 1L
)
expect_equal(as.vector(xval), rep_len(3.0, length(xval)))

# a warning raised inside a fit reaches the caller the same way at any
# n.threads, a worker's included: each distinct (class, message) once
warningLoss <- function(y.test, y.test.hat, weights) {
  warning(warningCondition("loss warned", class = "xbartTestWarning"))
  sqrt(mean((y.test - rowMeans(y.test.hat))^2))
}
runWarned <- function(n.threads) {
  warned <- character()
  loss <- withCallingHandlers(
    dbarts::xbart(
      x,
      y,
      n.samples = 6L,
      n.burn = c(5L, 3L),
      n.test = 3,
      n.reps = 2L,
      n.trees = 5L,
      loss = warningLoss,
      n.threads = n.threads,
      seed = 1L
    ),
    warning = function(w) {
      warned <<- c(warned, paste(class(w)[1L], conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  list(loss = loss, warned = warned)
}
warned.1 <- runWarned(1L)
warned.2 <- runWarned(2L)
expect_identical(warned.1$warned, "xbartTestWarning loss warned")
expect_identical(warned.2$warned, warned.1$warned)
expect_identical(warned.2$loss, warned.1$loss)

# in this process, an error still propagates with the warnings raised before
# it, and options(warn = 2) aborts at the first warning. xbart runs a loss
# function in the caller's environment, so the counter lives here.
lossCalls <- 0L
failingLoss <- function(y.test, y.test.hat, weights) {
  lossCalls <<- lossCalls + 1L
  if (lossCalls > 1L) {
    stop("loss failed")
  }
  warning(warningCondition("loss warned", class = "xbartTestWarning"))
  1
}
fitFailing <- function() {
  dbarts::xbart(
    x,
    y,
    n.samples = 6L,
    n.burn = c(5L, 3L),
    n.test = 3,
    n.reps = 1L,
    n.trees = 5L,
    loss = failingLoss,
    n.threads = 1L,
    seed = 1L
  )
}
runFailing <- function() {
  warned <- character()
  failed <- tryCatch(
    withCallingHandlers(
      fitFailing(),
      warning = function(w) {
        warned <<- c(warned, paste(class(w)[1L], conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = conditionMessage
  )
  list(failed = failed, warned = warned)
}
failing <- runFailing()
expect_identical(failing$failed, "loss failed")
expect_identical(failing$warned, "xbartTestWarning loss warned")

lossCalls <- -10L
oldOptions <- options(warn = 2L)
failed <- tryCatch(fitFailing(), error = conditionMessage)
options(oldOptions)
expect_true(grepl("loss warned", failed, fixed = TRUE))
expect_identical(lossCalls, -9L)

rm(oldOptions, failed, failing, runFailing, fitFailing, failingLoss, lossCalls)
rm(warned.2, warned.1, runWarned, warningLoss)
rm(xval, constantLoss, mad, base, power, k, n.trees, n.reps, y, x)


# A custom loss keeps its own environment: a closure from a factory reads what
# it captured, not a same-named variable in the caller's frame.
makeScaledLoss <- function(multiplier) {
  function(y.test, testSamples, weights) {
    multiplier * sqrt(mean((y.test - rowMeans(testSamples))^2))
  }
}
xbartArgs <- list(
  testData$x,
  testData$y,
  n.samples = 5L,
  n.burn = c(5L, 3L),
  n.reps = 1L,
  n.trees = 5L,
  k = 2,
  seed = 1L,
  verbose = FALSE,
  n.threads = 1L
)
plainLoss <- do.call(xbart, xbartArgs)
scaledLoss <- do.call(xbart, c(xbartArgs, loss = makeScaledLoss(10)))
expect_equal(scaledLoss, 10 * plainLoss)
shadowedCall <- function() {
  multiplier <- 1
  do.call(xbart, c(xbartArgs, loss = makeScaledLoss(10)))
}
expect_equal(shadowedCall(), 10 * plainLoss)
# the list form calls it from the given environment, closure intact
listLoss <- do.call(
  xbart,
  c(xbartArgs, loss = list(list(makeScaledLoss(10), globalenv())))
)
expect_equal(listLoss, 10 * plainLoss)
rm(makeScaledLoss, xbartArgs, plainLoss, scaledLoss, shadowedCall, listLoss)

# On a binary response the built-in rmse compares the response with the
# posterior mean probability, the link applied per draw and averaged last: the
# hand computation from the same fold's test draws (read through a supplied loss
# that records them, at the same seed) must match, for both links, with and
# without weights. The raw-latent value is a different number.
set.seed(31L)
nBin <- 40L
xBin <- matrix(runif(nBin * 2L), nBin, 2L)
yBin <- as.numeric(2 * xBin[, 1L] - 1 + rnorm(nBin, 0, 0.5) > 0)
wBin <- rep(c(1, 5), length.out = nBin)
binaryCaptured <- new.env(parent = emptyenv())
binaryCapture <- function(y.test, testSamples, weights) {
  binaryCaptured$calls <- c(
    binaryCaptured$calls,
    list(list(y = y.test, samples = testSamples, weights = weights))
  )
  0
}
binaryRun <- function(...) {
  dbarts::xbart(
    xBin,
    yBin,
    n.samples = 6L,
    n.burn = c(4L, 2L),
    method = "k-fold",
    n.test = 4,
    n.reps = 1L,
    n.trees = 5L,
    k = 2,
    seed = 5L,
    verbose = FALSE,
    n.threads = 1L,
    ...
  )
}
for (family in c("probit", "logistic")) {
  link <- if (family == "probit") pnorm else plogis
  # xbart takes weights on a binary fit only under the logistic family
  for (weighted in if (family == "probit") FALSE else c(FALSE, TRUE)) {
    args <- list(family = family)
    if (weighted) {
      args$weights <- wBin
    }
    binaryCaptured$calls <- list()
    do.call(binaryRun, c(args, loss = binaryCapture))
    folds <- binaryCaptured$calls
    expect_equal(length(folds), 4L)
    byFold <- function(transform) {
      vapply(
        folds,
        function(f) {
          resid2 <- (f$y - rowMeans(transform(f$samples)))^2
          if (is.null(f$weights)) {
            sqrt(mean(resid2))
          } else {
            sqrt(sum(f$weights * resid2) / sum(f$weights))
          }
        },
        0
      )
    }
    reported <- as.vector(do.call(binaryRun, c(args, loss = "rmse")))
    expect_equal(reported, mean(byFold(link)))
    # the mutant: rmse on the raw latent draws
    expect_false(isTRUE(all.equal(reported, mean(byFold(identity)))))
  }
}
rm(family, link, weighted, args, folds, byFold, reported)
rm(binaryRun, binaryCapture, binaryCaptured, wBin, yBin, xBin, nBin)

# A binomial response is scored as its rows written out one trial each: the
# built-in log and mcr losses on the binomial rows equal the same losses on the
# 0/1 expansion of each fold's held-out rows, each trial carrying its row's
# draws and weight, to rounding. A supplied loss is handed each row's observed
# proportion and its trials times its weight. The draws are the fold's own,
# recorded through a supplied loss at the seed the built-in runs take.
set.seed(37L)
nBnm <- 30L
bnmData <- data.frame(x1 = runif(nBnm), x2 = runif(nBnm))
bnmTrials <- rep_len(c(0L, 1L, 3L, 4L, 2L), nBnm)
bnmData$k <- rbinom(nBnm, bnmTrials, plogis(2 * bnmData$x1 - 1))
bnmData$f <- bnmTrials - bnmData$k
bnmWeights <- rep_len(c(1, 2, 1), nBnm)
bnmCaptured <- new.env(parent = emptyenv())
bnmCapture <- function(y.test, testSamples, weights) {
  bnmCaptured$calls <- c(
    bnmCaptured$calls,
    list(list(y = y.test, samples = testSamples, weights = weights))
  )
  0
}
bnmRun <- function(...) {
  suppressWarnings(dbarts::xbart(
    cbind(k, f) ~ x1 + x2,
    bnmData,
    family = binomial,
    weights = bnmWeights,
    n.samples = 6L,
    n.burn = c(4L, 2L),
    method = "k-fold",
    n.test = 3,
    n.reps = 1L,
    n.trees = 5L,
    k = 2,
    seed = 7L,
    verbose = FALSE,
    n.threads = 1L,
    ...
  ))
}
bnmCaptured$calls <- list()
bnmRun(loss = bnmCapture)
bnmFolds <- bnmCaptured$calls
expect_equal(length(bnmFolds), 3L)
for (fold in bnmFolds) {
  rows <- as.integer(names(fold$y))
  expect_equal(
    unname(fold$y),
    ifelse(bnmTrials[rows] > 0L, bnmData$k[rows] / bnmTrials[rows], 0)
  )
  expect_equal(fold$weights, bnmTrials[rows] * bnmWeights[rows])
}
bnmBinaryControl <- dbarts::dbartsControl()
bnmBinaryControl@binary <- TRUE
# the fold's rows written out: row i as k_i trials of 1 and n_i - k_i of 0,
# each with row i's draws and weight
bnmExpanded <- function(lossName, fold) {
  rows <- as.integer(names(fold$y))
  trials <- bnmTrials[rows]
  each <- rep(seq_along(rows), trials)
  y <- unlist(lapply(seq_along(rows), function(i) {
    rep(c(1, 0), c(bnmData$k[rows[i]], bnmData$f[rows[i]]))
  }))
  dbarts:::xbartLossFunction(lossName, bnmBinaryControl, "logistic")(
    y,
    fold$samples[each, , drop = FALSE],
    bnmWeights[rows][each]
  )
}
for (lossName in c("log", "mcr")) {
  reported <- as.vector(bnmRun(loss = lossName))
  expanded <- mean(vapply(bnmFolds, bnmExpanded, 0, lossName = lossName))
  expect_equal(reported, expanded, info = lossName)
}
# and each equal to its definition written by hand on the binomial rows, the
# fold probabilities on both sides of 0.5 so that mcr's direction is read
bnmProbabilities <- lapply(bnmFolds, function(fold) {
  rowMeans(plogis(fold$samples))
})
expect_true(any(unlist(bnmProbabilities) > 0.5))
expect_true(any(unlist(bnmProbabilities) < 0.5))
bnmByHand <- function(perRow) {
  mean(vapply(
    seq_along(bnmFolds),
    function(i) {
      fold <- bnmFolds[[i]]
      sum(fold$weights * perRow(fold$y, bnmProbabilities[[i]])) /
        sum(fold$weights)
    },
    0
  ))
}
expect_equal(
  as.vector(bnmRun(loss = "log")),
  bnmByHand(function(y, p) -(y * log(p) + (1 - y) * log(1 - p)))
)
expect_equal(
  as.vector(bnmRun(loss = "mcr")),
  bnmByHand(function(y, p) ifelse(p > 0.5, 1 - y, y))
)
rm(bnmProbabilities, bnmByHand)
rm(fold, rows, lossName, reported, expanded, bnmExpanded, bnmBinaryControl)
rm(bnmFolds, bnmRun, bnmCapture, bnmCaptured, bnmWeights, bnmTrials, bnmData)
rm(nBnm)

rm(testData)
