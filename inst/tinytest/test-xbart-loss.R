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

rm(testData)
