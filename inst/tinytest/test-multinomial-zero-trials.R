# A multinomial count row with zero trials: accepted at creation and through
# $setCounts, inert in the likelihood (the coupling composes it into its
# global active-row mask), still fitted, and warned about once per session.
#
# The leading oracle is a sampler whose empty rows instead carry other counts
# and are masked through $setActiveRows: the mask's own arms pin that path, so
# agreement there says a zero-trial row IS an inactive row. The fresh-x arm
# then ties the rows to the fit without them.

zeroTrialsKey <- "multinomialZeroTrials"
resetZeroTrialsKey <- function() {
  env <- dbarts:::onceWarnState
  env[[zeroTrialsKey]] <- NULL
}
# evaluates expr, returning its value and the zero-trial warnings it raised;
# every such warning is muffled, any other one propagates
countZeroTrialsWarnings <- function(expr) {
  numWarnings <- 0L
  value <- withCallingHandlers(
    expr,
    dbartsZeroTrialsWarning = function(w) {
      numWarnings <<- numWarnings + 1L
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, numWarnings = numWarnings)
}

set.seed(2609)
n <- 80L
numEmpty <- 30L
K <- 3L
x <- matrix(runif(n * 2L), n, 2L)
eta <- cbind(2 * (x[, 1L] - 0.5), x[, 2L] - 0.5, 0)
probs <- exp(eta) / rowSums(exp(eta))
counts <- t(vapply(
  seq_len(n),
  function(i) as.vector(rmultinom(1L, 1L + i %% 3L, probs[i, ])),
  integer(K)
))
storage.mode(counts) <- "integer"
# empty rows at fresh x strictly inside the data's range, so the default
# uniform cut grid is the same with or without them
xEmpty <- cbind(
  runif(numEmpty, min(x[, 1L]) + 0.01, max(x[, 1L]) - 0.01),
  runif(numEmpty, min(x[, 2L]) + 0.01, max(x[, 2L]) - 0.01)
)
xAll <- rbind(x, xEmpty)
dataRows <- seq_len(n)
emptyRows <- n + seq_len(numEmpty)
countsEmpty <- rbind(counts, matrix(0L, numEmpty, K))
# the same rows carrying counts, for the masked oracle
countsFilled <- rbind(counts, counts[seq_len(numEmpty), , drop = FALSE])
maskEmpty <- c(rep(1, n), rep(0, numEmpty))

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  updateState = FALSE,
  verbose = FALSE
)
build <- function(xs, cs, seed, mask = NULL) {
  set.seed(seed)
  sampler <- suppressWarnings(
    dbarts(xs, cs, family = "multinomial", control = control),
    classes = "dbartsZeroTrialsWarning"
  )
  if (!is.null(mask)) {
    sampler$setActiveRows(mask)
  }
  sampler
}
isSimplex <- function(p) {
  all(is.finite(p)) &&
    all(p >= 0) &&
    all(abs(apply(p, c(1L, 3L), sum) - 1) < 1e-12)
}

# --- the warning: once per session, only where a sampler takes the rows ------
resetZeroTrialsKey()
expect_silent(dbartsData(xAll, counts = countsEmpty))
created <- countZeroTrialsWarnings(
  dbarts(xAll, countsEmpty, family = "multinomial", control = control)
)
expect_identical(created$numWarnings, 1L)
expect_true(isTRUE(dbarts:::onceWarnState[[zeroTrialsKey]]))
# the class and the counts the message reports
resetZeroTrialsKey()
caught <- tryCatch(
  dbarts(xAll, countsEmpty, family = "multinomial", control = control),
  dbartsZeroTrialsWarning = function(w) w
)
expect_inherits(caught, c("dbartsZeroTrialsWarning", "dbartsWarning"))
expect_true(grepl(
  sprintf("zero trials (%d of %d)", numEmpty, n + numEmpty),
  conditionMessage(caught),
  fixed = TRUE
))
# a second creation is silent
second <- countZeroTrialsWarnings(
  dbarts(xAll, countsEmpty, family = "multinomial", control = control)
)
expect_identical(second$numWarnings, 0L)
# a creation refused after the counts validated does not spend the key
resetZeroTrialsKey()
expect_error(
  dbarts(
    xAll,
    countsEmpty,
    family = "multinomial",
    control = control,
    tree.prior = dart()
  ),
  "DART"
)
expect_null(dbarts:::onceWarnState[[zeroTrialsKey]])
# bart() reaches the same creation site
fitWarnings <- countZeroTrialsWarnings(
  bart(
    xAll,
    countsEmpty,
    family = "multinomial",
    n.trees = 20L,
    n.samples = 30L,
    n.burn = 20L,
    verbose = FALSE
  )
)
expect_identical(fitWarnings$numWarnings, 1L)
fit <- fitWarnings$value
# $setCounts after a creation warning is silent, and warns after a reset
sampler <- build(x, counts, 11L)
swapped <- countZeroTrialsWarnings(sampler$setCounts(countsEmpty[dataRows, ]))
expect_identical(swapped$numWarnings, 0L)
sampler <- build(xAll, countsFilled, 11L)
resetZeroTrialsKey()
swapped <- countZeroTrialsWarnings(sampler$setCounts(countsEmpty))
expect_identical(swapped$numWarnings, 1L)
expect_identical(sampler$data@counts, countsEmpty)
expect_identical(sampler$data@y, as.double(rowSums(countsEmpty)))

# --- empty rows are masked rows ------------------------------------------------
empty <- build(xAll, countsEmpty, 12L)
masked <- build(xAll, countsFilled, 12L, maskEmpty)
emptyRun <- empty$run(30L, 20L)
maskedRun <- masked$run(30L, 20L)
expect_identical(emptyRun$train[dataRows, , ], maskedRun$train[dataRows, , ])
expect_equal(
  emptyRun$train[emptyRows, , ],
  maskedRun$train[emptyRows, , ],
  tolerance = 1e-12
)
expect_true(isSimplex(emptyRun$train))
# composing a mask over rows that are already empty changes nothing
both <- build(xAll, countsEmpty, 12L, maskEmpty)
expect_identical(both$run(30L, 20L)$train, emptyRun$train)

# the veto: $sampleTreesFromPrior is the one reader of the veto precisions.
# The empty rows sit apart from the data rows, so a prior split between them
# leaves a leaf of only empty rows, which the veto must count absent
xApart <- rbind(x, cbind(runif(numEmpty, 1.2, 1.5), runif(numEmpty)))
empty <- build(xApart, countsEmpty, 13L)
masked <- build(xApart, countsFilled, 13L, maskEmpty)
empty$sampleTreesFromPrior()
masked$sampleTreesFromPrior()
expect_identical(
  empty$run(0L, 10L)$train[dataRows, , ],
  masked$run(0L, 10L)$train[dataRows, , ]
)

# --- against the fit without the rows --------------------------------------
# many more empty rows than data rows, so leaves holding only empty rows arise
numMany <- 300L
xMany <- rbind(
  x,
  cbind(
    runif(numMany, min(x[, 1L]) + 0.01, max(x[, 1L]) - 0.01),
    runif(numMany, min(x[, 2L]) + 0.01, max(x[, 2L]) - 0.01)
  )
)
countsMany <- rbind(counts, matrix(0L, numMany, K))
withRows <- build(xMany, countsMany, 14L)$run(30L, 20L)
withoutRows <- build(x, counts, 14L)$run(30L, 20L)
expect_equal(
  withRows$train[dataRows, , ],
  withoutRows$train,
  tolerance = 1e-12
)
# non-vacuity: the same rows with counts are a different posterior
filledRows <- build(
  xMany,
  rbind(counts, countsMany[n + seq_len(numMany), ] + 1L),
  14L
)
expect_false(isTRUE(all.equal(
  filledRows$run(30L, 20L)$train[dataRows, , ],
  withoutRows$train
)))

# --- mid-run: $setCounts emptying rows is $setActiveRows on them -------------
emptied <- build(xAll, countsFilled, 15L)
masked <- build(xAll, countsFilled, 15L)
invisible(emptied$run(20L, 5L))
invisible(masked$run(20L, 5L))
countZeroTrialsWarnings(emptied$setCounts(countsEmpty))
masked$setActiveRows(maskEmpty)
emptiedRun <- emptied$run(0L, 10L)
maskedRun <- masked$run(0L, 10L)
expect_identical(emptiedRun$train[dataRows, , ], maskedRun$train[dataRows, , ])
expect_equal(
  emptiedRun$train[emptyRows, , ],
  maskedRun$train[emptyRows, , ],
  tolerance = 1e-12
)
# restoring the counts is clearing the mask
emptied$setCounts(countsFilled)
masked$setActiveRows(NULL)
expect_equal(
  emptied$run(0L, 10L)$train,
  masked$run(0L, 10L)$train,
  tolerance = 1e-10
)
# a mask installed over empty rows and then cleared leaves them out
cleared <- build(xAll, countsEmpty, 16L, maskEmpty)
cleared$setActiveRows(NULL)
never <- build(xAll, countsEmpty, 16L)
expect_identical(cleared$run(20L, 10L)$train, never$run(20L, 10L)$train)

# --- readers on a bart() fit ------------------------------------------------
fitted <- fitted(fit)
expect_true(all(is.finite(fitted[emptyRows, ])))
expect_equal(rowSums(fitted[emptyRows, ]), rep(1, numEmpty), tolerance = 1e-12)
fitResiduals <- residuals(fit)
# glm's response residual at a zero-weight row: observed 0 minus the fit
expect_identical(fitResiduals[emptyRows, ], -fitted[emptyRows, ])
expect_true(all(is.finite(fitResiduals)))
loglik <- extract(fit, type = "loglik")
expect_true(all(loglik[, emptyRows] == 0))
expect_true(all(is.finite(loglik[, dataRows])))
expect_true(all(loglik[, dataRows] < 0))
ppd <- extract(fit, type = "ppd")
expect_true(all(ppd %in% seq_len(K)))
expect_silent(summary(fit))
# the single-trial panel, whose observed category an empty row does not have
oneHot <- matrix(0L, n, K)
oneHot[cbind(dataRows, max.col(counts, "first"))] <- 1L
fitOneHot <- suppressWarnings(
  bart(
    xAll,
    rbind(oneHot, matrix(0L, numEmpty, K)),
    family = "multinomial",
    n.chains = 1L,
    n.trees = 20L,
    n.samples = 10L,
    n.burn = 10L,
    verbose = FALSE
  ),
  classes = "dbartsZeroTrialsWarning"
)
pdf(NULL)
expect_silent(plot(fit))
expect_silent(plot(fitOneHot))
dev.off()

# --- all-zero counts ----------------------------------------------------------
allZero <- matrix(0L, n, K)
resetZeroTrialsKey()
created <- countZeroTrialsWarnings(
  dbarts(x, allZero, family = "multinomial", control = control)
)
expect_identical(created$numWarnings, 1L)
allZeroRun <- created$value$run(10L, 10L)
expect_true(isSimplex(allZeroRun$train))
sampler <- build(x, counts, 17L)
invisible(sampler$run(10L, 5L))
resetZeroTrialsKey()
swapped <- countZeroTrialsWarnings(sampler$setCounts(allZero))
expect_identical(swapped$numWarnings, 1L)
expect_true(isSimplex(sampler$run(0L, 10L)$train))
# plot has no observed panel to draw and says so, without warning
fitAllZero <- suppressWarnings(
  bart(
    x,
    allZero,
    family = "multinomial",
    n.chains = 1L,
    n.trees = 20L,
    n.samples = 10L,
    n.burn = 10L,
    verbose = FALSE
  ),
  classes = "dbartsZeroTrialsWarning"
)
plotWarnings <- 0L
pdf(NULL)
withCallingHandlers(
  plot(fitAllZero),
  warning = function(w) {
    plotWarnings <<- plotWarnings + 1L
    invokeRestart("muffleWarning")
  }
)
dev.off()
expect_identical(plotWarnings, 0L)

# --- loglik where a probability underflows to 0 -----------------------------
# a zero-count cell contributes 0, as in dmultinom, so the empty rows are
# exactly 0 and the data rows are dmultinom's value, never 0 * log(0) = NaN
countsNoFirst <- countsEmpty
countsNoFirst[, 2L] <- countsNoFirst[, 2L] + countsNoFirst[, 1L]
countsNoFirst[, 1L] <- 0L
extremeOffset <- matrix(0, n + numEmpty, K)
extremeOffset[, 1L] <- -1000
fitExtreme <- suppressWarnings(
  bart(
    xAll,
    countsNoFirst,
    family = "multinomial",
    offset = extremeOffset,
    n.chains = 1L,
    n.trees = 20L,
    n.samples = 10L,
    n.burn = 10L,
    verbose = FALSE
  ),
  classes = "dbartsZeroTrialsWarning"
)
extremeProbs <- extract(fitExtreme, type = "ev")
extremeLoglik <- extract(fitExtreme, type = "loglik")
expect_true(all(extremeProbs[,, 1L] == 0))
expect_true(all(extremeLoglik[, emptyRows] == 0))
expect_true(all(is.finite(extremeLoglik[, dataRows])))
expect_equal(
  unname(extremeLoglik[, dataRows]),
  t(vapply(
    seq_len(dim(extremeProbs)[1L]),
    function(s) {
      vapply(
        dataRows,
        function(i) {
          dmultinom(countsNoFirst[i, ], prob = extremeProbs[s, i, ], log = TRUE)
        },
        0
      )
    },
    numeric(n)
  )),
  tolerance = 1e-12
)

# --- re-creation keeps the rows inert and does not warn -------------------
# a re-created sampler is compared with a re-created one: re-creation is not
# bitwise the original even without empty rows
controlState <- control
controlState@updateState <- TRUE
buildState <- function(cs, seed, mask = NULL) {
  set.seed(seed)
  sampler <- suppressWarnings(
    dbarts(xAll, cs, family = "multinomial", control = controlState),
    classes = "dbartsZeroTrialsWarning"
  )
  if (!is.null(mask)) {
    sampler$setActiveRows(mask)
  }
  invisible(sampler$run(20L, 5L))
  sampler
}
empty <- buildState(countsEmpty, 18L)
masked <- buildState(countsFilled, 18L, maskEmpty)
# a copy is not a creation: silent even with the key unspent, and it leaves
# the key unspent
resetZeroTrialsKey()
copies <- countZeroTrialsWarnings(
  list(empty$copy(), masked$copy(), empty$copy(shallow = TRUE))
)
expect_identical(copies$numWarnings, 0L)
expect_null(dbarts:::onceWarnState[[zeroTrialsKey]])
set.seed(19L)
emptyCopyRun <- copies$value[[1L]]$run(0L, 10L)
set.seed(19L)
maskedCopyRun <- copies$value[[2L]]$run(0L, 10L)
expect_identical(
  emptyCopyRun$train[dataRows, , ],
  maskedCopyRun$train[dataRows, , ]
)

emptyFile <- tempfile(fileext = ".rds")
maskedFile <- tempfile(fileext = ".rds")
saveRDS(empty, emptyFile)
saveRDS(masked, maskedFile)
resetZeroTrialsKey()
restored <- countZeroTrialsWarnings({
  emptyRestored <- readRDS(emptyFile)
  maskedRestored <- readRDS(maskedFile)
  set.seed(20L)
  emptyRestoredRun <- emptyRestored$run(0L, 10L)
  set.seed(20L)
  maskedRestoredRun <- maskedRestored$run(0L, 10L)
  emptyRestored$copy()
  NULL
})
expect_identical(restored$numWarnings, 0L)
expect_identical(
  emptyRestoredRun$train[dataRows, , ],
  maskedRestoredRun$train[dataRows, , ]
)
unlink(c(emptyFile, maskedFile))
