# Predictor updates on a sampler with a monotone constraint that would leave a
# tree's leaf values out of order. An unforced one is refused and rolled back
# as one that would empty a leaf is; the calls that always complete set every
# leaf of such a tree to zero and say nothing. The engine's side, the column
# and single-sampler row forms included, is tests/cpp
# (testMonotoneMissingArrives).

source(
  system.file("common", "captureWarnings.R", package = "dbarts"),
  local = TRUE
)

set.seed(20261006L)
n <- 400L
x1 <- runif(n)
f <- factor(rep(letters[1:4], length.out = n))
y <- x1 + rnorm(n, sd = 0.1)
df <- data.frame(y, x1, f)
controlOf <- function(n.chains = 1L) {
  dbarts::dbartsControl(
    n.trees = 1L,
    n.chains = n.chains,
    n.threads = 1L,
    updateState = FALSE,
    verbose = FALSE
  )
}

# The tree, one per chain: x1 (increasing) cut once near 0.5, the low half
# split on f as {a, b} | {c, d} and the high half as {c, d} | {a, b}. Leaves
# in order: (low; a, b), (low; c, d), (high; c, d), (high; a, b). Without a
# missing value in f the order is each (low; S) below (high; S). A missing
# value goes left at both f rules and puts (low; a, b) below (high; c, d),
# which `breaks` violates and `holds` does not. A level mask is 64 bits in the
# machine's byte order.
breaks <- c(0.05, -0.05, -0.04, 0.06)
holds <- c(-0.05, -0.06, 0.04, 0.06)
maskBytes <- function(bits) {
  words <- c(as.integer(bits), 0L)
  writeBin(if (.Platform$endian == "big") rev(words) else words, raw())
}
handState <- function(sampler, values) {
  sampler$storeState()
  state <- sampler$state
  cuts <- attr(state, "cutPoints")[[1L]]
  cut <- cuts[which.max(cuts >= 0.5)]
  for (chain in seq_along(values)) {
    forest <- state[[chain]]$forests[[1L]]
    forest$tree.vars <- c(1L, 2L, -1L, -1L, 2L, -1L, -1L)
    forest$tree.values <- c(
      writeBin(cut, raw()),
      maskBytes(12L),
      writeBin(values[[chain]][1:2], raw()),
      maskBytes(3L),
      writeBin(values[[chain]][3:4], raw())
    )
    forest$tree.sizes <- 7L
    forest$tree.flags <- as.raw(c(2L, 4L, 0L, 0L, 4L, 0L, 0L))
    state[[chain]]$forests[[1L]] <- forest
  }
  state
}
make <- function(values = list(breaks), data = df) {
  sampler <- dbarts::dbarts(
    y ~ x1 + f,
    data,
    monotone = c(x1 = "increasing"),
    control = controlOf(length(values)),
    seed = 7L
  )
  stopifnot(isTRUE(sampler$setState(handState(sampler, values))))
  sampler
}
leaves <- function(sampler) {
  trees <- sampler$getTrees()
  trees$value[trees$var < 0L]
}
# a call's value, whether it was visible, and how many warnings it raised;
# `outcome` is that for a call that raised none
observe <- function(expr) {
  out <- NULL
  warnings <- captureWarnings(out <- withVisible(expr))
  c(out, numWarnings = length(warnings))
}
outcome <- function(value, visible = TRUE) {
  list(value = value, visible = visible, numWarnings = 0L)
}
# A refused update leaves the sampler its untouched twin: the stored state
# (trees, sigma, generator), the predictors and the predictions identical.
# The next draws are equal to rounding and not bit for bit, as after any
# refusal: a sweep sums a leaf's rows in the order the leaf holds them, which
# the refusal's two re-routes can change.
expectTwin <- function(sampler, twin) {
  sampler$storeState()
  twin$storeState()
  expect_identical(sampler$state, twin$state)
  expect_identical(sampler$data@x, twin$data@x)
  expect_identical(sampler$predict(df), twin$predict(df))
  draws <- sampler$run(0L, 5L)
  twinDraws <- twin$run(0L, 5L)
  expect_equal(draws$train, twinDraws$train, tolerance = 1e-12)
  expect_equal(draws$sigma, twinDraws$sigma, tolerance = 1e-12)
}
# the largest fall of the fit along x1 at any value of the grid's other
# column, x1 varying fastest over 101 points; 0 when monotone
maxDrop <- function(sampler, grid) {
  -min(apply(matrix(sampler$predict(grid), 101L), 2L, diff))
}
x1Grid <- seq(0, 1, length.out = 101L)

cut <- with(list(cuts = attr(make()$state, "cutPoints")[[1L]]), {
  cuts[which.max(cuts >= 0.5)]
})
# the rows given a missing value are rows it moves to another leaf, so a fit
# rebuilt from the new partition would show; `moved` changes level instead
naRows <- which(x1 > cut & f %in% c("a", "b"))[1:2]
moved <- which(x1 <= cut & f == "a")[1L]
codes <- cbind(x1 = x1, f = as.double(as.integer(f) - 1L))
xMissing <- codes
xMissing[naRows[1L], "f"] <- NA
fMoved <- codes[, "f"]
fMoved[moved] <- 2
fMissing <- fMoved
fMissing[naRows] <- NA

# ---- whole matrix of codes, unforced: refused, the sampler its twin ----

sampler <- make()
twin <- make()
xBefore <- sampler$data@x
expect_identical(leaves(sampler), breaks)
expect_identical(
  observe(sampler$setPredictor(xMissing, forceUpdate = FALSE)),
  outcome(FALSE)
)
expect_identical(sampler$data@x, xBefore)
expect_identical(leaves(sampler), breaks)
expectTwin(sampler, twin)

# with a cut refresh over a rescaled x1, which moves every cut point
xRescaled <- xMissing
xRescaled[, "x1"] <- 0.25 + 0.5 * x1
sampler <- make()
twin <- make()
expect_identical(
  observe(sampler$setPredictor(
    xRescaled,
    forceUpdate = FALSE,
    updateCutPoints = TRUE
  )),
  outcome(FALSE)
)
expectTwin(sampler, twin)

# two chains, the tree out of order in the second alone
sampler <- make(list(holds, breaks))
twin <- make(list(holds, breaks))
expect_identical(
  observe(sampler$setPredictor(xMissing, forceUpdate = FALSE)),
  outcome(FALSE)
)
expect_identical(leaves(sampler), c(holds, breaks))
expectTwin(sampler, twin)

# ---- the joint form: the rows bringing the value are declined ----

# the twin is given the column the call left, row by row too: the scan order
# is drawn either way
sampler <- make()
twin <- make()
installed <- observe(
  dbarts::updatePredictorPerObservationJointly(list(sampler), fMissing, "f")
)
expect_identical(installed, outcome(!seq_len(n) %in% naRows))
expect_identical(sampler$data@x$dense[[2L]], fMoved)
expect_true(all(
  dbarts::updatePredictorPerObservationJointly(list(twin), fMoved, "f")
))
expectTwin(sampler, twin)

# a plain sampler listed first declines them too, and holds the same column
plain <- dbarts::dbarts(y ~ x1 + f, df, control = controlOf(), seed = 3L)
sampler <- make()
installed <- dbarts::updatePredictorPerObservationJointly(
  list(plain, sampler),
  fMissing,
  "f"
)
expect_identical(which(!installed), naRows)
expect_identical(plain$data@x$dense[[2L]], fMoved)
expect_identical(sampler$data@x$dense[[2L]], fMoved)

# ---- proposing again: refused until the trees move ----

# other values with the same missing row are refused too; once the values are
# in order the first proposal is accepted
sampler <- make()
expect_false(sampler$setPredictor(xMissing, forceUpdate = FALSE))
xAgain <- xMissing
xAgain[moved, "f"] <- 2
expect_false(sampler$setPredictor(xAgain, forceUpdate = FALSE))
expect_true(sampler$setState(handState(sampler, list(holds))))
expect_true(sampler$setPredictor(xMissing, forceUpdate = FALSE))
expect_identical(leaves(sampler), holds)

# ---- values that stay in order: accepted and kept ----

sampler <- make(list(holds))
expect_identical(
  observe(sampler$setPredictor(xMissing, forceUpdate = FALSE)),
  outcome(TRUE)
)
expect_identical(leaves(sampler), holds)
expect_true(is.na(sampler$data@x[naRows[1L], "f"]))
sampler <- make(list(holds))
expect_true(all(
  dbarts::updatePredictorPerObservationJointly(list(sampler), fMissing, "f")
))
expect_identical(leaves(sampler), holds)
expect_identical(which(is.na(sampler$data@x$dense[[2L]])), naRows)

# a factor that holds a missing value from the start takes another
dfMissing <- df
dfMissing$f[moved] <- NA
xSecond <- xMissing
xSecond[moved, "f"] <- NA
sampler <- make(list(holds), dfMissing)
expect_true(sampler$setPredictor(xSecond, forceUpdate = FALSE))
expect_identical(leaves(sampler), holds)

# a first missing value in the constrained x1 relates no new leaves: accepted
# by column, whole and row by row, on the tree a missing f would break
x1Missing <- x1
x1Missing[naRows[1L]] <- NA
xFirst <- codes
xFirst[naRows[1L], "x1"] <- NA
sampler <- make()
expect_true(sampler$setPredictor(x1Missing, "x1"))
expect_identical(leaves(sampler), breaks)
sampler <- make()
expect_true(sampler$setPredictor(xFirst, forceUpdate = FALSE))
expect_identical(leaves(sampler), breaks)
sampler <- make()
expect_true(all(sampler$setPredictor(x1Missing, "x1", forceUpdate = "partial")))
expect_identical(leaves(sampler), breaks)

# ---- the calls that complete: taken in silence, the tree set to zero ----

dfArrived <- df
dfArrived$f[naRows[1L]] <- NA
levelGrid <- expand.grid(x1 = x1Grid, f = factor(c(letters[1:4], NA)))
completes <- function(call, value, onto = make()) {
  expect_identical(observe(call(onto)), outcome(value, visible = FALSE))
  expect_identical(leaves(onto), numeric(4L))
  invisible(onto$run(0L, 1L))
  expect_true(maxDrop(onto, levelGrid) <= 1e-8)
}
completes(function(s) s$setPredictor(xMissing, forceUpdate = TRUE), TRUE)
completes(function(s) s$setPredictor(xMissing), TRUE)
completes(
  function(s) s$setData(dbarts::dbartsData(y ~ x1 + f, dfArrived)),
  NULL
)
completes(
  function(s) s$installTrees(make()),
  NULL,
  dbarts::dbarts(
    y ~ x1 + f,
    dfArrived,
    monotone = c(x1 = "increasing"),
    control = controlOf(),
    seed = 7L
  )
)

# by column and "partial" the factor is given as a factor, and its first
# missing value is refused by name before the engine sees it, forced or not
fLabels <- f
fLabels[naRows[1L]] <- NA
sampler <- make()
twin <- make()
for (force in list(FALSE, TRUE, "partial")) {
  expect_error(
    sampler$setPredictor(fLabels, "f", forceUpdate = force),
    "column 'f' has missing values, which its training values do not"
  )
}
expectTwin(sampler, twin)

# ---- a merge that leaves a tree out of order, on numeric predictors ----

# x1 at 0.5; the low half splits x2 at 0.3 into J, K and the high half at 0.7
# into S1, S2, with J and K below S1 and K below S2. Emptying S1 merges the
# high half into one leaf with S2's value, which then sits below J.
x2 <- runif(n)
xNumeric <- cbind(x1 = x1, x2 = x2)
numericSampler <- function(x) {
  dbarts::dbarts(
    x,
    y,
    monotone = c(x1 = "increasing"),
    control = controlOf(),
    seed = 3L
  )
}
makeNumeric <- function() {
  sampler <- numericSampler(xNumeric)
  sampler$storeState()
  state <- sampler$state
  cuts <- attr(state, "cutPoints")
  at <- function(j, value) cuts[[j]][which.max(cuts[[j]] >= value)]
  forest <- state[[1L]]$forests[[1L]]
  forest$tree.vars <- c(1L, 2L, -1L, -1L, 2L, -1L, -1L)
  forest$tree.values <- writeBin(
    c(at(1L, 0.5), at(2L, 0.3), 0, -0.2, at(2L, 0.7), 0.01, -0.19),
    raw()
  )
  forest$tree.sizes <- 7L
  forest$tree.flags <- as.raw(c(2L, 2L, 0L, 0L, 2L, 0L, 0L))
  state[[1L]]$forests[[1L]] <- forest
  stopifnot(isTRUE(sampler$setState(state)))
  list(sampler = sampler, emptied = x1 > at(1L, 0.5) & x2 <= at(2L, 0.7))
}
built <- makeNumeric()
sampler <- built$sampler
xEmptied <- xNumeric
xEmptied[built$emptied, "x2"] <- 0.9

# unforced it is refused for the empty leaf; forced by column it resets
expect_false(sampler$setPredictor(xEmptied[, "x2"], "x2"))
expect_identical(
  observe(sampler$setPredictor(xEmptied[, "x2"], "x2", forceUpdate = TRUE)),
  outcome(TRUE, visible = FALSE)
)
expect_identical(leaves(sampler), numeric(3L))

# a setState that has to merge resets too, returning FALSE invisibly
donor <- makeNumeric()$sampler
donor$storeState()
sampler <- numericSampler(xEmptied)
expect_identical(
  observe(sampler$setState(donor$state)),
  outcome(FALSE, visible = FALSE)
)
expect_identical(leaves(sampler), numeric(3L))

# setCutPoints completes in silence and leaves the fit in order
sampler <- makeNumeric()$sampler
expect_identical(
  observe(sampler$setCutPoints(0.9, "x2")),
  outcome(NULL, visible = FALSE)
)
numericGrid <- as.matrix(expand.grid(x1 = x1Grid, x2 = c(0.1, 0.5, 0.95)))
expect_true(maxDrop(sampler, numericGrid) <= 1e-8)
