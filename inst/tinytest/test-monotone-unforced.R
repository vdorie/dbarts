# Predictor updates on a sampler with a monotone constraint that would leave a
# tree's leaf values out of order. An unforced one is refused and rolled back
# as one that would empty a leaf is, whichever tree of whichever chain it is
# and under either direction; the calls that always complete set every leaf
# of such a tree to zero and say nothing. The engine's side, the column and
# single-sampler row forms included, is tests/cpp
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
controlOf <- function(n.chains = 1L, n.trees = 1L) {
  dbarts::dbartsControl(
    n.trees = n.trees,
    n.chains = n.chains,
    n.threads = 1L,
    updateState = FALSE,
    verbose = FALSE
  )
}

# The hand tree: x1 (increasing) cut once near 0.5, the low half split on f
# as {a, b} | {c, d} and the high half as {c, d} | {a, b}. Leaves in order:
# (low; a, b), (low; c, d), (high; c, d), (high; a, b). Without a missing
# value in f the order is each (low; S) below (high; S). A missing value goes
# left at both f rules and puts (low; a, b) below (high; c, d), which `breaks`
# violates and `holds` does not. A decreasing constraint mirrors the order,
# and the values. A level mask is 64 bits in the machine's byte order.
breaks <- c(0.05, -0.05, -0.04, 0.06)
holds <- c(-0.05, -0.06, 0.04, 0.06)
maskBytes <- function(bits) {
  words <- c(as.integer(bits), 0L)
  writeBin(if (.Platform$endian == "big") rev(words) else words, raw())
}
# A state with the hand tree as tree `at` of each chain, values[[chain]] at
# its leaves, and the chain's other trees single leaves at zero.
handState <- function(sampler, values, at = 1L) {
  sampler$storeState()
  state <- sampler$state
  cuts <- attr(state, "cutPoints")[[1L]]
  cut <- cuts[which.max(cuts >= 0.5)]
  for (chain in seq_along(values)) {
    forest <- state[[chain]]$forests[[1L]]
    numBefore <- at - 1L
    numAfter <- length(forest$tree.sizes) - at
    forest$tree.vars <- c(
      rep(-1L, numBefore),
      c(1L, 2L, -1L, -1L, 2L, -1L, -1L),
      rep(-1L, numAfter)
    )
    forest$tree.values <- c(
      writeBin(numeric(numBefore), raw()),
      writeBin(cut, raw()),
      maskBytes(12L),
      writeBin(values[[chain]][1:2], raw()),
      maskBytes(3L),
      writeBin(values[[chain]][3:4], raw()),
      writeBin(numeric(numAfter), raw())
    )
    forest$tree.sizes <- c(rep(1L, numBefore), 7L, rep(1L, numAfter))
    forest$tree.flags <- as.raw(c(
      integer(numBefore),
      c(2L, 4L, 0L, 0L, 4L, 0L, 0L),
      integer(numAfter)
    ))
    state[[chain]]$forests[[1L]] <- forest
  }
  state
}
make <- function(
  values = list(breaks),
  data = df,
  n.trees = 1L,
  at = 1L,
  direction = "increasing"
) {
  sampler <- dbarts::dbarts(
    y ~ x1 + f,
    data,
    monotone = c(x1 = direction),
    control = controlOf(length(values), n.trees),
    seed = 7L
  )
  stopifnot(isTRUE(sampler$setState(handState(sampler, values, at))))
  sampler
}
# the leaf values of every tree that has split, in order
leaves <- function(sampler) {
  trees <- sampler$getTrees()
  trees$value[trees$var < 0L & trees$n < n]
}
# the factor column of data@x as codes from 0, however the sampler holds it
heldCodes <- function(sampler) {
  held <- sampler$data@x$dense[[2L]]
  if (is.factor(held)) as.double(as.integer(held) - 1L) else held
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
# (trees, sigma, generator), the predictors, the predictions and the cached
# fits identical. The next draws are equal to rounding and not bit for bit, as
# after any refusal: a sweep sums a leaf's rows in the order the leaf holds
# them, which the refusal's two re-routes can change.
expectTwin <- function(sampler, twin) {
  sampler$storeState()
  twin$storeState()
  expect_identical(sampler$state, twin$state)
  expect_identical(sampler$data@x, twin$data@x)
  expect_identical(sampler$predict(df), twin$predict(df))
  expect_identical(sampler$getForestFits(), twin$getForestFits())
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
# the rows given a missing value are rows it moves to another leaf, and
# twenty of them, so that a fit rebuilt from the new partition, or a
# partition left as the refused values route it, shows in the draws that
# follow; `moved` changes level instead
naRows <- which(x1 > cut & f %in% c("a", "b"))[1:20]
moved <- which(x1 <= cut & f == "a")[1L]
codes <- cbind(x1 = x1, f = as.double(as.integer(f) - 1L))
xMissing <- codes
xMissing[naRows, "f"] <- NA
fMoved <- codes[, "f"]
fMoved[moved] <- 2
fMissing <- fMoved
fMissing[naRows] <- NA
# the same columns as labels, which is what the joint form takes
asLabels <- function(codes, levels) factor(levels[codes + 1L], levels = levels)
fMovedLabels <- asLabels(fMoved, levels(f))
fMissingLabels <- asLabels(fMissing, levels(f))
missingRefusal <- "column 'f' has missing values, which its training values do not"

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

# three trees a chain with the hand tree second, out of order in the second
# chain alone: neither a chain's first tree nor the first chain shows it
several <- function() make(list(holds, breaks), n.trees = 3L, at = 2L)
sampler <- several()
twin <- several()
expect_identical(
  observe(sampler$setPredictor(xMissing, forceUpdate = FALSE)),
  outcome(FALSE)
)
expect_identical(leaves(sampler), c(holds, breaks))
expectTwin(sampler, twin)
expect_error(
  dbarts::updatePredictorPerObservationJointly(
    list(several()),
    fMissingLabels,
    "f"
  ),
  missingRefusal,
  fixed = TRUE
)

# a decreasing constraint, with the mirrored values
sampler <- make(list(-breaks), direction = "decreasing")
expect_false(sampler$setPredictor(xMissing, forceUpdate = FALSE))
expect_identical(leaves(sampler), -breaks)
sampler <- make(list(-holds), direction = "decreasing")
expect_true(sampler$setPredictor(xMissing, forceUpdate = FALSE))
expect_identical(leaves(sampler), -holds)

# ---- the joint form: a first missing value is refused before any row moves ----

# a column that holds no missing value takes none by label, so no row is
# declined for one; the refusal leaves the sampler its twin, and labels
# without one move every row
sampler <- make()
twin <- make()
expect_error(
  dbarts::updatePredictorPerObservationJointly(
    list(sampler),
    fMissingLabels,
    "f"
  ),
  missingRefusal,
  fixed = TRUE
)
expectTwin(sampler, twin)
expect_true(all(
  dbarts::updatePredictorPerObservationJointly(
    list(sampler),
    fMovedLabels,
    "f"
  )
))
expect_identical(heldCodes(sampler), fMoved)

# a plain sampler listed first refuses it too
plain <- dbarts::dbarts(y ~ x1 + f, df, control = controlOf(), seed = 3L)
sampler <- make()
expect_error(
  dbarts::updatePredictorPerObservationJointly(
    list(plain, sampler),
    fMissingLabels,
    "f"
  ),
  missingRefusal,
  fixed = TRUE
)

# ---- proposing again: refused until the trees move ----

# other values with the same missing row are refused too, whole and row by
# row, each refusal having left the column without a missing value; once the
# values are in order the first proposal is accepted
sampler <- make()
expect_false(sampler$setPredictor(xMissing, forceUpdate = FALSE))
xAgain <- xMissing
xAgain[moved, "f"] <- 2
expect_false(sampler$setPredictor(xAgain, forceUpdate = FALSE))
expect_error(
  dbarts::updatePredictorPerObservationJointly(
    list(sampler),
    fMissingLabels,
    "f"
  ),
  missingRefusal,
  fixed = TRUE
)
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
expect_identical(which(is.na(sampler$data@x[, "f"])), naRows)
sampler <- make(list(holds))
expect_true(all(
  dbarts::updatePredictorPerObservationJointly(
    list(sampler),
    fMovedLabels,
    "f"
  )
))
expect_identical(leaves(sampler), holds)

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
# `call` returns `value` invisibly and warns of nothing, the hand tree is
# left with `numLeaves` leaves at zero, and a sweep later the fit is monotone
# at every level, the missing one included where the factor holds one
completes <- function(call, value, onto = make(), numLeaves = 4L) {
  expect_identical(observe(call(onto)), outcome(value, visible = FALSE))
  expect_identical(leaves(onto), numeric(numLeaves))
  invisible(onto$run(0L, 1L))
  grid <- levelGrid[numLeaves == 4L | !is.na(levelGrid$f), ]
  expect_true(maxDrop(onto, grid) <= 1e-8)
}
completes(function(s) s$setPredictor(xMissing, forceUpdate = TRUE), TRUE)
completes(function(s) s$setPredictor(xMissing), TRUE)
completes(
  function(s) s$setData(dbarts::dbartsData(y ~ x1 + f, dfArrived)),
  NULL
)
# a warm start from a donor on another cut grid, onto a sampler whose factor
# holds a missing value: the route that maps the donor's rules onto the grid
dfShifted <- dfArrived
dfShifted$x1 <- 0.05 + 0.9 * x1
onto <- dbarts::dbarts(
  y ~ x1 + f,
  dfShifted,
  monotone = c(x1 = "increasing"),
  control = controlOf(),
  seed = 7L
)
onto$storeState()
expect_false(identical(
  attr(onto$state, "cutPoints"),
  attr(make()$state, "cutPoints")
))
completes(function(s) s$installTrees(make()), NULL, onto)

# a merge that leaves the tree out of order: with no row left in (high; a, b)
# the high half becomes one leaf holding the value of (high; c, d), which
# `breaks` puts below (low; a, b). Forced by column, and a state installed
# over such predictors, which returns FALSE for the merge.
fEmptied <- f
fEmptied[x1 > cut & f %in% c("a", "b")] <- "c"
completes(
  function(s) s$setPredictor(fEmptied, "f", forceUpdate = TRUE),
  TRUE,
  numLeaves = 3L
)
emptied <- function() {
  dbarts::dbarts(
    y ~ x1 + f,
    data.frame(y, x1, f = fEmptied),
    monotone = c(x1 = "increasing"),
    control = controlOf(),
    seed = 7L
  )
}
completes(
  function(s) s$setState(handState(s, list(breaks))),
  FALSE,
  emptied(),
  numLeaves = 3L
)
# where the merged tree stays in order its values are kept
sampler <- emptied()
expect_false(sampler$setState(handState(sampler, list(holds))))
expect_identical(leaves(sampler), holds[1:3])
# setCutPoints completes in silence too: a grid without the x1 cut leaves
# the tree a single leaf
sampler <- make()
expect_identical(
  observe(sampler$setCutPoints(c(0.25, 0.75), "x1")),
  outcome(NULL, visible = FALSE)
)
expect_identical(leaves(sampler), numeric())

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

# ---- a factor of more than 63 levels that regains a missing value ----

# Such a factor's rules keep their side for a missing value while the column
# has none, so a regained one can go right. Both rules here send it right:
# the low half splits levels 1-35 | 36-70 and missing, the high half 36-70 |
# 1-35 and missing. `crossed` is in order only while the two right-hand
# leaves share no position.
levels70 <- sprintf("L%02d", 1:70)
h <- factor(levels70[rep(1:70, length.out = n)], levels = levels70)
dfPooled <- data.frame(y, x1, h)
dfPooled$h[which(x1 > 0.6 & as.integer(h) <= 35L)[1:2]] <- NA
# the 128-bit level set of a rule as two words in the machine's byte order
poolBytes <- function(levelCodes) {
  set <- logical(128L)
  set[levelCodes + 1L] <- TRUE
  bytes <- packBits(set, "raw")
  if (.Platform$endian == "big") bytes[c(matrix(16:1, 8L)[, 2:1])] else bytes
}
installPooled <- function(sampler, values) {
  sampler$storeState()
  state <- sampler$state
  cuts <- attr(state, "cutPoints")[[1L]]
  forest <- state[[1L]]$forests[[1L]]
  forest$tree.vars <- c(1L, 2L, -1L, -1L, 2L, -1L, -1L)
  forest$tree.values <- c(
    writeBin(cuts[which.max(cuts >= 0.5)], raw()),
    maskBytes(0L),
    writeBin(values[1:2], raw()),
    maskBytes(2L),
    writeBin(values[3:4], raw())
  )
  forest$tree.flags <- as.raw(c(2L, 7L, 0L, 0L, 7L, 0L, 0L))
  forest$tree.sizes <- 7L
  forest$tree.masks <- c(poolBytes(35:69), poolBytes(0:34))
  state[[1L]]$forests[[1L]] <- forest
  sampler$setState(state)
}
pooledCodes <- cbind(x1 = x1, h = as.double(as.integer(h) - 1L))
sampler <- dbarts::dbarts(
  y ~ x1 + h,
  dfPooled,
  monotone = c(x1 = "increasing"),
  control = controlOf(),
  seed = 7L
)
expect_true(installPooled(sampler, holds))
# the column loses its missing values, and the values may then cross
expect_true(sampler$setPredictor(pooledCodes, forceUpdate = FALSE))
crossed <- c(-0.05, 0.05, 0.06, -0.04)
expect_true(installPooled(sampler, crossed))
xRegained <- pooledCodes
regained <- which(x1 <= cut & as.integer(h) <= 35L)[1L]
xRegained[regained, "h"] <- NA
expect_false(sampler$setPredictor(xRegained, forceUpdate = FALSE))
expect_identical(leaves(sampler), crossed)
expect_error(
  dbarts::updatePredictorPerObservationJointly(
    list(sampler),
    asLabels(xRegained[, "h"], levels70),
    "h"
  ),
  "column 'h' has missing values, which its training values do not",
  fixed = TRUE
)
