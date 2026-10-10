# Predictor updates on a sampler with a monotone constraint that would leave a
# tree's leaf values out of order. An unforced one is refused and rolled back
# as one that would empty a leaf is, whichever tree of whichever chain it is
# and under either direction; the calls that always complete set every leaf
# of such a tree to zero and say nothing. The only change that can break the
# order is an unordered factor's first missing value, which draws the side of
# each rule on the factor: whether the order then holds is the draw's, so each
# form of the update is held to both answers, on seeds searched for them. The
# coins themselves are tests/cpp (testMonotoneMissingArrives).

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
# value in f the order is each (low; S) below (high; S). A missing value sent
# left at both f rules puts (low; a, b) below (high; c, d), which `breaks`
# violates, and sent right at both puts (low; c, d) below (high; a, b), which
# `crossed` violates; `holds` is in order either way. A decreasing constraint
# mirrors the order, and the values. A level mask is 64 bits in the machine's
# byte order.
breaks <- c(0.05, -0.05, -0.04, 0.06)
crossed <- c(-0.05, 0.05, 0.06, -0.04)
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
  direction = "increasing",
  seed = 7L
) {
  sampler <- dbarts::dbarts(
    y ~ x1 + f,
    data,
    monotone = c(x1 = direction),
    control = controlOf(length(values), n.trees),
    seed = seed
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
  x <- sampler$data@x
  held <- if (is.matrix(x)) unname(x[, 2L]) else x$dense[[2L]]
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

# the whole predictors as a data frame, the codes read as labels; a missing
# code is a missing label
asFrame <- function(codes) {
  data.frame(x1 = codes[, "x1"], f = asLabels(codes[, "f"], levels(f)))
}
# a matrix of codes is refused on this design, by name
sampler <- make()
expect_error(
  sampler$setPredictor(xMissing, forceUpdate = FALSE),
  "the predictor 'f' is a factor",
  fixed = TRUE
)
expectTwin(sampler, make())

# ---- f's first missing value: the answer follows the sides it draws ----

# Each f rule draws the side a missing value takes when the first one
# arrives, the low rule first. Left at both leaves `breaks` out of order; any
# other sides keep it. A seed's sides are read off a sampler built the same
# way that holds `holds`, in order whatever they are, which draws the same
# coins and takes the update. Each case is met with sides that break the
# order and with sides that keep it: unforced it is refused, the sampler its
# twin, or taken with its values; forced it is taken, the tree set to zero or
# kept; row by row each row that brings the value is declined, the sampler
# the twin given the rows without them, or every row is taken.
xRescaled <- xMissing
xRescaled[, "x1"] <- 0.25 + 0.5 * x1
sidesOf <- function(sampler) {
  trees <- sampler$getTrees()
  matrix(trees$missing[trees$var == 2L], 2L)
}
jointly <- function(labels, plainFirst = FALSE) {
  function(s) {
    plain <- if (plainFirst) {
      list(dbarts::dbarts(y ~ x1 + f, df, control = controlOf(), seed = 3L))
    }
    dbarts::updatePredictorPerObservationJointly(c(plain, s), labels, "f")
  }
}
whole <- function(codes, ...) function(s) s$setPredictor(asFrame(codes), ...)
byColumn <- function(labels, force) {
  function(s) s$setPredictor(labels, "f", forceUpdate = force)
}
fArrived <- asLabels(xMissing[, "f"], levels(f))
rowsLeft <- outcome(!(seq_len(n) %in% naRows))
unforced <- list(taken = outcome(TRUE), declined = outcome(FALSE))
forced <- list(taken = outcome(NULL, visible = FALSE))
byRow <- list(taken = outcome(rep(TRUE, n)), declined = rowsLeft)
cases <- list(
  c(
    list(info = "whole", bring = whole(xMissing, forceUpdate = FALSE)),
    unforced
  ),
  c(
    list(
      info = "whole, cut refresh",
      bring = whole(
        xRescaled,
        forceUpdate = FALSE,
        updateCutPoints = "position"
      )
    ),
    unforced
  ),
  c(
    list(
      info = "the second chain's second tree",
      bring = whole(xMissing, forceUpdate = FALSE),
      values = list(holds, breaks),
      still = list(holds, holds),
      args = list(n.trees = 3L, at = 2L)
    ),
    unforced
  ),
  c(
    list(
      info = "decreasing",
      bring = whole(xMissing, forceUpdate = FALSE),
      values = list(-breaks),
      still = list(-holds),
      args = list(direction = "decreasing")
    ),
    unforced
  ),
  c(list(info = "by column", bring = byColumn(fArrived, FALSE)), unforced),
  c(list(info = "whole, forced", bring = whole(xMissing)), forced),
  c(list(info = "by column, forced", bring = byColumn(fArrived, TRUE)), forced),
  c(
    list(
      info = "row by row",
      bring = byColumn(fMissingLabels, "partial"),
      leave = byColumn(fMovedLabels, "partial")
    ),
    byRow
  ),
  c(
    list(
      info = "jointly",
      bring = jointly(fMissingLabels),
      leave = jointly(fMovedLabels)
    ),
    byRow
  ),
  c(
    list(
      info = "jointly, a plain sampler first",
      bring = jointly(fMissingLabels, TRUE),
      leave = jointly(fMovedLabels, TRUE)
    ),
    byRow
  )
)
brokenSeed <- NULL
for (case in cases) {
  values <- if (is.null(case$values)) list(breaks) else case$values
  still <- if (is.null(case$still)) list(holds) else case$still
  build <- function(values, seed) {
    do.call(make, c(list(values, seed = seed), case$args))
  }
  met <- c(broken = FALSE, kept = FALSE)
  for (seed in 7:60) {
    probe <- build(still, seed)
    case$bring(probe)
    sides <- sidesOf(probe)
    arm <- if (all(sides[, ncol(sides)] == "L")) "broken" else "kept"
    if (met[[arm]]) {
      next
    }
    met[[arm]] <- TRUE
    info <- paste(case$info, arm, sep = ", ")
    sampler <- build(values, seed)
    result <- observe(case$bring(sampler))
    if (arm == "kept") {
      expect_identical(result, case$taken, info = info)
      expect_identical(leaves(sampler), unlist(values), info = info)
      expect_identical(sidesOf(sampler), sides, info = info)
      expect_identical(which(is.na(heldCodes(sampler))), naRows, info = info)
    } else if (is.null(case$declined)) {
      expect_identical(result, case$taken, info = info)
      expect_identical(leaves(sampler), numeric(4L), info = info)
    } else {
      expect_identical(result, case$declined, info = info)
      expect_identical(leaves(sampler), unlist(values), info = info)
      if (is.null(brokenSeed)) {
        # proposed again with another row changed, it is refused again: the
        # refusal put the generator back
        brokenSeed <- seed
        xAgain <- xMissing
        xAgain[moved, "f"] <- 2
        expect_false(sampler$setPredictor(asFrame(xAgain), forceUpdate = FALSE))
      }
      twin <- build(values, seed)
      if (!is.null(case$leave)) {
        expect_identical(observe(case$leave(twin)), case$taken, info = info)
      }
      expect_null(dbarts:::dataMissingSeen(sampler$data), info = info)
      expectTwin(sampler, twin)
    }
    if (all(met)) {
      break
    }
  }
  expect_identical(met, c(broken = TRUE, kept = TRUE), info = case$info)
}
# once the values are in order the proposal a seed refused is taken
sampler <- make(seed = brokenSeed)
expect_false(sampler$setPredictor(asFrame(xMissing), forceUpdate = FALSE))
expect_true(sampler$setState(handState(sampler, list(holds))))
expect_true(sampler$setPredictor(asFrame(xMissing), forceUpdate = FALSE))
expect_identical(leaves(sampler), holds)

# ---- values that stay in order: accepted and kept ----

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
expect_true(sampler$setPredictor(asFrame(xSecond), forceUpdate = FALSE))
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
expect_true(sampler$setPredictor(asFrame(xFirst), forceUpdate = FALSE))
expect_identical(leaves(sampler), breaks)
sampler <- make()
expect_true(all(sampler$setPredictor(x1Missing, "x1", forceUpdate = "partial")))
expect_identical(leaves(sampler), breaks)

sampler <- make()
expect_true(all(
  dbarts::updatePredictorPerObservationJointly(list(sampler), x1Missing, "x1")
))
expect_identical(leaves(sampler), breaks)
expect_true(is.na(sampler$data@x[naRows[1L], "x1"]))

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
# setData is the call that brings a factor its first missing value. Each f
# rule then draws the side the value goes to, the low one first, and only two
# like sides relate leaves the order did not. The draw is the first the
# sampler's generator makes, so a seed's sides are read off a sampler holding
# `holds`. Left at both is where a rule that drew nothing sends the value, so
# the first seed that sends it right at both gets `crossed`, which no other
# sides break.
arrive <- function(s) s$setData(dbarts::dbartsData(y ~ x1 + f, dfArrived))
sidesDrawn <- function(seed) {
  probe <- make(list(holds), seed = seed)
  arrive(probe)
  trees <- probe$getTrees()
  trees$missing[trees$var == 2L]
}
arrivalSeed <- 7L
while (any((sides <- sidesDrawn(arrivalSeed)) != "R") && arrivalSeed < 60L) {
  arrivalSeed <- arrivalSeed + 1L
}
expect_identical(sides, c("R", "R"))
completes(arrive, NULL, make(list(crossed), seed = arrivalSeed))
# A state stored while f could hold no missing value, installed where it can:
# the install draws the side of each f rule from the state's own generator
# before the order is judged. Left at both leaves `breaks` out of order: the
# state is declined, the sampler its twin, and forced it goes in with the
# tree set to zero. Any other sides keep the order, and the state is clean.
# Seeds are tried until each has been met.
flagged <- function(seed) make(list(holds), dfArrived, seed = seed)
met <- c(declined = FALSE, clean = FALSE)
for (seed in 7:60) {
  earlier <- make(seed = seed)$state
  sampler <- flagged(seed)
  status <- sampler$setState(earlier)
  arm <- if (status) "clean" else "declined"
  if (met[[arm]]) {
    next
  }
  met[[arm]] <- TRUE
  if (status) {
    trees <- sampler$getTrees()
    expect_true(any(trees$missing[trees$var == 2L] == "R"))
    expect_identical(leaves(sampler), breaks)
  } else {
    expectTwin(sampler, flagged(seed))
    forced <- flagged(seed)
    expect_identical(
      observe(forced$setState(earlier, forceUpdate = TRUE)),
      outcome(NULL, visible = FALSE)
    )
    expect_identical(leaves(forced), numeric(4L))
    trees <- forced$getTrees()
    expect_identical(trees$missing[trees$var == 2L], c("L", "L"))
  }
  if (all(met)) {
    break
  }
}
expect_identical(met, c(declined = TRUE, clean = TRUE))
# a warm start from a donor on another cut grid, onto a sampler whose factor
# holds a missing value: the route that maps the donor's rules onto the grid.
# The donor's f could hold none, so each f rule draws its side from the
# receiving sampler's generator before the trees are judged. Left at both
# leaves `breaks` out of order and the tree is set to zero; any other sides
# keep it, and its values. The receiving seed is searched for each, its sides
# read off a start from a donor that holds `holds`.
dfShifted <- dfArrived
dfShifted$x1 <- 0.05 + 0.9 * x1
shifted <- function(seed) {
  dbarts::dbarts(
    y ~ x1 + f,
    dfShifted,
    monotone = c(x1 = "increasing"),
    control = controlOf(),
    seed = seed
  )
}
onto <- shifted(7L)
onto$storeState()
expect_false(identical(
  attr(onto$state, "cutPoints"),
  attr(make()$state, "cutPoints")
))
met <- c(broken = FALSE, kept = FALSE)
for (seed in 7:60) {
  probe <- shifted(seed)
  probe$installTrees(make(list(holds)))
  sides <- sidesOf(probe)
  arm <- if (all(sides == "L")) "broken" else "kept"
  if (met[[arm]]) {
    next
  }
  met[[arm]] <- TRUE
  onto <- shifted(seed)
  if (arm == "broken") {
    completes(function(s) s$installTrees(make()), NULL, onto)
  } else {
    expect_identical(
      observe(onto$installTrees(make())),
      outcome(NULL, visible = FALSE)
    )
    expect_identical(leaves(onto), breaks)
    expect_identical(sidesOf(onto), sides)
  }
  if (all(met)) {
    break
  }
}
expect_identical(met, c(broken = TRUE, kept = TRUE))

# a merge that leaves the tree out of order: with no row left in (high; a, b)
# the high half becomes one leaf holding the value of (high; c, d), which
# `breaks` puts below (low; a, b). Forced by column, and a state forced in
# over such predictors, which without force is declined for the merge.
fEmptied <- f
fEmptied[x1 > cut & f %in% c("a", "b")] <- "c"
completes(
  function(s) s$setPredictor(fEmptied, "f", forceUpdate = TRUE),
  NULL,
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
sampler <- emptied()
expect_identical(sampler$setState(handState(sampler, list(breaks))), FALSE)
expect_identical(leaves(sampler), numeric(0L))
completes(
  function(s) s$setState(handState(s, list(breaks)), forceUpdate = TRUE),
  NULL,
  emptied(),
  numLeaves = 3L
)
# where the merged tree stays in order its values are kept
sampler <- emptied()
merging <- handState(sampler, list(holds))
expect_identical(sampler$setState(merging), FALSE)
expect_null(sampler$setState(merging, forceUpdate = TRUE))
expect_identical(leaves(sampler), holds[1:3])
# setCutPoints completes in silence too: on a shorter grid the x1 split
# keeps its place, rescaled, and the tree its leaves; a grid of one point
# below every x1 value leaves the split a side no row reaches, and the tree
# a single leaf
sampler <- make()
expect_identical(
  observe(sampler$setCutPoints(c(0.25, 0.75), "x1")),
  outcome(NULL, visible = FALSE)
)
expect_identical(length(leaves(sampler)), 4L)
sampler <- make()
expect_identical(
  observe(sampler$setCutPoints(min(x1) - 1, "x1")),
  outcome(NULL, visible = FALSE)
)
expect_identical(leaves(sampler), numeric())

# ---- a factor of more than 63 levels that regains a missing value ----

# A column that has held a missing value can hold one again: its rules keep
# their side for it, the order is judged with it, and a regained one is taken
# as a label, where a first one is refused. Both rules here send it right:
# the low half splits levels 1-35 | 36-70 and missing, the high half 36-70 |
# 1-35 and missing. `crossed` is in order only if the two right-hand leaves
# could share no position.
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
installPooled <- function(sampler, values, ...) {
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
  sampler$setState(state, ...)
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
# the column loses its missing values and may still hold one, so values that
# cross with one are still out of order: declined, and reseeded when forced
asFrameH <- function(codes) {
  data.frame(x1 = codes[, "x1"], h = asLabels(codes[, "h"], levels70))
}
expect_true(sampler$setPredictor(asFrameH(pooledCodes), forceUpdate = FALSE))
expect_identical(dbarts:::dataMissingSeen(sampler$data), c(FALSE, TRUE))
expect_identical(installPooled(sampler, crossed), FALSE)
expect_identical(leaves(sampler), holds)
reseeded <- sampler$copy()
expect_null(installPooled(reseeded, crossed, forceUpdate = TRUE))
expect_identical(leaves(reseeded), numeric(4L))
rm(reseeded)
xRegained <- pooledCodes
regained <- which(x1 <= cut & as.integer(h) <= 35L)[1L]
xRegained[regained, "h"] <- NA
taken <- NULL
expect_silent(
  taken <- sampler$setPredictor(asFrameH(xRegained), forceUpdate = FALSE)
)
expect_true(taken)
expect_identical(leaves(sampler), holds)
expect_true(dbarts:::sourceHasNA(sampler$data@x))
expect_silent(
  taken <- dbarts::updatePredictorPerObservationJointly(
    list(sampler),
    asLabels(xRegained[, "h"], levels70),
    "h"
  )
)
expect_true(all(taken))
