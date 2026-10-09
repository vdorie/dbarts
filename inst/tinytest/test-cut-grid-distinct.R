# No cut grid repeats a point. Every grid a sampler derives - at creation, at
# setData, at a refresh through setPredictor - holds each point once, so a
# column whose values supply fewer points than n.cuts holds fewer and a
# constant column one. Under either rule a column with fewer distinct values
# than n.cuts holds one point in each gap between them, which the default
# rule weighs by the gap's width. A refresh derives the grid a creation would for the
# values, whatever number of points the column held, and says where the
# splits on the column go: on their position ("position") or to the point
# nearest their threshold ("value"); setCutPoints takes the same choice as
# splits. A state whose grid repeats a point is refused by name, and a stored
# tree that stacks two splits on one value installs with the lower one merged.

cutPointsOf <- function(sampler) {
  sampler$storeState()
  attr(sampler$state, "cutPoints")
}
strictly <- function(grids) {
  all(vapply(
    grids,
    function(grid) length(grid) >= 1L && !is.unsorted(grid, strictly = TRUE),
    NA
  ))
}
controlWith <- function(n.cuts = 100L, useQuantiles = FALSE, n.trees = 10L) {
  dbartsControl(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = n.trees,
    n.cuts = n.cuts,
    useQuantiles = useQuantiles,
    updateState = FALSE,
    seed = 11L
  )
}
printedTrees <- function(sampler) {
  capture.output(sampler$printTrees())
}
# the splits on a column, one row a node, in the order getTrees gives them
splitsOn <- function(sampler, column) {
  trees <- sampler$getTrees()
  trees[trees$var == column, c("tree", "n", "value")]
}
warningsOf <- function(expr) {
  seen <- list()
  withCallingHandlers(
    expr,
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  seen
}
refreshWords <- c("position", "value")

set.seed(9173L)
n <- 200L
eps <- .Machine$double.eps
z <- rnorm(n)
w <- rnorm(n)
y <- z + sin(2 * w) + rnorm(n, 0, 0.3)
narrow <- 1 + eps * (seq_len(n) %% 5L)
const <- rep(3, n)
constNA <- replace(const, 1:40, NA)
binary <- as.double(seq_len(n) %% 2L)
six <- as.double(seq_len(n) %% 6L)

# ---- the rules, by their examples -------------------------------------

for (useQuantiles in c(FALSE, TRUE)) {
  rule <- if (useQuantiles) "quantile" else "uniform"
  # five adjacent doubles: a point in each gap, on its lower value
  numNarrow <- 4L

  # a grid is a set of distinct thresholds
  ctl <- controlWith(100L, useQuantiles)
  sampler <- dbarts(cbind(z, narrow), y, control = ctl)
  expect_identical(
    lengths(cutPointsOf(sampler)),
    c(100L, numNarrow),
    info = rule
  )
  sampler <- dbarts(cbind(z, constNA), y, control = ctl)
  expect_identical(cutPointsOf(sampler)[[2L]], 3, info = rule)

  # a refresh derives up to n.cuts points from the new values, whatever the
  # column held; setData already did
  ctl <- controlWith(20L, useQuantiles)
  sampler <- dbarts(cbind(z, latent = 0), y, control = ctl)
  expect_identical(cutPointsOf(sampler)[[2L]], 0, info = rule)
  expect_null(
    sampler$setPredictor(
      w,
      2L,
      forceUpdate = TRUE,
      updateCutPoints = "position"
    ),
    info = rule
  )
  expect_identical(lengths(cutPointsOf(sampler)), c(20L, 20L), info = rule)
  sampler$setCutPoints(seq(-2, 2, length.out = 60L), 2L)
  expect_identical(lengths(cutPointsOf(sampler)), c(20L, 60L), info = rule)
  sampler$setData(dbartsData(cbind(z, latent = w), y))
  expect_identical(lengths(cutPointsOf(sampler)), c(20L, 20L), info = rule)

  # a refresh never fails for want of distinct points: the grid shrinks
  expect_null(
    sampler$setPredictor(
      narrow,
      2L,
      forceUpdate = TRUE,
      updateCutPoints = "position"
    ),
    info = rule
  )
  expect_identical(length(cutPointsOf(sampler)[[2L]]), 4L, info = rule)
  # and onto a single value it is one point
  expect_null(
    sampler$setPredictor(
      const,
      2L,
      forceUpdate = TRUE,
      updateCutPoints = "value"
    ),
    info = rule
  )
  expect_identical(cutPointsOf(sampler)[[2L]], 3, info = rule)
}

# ---- no grid repeats, on any path ---------------------------------------

columns <- cbind(z, narrow, const, constNA, binary, six)
others <- matrix(rnorm(n * ncol(columns)), n, dimnames = dimnames(columns))
for (useQuantiles in c(FALSE, TRUE)) {
  rule <- if (useQuantiles) "quantile" else "uniform"
  ctl <- controlWith(100L, useQuantiles)
  counts <- c(100L, 4L, 1L, 1L, 1L, 5L)

  sampler <- dbarts(columns, y, control = ctl)
  created <- cutPointsOf(sampler)
  expect_true(strictly(created), info = rule)
  expect_identical(lengths(created), counts, info = rule)
  expect_identical(created[3L:4L], list(3, 3), info = rule)

  # setData and a refresh from a sampler created on other values
  replaced <- dbarts(others, y, control = ctl)
  replaced$setData(dbartsData(columns, y))
  expect_identical(cutPointsOf(replaced), created, info = rule)
  for (word in refreshWords) {
    refreshed <- dbarts(others, y, control = ctl)
    expect_null(
      refreshed$setPredictor(
        columns,
        forceUpdate = TRUE,
        updateCutPoints = word
      ),
      info = paste(rule, word)
    )
    expect_identical(cutPointsOf(refreshed), created, info = paste(rule, word))
  }

  # a state, a copy and a reload carry the grids as they are
  invisible(sampler$run(0L, 5L))
  sampler$storeState()
  expect_true(sampler$setState(sampler$state), info = rule)
  expect_identical(cutPointsOf(sampler), created, info = rule)
  expect_identical(cutPointsOf(sampler$copy()), created, info = rule)
  path <- tempfile(fileext = ".rds")
  saveRDS(sampler, path)
  expect_identical(cutPointsOf(readRDS(path)), created, info = rule)
  unlink(path)
  # and the grids a sampler reports go back into it
  expect_silent(sampler$setCutPoints(created))
  expect_identical(cutPointsOf(sampler), created, info = rule)

  # a sparse design derives the dense design's grids, created and refreshed
  if (requireNamespace("Matrix", quietly = TRUE)) {
    plain <- cbind(binary, narrow, zero = 0, six)
    start <- cbind(binary, narrow = w, zero = z, six = w^2)
    dense <- cutPointsOf(dbarts(plain, y, control = ctl))
    sparse <- dbarts(
      Matrix::Matrix(plain, sparse = TRUE),
      y,
      control = ctl,
      sigest = sd(y)
    )
    expect_identical(cutPointsOf(sparse), dense, info = rule)
    sparse <- dbarts(
      Matrix::Matrix(start, sparse = TRUE),
      y,
      control = ctl,
      sigest = sd(y)
    )
    expect_null(
      sparse$setPredictor(
        Matrix::Matrix(plain, sparse = TRUE),
        forceUpdate = TRUE,
        updateCutPoints = "position"
      ),
      info = rule
    )
    expect_identical(cutPointsOf(sparse), dense, info = rule)
    expect_true(strictly(dense), info = rule)
  }
}

# ---- a refresh derives a creation's grid --------------------------------

for (useQuantiles in c(FALSE, TRUE)) {
  rule <- if (useQuantiles) "quantile" else "uniform"
  ctl <- controlWith(20L, useQuantiles)
  createdOn <- function(values) {
    cutPointsOf(dbarts(cbind(z, v = values), y, control = ctl))[[2L]]
  }
  cases <- list(
    placeholder = list(rep(0, n), w),
    fewToMany = list(six, w),
    manyToFew = list(w, six),
    manyToMany = list(w, exp(w))
  )
  for (case in names(cases)) {
    from <- cases[[case]][[1L]]
    onto <- cases[[case]][[2L]]
    for (word in refreshWords) {
      info <- paste(rule, case, word)
      sampler <- dbarts(cbind(z, v = from), y, control = ctl)
      expect_null(
        sampler$setPredictor(
          onto,
          2L,
          forceUpdate = TRUE,
          updateCutPoints = word
        ),
        info = info
      )
      expect_identical(cutPointsOf(sampler)[[2L]], createdOn(onto), info = info)
    }
    # without a refresh the grid is kept
    sampler <- dbarts(cbind(z, v = from), y, control = ctl)
    expect_null(sampler$setPredictor(onto, 2L, forceUpdate = TRUE), info = case)
    expect_identical(cutPointsOf(sampler)[[2L]], createdOn(from), info = case)
    expect_null(
      sampler$setPredictor(
        onto,
        2L,
        forceUpdate = TRUE,
        updateCutPoints = "none"
      ),
      info = case
    )
    expect_identical(cutPointsOf(sampler)[[2L]], createdOn(from), info = case)
  }
  # over a grid set at 60 points
  sampler <- dbarts(cbind(z, v = w), y, control = ctl)
  sampler$setCutPoints(seq(-2, 2, length.out = 60L), 2L)
  expect_null(
    sampler$setPredictor(
      exp(w),
      2L,
      forceUpdate = TRUE,
      updateCutPoints = "position"
    ),
    info = rule
  )
  expect_identical(cutPointsOf(sampler)[[2L]], createdOn(exp(w)), info = rule)
}

# ---- a shrink -------------------------------------------------------------

warmed <- function(x, response = y, control = controlWith(n.trees = 20L), ...) {
  sampler <- dbarts(x, response, control = control, ...)
  invisible(sampler$run(20L, 5L))
  sampler$storeState()
  sampler
}
for (word in refreshWords) {
  # unforced, the trees cannot hold four points where they held a hundred:
  # the call declines and the sampler is its untouched twin's equal
  sampler <- warmed(cbind(z, w))
  twin <- warmed(cbind(z, w))
  expect_true(nrow(splitsOn(sampler, 1L)) > 3L, info = word)
  expect_identical(
    sampler$setPredictor(narrow, 1L, updateCutPoints = word),
    FALSE,
    info = word
  )
  expect_identical(sampler$data@x, twin$data@x, info = word)
  expect_identical(cutPointsOf(sampler), cutPointsOf(twin), info = word)
  expect_identical(sampler$getTrees(), twin$getTrees(), info = word)
  expect_identical(printedTrees(sampler), printedTrees(twin), info = word)
  # a rollback restores the state exactly and the next draws to rounding, as
  # on the build before this rule: a sweep sums a leaf's rows in the order
  # the leaf holds them
  expect_equal(sampler$run(0L, 3L)$train, twin$run(0L, 3L)$train, info = word)

  # forced, the grid is the four points, the sampler runs on, and its state
  # goes back in as stored, into itself and into a copy and a reload
  sampler <- warmed(cbind(z, w))
  expect_null(
    sampler$setPredictor(
      narrow,
      1L,
      forceUpdate = TRUE,
      updateCutPoints = word
    ),
    info = word
  )
  expect_identical(lengths(cutPointsOf(sampler)), c(4L, 100L), info = word)
  expect_true(all(is.finite(sampler$run(0L, 3L)$train)), info = word)
  sampler$storeState()
  stored <- sampler$state
  printed <- printedTrees(sampler)
  expect_true(sampler$setState(stored), info = word)
  expect_identical(printedTrees(sampler), printed, info = word)
  expect_identical(printedTrees(sampler$copy()), printed, info = word)
}

# ---- restores are exact where grids once repeated -----------------------

narrowNA <- replace(narrow, 1:40, NA)
restoreCases <- list(
  narrow = narrow,
  narrowNA = narrowNA,
  constNA = constNA
)
for (case in names(restoreCases)) {
  column <- restoreCases[[case]]
  # a response the column explains, so trees split on it
  signal <- if (case == "constNA") {
    2 * is.na(column)
  } else {
    2 * ((column - 1) / eps >= 2)
  }
  signal[is.na(signal)] <- -2
  response <- signal + rnorm(n, 0, 0.3)
  for (useQuantiles in c(FALSE, TRUE)) {
    info <- paste(case, if (useQuantiles) "quantile" else "uniform")
    sampler <- warmed(
      cbind(z, column),
      response,
      controlWith(100L, useQuantiles, 20L)
    )
    expect_true(strictly(cutPointsOf(sampler)), info = info)
    expect_true(nrow(splitsOn(sampler, 2L)) > 0L, info = info)
    stored <- sampler$state
    printed <- printedTrees(sampler)

    expect_true(sampler$setState(stored), info = info)
    expect_identical(printedTrees(sampler), printed, info = info)
    sampler$storeState()
    expect_identical(sampler$state, stored, info = info)

    copied <- sampler$copy()
    expect_identical(printedTrees(copied), printed, info = info)
    path <- tempfile(fileext = ".rds")
    saveRDS(sampler, path)
    reloaded <- readRDS(path)
    unlink(path)
    expect_identical(printedTrees(reloaded), printed, info = info)
    expect_identical(
      copied$run(0L, 3L)$train,
      reloaded$run(0L, 3L)$train,
      info = info
    )
  }
}

# ---- a grid that repeats a point is refused -----------------------------

stateRefusal <- paste(
  "cut points of column 2 in bartcore state repeat a value",
  "(0.5): a cut grid holds each point once"
)
donorRefusal <- paste(
  "cut points of column 2 in warm-start donor repeat a value",
  "(0.5): a cut grid holds each point once"
)
cutsRefusal <- "a cut point may appear only once in 'cuts'"
ctl <- controlWith(20L, n.trees = 10L)
sampler <- warmed(cbind(z, w), control = ctl)
twin <- warmed(cbind(z, w), control = ctl)
good <- sampler$state
sampler$setCutPoints(c(0.25, 0.5, 0.75), 2L)
twin$setCutPoints(c(0.25, 0.5, 0.75), 2L)
sampler$storeState()
good <- sampler$state
doubled <- good
attr(doubled, "cutPoints")[[2L]] <- c(0.25, 0.5, 0.5, 0.75)

expect_error(sampler$setState(doubled), pattern = stateRefusal, fixed = TRUE)
# a copy and a reload install the state the field holds, and refuse it
sampler$state <- doubled
expect_error(sampler$copy(), pattern = stateRefusal, fixed = TRUE)
path <- tempfile(fileext = ".rds")
saveRDS(sampler, path)
reloaded <- readRDS(path)
unlink(path)
expect_error(reloaded$run(0L, 1L), pattern = stateRefusal, fixed = TRUE)
# a warm start from it is refused in the donor's own words
recipient <- dbarts(cbind(z, w), y, control = ctl)
recipientTwin <- dbarts(cbind(z, w), y, control = ctl)
expect_error(
  recipient$installTrees(doubled),
  pattern = donorRefusal,
  fixed = TRUE
)
expect_identical(recipient$run(0L, 3L), recipientTwin$run(0L, 3L))
expect_error(
  bart(
    cbind(z, w),
    y,
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    n.samples = 5L,
    n.burn = 0L,
    verbose = FALSE,
    warm.start = doubled
  ),
  pattern = donorRefusal,
  fixed = TRUE
)
sampler$state <- good

# setCutPoints refuses a repeat as it did, and takes the grid the column holds
expect_error(
  sampler$setCutPoints(c(0, 0.5, 0.5, 1), 2L),
  pattern = cutsRefusal,
  fixed = TRUE
)
expect_identical(cutPointsOf(sampler)[[2L]], c(0.25, 0.5, 0.75))
# every refusal left the sampler as it was: it draws what its twin draws
expect_identical(sampler$run(0L, 3L), twin$run(0L, 3L))
expect_silent(sampler$setCutPoints(cutPointsOf(sampler)[[2L]], 2L))
expect_identical(cutPointsOf(sampler)[[2L]], c(0.25, 0.5, 0.75))

# ---- a stored tree that stacks two splits on one value ----------------

# tree 1 of the state: a root on column 2 sending missing values left, its
# left child on the same value sending them right, three leaves. Beside
# missing values both sides of the child hold rows, so no side is empty.
stackedState <- function(sampler, childCut) {
  sampler$storeState()
  state <- sampler$state
  cuts <- attr(state, "cutPoints")[[2L]]
  cut <- cuts[length(cuts) %/% 2L]
  forest <- state[[1L]]$forests[[1L]]
  numAfter <- length(forest$tree.sizes) - 1L
  forest$tree.vars <- c(c(2L, 2L, -1L, -1L, -1L), rep(-1L, numAfter))
  forest$tree.values <- c(
    writeBin(c(cut, childCut(cuts), 0.01, -0.01, 0.02), raw()),
    writeBin(numeric(numAfter), raw())
  )
  forest$tree.sizes <- c(5L, rep(1L, numAfter))
  # 2 tags a threshold split, 1 sends missing values right
  forest$tree.flags <- as.raw(c(c(2L, 3L, 0L, 0L, 0L), integer(numAfter)))
  state[[1L]]$forests[[1L]] <- forest
  state
}
wNA <- replace(w, seq(1L, n, by = 4L), NA)
sampler <- dbarts(cbind(z, wNA), y, control = controlWith(20L, n.trees = 5L))
ordered <- stackedState(sampler, function(cuts) cuts[length(cuts) %/% 4L])
expect_true(sampler$setState(ordered))
expect_identical(sum(sampler$getTrees()$var == 2L), 2L)
stacked <- stackedState(sampler, function(cuts) cuts[length(cuts) %/% 2L])
expect_identical(sampler$setState(stacked), FALSE)
trees <- sampler$getTrees()
expect_identical(sum(trees$var == 2L), 1L)
expect_identical(nrow(trees[trees$tree == 1L, ]), 3L)
expect_true(all(is.finite(sampler$run(0L, 5L)$train)))
sampler$storeState()
expect_true(sampler$setState(sampler$state))

# ---- updateCutPoints: three words, and a logical for one more release ----

wordRefusal <- paste(
  "'updateCutPoints' must be one of",
  "\"none\", \"position\", \"value\""
)
sampler <- warmed(cbind(z, w))
for (bad in list(
  NA,
  NA_character_,
  1,
  c("position", "value"),
  "refresh",
  "",
  NULL
)) {
  expect_error(
    sampler$setPredictor(w, 2L, updateCutPoints = bad),
    pattern = wordRefusal,
    fixed = TRUE
  )
}
# a unique abbreviation is taken, as match.arg takes one
abbreviated <- warmed(cbind(z, w))
spelled <- warmed(cbind(z, w))
for (pair in list(c("p", "position"), c("val", "value"), c("no", "none"))) {
  expect_null(
    abbreviated$setPredictor(
      exp(w),
      2L,
      forceUpdate = TRUE,
      updateCutPoints = pair[1L]
    ),
    info = pair[2L]
  )
  spelled$setPredictor(
    exp(w),
    2L,
    forceUpdate = TRUE,
    updateCutPoints = pair[2L]
  )
  expect_identical(
    cutPointsOf(abbreviated),
    cutPointsOf(spelled),
    info = pair[2L]
  )
  expect_identical(abbreviated$getTrees(), spelled$getTrees(), info = pair[2L])
}
# a row-by-row update keeps the grid, so it takes "none" alone
for (word in refreshWords) {
  expect_error(
    sampler$setPredictor(
      w,
      2L,
      forceUpdate = "partial",
      updateCutPoints = word
    ),
    pattern = "partial updates cannot also update cut points"
  )
}
expect_true(is.logical(
  sampler$setPredictor(w, 2L, forceUpdate = "partial", updateCutPoints = "none")
))

# a logical warns once in a session, whichever value comes first, and is
# read as the word it stood for: TRUE draws what "position" draws
warnEnv <- dbarts:::onceWarnState
warnEnv[["tombstone.updateCutPoints.logical"]] <- NULL
logical <- warmed(cbind(z, w))
worded <- warmed(cbind(z, w))
seen <- warningsOf({
  viaTrue <- logical$setPredictor(
    exp(w),
    2L,
    forceUpdate = TRUE,
    updateCutPoints = TRUE
  )
  viaFalse <- logical$setPredictor(
    w,
    2L,
    forceUpdate = TRUE,
    updateCutPoints = FALSE
  )
  again <- logical$setPredictor(
    exp(w),
    2L,
    forceUpdate = TRUE,
    updateCutPoints = TRUE
  )
})
expect_identical(length(seen), 1L)
expect_true(inherits(seen[[1L]], "dbartsDeprecatedWarning"))
for (fragment in c(
  "\"none\"",
  "\"position\"",
  "\"value\"",
  "1.1-0",
  "updateCutPoints"
)) {
  expect_true(
    grepl(fragment, conditionMessage(seen[[1L]]), fixed = TRUE),
    info = fragment
  )
}
# forced, each returns NULL
expect_identical(list(viaTrue, viaFalse, again), list(NULL, NULL, NULL))
worded$setPredictor(
  exp(w),
  2L,
  forceUpdate = TRUE,
  updateCutPoints = "position"
)
worded$setPredictor(w, 2L, forceUpdate = TRUE, updateCutPoints = "none")
worded$setPredictor(
  exp(w),
  2L,
  forceUpdate = TRUE,
  updateCutPoints = "position"
)
expect_identical(cutPointsOf(logical), cutPointsOf(worded))
expect_identical(logical$run(0L, 3L), worded$run(0L, 3L))
# the words never warn
warnEnv[["tombstone.updateCutPoints.logical"]] <- NULL
expect_identical(
  length(warningsOf(
    for (word in c("none", refreshWords)) {
      worded$setPredictor(w, 2L, forceUpdate = TRUE, updateCutPoints = word)
    }
  )),
  0L
)

# ---- where the splits go ---------------------------------------------------

# a column moved up by three and a half points' spacing. The grid moves with
# it, so under "position" every threshold moves with the column; under
# "value" every threshold stays as near where it was as the grid allows
positioned <- warmed(cbind(z, w))
valued <- warmed(cbind(z, w))
before <- splitsOn(positioned, 2L)
oldGrid <- cutPointsOf(positioned)[[2L]]
spacing <- oldGrid[2L] - oldGrid[1L]
shifted <- w + 3.5 * spacing
expect_true(nrow(before) > 3L)
expect_null(positioned$setPredictor(
  shifted,
  2L,
  forceUpdate = TRUE,
  updateCutPoints = "position"
))
afterPosition <- splitsOn(positioned, 2L)
newGrid <- cutPointsOf(positioned)[[2L]]
expect_identical(afterPosition[c("tree", "n")], before[c("tree", "n")])
expect_identical(
  match(afterPosition$value, newGrid),
  match(before$value, oldGrid)
)
expect_null(valued$setPredictor(
  shifted,
  2L,
  forceUpdate = TRUE,
  updateCutPoints = "value"
))
expect_identical(cutPointsOf(valued)[[2L]], newGrid)
afterValue <- splitsOn(valued, 2L)
expect_true(all(afterValue$value %in% newGrid))
# the lowest threshold on the column cannot have been held back by a split
# above it in its tree: it sits on the new point nearest where it was
lowest <- min(before$value)
expect_identical(
  min(afterValue$value),
  newGrid[which.min(abs(newGrid - lowest))]
)
expect_true(abs(min(afterValue$value) - lowest) <= spacing / 2)
expect_true(abs(min(afterPosition$value) - lowest) > 3 * spacing)

# ---- setCutPoints takes the same choice ---------------------------------

expect_error(
  sampler$setCutPoints(oldGrid, 2L, splits = "none"),
  pattern = "'splits' must be one of \"position\", \"value\"",
  fixed = TRUE
)
# a grid holding every old point and one more above each: by value each
# threshold stays on its point, no row moves and nothing changes but the
# grid; by position each split moves up to the added point above it
finer <- sort(c(oldGrid, oldGrid + spacing / 2))
byValue <- warmed(cbind(z, w))
byPosition <- warmed(cbind(z, w))
untouched <- warmed(cbind(z, w))
before <- splitsOn(untouched, 2L)
byValue$setCutPoints(finer, 2L, splits = "value")
expect_identical(cutPointsOf(byValue)[[2L]], finer)
expect_identical(splitsOn(byValue, 2L), before)
expect_identical(byValue$getTrees(), untouched$getTrees())
byPosition$setCutPoints(finer, 2L)
moved <- splitsOn(byPosition, 2L)
expect_identical(moved$tree, before$tree)
expect_identical(match(moved$value, finer), 2L * match(before$value, oldGrid))
# by position, a grid of the held count keeps every split where it is on the
# grid, whatever the points: the default, and what the column did before
sameCount <- warmed(cbind(z, w))
sameCount$setCutPoints(oldGrid + spacing / 4, 2L, splits = "position")
kept <- splitsOn(sameCount, 2L)
expect_identical(kept$tree, before$tree)
expect_identical(
  match(kept$value, oldGrid + spacing / 4),
  match(before$value, oldGrid)
)
# and a shorter grid rescales: on half the points no split is past the end
shorter <- warmed(cbind(z, w))
shorter$setCutPoints(oldGrid[seq(2L, 100L, by = 2L)], 2L)
expect_true(all(
  splitsOn(shorter, 2L)$value %in% oldGrid[seq(2L, 100L, by = 2L)]
))
expect_true(all(is.finite(shorter$run(0L, 3L)$train)))
shorter$storeState()
expect_true(shorter$setState(shorter$state))

# ---- an unforced refresh moves the splits of every forest -----------------

# Two forests on one design. Column 2 is set to every other point of its grid
# and then refreshed, unforced, onto the values it holds: the grid is the
# hundred points again, and place i of fifty goes to place 2 i of the
# hundred, the point the split was on. So every tree of both forests is what
# it was, which it is not where a forest's splits were left on their old
# places.
arm <- seq_len(n) %% 2L
twoForests <- dbarts(
  cbind(z, w),
  y,
  forests = list(forest(), forest(basis = ~ factor(arm))),
  control = controlWith(n.trees = 20L)
)
invisible(twoForests$run(20L, 5L))
created <- cutPointsOf(twoForests)[[2L]]
twoForests$setCutPoints(created[seq(2L, 100L, by = 2L)], 2L, splits = "value")
before <- twoForests$getTrees()
expect_true(all(c(1L, 2L) %in% before$forest[before$var == 2L]))
expect_true(twoForests$setPredictor(w, 2L, updateCutPoints = "position"))
expect_identical(cutPointsOf(twoForests)[[2L]], created)
expect_identical(twoForests$getTrees(), before)

# ---- a forced merge weighs each leaf by the rows it holds -----------------

# Given one point, a column keeps the first split on it along a path, moved
# to the point, and a split on it beneath another has no point left: it is
# merged with everything under it. The leaf that results is the mean of the
# leaves it replaces, each weighed by the case weights of the rows it held.
# Worked out here from the trees as they were, for a sampler and for its
# copy, whose node statistics no sweep has filled.
caseWeights <- rep(c(0.2, 1, 3, 1), n / 4L)
design <- cbind(z, w)
nested <- function(nodes) {
  at <- 0L
  build <- function() {
    at <<- at + 1L
    node <- list(var = nodes$var[at], value = nodes$value[at])
    if (node$var > 0L) {
      node$left <- build()
      node$right <- build()
    }
    node
  }
  build()
}
flattened <- function(node) {
  here <- data.frame(var = node$var, value = node$value)
  if (node$var < 0L) {
    return(here)
  }
  rbind(here, flattened(node$left), flattened(node$right))
}
leavesUnder <- function(node, rows) {
  if (node$var < 0L) {
    return(data.frame(value = node$value, weight = sum(caseWeights[rows])))
  }
  goesLeft <- design[rows, node$var] <= node$value
  rbind(
    leavesUnder(node$left, rows[goesLeft]),
    leavesUnder(node$right, rows[!goesLeft])
  )
}
onOnePoint <- function(node, rows, point, beneath = FALSE) {
  if (node$var < 0L) {
    return(node)
  }
  if (node$var == 1L && beneath) {
    leaves <- leavesUnder(node, rows)
    return(list(
      var = -1L,
      value = sum(leaves$weight * leaves$value) / sum(leaves$weight)
    ))
  }
  goesLeft <- design[rows, node$var] <= node$value
  if (node$var == 1L) {
    node$value <- point
    beneath <- TRUE
  }
  node$left <- onOnePoint(node$left, rows[goesLeft], point, beneath)
  node$right <- onOnePoint(node$right, rows[!goesLeft], point, beneath)
  node
}
for (door in c("sampler", "copy")) {
  sampler <- warmed(design, weights = caseWeights)
  if (door == "copy") {
    sampler <- sampler$copy()
  }
  before <- sampler$getTrees()
  sampler$setCutPoints(0, 1L)
  after <- sampler$getTrees()
  numMerged <- 0L
  for (tree in unique(before$tree)) {
    held <- before[before$tree == tree, ]
    byHand <- flattened(onOnePoint(nested(held), seq_len(n), 0))
    left <- after[after$tree == tree, ]
    expect_identical(left$var, byHand$var, info = paste(door, tree))
    expect_equal(
      left$value,
      byHand$value,
      tolerance = 1e-12,
      info = paste(door, tree)
    )
    numMerged <- numMerged + (nrow(byHand) < nrow(held))
  }
  expect_true(numMerged >= 3L, info = door)
}

# ---- the default rule weighs one cut per gap ------------------------------

# A column with fewer distinct values than n.cuts takes the quantile rule's
# grid, a point halfway along each gap, and each gap is chosen with
# probability in proportion to its width. The weights ride the state as the
# column's sorted values; a column whose widths are equal carries none.
massOf <- function(sampler) {
  sampler$storeState()
  attr(sampler$state, "cutMass")
}
squares <- as.double((seq_len(n) %% 6L)^2) # widths 1, 3, 5, 7, 9
sampler <- warmed(cbind(z, squares, binary), y)
grid <- cutPointsOf(sampler)
mass <- massOf(sampler)
expect_identical(grid[[2L]], c(0.5, 2.5, 6.5, 12.5, 20.5))
expect_identical(grid[[3L]], 0.5)
expect_identical(mass, list(NULL, c(0, 1, 4, 9, 16, 25), NULL))
quantiled <- dbarts(
  cbind(z, squares, binary),
  y,
  control = controlWith(
    useQuantiles = TRUE
  )
)
expect_identical(cutPointsOf(quantiled)[2L:3L], grid[2L:3L])
expect_null(massOf(quantiled))

# a 0/1 column holds its one point whatever n.cuts asks: the draws do not move
onBinary <- function(n.cuts) {
  dbarts(cbind(binary), y, control = controlWith(n.cuts))$run(5L, 5L)
}
expect_identical(onBinary(3L), onBinary(100L))

# a state, a copy and a reload carry the weights and continue the chain
stored <- sampler$state
copied <- sampler$copy()
path <- tempfile(fileext = ".rds")
saveRDS(sampler, path)
reloaded <- readRDS(path)
unlink(path)
expect_identical(massOf(copied), mass)
expect_identical(massOf(reloaded), mass)
expect_identical(copied$run(0L, 3L)$train, reloaded$run(0L, 3L)$train)
expect_true(sampler$setState(stored))

# a state without them installs the grids unweighted, and a malformed one is
# refused by column
bare <- stored
attr(bare, "cutMass") <- NULL
expect_true(sampler$setState(bare))
expect_null(massOf(sampler))
expect_true(sampler$setState(stored))
expect_identical(massOf(sampler), mass)
malformed <- stored
attr(malformed, "cutMass")[[2L]] <- c(0, 1, 1, 9, 16, 25)
expect_error(
  sampler$setState(malformed),
  pattern = paste(
    "cut weights of column 2 in bartcore state are not 6 finite values",
    "strictly increasing"
  ),
  fixed = TRUE
)
expect_identical(massOf(sampler), mass)

# setCutPoints keeps them for the grid the column holds and drops them for
# another; a refresh derives them again, and "none" keeps the column's
sampler$setCutPoints(grid[[2L]], 2L)
expect_identical(massOf(sampler), mass)
sampler$setCutPoints(grid[[2L]] + 0.25, 2L)
expect_null(massOf(sampler))
expect_null(sampler$setPredictor(
  squares,
  2L,
  forceUpdate = TRUE,
  updateCutPoints = "value"
))
expect_identical(massOf(sampler), mass)
expect_null(sampler$setPredictor(
  w,
  2L,
  forceUpdate = TRUE,
  updateCutPoints = "none"
))
expect_identical(cutPointsOf(sampler), grid)
expect_identical(massOf(sampler), mass)
expect_null(sampler$setPredictor(
  w,
  2L,
  forceUpdate = TRUE,
  updateCutPoints = "position"
))
expect_null(massOf(sampler))
expect_identical(length(cutPointsOf(sampler)[[2L]]), 100L)
