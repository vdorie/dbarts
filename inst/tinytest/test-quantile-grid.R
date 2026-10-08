# With useQuantiles = TRUE a numeric column's split points are midpoints
# between its sorted distinct values. A column with no more distinct values
# than the cut count plus one takes every midpoint; a longer one takes the cut
# count of them, spread evenly over all the midpoints, so both ends of the
# column are reached. Every entry point builds that one grid, and a refresh
# through setPredictor derives it again from the new values, counting from
# n.cuts whatever number of points the column held.

cutPointsOf <- function(sampler) {
  sampler$storeState()
  attr(sampler$state, "cutPoints")
}
allMidpoints <- function(values) {
  u <- sort(unique(values))
  (u[-1L] + u[-length(u)]) / 2
}
# cut k of count, k from zero, is midpoint floor((2k + 1) M / (2 count)) of M
spreadMidpoints <- function(values, count) {
  mid <- allMidpoints(values)
  count <- min(count, length(mid))
  mid[((2 * seq_len(count) - 1) * length(mid)) %/% (2 * count) + 1L]
}
warningsOf <- function(expr) {
  seen <- character(0L)
  withCallingHandlers(
    expr,
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seen
}

set.seed(8125L)
n <- 199L
few <- rnorm(150L)
x <- cbind(
  x1 = rnorm(n), # 199 distinct values
  x2 = c(few, sample(few, n - 150L, replace = TRUE)), # 150
  x3 = rep(seq_len(101L) / 101, length.out = n), # the cut count plus one
  x4 = rep(seq_len(11L), length.out = n)
)
y <- x[, 1L] + sin(3 * x[, 2L]) + rnorm(n, 0, 0.2)
expect_identical(
  unname(apply(x, 2L, function(v) length(unique(v)))),
  c(199L, 150L, 101L, 11L)
)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 5L,
  useQuantiles = TRUE,
  updateState = FALSE
)
expect_identical(control@n.cuts, 100L)
sampler <- dbarts(x, y, control = control)
cuts <- cutPointsOf(sampler)
expect_identical(lengths(cuts), c(100L, 100L, 100L, 10L))

# the long columns: no part of the range is left without a split point
for (j in 1L:2L) {
  expect_true(max(cuts[[j]]) > quantile(x[, j], 0.98, names = FALSE))
  expect_true(min(cuts[[j]]) < quantile(x[, j], 0.02, names = FALSE))
  expect_identical(cuts[[j]], spreadMidpoints(x[, j], 100L))
  expect_true(all(diff(cuts[[j]]) > 0))
  expect_false(any(cuts[[j]] %in% x[, j]))
}
# the short ones: every midpoint, as before
expect_identical(cuts[[3L]], allMidpoints(x[, 3L]))
expect_identical(cuts[[4L]], seq_len(10L) + 0.5)

# the doors that return a sampler build the same grid: the modern one, the
# BayesTree-style one, a BayesTree-spelled call to the modern one, which is
# forwarded, the grouped one and the two partial-dependence ones
onceState <- dbarts:::onceWarnState
# evaluates expr with the named once-per-session notices set to value, then
# puts them back; the notices themselves are pinned in their own files
withNotices <- function(keys, value, expr) {
  saved <- mget(keys, envir = onceState, ifnotfound = list(NULL))
  on.exit(
    for (key in keys) {
      onceState[[key]] <- saved[[key]]
    }
  )
  for (key in keys) {
    onceState[[key]] <- value
  }
  expr
}
frame <- data.frame(x, y = y, g = rep(seq_len(4L), length.out = n))
modernArgs <- list(
  y ~ x1 + x2 + x3 + x4,
  frame,
  useQuantiles = TRUE,
  n.trees = 5L,
  n.samples = 4L,
  n.burn = 2L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
legacyArgs <- list(
  x.train = x,
  y.train = y,
  usequants = TRUE,
  ntree = 5L,
  ndpost = 4L,
  nskip = 2L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  keepsampler = TRUE
)
fitModern <- do.call(dbarts::bart, c(modernArgs, list(keepSampler = TRUE)))
fitLegacy <- do.call(dbarts::bartBT, legacyArgs)
shimWarnings <- warningsOf(withNotices(
  "tombstone.bartShim",
  NULL,
  fitShim <- do.call(dbarts::bart, legacyArgs)
))
expect_identical(length(shimWarnings), 1L)
expect_true(grepl("'bartBT'", shimWarnings[1L], fixed = TRUE))
expect_identical(dbarts:::callName(fitShim$call), "bartBT")
otherWarnings <- warningsOf(withNotices(
  c("tombstone.rbart_vi", "tombstone.pdbartDefaultsMessage"),
  TRUE,
  {
    fitGrouped <- do.call(
      dbarts::rbart_vi,
      c(modernArgs, list(group.by = quote(g), n.thin = 1L, keepSampler = TRUE))
    )
    pd <- do.call(dbarts::pdbart, c(modernArgs, list(xind = "x1", pl = FALSE)))
    pd2 <- do.call(
      dbarts::pd2bart,
      c(modernArgs, list(xind = c("x1", "x2"), pl = FALSE))
    )
  }
))
expect_identical(length(otherWarnings), 0L)
doors <- list(
  bart = fitModern$fit,
  bartBT = fitLegacy$fit,
  forwarded = fitShim$fit,
  rbart_vi = fitGrouped$fit[[1L]],
  pdbart = pd$fit,
  pd2bart = pd2$fit
)
for (door in names(doors)) {
  expect_true(doors[[door]]$control@useQuantiles, info = door)
  expect_identical(
    unname(cutPointsOf(doors[[door]])),
    unname(cuts),
    info = door
  )
}

# a refresh spreads n.cuts over the new values' midpoints, whatever count the
# column holds, and raises nothing: 10 cuts from 11 distinct values, refreshed
# onto 60 under n.cuts = 10 and, below, under n.cuts = 100
xRefresh <- cbind(rep(seq_len(11L), length.out = n), x[, 1L])
controlTen <- control
controlTen@n.cuts <- 10L
sampler <- dbarts(xRefresh, y, control = controlTen)
expect_identical(cutPointsOf(sampler)[[1L]], seq_len(10L) + 0.5)
replacement <- rep(seq_len(60L), length.out = n)
refreshWarnings <- warningsOf(
  sampler$setPredictor(
    replacement,
    1L,
    forceUpdate = TRUE,
    updateCutPoints = "position"
  )
)
expect_identical(length(refreshWarnings), 0L)
refreshed <- cutPointsOf(sampler)
expect_identical(refreshed[[1L]], spreadMidpoints(replacement, 10L))
expect_identical(refreshed[[1L]], c(3, 9, 15, 21, 27, 33, 39, 45, 51, 57) + 0.5)
expect_identical(refreshed[[2L]], spreadMidpoints(x[, 1L], 10L))
# under n.cuts = 100 the same refresh takes every midpoint of the 60 values,
# where the column held 10 cuts
sampler <- dbarts(xRefresh, y, control = control)
expect_identical(cutPointsOf(sampler)[[1L]], seq_len(10L) + 0.5)
expect_true(sampler$setPredictor(
  replacement,
  1L,
  forceUpdate = TRUE,
  updateCutPoints = "position"
))
expect_identical(cutPointsOf(sampler)[[1L]], seq_len(59L) + 0.5)
# fewer distinct values than cuts held shrink the grid to their midpoints
expect_true(sampler$setPredictor(
  rep(seq_len(5L), length.out = n),
  1L,
  forceUpdate = TRUE,
  updateCutPoints = "position"
))
expect_identical(cutPointsOf(sampler)[[1L]], seq_len(4L) + 0.5)

rm(sampler, fitModern, fitLegacy, fitShim, fitGrouped, pd, pd2, doors)
