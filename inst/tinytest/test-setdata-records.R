# setData and setPredictor under the record a data object carries. setData
# follows updateCutPoints, updateScale and updateStandardization where given,
# deriving the value from the new data or keeping the sampler's; left NULL it
# uses a value the data object records, whoever wrote it (dec-B458), and
# derives one it does not, so setData(sampler$data) after an edit holds the
# grid, the range and the standardization, and a fresh dbartsData() derives
# them as before. setPredictor(updateStandardization = TRUE) derives a linear
# or gp leaf's standardization over the values the update leaves, a linear
# leaf's coefficients rewritten into it. A gp sampler holding kept draws is
# refused only a call that would change its standardization or its range.

set.seed(71L, sample.kind = "Rejection")
n <- 120L
x <- cbind(a = runif(n), b = runif(n))
y <- sin(4 * x[, 1L]) + x[, 2L] + rnorm(n, 0, 0.2)
x.probe <- cbind(a = runif(30L), b = runif(30L))
# other values for the same predictors
x.other <- cbind(a = 2 * x[, 1L], b = 3 * x[, 2L] + 5)
y.other <- 3 * y + 10
controlWith <- function(...) {
  settings <- list(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 15L,
    n.samples = 4L,
    n.burn = 0L,
    keepTrees = TRUE,
    updateState = FALSE,
    seed = 7L
  )
  do.call(dbarts::dbartsControl, utils::modifyList(settings, list(...)))
}
# a sampler of linear leaves on column b unless another leaf prior is named,
# as written or quoted in a variable, run so that it holds kept draws
make <- function(leaf.prior, data = x, response = y, ...) {
  prior <- if (missing(leaf.prior)) {
    quote(linear(columns = "b"))
  } else {
    substitute(leaf.prior)
  }
  if (is.name(prior)) {
    prior <- leaf.prior
  }
  sampler <- eval(bquote(dbarts::dbarts(
    data,
    response,
    leaf.prior = .(prior),
    control = controlWith(...)
  )))
  invisible(sampler$run(10L, 4L))
  sampler
}
readerOf <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_derivedRecords, sampler$getPointer(), TRUE)
}
# the record, the covariates' names left off as the reader leaves them
recordOf <- function(sampler) {
  standardization <- sampler$data@standardization
  if (!is.null(standardization)) {
    standardization[[1L]] <- lapply(standardization[[1L]], unname)
  }
  list(
    cutPoints = sampler$data@cutPoints,
    response.range = sampler$data@response.range,
    forest.calibration = sampler$data@forest.calibration,
    standardization = standardization
  )
}
expectRecordIsReader <- function(sampler, info) {
  expect_identical(recordOf(sampler), readerOf(sampler), info = info)
}
stateOf <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_storeState, sampler$getPointer())
}
kept <- function(sampler) sampler$predict(x.probe)
liveFits <- function(sampler) as.vector(sampler$getFitsWithoutOffset())
worstGap <- function(a, b) {
  max(abs(unlist(a) - unlist(b)) / pmax(1, abs(unlist(a))))
}
unnamedStandardization <- function(record) {
  lapply(record$standardization[[1L]], unname)
}

# ---- fresh object neutral --------------------------------------------------

# a data object that carries no record derives every value, exactly as the
# three arguments asking for the derivation do: the state, the record and the
# next draws
for (prior in c("linear", "normal")) {
  leaf <- if (prior == "linear") {
    quote(linear(columns = "b"))
  } else {
    quote(normal())
  }
  absent <- make(leaf)
  explicit <- make(leaf)
  absent$setData(dbarts::dbartsData(x.other, y.other))
  explicit$setData(
    dbarts::dbartsData(x.other, y.other),
    updateCutPoints = "value",
    updateScale = TRUE,
    updateStandardization = TRUE
  )
  expect_identical(stateOf(absent), stateOf(explicit), info = prior)
  expect_identical(recordOf(absent), recordOf(explicit), info = prior)
  expect_identical(absent$run(0L, 5L), explicit$run(0L, 5L), info = prior)
  # and the derivation is the one a sampler made on those data holds
  made <- make(leaf, data = x.other, response = y.other)
  expect_identical(recordOf(absent), recordOf(made), info = prior)
}

# ---- carried held ----------------------------------------------------------

# the sampler's own data object, edited: the grid and the range are held, by
# the reader, and the kept draws are not rewritten
sampler <- make()
before <- readerOf(sampler)
keptBefore <- kept(sampler)
edited <- sampler$data
edited@x[, 1L] <- 2 * x[, 1L]
sampler$setData(edited)
expect_identical(readerOf(sampler)$cutPoints, before$cutPoints)
expect_identical(readerOf(sampler)$standardization, before$standardization)
expect_identical(kept(sampler), keptBefore)
edited <- sampler$data
edited@y <- 3 * y + 10
sampler$setData(edited)
expect_identical(readerOf(sampler)$response.range, before$response.range)
expect_identical(kept(sampler), keptBefore)
expectRecordIsReader(sampler, "carried held")
# a count sampler's center is held as its range, and a mask its 0/1 weights
# install does not read it again
counts <- rpois(n, exp(1 + x[, 1L]))
countSampler <- dbarts::dbarts(
  x,
  counts,
  family = "nbinom",
  control = controlWith()
)
center <- countSampler$data@response.range
edited <- countSampler$data
edited@y <- 2 * counts
countSampler$setData(edited)
expect_identical(countSampler$data@response.range, center)
edited@weights <- as.double(x[, 1L] > 0.5)
countSampler$setData(edited)
expect_identical(countSampler$data@response.range, center)
expect_identical(countSampler$activeRows, edited@weights)
countSampler$setData(edited, updateScale = TRUE)
expect_equal(
  countSampler$data@response.range[[1L]],
  log(mean(2 * counts[x[, 1L] > 0.5])),
  tolerance = 1e-14
)
# a count record that is no center's is refused naming the slot; one of
# another width is read by its first number, the center, and recorded as
# the window the family holds
edited <- countSampler$data
edited@response.range <- c(2, 2)
expect_error(
  countSampler$setData(edited),
  paste0(
    "the data's 'response.range' record must be two finite numbers, the ",
    "second above the first; set it to NULL"
  ),
  fixed = TRUE
)
expect_error(
  dbarts::dbarts(edited, family = "nbinom", control = controlWith()),
  "'response.range' record must be two finite numbers, the second above",
  fixed = TRUE
)
edited@response.range <- c(2, 5)
countSampler$setData(edited)
expect_identical(countSampler$data@response.range, c(2, 3))
# a Student-t sampler holds its range as a gaussian one does
student <- dbarts::dbarts(x, y, family = "student", control = controlWith())
invisible(student$run(5L, 4L))
before <- readerOf(student)$response.range
keptBefore <- kept(student)
edited <- student$data
edited@y <- 3 * y + 10
student$setData(edited)
expect_identical(readerOf(student)$response.range, before)
expect_identical(student$data@response.range, before)
expect_identical(kept(student), keptBefore)
# the sampler's own data object handed back moves no stored number: a held
# transform that does not move leaves sigma as it was, which a round trip
# through the range would not on every chain
moved <- 0L
for (seed in 1:20) {
  sampler <- dbarts::dbarts(
    x,
    y,
    control = controlWith(seed = seed, keepTrees = FALSE)
  )
  invisible(sampler$run(10L, 0L))
  before <- stateOf(sampler)
  sampler$setData(sampler$data)
  moved <- moved + !identical(stateOf(sampler), before)
}
expect_identical(moved, 0L)
# an ordered factor's cut points are its levels', so a record giving it
# others is refused naming the slot, as at creation, the sampler untouched
ordered <- data.frame(
  a = x[, 1L],
  o = factor(rep_len(c("l", "m", "h"), n), c("l", "m", "h"), ordered = TRUE)
)
sampler <- dbarts::dbarts(ordered, y, control = controlWith())
before <- stateOf(sampler)
edited <- sampler$data
edited@cutPoints[[2L]] <- 1.5
for (case in c("setData", "setData, none", "creation")) {
  expect_error(
    switch(
      case,
      "setData" = sampler$setData(edited),
      "setData, none" = sampler$setData(edited, updateCutPoints = "none"),
      "creation" = dbarts::dbarts(edited, control = controlWith())
    ),
    paste0(
      "the data's 'cutPoints' record gives predictor 2, an ordered factor, ",
      "other cut points than those of its levels"
    ),
    fixed = TRUE,
    info = case
  )
}
expect_identical(stateOf(sampler), before)
expect_identical(sampler$data@cutPoints[[2L]], c(0.5, 1.5))

# ---- arguments win ---------------------------------------------------------

# each argument at each value, over a data object that carries a record other
# than the sampler's and one that carries none: an argument given is
# followed, deriving or keeping the sampler's, and NULL takes the record or
# derives
derivedOn <- readerOf(make(data = x.other, response = y.other))
carried <- dbarts::dbartsData(x.other, y.other)
carried@cutPoints <- derivedOn$cutPoints
carried@cutPoints[[1L]] <- c(0.5, 1, 1.5)
carried@response.range <- c(-5, 25)
carried@standardization <- list(list(
  center = c(b = 6),
  scale = c(b = 0.5)
))
fresh <- dbarts::dbartsData(x.other, y.other)
arguments <- list(
  cutPoints = list(
    name = "updateCutPoints",
    values = list(NULL, "none", "value")
  ),
  response.range = list(name = "updateScale", values = list(NULL, FALSE, TRUE)),
  standardization = list(
    name = "updateStandardization",
    values = list(NULL, FALSE, TRUE)
  )
)
recordValue <- function(record, slot) {
  if (slot == "standardization") {
    unnamedStandardization(record)
  } else {
    record[[slot]]
  }
}
for (objectName in c("carried", "fresh")) {
  object <- if (objectName == "carried") carried else fresh
  for (slot in names(arguments)) {
    for (k in seq_along(arguments[[slot]]$values)) {
      value <- arguments[[slot]]$values[[k]]
      info <- paste(objectName, arguments[[slot]]$name, deparse(value))
      sampler <- make()
      own <- readerOf(sampler)
      call <- list(object)
      call[arguments[[slot]]$name] <- list(value)
      do.call(sampler$setData, call)
      expected <- if (k == 2L) {
        own
      } else if (k == 3L || objectName == "fresh") {
        derivedOn
      } else {
        list(
          cutPoints = carried@cutPoints,
          response.range = carried@response.range,
          standardization = carried@standardization
        )
      }
      expect_equal(
        recordValue(readerOf(sampler), slot),
        recordValue(expected, slot),
        tolerance = 0,
        info = info
      )
      expectRecordIsReader(sampler, info)
    }
  }
}

# ---- keep is the sampler's -------------------------------------------------

# "none" and FALSE keep the sampler's own values, never the carried ones
sampler <- make()
own <- readerOf(sampler)
sampler$setData(
  carried,
  updateCutPoints = "none",
  updateScale = FALSE,
  updateStandardization = FALSE
)
expect_identical(readerOf(sampler)[-3L], own[-3L])

# ---- placement -------------------------------------------------------------

# a derived grid: "position" keeps each split's place on its column's grid,
# "value" moves it to the point nearest its threshold; a carried grid other
# than the sampler's moves the splits by value
splitsOn <- function(state, column) {
  forest <- state[[1L]]$forests[[1L]]
  values <- readBin(forest$tree.values, "double", length(forest$tree.vars))
  values[forest$tree.vars == column]
}
doubled <- dbarts::dbartsData(cbind(a = 2 * x[, 1L], b = x[, 2L]), y)
for (rule in c("position", "value")) {
  sampler <- make(normal())
  oldGrid <- sampler$data@cutPoints[[1L]]
  oldSplits <- splitsOn(stateOf(sampler), 1L)
  sampler$setData(doubled, updateCutPoints = rule)
  newGrid <- sampler$data@cutPoints[[1L]]
  newSplits <- splitsOn(stateOf(sampler), 1L)
  expect_identical(length(newGrid), length(oldGrid), info = rule)
  expect_true(length(oldSplits) > 0L, info = rule)
  if (rule == "position") {
    # the same rows doubled: no leaf empties, so every split survives
    expect_identical(newSplits, newGrid[match(oldSplits, oldGrid)])
  } else {
    # each split near its old threshold, inside the old column's range
    expect_true(max(newSplits) <= max(oldGrid) + diff(newGrid[1:2]))
  }
}
# a carried grid of the column's own count with every point moved: by value
# each split takes the point nearest its threshold, where by position it
# would keep its index
sampler <- make(normal(), n.trees = 5L, keepTrees = FALSE)
invisible(sampler$run(10L, 0L))
own <- sampler$data@cutPoints[[1L]]
oldSplits <- splitsOn(stateOf(sampler), 1L)
given <- sampler$data
given@cutPoints[[1L]] <- 0.9 * own + 0.05
grid <- given@cutPoints[[1L]]
nearest <- vapply(oldSplits, function(v) grid[which.min(abs(grid - v))], 0)
expect_true(length(oldSplits) > 0L)
expect_false(identical(nearest, grid[match(oldSplits, own)]))
sampler$setData(given)
expect_identical(splitsOn(stateOf(sampler), 1L), nearest)
expect_identical(sampler$data@cutPoints[[1L]], grid)

# ---- new row count ---------------------------------------------------------

# kept with 40 rows added, the first 40 rows again: every row is read on the
# grid and under the constants in force, so the old rows fit as they did and
# the added ones as the rows they repeat
sampler <- make()
before <- readerOf(sampler)
fitsBefore <- liveFits(sampler)
longer <- dbarts::dbartsData(rbind(x, x[1:40, ]), c(y, y[1:40]))
sampler$setData(
  longer,
  updateCutPoints = "none",
  updateScale = FALSE,
  updateStandardization = FALSE
)
expect_identical(readerOf(sampler)[-3L], before[-3L])
fitsAfter <- liveFits(sampler)
expect_equal(fitsAfter[seq_len(n)], fitsBefore, tolerance = 1e-12)
expect_identical(fitsAfter[n + 1:40], fitsAfter[1:40])
expect_true(all(is.finite(sampler$run(0L, 3L)$train)))
# a gp leaf's covariates are gathered for every row as well: kept draws as
# they were, and the next draws those of a twin handed its own record
gpSampler <- make(gp(columns = "b"))
gpTwin <- make(gp(columns = "b"))
gpKept <- kept(gpSampler)
ownRecord <- longer
ownRecord@cutPoints <- gpTwin$data@cutPoints
ownRecord@response.range <- gpTwin$data@response.range
ownRecord@standardization <- gpTwin$data@standardization
gpSampler$setData(
  longer,
  updateCutPoints = "none",
  updateScale = FALSE,
  updateStandardization = FALSE
)
gpTwin$setData(ownRecord)
expect_identical(kept(gpSampler), gpKept)
expect_identical(length(liveFits(gpSampler)), n + 40L)
expect_identical(gpSampler$run(0L, 3L), gpTwin$run(0L, 3L))

# ---- linear converts -------------------------------------------------------

# updateStandardization = TRUE at both setters: the standardization derived
# over the new values, the kept draws and the live fit the functions they
# were, to rounding
stretched <- 3 * x[, 2L] + 5
for (setter in c("matrix", "column")) {
  sampler <- make()
  twin <- make()
  keptBefore <- kept(sampler)
  if (setter == "matrix") {
    sampler$setPredictor(
      cbind(a = x[, 1L], b = stretched),
      forceUpdate = TRUE,
      updateStandardization = TRUE
    )
    twin$setPredictor(cbind(a = x[, 1L], b = stretched), forceUpdate = TRUE)
  } else {
    sampler$setPredictor(
      stretched,
      "b",
      forceUpdate = TRUE,
      updateStandardization = TRUE
    )
    twin$setPredictor(stretched, "b", forceUpdate = TRUE)
  }
  record <- sampler$data@standardization[[1L]]
  expect_equal(record$center[["b"]], mean(stretched), tolerance = 1e-14)
  expect_equal(record$scale[["b"]], sd(stretched), tolerance = 1e-14)
  expectRecordIsReader(sampler, setter)
  expect_true(worstGap(keptBefore, kept(sampler)) < 1e-12, info = setter)
  expect_true(worstGap(liveFits(twin), liveFits(sampler)) < 1e-12)
  # FALSE, the default, keeps it
  expect_identical(
    twin$data@standardization,
    make()$data@standardization,
    info = setter
  )
}
# TRUE derives the replaced columns' alone: b, replaced and held, keeps the
# center and scale it holds when a is replaced and derived
sampler <- make(linear(columns = c("a", "b")))
keptBefore <- kept(sampler)
sampler$setPredictor(stretched, "b", forceUpdate = TRUE)
held <- sampler$data@standardization[[1L]]
sampler$setPredictor(
  x[, 1L] / 2,
  "a",
  forceUpdate = TRUE,
  updateStandardization = TRUE
)
record <- sampler$data@standardization[[1L]]
expect_identical(record$center[["b"]], held$center[["b"]])
expect_identical(record$scale[["b"]], held$scale[["b"]])
expect_equal(record$center[["a"]], mean(x[, 1L] / 2), tolerance = 1e-14)
expect_equal(record$scale[["a"]], sd(x[, 1L] / 2), tolerance = 1e-14)
expectRecordIsReader(sampler, "replaced columns")
expect_true(worstGap(keptBefore, kept(sampler)) < 1e-12)
# an unforced update rolled back derives nothing: a standardization held
# over values it was not derived from stays
sampler <- make()
sampler$setPredictor(
  stretched,
  "b",
  forceUpdate = TRUE,
  updateCutPoints = "position"
)
stale <- recordOf(sampler)$standardization
expect_true(stale[[1L]]$center < 1)
result <- sampler$setPredictor(
  rep(0.5, n),
  "b",
  forceUpdate = FALSE,
  updateStandardization = TRUE
)
expect_false(result)
expect_identical(recordOf(sampler)$standardization, stale)
expect_identical(readerOf(sampler)$standardization, stale)
expect_error(
  sampler$setPredictor(x[, 2L], "b", "partial", updateStandardization = TRUE),
  "cannot also update the standardization"
)
expect_error(
  sampler$setPredictor(x[, 2L], "b", TRUE, updateStandardization = NA),
  "'updateStandardization' must be TRUE or FALSE"
)

# ---- gp narrowed -----------------------------------------------------------

# kept gp draws: a call that keeps the standardization and the range is
# served, one that derives either or carries another is refused, whatever
# argument says so, and a refusal touches nothing
refusedData <- paste0(
  "$setData: 'newData' cannot change the covariate standardization or the ",
  "response range of a sampler with gp leaves that holds saved draws"
)
gpSampler <- make(gp(columns = "b"))
gpKept <- kept(gpSampler)
expect_silent(gpSampler$setData(gpSampler$data))
expect_identical(kept(gpSampler), gpKept)
reordered <- gpSampler$data
reordered@x <- reordered@x[n:1, ]
reordered@y <- reordered@y[n:1]
expect_silent(gpSampler$setData(reordered))
expect_silent(gpSampler$setData(
  dbarts::dbartsData(x.other, y.other),
  updateScale = FALSE,
  updateStandardization = FALSE
))
# the grid never refuses, kept draws holding their thresholds by value
regridded <- gpSampler$data
regridded@cutPoints[[1L]] <- c(0.25, 0.5, 0.75)
expect_silent(gpSampler$setData(regridded))
expect_identical(kept(gpSampler), gpKept)
refusing <- list(
  "a fresh object" = list(dbarts::dbartsData(x, y)),
  "updateScale = TRUE" = list(gpSampler$data, updateScale = TRUE),
  "updateStandardization = TRUE" = list(
    gpSampler$data,
    updateStandardization = TRUE
  )
)
editEntry <- function(name, value) {
  d <- gpSampler$data
  entry <- d@standardization[[1L]]
  entry[name] <- list(value)
  d@standardization[[1L]] <- entry
  d
}
own <- gpSampler$data@standardization[[1L]]
range.moved <- gpSampler$data
range.moved@response.range[[2L]] <- range.moved@response.range[[2L]] + 1e-9
oneNumber <- list(
  center = editEntry("center", own$center + 1e-9),
  scale = editEntry("scale", 2 * own$scale),
  lengthscale = editEntry("lengthscale", 0.1),
  "no lengthscale" = editEntry("lengthscale", NULL),
  range = range.moved
)
for (edit in names(oneNumber)) {
  case <- paste("a record differing in its", edit)
  refusing[[case]] <- list(oneNumber[[edit]])
}
gpSampler$storeState()
heldBefore <- list(
  state = stateOf(gpSampler),
  data = gpSampler$data,
  kept = kept(gpSampler)
)
for (case in names(refusing)) {
  expect_error(
    do.call(gpSampler$setData, refusing[[case]]),
    refusedData,
    fixed = TRUE,
    info = case
  )
  expect_identical(
    list(
      state = stateOf(gpSampler),
      data = gpSampler$data,
      kept = kept(gpSampler)
    ),
    heldBefore,
    info = case
  )
}
expect_error(
  gpSampler$setPredictor(x[, 2L], "b", TRUE, updateStandardization = TRUE),
  paste0(
    "$setPredictor: 'updateStandardization' cannot be TRUE for a sampler ",
    "with gp leaves that holds saved draws"
  ),
  fixed = TRUE
)
expect_identical(stateOf(gpSampler), heldBefore$state)
# a probit gp sampler has no range to change: its own data object is served
probit <- dbarts::dbarts(
  x,
  as.double(y > median(y)),
  family = "probit",
  leaf.prior = gp(columns = "b"),
  control = controlWith()
)
invisible(probit$run(5L, 4L))
expect_silent(probit$setData(probit$data))
# with nothing kept, a derivation is served
unkept <- dbarts::dbarts(
  x,
  y,
  leaf.prior = gp(columns = "b"),
  control = controlWith(keepTrees = FALSE)
)
invisible(unkept$run(5L, 2L))
expect_silent(unkept$setData(dbarts::dbartsData(x.other, y.other)))
expect_silent(unkept$setPredictor(
  stretched,
  "b",
  TRUE,
  updateStandardization = TRUE
))
expectRecordIsReader(unkept, "gp unkept")
record <- unkept$data@standardization[[1L]]
expect_equal(record$center[["b"]], mean(stretched), tolerance = 1e-14)
expect_equal(record$scale[["b"]], sd(stretched), tolerance = 1e-14)
# and for the replaced columns alone, a covariate not replaced keeping its
# center, scale and lengthscale
unkept <- make(gp(columns = c("a", "b")), keepTrees = FALSE)
unkept$setPredictor(stretched, "b", forceUpdate = TRUE)
held <- unkept$data@standardization[[1L]]
unkept$setPredictor(x[, 1L] / 2, "a", TRUE, updateStandardization = TRUE)
record <- unkept$data@standardization[[1L]]
expect_identical(lapply(record, `[[`, "b"), lapply(held, `[[`, "b"))
expect_equal(record$center[["a"]], mean(x[, 1L] / 2), tolerance = 1e-14)
# a covariate without spread records an NA scale, which the leaf divides by
# 1: compared as the leaf divides, its own data object is served
flat <- dbarts::dbarts(
  cbind(a = x[, 1L], b = rep(0.3, n)),
  y,
  leaf.prior = gp(columns = "b"),
  control = controlWith()
)
invisible(flat$run(5L, 4L))
expect_true(is.na(flat$data@standardization[[1L]]$scale[[1L]]))
expect_silent(flat$setData(flat$data))

# ---- record written --------------------------------------------------------

# after every call the slots are the reader's
sampler <- make()
calls <- list(
  "setData, record" = function(s) s$setData(carried),
  "setData, keep" = function(s) {
    s$setData(fresh, updateCutPoints = "none", updateScale = FALSE)
  },
  "setData, derive" = function(s) s$setData(fresh),
  "setData, position" = function(s) {
    s$setData(carried, updateCutPoints = "position")
  },
  "setPredictor" = function(s) {
    s$setPredictor(x[, 2L], "b", TRUE, updateStandardization = TRUE)
  }
)
for (name in names(calls)) {
  calls[[name]](sampler)
  expectRecordIsReader(sampler, name)
}

# ---- fallback where derived ------------------------------------------------

# rows in that hold one value of the response: a derived range reads every
# row and says so, a held one reads none and is silent
oneIn <- c(1, rep(0, n - 1L))
fallbacks <- function(expr) {
  caught <- list()
  withCallingHandlers(
    expr,
    dbartsScaleFallbackWarning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  length(caught)
}
sampler <- make(normal())
replacement <- sampler$data
replacement@weights <- oneIn
before <- sampler$data@response.range
expect_identical(fallbacks(sampler$setData(replacement)), 0L)
expect_identical(sampler$data@response.range, before)
expect_identical(
  fallbacks(sampler$setData(replacement, updateScale = TRUE)),
  1L
)
expect_identical(sampler$data@response.range, range(y))

# ---- table sequences -------------------------------------------------------

# each row of the help's table: the change, the value put back by the call
# the table names, then setState of a state stored before the change installs
# and the next draws are those of a twin that only restored
sequences <- list(
  "cut points, setCutPoints" = list(
    change = function(s) s$setCutPoints(c(0.25, 0.5, 0.75), 1L),
    back = function(s, d) s$setCutPoints(d@cutPoints)
  ),
  "cut points, setData" = list(
    change = function(s) s$setCutPoints(c(0.25, 0.5, 0.75), 1L),
    back = function(s, d) s$setData(d)
  ),
  "weights, setData" = list(
    change = function(s) s$setCutPoints(c(25, 50, 75), 1L),
    back = function(s, d) s$setData(d)
  ),
  "weights, setPredictor" = list(
    change = function(s) s$setCutPoints(c(25, 50, 75), 1L),
    back = function(s, d) {
      s$setPredictor(weighted[, 1L], 1L, TRUE, updateCutPoints = "value")
    }
  ),
  "standardization, setPredictor" = list(
    change = function(s) {
      s$setPredictor(stretched, "b", TRUE, updateStandardization = TRUE)
    },
    back = function(s, d) {
      s$setPredictor(x[, 2L], "b", TRUE, updateStandardization = TRUE)
    }
  ),
  "standardization, setData" = list(
    change = function(s) {
      s$setPredictor(stretched, "b", TRUE, updateStandardization = TRUE)
    },
    back = function(s, d) s$setData(d)
  ),
  "range, setResponse" = list(
    change = function(s) s$setResponse(y.other, updateScale = TRUE),
    back = function(s, d) s$setResponse(y, updateScale = TRUE)
  ),
  "range, setData" = list(
    change = function(s) s$setResponse(y.other, updateScale = TRUE),
    back = function(s, d) s$setData(d)
  )
)
# a column of few values unevenly spaced, whose points carry weights
weighted <- x
weighted[, 1L] <- round(10 * x[, 1L])^2
for (name in names(sequences)) {
  data <- if (startsWith(name, "weights")) weighted else x
  sampler <- make(data = data, keepTrees = FALSE)
  twin <- make(data = data, keepTrees = FALSE)
  sampler$storeState()
  st <- sampler$state
  d <- sampler$data
  if (startsWith(name, "weights")) {
    expect_false(is.null(attr(d@cutPoints[[1L]], "mass")), info = name)
  }
  invisible(sampler$run(0L, 3L))
  sequences[[name]]$change(sampler)
  sequences[[name]]$back(sampler, d)
  expect_identical(readerOf(sampler), readerOf(twin), info = name)
  expect_true(sampler$setState(st), info = name)
  expect_true(twin$setState(st), info = name)
  expect_identical(sampler$run(0L, 5L), twin$run(0L, 5L), info = name)
}
