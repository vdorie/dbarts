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
sampler <- make(normal())
given <- sampler$data
given@cutPoints[[1L]] <- c(0.25, 0.5, 0.75)
sampler$setData(given)
expect_true(all(splitsOn(stateOf(sampler), 1L) %in% c(0.25, 0.5, 0.75)))
expect_identical(sampler$data@cutPoints[[1L]], c(0.25, 0.5, 0.75))

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
# an unforced update rolled back derives nothing
sampler <- make()
own <- recordOf(sampler)$standardization
result <- sampler$setPredictor(
  rep(0.5, n),
  "b",
  forceUpdate = FALSE,
  updateStandardization = TRUE
)
if (isFALSE(result)) {
  expect_identical(recordOf(sampler)$standardization, own)
  expect_identical(readerOf(sampler)$standardization, own)
}
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
  x[, 2L],
  "b",
  TRUE,
  updateStandardization = TRUE
))
expectRecordIsReader(unkept, "gp unkept")

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
