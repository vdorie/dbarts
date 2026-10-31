# The tree count and the cut settings have their own arguments at every door
# that can state them. A value written on a control stands for that argument,
# with one deprecation warning a session per name; a value a fit carried states
# nothing; a fresh control's is never read. Where a control's value differs from
# what the call takes the not-used warning says so on every call. Each row is
# read back from the engine or the objects, with the warnings counted.

forest <- dbartsForests$forest

set.seed(31)
n <- 60L
frame <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))
frame$z <- rbinom(n, 1L, 0.5)
frame$zn <- rnorm(n)
z <- frame$z
zn <- frame$zn
frame$y <- with(frame, x1 + z * (1 + x2) + rnorm(n, sd = 0.3))
x <- as.matrix(frame[c("x1", "x2", "x3")])
y <- frame$y
yb <- as.numeric(y > median(y))
breaks <- quantile(y, c(0, 0.3, 0.7, 1))
yo <- cut(y, breaks, include.lowest = TRUE, ordered_result = TRUE)
ycount <- rpois(n, 2)
ycat <- factor(sample(c("a", "b", "c"), n, TRUE))
time <- rexp(n)
status <- rbinom(n, 1L, 0.7)
yh <- abs(y)
yh[y < 0.5] <- 0

run <- function(...) {
  dbartsControl(
    n.samples = 2L,
    n.burn = 1L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 9L,
    ...
  )
}
# a control made with its constructor's warning set aside and its key cleared, so
# that a door is counted on its own
wrote <- function(...) {
  control <- suppressWarnings(run(...), classes = "dbartsDeprecatedWarning")
  resetKeys()
  control
}
gauss <- function(control = run(), ...) dbarts(x, y, ..., control = control)
terms <- function(formula, control = run(), ...) {
  dbarts(formula, frame, ..., control = control)
}
barts <- function(...) {
  bart(
    ...,
    verbose = FALSE,
    n.samples = 2L,
    n.burn = 1L,
    n.chains = 1L,
    n.threads = 1L,
    samplerOnly = TRUE
  )
}
legacy <- function(...) {
  bartBT(..., verbose = FALSE, ndpost = 2L, nskip = 1L, sampleronly = TRUE)
}
xbarts <- function(...) {
  xbart(
    ...,
    n.samples = 2L,
    method = "k-fold",
    n.test = 2L,
    n.reps = 1L,
    n.burn = c(1L, 1L),
    n.threads = 1L,
    seed = 3L
  )
}
specSampler <- function(...) {
  spec <- dbartsSpec(...)
  new("dbartsSampler", run(), spec$model, spec$data)
}
direct <- function(control, from) {
  new("dbartsSampler", control, from$model, from$data)
}

onceState <- dbarts:::onceWarnState
copyNames <- c("n.trees", "n.cuts", "useQuantiles")
resetKeys <- function() {
  for (name in copyNames) {
    onceState[[paste0("tombstone.control.", name)]] <- NULL
  }
}
# a call with its warnings counted: its value or its error, the names deprecated,
# the names not used and any other warning
counted <- function(expr) {
  resetKeys()
  deprecated <- character()
  notUsed <- character()
  other <- character()
  note <- function(w) {
    message <- conditionMessage(w)
    if (inherits(w, "dbartsControlNotUsedWarning")) {
      notUsed <<- c(notUsed, sub("^the control's ([^ ]+) .*$", "\\1", message))
    } else if (startsWith(message, "dbartsControl(")) {
      deprecated <<- c(
        deprecated,
        sub("^dbartsControl\\(([^ ]+) .*$", "\\1", message)
      )
    } else {
      other <<- c(other, message)
    }
    invokeRestart("muffleWarning")
  }
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = note
  )
  list(
    value = value,
    deprecated = sort(deprecated),
    notUsed = sort(notUsed),
    other = other
  )
}
forestCounts <- function(sampler) {
  ptr <- sampler$getPointer()
  vapply(
    seq_len(dbarts:::bartcoreNumForests(ptr)),
    function(index) dbarts:::bartcoreForestTreeCount(ptr, index - 1L),
    0L
  )
}
statedOf <- function(sampler) {
  if (!inherits(sampler, "error")) sampler$model@forests[[1L]]@stated
}
gridOf <- function(sampler) {
  if (!inherits(sampler, "error")) sampler$data@cutPoints[[1L]]
}
isQuantile <- function(sampler) {
  gaps <- diff(gridOf(sampler))
  !isTRUE(all.equal(gaps, rep(mean(gaps), length(gaps))))
}
# one door's reading: the engine's counts, the cut count, the rule held by the
# data object (the mirror then its isTRUE), and the warnings
expectRow <- function(
  row,
  trees = NULL,
  cuts = NULL,
  rule = NULL,
  d = character(),
  u = character(),
  info = ""
) {
  expect_false(inherits(row$value, "error"), info = info)
  if (inherits(row$value, "error")) {
    return(invisible(conditionMessage(row$value)))
  }
  if (!is.null(trees)) {
    expect_identical(forestCounts(row$value), as.integer(trees), info = info)
  }
  if (!is.null(cuts)) {
    held <- unique(row$value$data@n.cuts)
    expect_identical(held, as.integer(cuts), info = info)
  }
  # a rule of NA is one no caller stated, which a sampler's data object reads
  # as FALSE
  if (!is.null(rule)) {
    control <- row$value$control
    expect_identical(row$value$data@useQuantiles, isTRUE(rule), info = info)
    expect_identical(
      dbarts:::dataStatesRule(row$value$data),
      !is.na(rule),
      info = info
    )
    expect_identical(control@useQuantiles, isTRUE(rule), info = info)
  }
  expect_identical(row$deprecated, sort(d), info = info)
  expect_identical(row$notUsed, sort(u), info = info)
  expect_identical(row$other, character(), info = info)
}
# a slot of the control a row's sampler holds; a refused call has none
controlSlot <- function(row, name) {
  if (!inherits(row$value, "error")) slot(row$value$control, name)
}
# a refused call: by name, and a door's warnings unspent
expectRefused <- function(row, pattern, info = "") {
  expect_inherits(row$value, "error", info = info)
  if (inherits(row$value, "error")) {
    expect_true(grepl(pattern, conditionMessage(row$value)), info = info)
  }
  expect_identical(row$notUsed, character(), info = info)
}

# --- dbartsControl(): three arguments that default to NULL and warn ---

for (name in copyNames) {
  expect_null(formals(dbartsControl)[[name]], info = name)
}
row <- counted(dbartsControl())
expect_identical(row$deprecated, character())
expect_identical(row$value@n.trees, 75L)
expect_false(row$value@useQuantiles)
# a fresh control holds no cut count, where 0.9-34's read 100
expect_error(row$value@n.cuts, "no slot")
expect_false("n.cuts" %in% slotNames("dbartsControl"))

row <- counted(dbartsControl(n.trees = 50L, n.cuts = 20L, useQuantiles = TRUE))
expect_identical(row$deprecated, sort(copyNames))
expect_identical(row$value@n.trees, 50L)
expect_identical(row$value@n.cuts, 20L)
expect_true(row$value@useQuantiles)
# one message a session per name, the constructor's and a door's being one key
resetKeys()
expect_warning(dbartsControl(n.trees = 50L), "dbarts\\(n.trees = \\)")
expect_silent(dbartsControl(n.trees = 51L))
expect_warning(dbartsControl(n.cuts = 50L), "dbartsData\\(n.cuts = \\)")
expect_silent(dbartsControl(n.cuts = 51L))
expect_warning(dbartsControl(useQuantiles = TRUE), "useQuantiles = \\) or")
expect_silent(dbartsControl(useQuantiles = TRUE))
edited <- run()
edited@n.trees <- 33L
resetKeys()
expect_warning(gauss(edited), "dbarts\\(n.trees = \\)")
expect_silent(gauss(edited))
resetKeys()
expect_warning(dbartsControl(n.trees = 50L), "dbarts\\(n.trees = \\)")
expect_silent(gauss(edited))

# a written NULL is the name left out, by value, however it is written
nullName <- NULL
wrapped <- function(...) dbartsControl(...)
for (row in list(
  counted(dbartsControl(n.trees = NULL, n.cuts = NULL, useQuantiles = NULL)),
  counted(dbartsControl(n.trees = nullName, n.cuts = nullName)),
  counted(wrapped(n.trees = NULL, useQuantiles = nullName))
)) {
  expect_identical(row$deprecated, character())
  expect_identical(row$value@n.trees, 75L)
  expect_false(row$value@useQuantiles)
  expect_null(attr(row$value, "n.cuts"))
  named <- dbarts:::controlSuppliedSlots(row$value)
  expect_false(any(copyNames %in% named))
}
# what it refuses it refuses by name
for (bad in list(0L, NA, 1:2)) {
  expect_error(suppressWarnings(dbartsControl(n.trees = bad)), "'n.trees'")
}
for (bad in list(0L, NA)) {
  expect_error(suppressWarnings(dbartsControl(n.cuts = bad)), "'n.cuts'")
}
expect_error(dbartsControl(useQuantiles = NA), "'useQuantiles'")

# the registry's rows
registry <- dbarts:::dbartsTombstones
rowNames <- vapply(registry, `[[`, "", "name")
owners <- vapply(registry, `[[`, "", "owner")
for (name in copyNames) {
  expect_true(any(rowNames == name & owners == "dbartsControl"), info = name)
}
rowName <- "a count, cut count or quantile rule read from a control"
expect_true(rowName %in% rowNames)
for (name in c("n.trees slot", "useQuantiles slot", "binary slot")) {
  condition <- registry[[which(rowNames == name)]]$condition
  expect_true(grepl("no CRAN package", condition, fixed = TRUE), info = name)
}

# --- the flat arguments at dbarts() ---

row <- counted(gauss(n.trees = 30L, n.cuts = 20L, useQuantiles = TRUE))
expectRow(row, 30L, 20L, TRUE)
expect_identical(statedOf(row$value), "n.trees")
expect_true(isQuantile(row$value))
expect_identical(controlSlot(row, "n.trees"), 30L)
expect_true(length(gridOf(row$value)) <= 20L)
row <- counted(gauss(n.cuts = 20L, useQuantiles = FALSE))
expectRow(row, 75L, 20L, FALSE)
expect_false(isQuantile(row$value))
expect_identical(statedOf(row$value), character())
row <- counted(gauss())
expectRow(row, 75L, 100L, NA)
# per column, recycled
expect_identical(gauss(n.cuts = c(5L, 9L))$data@n.cuts, c(5L, 9L, 5L))
expect_error(gauss(n.cuts = 1:4), "'n.cuts' has 4")
# bad values, by name
for (bad in list(0L, NA, 1:2)) {
  expect_error(gauss(n.trees = bad), "'n.trees'")
}
for (bad in list(0L, NA)) {
  expect_error(gauss(n.cuts = bad), "'n.cuts'")
}
expect_error(gauss(useQuantiles = NA), "'useQuantiles'")
expect_error(gauss(useQuantiles = c(TRUE, TRUE)), "'useQuantiles'")
# a written NULL is the name left out, by value: a wrapper's variable and a NULL
# passed through dots read as one written in place
viaDots <- function(...) gauss(...)
viaName <- function(count = NULL) dbarts(x, y, n.trees = count, control = run())
for (row in list(
  counted(gauss(n.trees = NULL, n.cuts = NULL, useQuantiles = NULL)),
  counted(gauss(n.trees = nullName, n.cuts = nullName)),
  counted(viaDots(n.trees = NULL, useQuantiles = NULL)),
  counted(viaName())
)) {
  expectRow(row, 75L, 100L, NA)
  expect_identical(statedOf(row$value), character())
}
# nothing stated, so nothing refused beside forests that all have a basis
forestsAll <- list(forest(basis = ~ factor(z)), forest(basis = ~zn))
row <- counted(terms(y ~ x1 + x2, forests = forestsAll, n.trees = NULL))
expectRow(row, c(50L, 50L))
# each written argument is read once
reads <- 0L
twelve <- function() {
  reads <<- reads + 1L
  12L
}
expect_identical(forestCounts(gauss(n.trees = twelve())), 12L)
expect_identical(reads, 1L)
sampler <- gauss(n.cuts = twelve())
expect_identical(reads, 2L)
reads <- 0L
sampler <- barts(x, y, n.trees = twelve(), n.cuts = twelve())
expect_identical(reads, 2L)
expect_identical(forestCounts(sampler), 12L)

# --- each way a value is written, carried or silent, under every family ---

fitted <- gauss(n.trees = 50L, n.cuts = 25L, useQuantiles = TRUE)
slotEdited <- run()
slotEdited@n.trees <- 33L
slotEdited@useQuantiles <- TRUE
copyEdited <- fitted$control
copyEdited@n.trees <- 44L
copyEditedDefault <- fitted$control
copyEditedDefault@n.trees <- 75L
allWritten <- wrote(n.trees = 30L, n.cuts = 20L, useQuantiles = TRUE)
beat <- list(n.trees = 40L, n.cuts = 10L, useQuantiles = FALSE)
noneGiven <- list(n.trees = NULL, n.cuts = NULL, useQuantiles = NULL)
carriedBack <- list(n.trees = 50L, useQuantiles = TRUE)
quantile <- "useQuantiles"
ways <- list(
  silent = list(run(), list(), 75L, 100L, NA, character(), character()),
  flat = list(
    run(),
    list(n.trees = 30L, n.cuts = 20L, useQuantiles = TRUE),
    30L,
    20L,
    TRUE,
    character(),
    character()
  ),
  flatNull = list(run(), noneGiven, 75L, 100L, NA, character(), character()),
  written = list(allWritten, list(), 30L, 20L, TRUE, copyNames, character()),
  beaten = list(allWritten, beat, 40L, 10L, FALSE, copyNames, copyNames),
  carried = list(
    fitted$control,
    list(),
    75L,
    100L,
    NA,
    character(),
    copyNames
  ),
  carriedAgreeing = list(
    fitted$control,
    carriedBack,
    50L,
    100L,
    TRUE,
    character(),
    "n.cuts"
  ),
  slotEdited = list(
    slotEdited,
    list(),
    33L,
    100L,
    TRUE,
    c("n.trees", quantile),
    character()
  ),
  copyEdited = list(
    copyEdited,
    list(),
    44L,
    100L,
    NA,
    "n.trees",
    c("n.cuts", quantile)
  ),
  copyEditedDefault = list(
    copyEditedDefault,
    list(),
    75L,
    100L,
    NA,
    "n.trees",
    c("n.cuts", quantile)
  )
)
families <- list(
  gaussian = function(...) dbarts(x, y, ...),
  student = function(...) dbarts(x, y, family = dbartsFamilies$student(3), ...),
  probit = function(...) dbarts(x, yb, ...),
  logistic = function(...) dbarts(x, yb, family = "logistic", ...),
  aft = function(...) dbarts(x, cbind(time, status), family = "aft", ...),
  hazard = function(...) {
    dbarts(x, cbind(ceiling(time * 3), status), family = "hazard", ...)
  },
  ordinal = function(...) dbarts(x, yo, family = "ordinal", ...),
  nbinom = function(...) dbarts(x, ycount, family = "nbinom", ...),
  multinomial = function(...) dbarts(x, ycat, family = "multinomial", ...),
  variance = function(...) dbarts(x, y, variance = TRUE, ...)
)
binaryFamilies <- c("probit", "logistic", "hazard")
for (family in names(families)) {
  for (way in names(ways)) {
    info <- paste(family, way)
    w <- ways[[way]]
    row <- counted(do.call(
      families[[family]],
      c(list(control = w[[1L]]), w[[2L]])
    ))
    numForests <- if (identical(family, "multinomial")) 3L else 1L
    trees <- rep_len(w[[3L]], numForests)
    expectRow(row, trees, w[[4L]], w[[5L]], w[[6L]], w[[7L]], info = info)
    if (inherits(row$value, "error")) {
      next
    }
    # the mirrors, written from the model and data
    control <- row$value$control
    expect_identical(control@n.trees, w[[3L]], info = info)
    expect_identical(control@binary, family %in% binaryFamilies, info = info)
    expect_null(attr(control, "n.cuts"), info = info)
    # and the stored control speaks for nothing the next time
    again <- counted(row$value$copy())
    expectRow(again, trees, info = info)
  }
}

# --- a flat value beats a written one, the not-used warning naming it ---

notUsedMessage <- function(expr) {
  tryCatch(
    {
      suppressWarnings(expr, classes = "dbartsDeprecatedWarning")
      ""
    },
    dbartsControlNotUsedWarning = conditionMessage
  )
}
message <- notUsedMessage(gauss(wrote(n.trees = 9L), n.trees = 12L))
expect_true(grepl("n.trees (9) is not used", message, fixed = TRUE))
expect_true(grepl("this fit has 12 trees", message, fixed = TRUE))
expect_true(grepl("dbarts(n.trees = )", message, fixed = TRUE))
message <- notUsedMessage(gauss(wrote(n.cuts = 9L), n.cuts = 12L))
expect_true(grepl("n.cuts (9) is not used", message, fixed = TRUE))
# on every call, where the deprecation is once
nineTrees <- wrote(n.trees = 9L)
tally <- function() {
  deprecated <- 0L
  notUsed <- 0L
  withCallingHandlers(
    gauss(nineTrees, n.trees = 12L),
    dbartsControlNotUsedWarning = function(w) {
      notUsed <<- notUsed + 1L
      invokeRestart("muffleWarning")
    },
    dbartsDeprecatedWarning = function(w) {
      deprecated <<- deprecated + 1L
      invokeRestart("muffleWarning")
    }
  )
  c(deprecated, notUsed)
}
expect_identical(tally(), c(1L, 1L))
expect_identical(tally(), c(0L, 1L))

# --- a count beside a forest's own ---

plainThirty <- y ~ forest(x1 + x2, n.trees = 30L)
beside <- y ~ forest(x1, n.trees = 30L) + forest(x2, basis = zn)
onControl <- "'n.trees' is given to the control and to the forest with no basis"
onFitting <- "'n.trees' is given to the fitting function and to the forest"
for (row in list(
  counted(terms(plainThirty, wrote(n.trees = 50L))),
  counted(terms(y ~ forest(x1, n.trees = 50L), wrote(n.trees = 50L))),
  counted(terms(y ~ x1, wrote(n.trees = 50L), forests = forest(n.trees = 30L))),
  counted(terms(
    y ~ x1,
    wrote(n.trees = 50L),
    forests = list(forest(n.trees = 30L))
  )),
  counted(terms(beside, wrote(n.trees = 50L))),
  counted(barts(plainThirty, frame, control = wrote(n.trees = 50L)))
)) {
  expectRefused(row, onControl)
  expect_identical(row$deprecated, "n.trees")
  expectRefused(row, "leave it out of dbartsControl")
}
for (row in list(
  counted(terms(plainThirty, n.trees = 50L)),
  counted(terms(y ~ x1, forests = forest(n.trees = 30L), n.trees = 50L)),
  counted(barts(plainThirty, frame, n.trees = 50L))
)) {
  expectRefused(row, onFitting)
}
# a forest with a basis stating its own stands beside the plain forest's
twoForests <- y ~ forest(x1) + forest(x2, basis = zn, n.trees = 5L)
row <- counted(terms(twoForests, wrote(n.trees = 30L)))
expectRow(row, c(30L, 5L), d = "n.trees")
row <- counted(terms(twoForests, n.trees = 30L))
expectRow(row, c(30L, 5L))
# a carried count is not a statement, and a written NULL states nothing
row <- counted(terms(plainThirty, fitted$control))
expectRow(row, 30L, u = copyNames)
row <- counted(terms(y ~ forest(x1, n.trees = 50L), fitted$control))
expectRow(row, 50L, u = c("n.cuts", quantile))
row <- counted(terms(plainThirty, n.trees = NULL))
expectRow(row, 30L)

# --- a count where every forest has a basis ---

allTerms <- y ~ forest(x1, basis = factor(z)) + forest(x2, basis = zn)
flatRefusal <- "'n.trees' given to the fitting function is the tree count"
writtenRefusal <- "the control's 'n.trees' is the tree count"
plainData <- dbartsData(y ~ x1 + x2 + x3, frame)
for (row in list(
  counted(terms(y ~ x1, forests = forestsAll, n.trees = 20L)),
  counted(dbartsSpec(plainData, run(), forests = forestsAll, n.trees = 20L)),
  counted(barts(allTerms, frame, n.trees = 20L))
)) {
  expectRefused(row, flatRefusal)
}
for (row in list(
  counted(terms(y ~ x1, wrote(n.trees = 20L), forests = forestsAll)),
  counted(barts(allTerms, frame, control = wrote(n.trees = 20L)))
)) {
  expectRefused(row, writtenRefusal)
  expect_identical(row$deprecated, "n.trees")
}
row <- counted(terms(y ~ x1, fitted$control, forests = forestsAll))
expectRow(row, c(50L, 50L), u = c("n.cuts", quantile))
# a wrapper's unstated default is no statement
row <- counted(barts(allTerms, frame))
expectRow(row, c(50L, 50L))
expect_identical(statedOf(row$value), character())
expectRow(counted(barts(plainThirty, frame)), 30L)
# a bart() fit's control into a model whose every forest has a basis is taken
row <- counted(terms(
  y ~ x1,
  barts(x, y, n.trees = 12L)$control,
  forests = forestsAll
))
expectRow(row, c(50L, 50L), u = "n.trees")

# --- cut settings beside a data object ---

expect_true(is.na(plainData@useQuantiles))
expect_true(anyNA(plainData@n.cuts))
ownData <- dbartsData(
  y ~ x1 + x2 + x3,
  frame,
  n.cuts = 10L,
  useQuantiles = TRUE
)
expect_identical(ownData@n.cuts, rep(10L, 3L))
expect_true(ownData@useQuantiles)
# no grid: what the data object states is used, at every door that takes one
for (row in list(
  counted(dbarts(ownData, control = run())),
  counted(specSampler(ownData, run())),
  counted(barts(ownData))
)) {
  expectRow(row, cuts = 10L, rule = TRUE)
  expect_true(length(gridOf(row$value)) <= 10L)
  expect_true(isQuantile(row$value))
}
# neither: the defaults
for (row in list(
  counted(dbarts(plainData, control = run())),
  counted(specSampler(plainData, run())),
  counted(barts(plainData))
)) {
  expectRow(row, cuts = 100L, rule = NA)
}
# flat beats what it holds, a written value likewise
row <- counted(dbarts(
  ownData,
  n.cuts = 5L,
  useQuantiles = FALSE,
  control = run()
))
expectRow(row, cuts = 5L, rule = FALSE)
expect_false(isQuantile(row$value))
writtenBoth <- wrote(n.cuts = 5L, useQuantiles = FALSE)
for (row in list(
  counted(dbarts(ownData, control = writtenBoth)),
  counted(specSampler(ownData, writtenBoth))
)) {
  expectRow(row, cuts = 5L, rule = FALSE, d = c("n.cuts", quantile))
}
expectRow(counted(barts(ownData, n.cuts = 5L)), cuts = 5L, rule = TRUE)
# what a fit carried leaves what it holds, the not-used warning naming each
row <- counted(dbarts(ownData, control = fitted$control))
expectRow(row, cuts = 10L, rule = TRUE, u = c("n.cuts", "n.trees"))
# a grid: kept where nothing is stated, derived again at a stated setting
fit20 <- gauss(n.cuts = 20L)
expect_false(is.null(fit20$data@cutPoints))
row <- counted(dbarts(fit20$data, control = run()))
expectRow(row, cuts = 20L)
expect_identical(gridOf(row$value), gridOf(fit20))
gridSet <- fit20$copy()
gridSet$setCutPoints(c(0.25, 0.5, 0.75), 1L)
points <- function(row) length(gridOf(row$value))
for (make in list(
  function(...) dbarts(gridSet$data, ..., control = run()),
  function(...) barts(gridSet$data, ...)
)) {
  expect_identical(points(counted(make())), 3L)
  expect_true(points(counted(make(n.cuts = 10L))) > 3L)
  row <- counted(make(useQuantiles = TRUE))
  expect_true(points(row) > 3L)
  expect_true(!inherits(row$value, "error") && row$value$data@useQuantiles)
}
row <- counted(dbarts(gridSet$data, control = wrote(n.cuts = 10L)))
expectRow(row, cuts = 10L, d = "n.cuts")
expect_true(points(row) > 3L)
row <- counted(dbarts(gridSet$data, control = wrote(useQuantiles = TRUE)))
expect_true(points(row) > 3L)
row <- counted(specSampler(gridSet$data, wrote(n.cuts = 10L)))
expect_true(points(row) > 3L)
row <- counted(specSampler(gridSet$data, run()))
expect_identical(points(row), 3L)
# a written NULL states none, so the grid stands
row <- counted(dbarts(
  gridSet$data,
  n.cuts = NULL,
  useQuantiles = NULL,
  control = run()
))
expect_identical(points(row), 3L)

# the data door
for (bad in list(0L, NA, -1L)) {
  expect_error(dbartsData(y ~ x1, frame, n.cuts = bad), "'n.cuts'")
}
expect_error(dbartsData(y ~ x1, frame, useQuantiles = NA), "'useQuantiles'")
expect_error(dbartsData(y ~ x1, frame, n.cuts = 1:3), "'n.cuts' has 3")
expect_silent(nulls <- dbartsData(y ~ x1, frame, n.cuts = NULL))
expect_silent(nulls <- dbartsData(y ~ x1, frame, useQuantiles = nullName))
expect_true(is.na(nulls@useQuantiles))
expect_true(anyNA(nulls@n.cuts))
matrixData <- dbartsData(x, y, n.cuts = c(4L, 8L), useQuantiles = FALSE)
expect_identical(matrixData@n.cuts, c(4L, 8L, 4L))
expect_false(matrixData@useQuantiles)
# a data object takes no more arguments, as before, and a NULL is not one
expect_warning(dbartsData(ownData, n.cuts = 3L), "are ignored: 'n.cuts'")
expect_silent(same <- dbartsData(ownData, n.cuts = NULL, useQuantiles = NULL))
expect_identical(same@n.cuts, ownData@n.cuts)
# a hazard fit's appended column takes the data object's count
hazardFrame <- cbind(frame, time = ceiling(time * 3), status = status)
hazardData <- dbartsData(
  survival::Surv(time, status) ~ x1 + x2,
  hazardFrame,
  n.cuts = 7L
)
row <- counted(dbarts(hazardData, family = "hazard", control = run()))
expect_false(inherits(row$value, "error"))
if (!inherits(row$value, "error")) {
  expect_identical(row$value$data@n.cuts, rep(7L, 3L))
}

# xbart's comparison is with what it takes: a data object's count and rule
agreeing <- dbartsData(x, y, n.cuts = 25L, useQuantiles = TRUE)
quantileFit <- gauss(n.trees = 50L, n.cuts = 25L, useQuantiles = TRUE)
row <- counted(xbarts(agreeing, n.trees = 50L, control = quantileFit$control))
expect_identical(row$notUsed, character())
own <- dbartsData(x, y, n.cuts = 10L, useQuantiles = TRUE)
row <- counted(xbarts(own, n.trees = 5L, control = quantileFit$control))
expect_identical(row$notUsed, sort(c("n.trees", "n.cuts")))
# a per-column count is not the first of them
perColumn <- gauss(n.cuts = c(5L, 7L, 9L))
row <- counted(xbarts(x, y, n.cuts = 5L, control = perColumn$control))
expect_identical(row$notUsed, "n.cuts")
row <- counted(gauss(perColumn$control, n.cuts = 5L))
expect_identical(row$notUsed, "n.cuts")
# one per column of another count of columns, all alike, is the same count
sameCuts <- gauss(n.cuts = 8L)
row <- counted(terms(y ~ x1, sameCuts$control, n.cuts = 8L))
expect_identical(row$notUsed, character())
# and one that differs among them is not, however few columns this fit has
mixedCuts <- gauss(n.cuts = c(8L, 8L, 5L))
row <- counted(terms(y ~ x1, mixedCuts$control, n.cuts = 8L))
expect_identical(row$notUsed, "n.cuts")

# what each door says, whatever the comparison
said <- function(expr) {
  told <- character()
  withCallingHandlers(
    tryCatch(expr, error = function(e) stop(e)),
    dbartsControlNotUsedWarning = function(w) {
      told <<- c(told, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  told
}
expect_identical(
  said(gauss(quantileFit$control)),
  c(
    "the control's n.trees (50) is not used: this fit has 75 trees; state the count as dbarts(n.trees = ) or forest(n.trees = )",
    "the control's useQuantiles (TRUE) is not used: this fit takes FALSE; state the rule as dbarts(useQuantiles = ) or dbartsData(useQuantiles = )",
    "the control's n.cuts (25) is not used: this fit takes 100; state the count as dbarts(n.cuts = ) or dbartsData(n.cuts = )"
  )
)
expect_identical(
  said(dbartsSpec(dbartsData(x, y), quantileFit$control)),
  c(
    "the control's n.trees (50) is not used: this fit has 75 trees; state the count as dbartsSpec(n.trees = ) or forest(n.trees = )",
    "the control's useQuantiles (TRUE) is not used: this fit takes FALSE; state the rule as dbartsData(useQuantiles = )",
    "the control's n.cuts (25) is not used: this fit takes 100; state the count as dbartsData(n.cuts = )"
  )
)
expect_identical(
  said(barts(x, y, control = quantileFit$control))[1L],
  "the control's n.trees (50) is not used: this fit has 75 trees; state the count as bart(n.trees = ) or forest(n.trees = )"
)
# bart() given a data object takes no forest
expect_identical(
  said(barts(quantileFit$data, n.trees = 75L, control = quantileFit$control)),
  "the control's n.trees (50) is not used: this fit has 75 trees; state the count as bart(n.trees = )"
)
expect_identical(
  said(xbarts(x, y, n.trees = 5L, control = quantileFit$control)),
  c(
    "the control's n.trees (50) is not used: this cross-validation takes 5; state it as xbart(n.trees = )",
    "the control's useQuantiles (TRUE) is not used: this cross-validation takes FALSE; state it as xbart(useQuantiles = )",
    "the control's n.cuts (25) is not used: this cross-validation takes 100; state it as xbart(n.cuts = )"
  )
)
expect_identical(
  said(terms(
    y ~ x1 + x2 + x3,
    quantileFit$control,
    forests = list(forest(), forest(basis = ~zn))
  ))[1L],
  "the control's n.trees (50) is not used: this fit has 75 trees in its forest with no basis; state the count as dbarts(n.trees = ) or forest(n.trees = )"
)
expect_identical(
  said(terms(
    y ~ x1 + x2 + x3,
    quantileFit$control,
    forests = list(
      forest(x1, basis = ~zn, n.trees = 7L),
      forest(x2, basis = ~z)
    )
  ))[1L],
  "the control's n.trees (50) is not used: this fit's first forest has 7 trees; state the count as forest(n.trees = )"
)
# the same at bart() given a data object whose every forest has a basis
allBasisData <- dbartsData(
  y ~ x1 + x2 + x3,
  frame,
  bases = list(as.matrix(zn), as.matrix(z))
)
expect_identical(
  said(barts(allBasisData, control = gauss(n.trees = 60L)$control))[1L],
  "the control's n.trees (60) is not used: this fit's first forest has 50 trees; fit with dbarts() to state the count"
)

# --- dbartsSpec() ---

row <- counted(specSampler(plainData, run(), n.trees = 30L))
expectRow(row, 30L)
expect_identical(statedOf(row$value), "n.trees")
for (row in list(
  counted(specSampler(plainData, run(), n.trees = NULL)),
  counted(specSampler(plainData, run(), n.trees = nullName))
)) {
  expectRow(row, 75L)
  expect_identical(statedOf(row$value), character())
}
row <- counted(specSampler(plainData, wrote(n.trees = 40L)))
expectRow(row, 40L, d = "n.trees")
expect_identical(statedOf(row$value), "n.trees")
row <- counted(specSampler(plainData, wrote(n.trees = 40L), n.trees = 30L))
expectRow(row, 30L, d = "n.trees", u = "n.trees")
row <- counted(specSampler(plainData, fitted$control))
expectRow(row, 75L, u = copyNames)
row <- counted(specSampler(plainData, run(), forests = forest(n.trees = 12L)))
expectRow(row, 12L)
row <- counted(specSampler(
  plainData,
  run(),
  forests = forest(n.trees = 12L),
  n.trees = 13L
))
expectRefused(row, onFitting)
row <- counted(specSampler(
  plainData,
  wrote(n.trees = 13L),
  forests = forest(n.trees = 12L)
))
expectRefused(row, onControl)
expect_error(dbartsSpec(plainData, run(), n.trees = 0L), "'n.trees'")
# the rule rides the data object, so a fresh control at the direct door keeps it
row <- counted(specSampler(plainData, wrote(useQuantiles = TRUE)))
expectRow(row, rule = TRUE, d = quantile)
expect_true(isQuantile(row$value))

# --- the bare forest() ---

expectRow(counted(terms(y ~ x1, forests = forest(n.trees = 50L))), 50L)
held <- forest(n.trees = 50L, vars = ~x1)
expectRow(counted(terms(y ~ x1 + x2, forests = held)), 50L)
withList <- terms(y ~ x1 + x2, forests = list(forest(n.trees = 50L)))
withBare <- terms(y ~ x1 + x2, forests = forest(n.trees = 50L))
expect_identical(withBare$run()$train, withList$run()$train)
# one forest with a basis is still refused, in either spelling
expect_error(
  terms(y ~ x1 + x2, forests = forest(basis = ~zn)),
  "at least two forests"
)

# --- the wrappers: flat on every call, and a control naming none ---

row <- counted(barts(x, y))
expectRow(row, 75L, 100L, NA)
expect_identical(statedOf(row$value), character())
row <- counted(barts(x, y, n.trees = 30L, n.cuts = 20L, useQuantiles = TRUE))
expectRow(row, 30L, 20L, TRUE)
expect_identical(statedOf(row$value), "n.trees")
row <- counted(barts(x, y, control = allWritten))
expectRow(row, 30L, 20L, TRUE, d = copyNames)
row <- counted(barts(x, y, n.trees = 12L, control = wrote(n.trees = 9L)))
expectRow(row, 12L, d = "n.trees", u = "n.trees")
row <- counted(barts(x, y, n.cuts = 12L, control = wrote(n.cuts = 9L)))
expectRow(row, cuts = 12L, d = "n.cuts", u = "n.cuts")
row <- counted(barts(
  x,
  y,
  useQuantiles = TRUE,
  control = wrote(useQuantiles = FALSE)
))
expectRow(row, rule = TRUE, d = quantile, u = quantile)
row <- counted(barts(x, y, control = fitted$control))
expectRow(row, 75L, 100L, NA, u = copyNames)
row <- counted(barts(x, y, control = slotEdited))
expectRow(row, 33L, rule = TRUE, d = c("n.trees", quantile))
# a signature whose default is not NULL refuses a written NULL
for (name in copyNames) {
  arguments <- list(x, y)
  arguments[name] <- list(NULL)
  expect_error(do.call(barts, arguments), paste0("'", name, "' cannot be NULL"))
}
expect_error(barts(x, y, n.trees = 0L), "'n.trees'")
expect_error(barts(x, y, n.cuts = 0L), "'n.cuts'")
expect_error(barts(x, y, useQuantiles = NA), "'useQuantiles'")
# a fit's reported count is the fit's, on the families packaged apart
for (make in list(
  function(...) bart(x, yo, family = "ordinal", ...),
  function(...) bart(x, ycount, family = "nbinom", ...),
  function(...) bart(x, ycat, family = "multinomial", ...)
)) {
  fit <- make(
    n.trees = 7L,
    verbose = FALSE,
    n.samples = 2L,
    n.burn = 1L,
    n.chains = 1L,
    n.threads = 1L
  )
  expect_identical(fit$n.trees, 7L)
}
# the two parts of a hurdle fit, one not-used warning between them
eightTrees <- wrote(n.trees = 8L)
hurdle <- counted(bart(
  x,
  yh,
  family = "hurdle.lognormal",
  n.trees = 7L,
  n.cuts = 9L,
  control = eightTrees,
  verbose = FALSE,
  keepTrees = TRUE,
  n.burn = 1L,
  n.samples = 2L,
  n.chains = 1L,
  n.threads = 1L
))
expect_identical(hurdle$notUsed, "n.trees")
expect_identical(hurdle$deprecated, "n.trees")
expect_identical(hurdle$value$zero$fit$control@n.trees, 7L)
expect_identical(hurdle$value$positive$fit$control@n.trees, 7L)

# bartBT keeps BayesTree's 200, passed flat and not stated
row <- counted(legacy(x, y))
expectRow(row, 200L, 100L, NA)
expect_identical(statedOf(row$value), character())
row <- counted(legacy(x, y, ntree = 30L, numcut = 20L, usequants = TRUE))
expectRow(row, 30L, 20L, TRUE)
expect_identical(statedOf(row$value), "n.trees")
expect_error(legacy(x, y, ntree = NULL), "'ntree' cannot be NULL")
expect_error(legacy(x, y, numcut = NULL), "'numcut' cannot be NULL")
expect_error(legacy(x, y, usequants = NULL), "'usequants' cannot be NULL")

# rbart_vi fits through dbarts()
rbartFit <- counted(rbart_vi(
  y ~ x1 + x2,
  frame,
  group.by = rep(1:3, 20L),
  n.trees = 5L,
  n.cuts = 9L,
  useQuantiles = TRUE,
  n.samples = 2L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
))
rbartSampler <- rbartFit$value$fit[[1L]]
expect_identical(rbartSampler$control@n.trees, 5L)
expect_identical(rbartSampler$data@n.cuts, rep(9L, 2L))
expect_true(rbartSampler$data@useQuantiles)
expect_identical(rbartFit$deprecated, character())
expect_identical(rbartFit$notUsed, character())
expect_error(
  rbart_vi(y ~ x1, frame, group.by = rep(1:3, 20L), n.cuts = NULL),
  "'n.cuts' cannot be NULL"
)

# xbart: a written count is a grid of one
xbartGrid <- function(...) xbarts(x, y, drop = FALSE, ...)
expectGrid <- function(row, grid, d = character(), u = character()) {
  expect_identical(dimnames(row$value)$n.trees, grid)
  expect_identical(row$deprecated, sort(d))
  expect_identical(row$notUsed, sort(u))
}
expectGrid(counted(xbartGrid()), "75")
expectGrid(counted(xbartGrid(n.trees = c(5L, 7L))), c("5", "7"))
expectGrid(counted(xbartGrid(control = wrote(n.trees = 9L))), "9", "n.trees")
row <- counted(xbartGrid(n.trees = 5L, control = wrote(n.trees = 9L)))
expectGrid(row, "5", "n.trees", "n.trees")
row <- counted(xbartGrid(
  n.cuts = 5L,
  useQuantiles = TRUE,
  control = wrote(n.cuts = 9L, useQuantiles = FALSE)
))
expectGrid(row, "75", c("n.cuts", quantile), c("n.cuts", quantile))
row <- counted(xbartGrid(control = fitted$control))
expectGrid(row, "75", u = copyNames)
for (name in copyNames) {
  arguments <- list()
  arguments[name] <- list(NULL)
  expect_error(do.call(xbartGrid, arguments), paste0("'", name, "' cannot be"))
}
# its cut settings reach its draws: flat, written and on a data object agree
scoresOf <- function(...) suppressWarnings(xbarts(...))
flatCuts <- scoresOf(x, y, n.cuts = 6L, useQuantiles = TRUE)
writtenCuts <- scoresOf(x, y, control = wrote(n.cuts = 6L, useQuantiles = TRUE))
expect_identical(writtenCuts, flatCuts)
dataCuts <- scoresOf(dbartsData(x, y, n.cuts = 6L, useQuantiles = TRUE))
expect_identical(dataCuts, flatCuts)
expect_false(identical(scoresOf(x, y, n.cuts = 6L), flatCuts))
expect_false(identical(scoresOf(x, y), flatCuts))
writtenCount <- scoresOf(x, y, control = wrote(n.trees = 6L))
expect_identical(writtenCount, scoresOf(x, y, n.trees = 6L))
expect_false(identical(writtenCount, scoresOf(x, y)))

# --- the direct door ---

fit75 <- gauss()
fit200 <- gauss(n.trees = 200L)
fitRule <- gauss(useQuantiles = TRUE)
fitUniform <- gauss(useQuantiles = FALSE)
# a fresh control is silent; the model's count and the data object's rule stand
expectRow(counted(direct(run(), fit75)), 75L, 100L)
expectRow(counted(direct(run(), fit200)), 200L)
expectRow(counted(direct(run(), fitRule)), 75L, rule = TRUE)
# a copied control: the model's, with the not-used warning where it differs
expectRow(counted(direct(fit200$control, fit75)), 75L, u = "n.trees")
expectRow(counted(direct(gauss(n.cuts = 9L)$control, fit75)), 75L, u = "n.cuts")
expectRow(counted(direct(fit75$control, fit200)), 200L, u = "n.trees")
expectRow(counted(direct(fit75$control, fit75)), 75L)
expectRow(counted(direct(fitRule$control, fit75)), 75L, u = quantile)
expectRow(
  counted(direct(fit75$control, fitRule)),
  75L,
  rule = TRUE,
  u = quantile
)
# written, the model stating none: used, marked, warned once, the model handed
# in left as it was
row <- counted(direct(wrote(n.trees = 20L), fit75))
expectRow(row, 20L, d = "n.trees")
expect_identical(statedOf(row$value), "n.trees")
expect_identical(controlSlot(row, "n.trees"), 20L)
expect_identical(fit75$model@forests[[1L]]@n.trees, 75L)
# written by an edit of a fit's control
copyWritten <- fit75$control
copyWritten@n.trees <- 21L
expectRow(counted(direct(copyWritten, fit75)), 21L, d = "n.trees")
# written beside a count the model states: refused where it differs, by name
row <- counted(direct(wrote(n.trees = 20L), fit200))
expectRefused(
  row,
  "'n.trees' \\(20\\) differs from the count the model states \\(200\\)"
)
expectRow(counted(direct(wrote(n.trees = 200L), fit200)), 200L, d = "n.trees")
stating <- terms(y ~ x1 + x2, forests = forest(n.trees = 25L))
row <- counted(direct(wrote(n.trees = 20L), stating))
expectRefused(row, "differs from the count the model states \\(25\\)")
# several forests: the one with no basis takes it, first or not
bcf <- terms(y ~ forest(x1 + x2) + forest(x1, basis = factor(z)))
row <- counted(direct(wrote(n.trees = 20L), bcf))
expectRow(row, c(20L, 50L), d = "n.trees")
firstBasis <- terms(y ~ x1, forests = list(forest(basis = ~zn), forest()))
row <- counted(direct(wrote(n.trees = 20L), firstBasis))
expectRow(row, c(50L, 20L), d = "n.trees")
expectRow(counted(direct(run(), firstBasis)), c(50L, 75L))
# every forest with a basis: no forest for it
everyBasis <- terms(y ~ x1 + x2, forests = forestsAll)
row <- counted(direct(wrote(n.trees = 20L), everyBasis))
expectRefused(
  row,
  "'n.trees' \\(20\\) is the tree count of the forest with no basis"
)
expectRow(counted(direct(run(), everyBasis)), c(50L, 50L))
# a multinomial model's one record
multinomial <- dbarts(x, ycat, family = "multinomial", control = run())
row <- counted(direct(wrote(n.trees = 20L), multinomial))
expectRow(row, rep(20L, 3L), d = "n.trees")
# a written rule, the data object stating none: used, the grid derived again
gridSet <- fit75$copy()
gridSet$setCutPoints(c(0.25, 0.5, 0.75), 1L)
row <- counted(direct(wrote(useQuantiles = TRUE), gridSet))
expectRow(row, rule = TRUE, d = quantile)
expect_true(points(row) > 3L)
expect_true(isQuantile(row$value))
expect_false(gridSet$data@useQuantiles)
expect_identical(length(gridOf(gridSet)), 3L)
row <- counted(direct(wrote(useQuantiles = FALSE), gridSet))
expectRow(row, rule = FALSE, d = quantile)
expect_true(points(row) > 3L)
# a sampler's data object reads its rule plain, whoever stated it, and what it
# reads may be stated
reloadedFit <- {
  saved <- tempfile(fileext = ".rds")
  saveRDS(fit75, saved)
  readRDS(saved)
}
unlink(saved)
for (held in list(fit75, fit75$copy(), reloadedFit)) {
  expect_identical(held$data@useQuantiles, FALSE)
  expect_false(dbarts:::dataStatesRule(held$data))
  row <- counted(direct(wrote(useQuantiles = TRUE), held))
  expectRow(row, rule = TRUE, d = "useQuantiles")
}
expect_true(is.na(dbartsData(x, y)@useQuantiles))
restated <- dbartsData(x, y, useQuantiles = fit75$data@useQuantiles)
expect_identical(restated@useQuantiles, FALSE)
expect_true(dbarts:::dataStatesRule(restated))
expect_true(dbarts:::dataStatesRule(gauss(useQuantiles = FALSE)$data))

# a written rule beside one the data object states: refused where it differs
row <- counted(direct(wrote(useQuantiles = FALSE), fitRule))
expectRefused(
  row,
  "'useQuantiles' \\(FALSE\\) differs from the rule the data object states \\(TRUE\\)"
)
row <- counted(direct(wrote(useQuantiles = TRUE), fitUniform))
expectRefused(row, "differs from the rule the data object states \\(FALSE\\)")
row <- counted(direct(wrote(useQuantiles = TRUE), fitRule))
expectRow(row, rule = TRUE, d = quantile)
expect_identical(gridOf(row$value), gridOf(fitRule))
# a written cut count: the data object's, told where it differs
row <- counted(direct(wrote(n.cuts = 20L), fit75))
expectRow(row, cuts = 100L, d = "n.cuts", u = "n.cuts")
expect_identical(gridOf(row$value), gridOf(fit75))
expect_null(attr(
  if (!inherits(row$value, "error")) row$value$control,
  "n.cuts"
))
row <- counted(direct(wrote(n.cuts = 100L), fit75))
expectRow(row, cuts = 100L, d = "n.cuts")
# the spec's data object, with no grid yet
spec <- dbartsSpec(plainData, run())
row <- counted(direct(wrote(useQuantiles = TRUE), spec))
expectRow(row, rule = TRUE, d = quantile)
expect_true(isQuantile(row$value))
# a refusal leaves the sampler its objects came from as it was
before <- fit200$run()$train
row <- counted(direct(wrote(n.trees = 20L), fit200))
expect_inherits(row$value, "error")
expect_identical(forestCounts(fit200), 200L)
expect_identical(fit200$model@forests[[1L]]@n.trees, 200L)
expect_identical(fit200$control@n.trees, 200L)
expect_identical(dim(fit200$run()$train), dim(before))
# re-creation: the stored model's, the mirror agreeing, no warning
written20 <- counted(direct(wrote(n.trees = 20L), fit75))$value
for (made in list(fit75, fit200, fitRule, written20)) {
  made$storeState()
  file <- tempfile(fileext = ".rds")
  saveRDS(made, file)
  trees <- forestCounts(made)
  expectRow(counted(readRDS(file)), trees)
  unlink(file)
  expectRow(counted(made$copy()), trees)
  expectRow(counted(direct(made$control, made)), trees)
}

# a model made by hand records no forest: a written count is its one forest's,
# and comes back through $setControl, $copy and a reload
strippedModel <- fit75$model
strippedModel@forests <- list()
for (handModel in list(
  strippedModel,
  new(
    "dbartsModel",
    fit75$model@tree.prior,
    fit75$model@leaf.prior,
    fit75$model@leaf.hyperprior,
    fit75$model@resid.prior,
    leaf.scale = fit75$model@leaf.scale
  )
)) {
  row <- counted(direct(
    wrote(n.trees = 20L),
    list(model = handModel, data = fit75$data)
  ))
  expectRow(row, 20L, d = "n.trees")
  handSampler <- row$value
  expect_identical(controlSlot(row, "n.trees"), 20L)
  expect_identical(statedOf(handSampler), "n.trees")
  if (!inherits(handSampler, "error")) {
    expectRow(
      counted(handSampler$setControl(run(n.thin = 2L))),
      d = character()
    )
    expect_identical(forestCounts(handSampler), 20L)
    expect_identical(handSampler$control@n.trees, 20L)
    expect_identical(handSampler$control@n.thin, 2L)
    expectRow(counted(handSampler$copy()), 20L)
    handSampler$storeState()
    file <- tempfile(fileext = ".rds")
    saveRDS(handSampler, file)
    expectRow(counted(readRDS(file)), 20L)
    unlink(file)
  }
}

# --- a bad value written on a slot ---

badCount <- run()
badCount@n.trees <- 0L
badTwoCounts <- run()
badTwoCounts@n.trees <- 1:2
badRule <- run()
badRule@useQuantiles <- NA
badTwoRules <- run()
badTwoRules@useQuantiles <- c(TRUE, FALSE)
badCuts <- run()
attr(badCuts, "n.cuts") <- -1L
badRows <- list(
  list(badCount, "'n.trees' must be a positive integer"),
  list(badTwoCounts, "'n.trees' must be of length 1"),
  list(badRule, "'useQuantiles' must be TRUE/FALSE"),
  list(badTwoRules, "'useQuantiles' must be of length 1"),
  list(badCuts, "'n.cuts' must contain only positive integers")
)
for (bad in badRows) {
  expectRefused(counted(gauss(bad[[1L]])), bad[[2L]], info = bad[[2L]])
  expectRefused(counted(direct(bad[[1L]], fit75)), bad[[2L]], info = bad[[2L]])
  expectRefused(
    counted(barts(x, y, control = bad[[1L]])),
    bad[[2L]],
    info = bad[[2L]]
  )
  expectRefused(
    counted(xbarts(x, y, control = bad[[1L]])),
    bad[[2L]],
    info = bad[[2L]]
  )
  expectRefused(
    counted(dbartsSpec(fit75$data, bad[[1L]])),
    bad[[2L]],
    info = bad[[2L]]
  )
  # at $setControl it is a change, which is refused first
  expectRefused(
    counted(gauss()$setControl(bad[[1L]])),
    "is not available on an existing",
    info = bad[[2L]]
  )
}

# what this door says
probitFit <- dbarts(x, yb, family = "probit", control = run())
expect_identical(
  said(direct(quantileFit$control, fit75)),
  c(
    "the control's n.trees (50) is not used: a sampler assembled by hand takes the model's count, 75",
    "the control's useQuantiles (TRUE) is not used: a sampler assembled by hand takes the data object's rule, FALSE",
    "the control's n.cuts (25) is not used: a sampler assembled by hand takes the data object's count, 100"
  )
)
expect_identical(
  said(direct(probitFit$control, fit75)),
  "the control's binary (TRUE) is not used: a sampler assembled by hand takes the model's family, \"gaussian\""
)
expect_identical(
  said(fit75$copy()$setControl(quantileFit$control)),
  c(
    "the control's n.trees (50) is not used: this sampler keeps its 75, fixed when it was created",
    "the control's useQuantiles (TRUE) is not used: this sampler keeps its FALSE, fixed when it was created",
    "the control's n.cuts (25) is not used: this sampler keeps its 100, fixed when it was created"
  )
)

# --- the scripts ---

# twice with one control: one warning a session
resetKeys()
warnings <- 0L
oneControl <- withCallingHandlers(
  {
    ctl <- run(n.trees = 50L)
    list(gauss(ctl), gauss(ctl))
  },
  dbartsDeprecatedWarning = function(w) {
    warnings <<- warnings + 1L
    invokeRestart("muffleWarning")
  }
)
expect_identical(warnings, 1L)
expect_identical(forestCounts(oneControl[[1L]]), 50L)
expect_identical(forestCounts(oneControl[[2L]]), 50L)
# the rule written on a control beside a fit's data object
row <- counted(dbarts(fit75$data, control = wrote(useQuantiles = TRUE)))
expectRow(row, rule = TRUE, d = quantile)
expect_true(isQuantile(row$value))
# a survival fit's data object under another family: its status stays behind
aftFit <- dbarts(x, cbind(time, status), family = "aft", control = run())
expect_false(is.null(aftFit$data@status))
aftSpec <- function() {
  spec <- dbartsSpec(aftFit$data, run())
  new("dbartsSampler", run(), spec$model, spec$data)
}
for (row in list(
  counted(dbarts(aftFit$data, control = run())),
  counted(barts(aftFit$data)),
  counted(aftSpec()),
  counted(direct(run(), list(model = fit75$model, data = aftFit$data)))
)) {
  expectRow(row, 75L)
  if (!inherits(row$value, "error")) {
    expect_identical(row$value$model@family, "gaussian")
    expect_null(row$value$data@status)
  }
}
row <- counted(dbarts(aftFit$data, control = run()))
expect_false(inherits(row$value, "error"))
if (!inherits(row$value, "error")) {
  expect_identical(
    row$value$run()$train,
    dbarts(aftFit$data@x, aftFit$data@y, control = run())$run()$train
  )
}
# and the sampler that has a status keeps it, and one that has none refuses one
expect_false(is.null(aftFit$copy()$data@status))
expect_error(gauss()$setResponse(y, status = rep(1, n)), "non-aft")
# a variance fit's control into a plain fit
varianceFit <- gauss(variance = TRUE)
for (row in list(
  counted(gauss(varianceFit$control)),
  counted(barts(x, y, control = varianceFit$control))
)) {
  expectRow(row, 75L)
  expect_true(
    !inherits(row$value, "error") && is.null(row$value$model@variance)
  )
}
# a probit fit's control into a gaussian fit, with no word of binary
row <- counted(gauss(dbarts(x, yb, control = run())$control))
expectRow(row, 75L)
expect_false(controlSlot(row, "binary"))
# a stored fit's readers
fit <- bart(
  x,
  y,
  n.trees = 11L,
  verbose = FALSE,
  n.samples = 2L,
  n.burn = 1L,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE
)
expect_identical(fit$fit$control@n.trees, 11L)
expect_false(fit$fit$control@binary)
expect_false(fit$fit$control@useQuantiles)
expect_error(fit$fit$control@n.cuts, "no slot")
