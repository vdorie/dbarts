# What a sampler derives from its data and holds while the data move - the
# cut grid and its weights, the response range and a calibration map's anchor
# and unit - is recorded on its data object, cutPoints, response.range and
# forest.calibration, written from the engine after creation and after every
# call that changes one. A creation from a data object takes each value it
# records and derives the rest, so a copy and a reload are built on the
# record and none from the state; a fresh data object records nothing and
# derives as before. A copy of a sampler whose engine exists in this session
# takes its engine's current chain; a reload takes the stored state, its
# splits moved onto the recorded grid by value. xbart's folds derive their
# own. The object a caller handed in is never written.

set.seed(23L, sample.kind = "Rejection")
n <- 120L
x <- cbind(a = runif(n), b = runif(n))
y <- sin(4 * x[, 1L]) + x[, 2L] + rnorm(n, 0, 0.2)
controlWith <- function(...) {
  settings <- list(
    n.chains = 1L,
    n.threads = 1L,
    n.trees = 15L,
    n.samples = 4L,
    n.burn = 0L,
    updateState = FALSE,
    seed = 7L
  )
  do.call(dbarts::dbartsControl, utils::modifyList(settings, list(...)))
}
make <- function(..., data = x, response = y, run = TRUE) {
  sampler <- dbarts::dbarts(data, response, control = controlWith(...))
  if (run) {
    invisible(sampler$run(10L, 1L))
  }
  sampler
}
readerOf <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_derivedRecords, sampler$getPointer(), TRUE)
}
recordOf <- function(sampler) {
  list(
    cutPoints = sampler$data@cutPoints,
    response.range = sampler$data@response.range,
    forest.calibration = sampler$data@forest.calibration,
    standardization = sampler$data@standardization
  )
}
stateOf <- function(sampler) {
  .Call(dbarts:::C_dbarts_bartcore_storeState, sampler$getPointer())
}
cutPointsOf <- function(sampler) attr(stateOf(sampler), "cutPoints")
treeFields <- c("tree.vars", "tree.values", "tree.sizes", "tree.flags")
liveOf <- function(state) {
  lapply(state, function(chain) lapply(chain$forests, `[`, treeFields))
}
reloadOf <- function(sampler) {
  path <- tempfile(fileext = ".rds")
  saveRDS(sampler, path)
  reloaded <- readRDS(path)
  unlink(path)
  reloaded
}
expectRecordIsReader <- function(sampler, info) {
  expect_identical(recordOf(sampler), readerOf(sampler), info = info)
}
three <- c(0.25, 0.5, 0.75)

# ---- written at creation ---------------------------------------------------

# at every door the record is the engine's, and the object handed in carries
# none
handed <- dbarts::dbartsData(x, y)
sampler <- dbarts::dbarts(handed, control = controlWith())
expectRecordIsReader(sampler, "dbarts")
expect_identical(lengths(sampler$data@cutPoints), c(100L, 100L))
expect_identical(sampler$data@response.range, range(y))
expect_null(sampler$data@forest.calibration)
expect_null(handed@cutPoints)
expect_null(handed@response.range)
expect_null(attr(sampler$model, "response.range", exact = TRUE))
expectRecordIsReader(sampler$copy(), "copy")
sampler$storeState()
expectRecordIsReader(reloadOf(sampler), "reload")
resolved <- sampler$data
resolved@cutPoints <- NULL
resolved@response.range <- NULL
viaNew <- new("dbartsSampler", sampler$control, sampler$model, resolved)
expectRecordIsReader(viaNew, "new")
expect_null(resolved@cutPoints)
expect_null(resolved@response.range)
# a family that reads no range records none, where its engine reports (0, 0)
yBinary <- as.integer(y > median(y))
yOrdinal <- cut(y, 3L, labels = FALSE)
noRange <- list(
  probit = dbarts::dbarts(x, yBinary, control = controlWith()),
  logistic = dbarts::dbarts(
    x,
    yBinary,
    family = "logistic",
    control = controlWith()
  ),
  ordinal = dbarts::dbarts(
    x,
    yOrdinal,
    family = "ordinal",
    control = controlWith()
  ),
  multinomial = dbarts::dbarts(
    x,
    factor(yOrdinal),
    family = "multinomial",
    control = controlWith()
  )
)
for (family in names(noRange)) {
  expect_null(noRange[[family]]$data@response.range, info = family)
  expect_identical(
    .Call(dbarts:::C_dbarts_bartcore_anchor, noRange[[family]]$pointer, NULL),
    c(0, 0),
    info = family
  )
  expectRecordIsReader(noRange[[family]], family)
}
expect_null(noRange$multinomial$data@forest.calibration)

# ---- copy and reload bitwise -----------------------------------------------

# after storeState the stored chain is the live one: a copy, a reload and a
# twin handed the state hold the same state and record and draw alike
caseWeights <- runif(n, 0.5, 2)
fixtures <- list(
  plain = function() make(),
  weights = function() {
    sampler <- dbarts::dbarts(
      x,
      y,
      weights = caseWeights,
      control = controlWith()
    )
    invisible(sampler$run(10L, 1L))
    sampler
  },
  factor = function() {
    frame <- data.frame(
      a = x[, 1L],
      f = factor(rep_len(c("p", "q", "r"), n)),
      o = ordered(rep_len(1:4, n))
    )
    sampler <- dbarts::dbarts(frame, y, control = controlWith())
    invisible(sampler$run(10L, 1L))
    sampler
  },
  probit = function() {
    sampler <- dbarts::dbarts(x, yBinary, control = controlWith())
    invisible(sampler$run(10L, 1L))
    sampler
  },
  forests = function() {
    z <- rep_len(0:1, n)
    sampler <- dbarts::dbarts(
      x,
      y,
      forests = list(
        dbarts::dbartsForests$forest(),
        dbarts::dbartsForests$forest(basis = factor(z))
      ),
      control = controlWith()
    )
    invisible(sampler$run(10L, 1L))
    sampler
  },
  multinomial = function() {
    sampler <- dbarts::dbarts(
      x,
      factor(yOrdinal),
      family = "multinomial",
      control = controlWith()
    )
    invisible(sampler$run(10L, 1L))
    sampler
  }
)
for (name in names(fixtures)) {
  sampler <- fixtures[[name]]()
  sampler$storeState()
  stored <- sampler$state
  copied <- sampler$copy()
  reloaded <- reloadOf(sampler)
  twin <- fixtures[[name]]()
  expect_true(twin$setState(stored), info = name)
  for (door in list(copied, reloaded)) {
    expect_identical(stateOf(door), stored, info = name)
    expect_identical(recordOf(door), recordOf(sampler), info = name)
    expectRecordIsReader(door, name)
  }
  expect_identical(copied$state, stored, info = name)
  expect_identical(reloaded$state, stored, info = name)
  draws <- twin$run(0L, 20L)
  expect_identical(copied$run(0L, 20L), draws, info = name)
  expect_identical(reloaded$run(0L, 20L), draws, info = name)
}
expect_true(!is.null(fixtures$forests()$data@forest.calibration))

# ---- copy is the live sampler ----------------------------------------------

# with updateState = FALSE the field is older than the chain: a copy has the
# chain as it stands - the three points, the current trees and generator -
# and its field is that chain; the source's field is left as it was
sampler <- make()
sampler$storeState()
older <- sampler$state
sampler$setCutPoints(three, 1L)
invisible(sampler$run(0L, 20L))
current <- stateOf(sampler)
expect_false(identical(liveOf(current), liveOf(older)))
copied <- sampler$copy()
expect_identical(lengths(cutPointsOf(copied)), c(3L, 100L))
expect_identical(stateOf(copied), current)
expect_identical(copied$state, current)
expect_identical(sampler$state, older)
expect_identical(recordOf(copied), recordOf(sampler))
expect_identical(sampler$data@cutPoints[[1L]], three)
# the copy draws as a second copy of the same chain, and continues the
# source, an install agreeing with a chain run on to the last few digits
expect_identical(copied$run(0L, 20L), sampler$copy()$run(0L, 20L))
continued <- sampler$copy()
expect_equal(continued$run(0L, 20L), sampler$run(0L, 20L), tolerance = 1e-12)
# and so for a multinomial sampler
sampler <- fixtures$multinomial()
sampler$storeState()
older <- sampler$state
sampler$setCutPoints(three, 1L)
invisible(sampler$run(0L, 5L))
copied <- sampler$copy()
expect_identical(stateOf(copied), stateOf(sampler))
expect_identical(lengths(cutPointsOf(copied)), c(3L, 100L))
expect_identical(sampler$state, older)

# ---- reload moves the stored trees -----------------------------------------

# the same object saved with the older field: a reload is on the three points
# its record holds, the stored splits moved onto them by value, as a live
# sampler forced to take that state has them; its field is then its engine's
sampler <- make()
sampler$storeState()
older <- sampler$state
sampler$setCutPoints(three, 1L)
invisible(sampler$run(0L, 5L))
reloaded <- reloadOf(sampler)
expect_identical(reloaded$state, older)
invisible(reloaded$getPointer())
expect_identical(lengths(cutPointsOf(reloaded)), c(3L, 100L))
twin <- make()
twin$setCutPoints(three, 1L)
expect_null(twin$setState(older, forceUpdate = TRUE))
expect_identical(liveOf(stateOf(reloaded)), liveOf(stateOf(twin)))
expect_identical(reloaded$state, stateOf(reloaded))
expect_false(identical(liveOf(reloaded$state), liveOf(older)))
# a reload handed its own field at its first call takes it, the field then
# its engine's; a copy made before the engine is rebuilt holds its engine's
# state as its field and takes it back
first <- reloadOf(sampler)
expect_true(first$setState(first$state))
expect_identical(first$state, stateOf(first))
deadCopy <- reloadOf(sampler)$copy()
expect_identical(deadCopy$state, stateOf(deadCopy))
expect_true(deadCopy$setState(deadCopy$state))
# a hand-edited grid out of order is sorted at the reload, the stored splits
# on it already
sampler$storeState()
edited <- reloadOf(sampler)
edited$data@cutPoints[[1L]] <- rev(three)
expect_identical(liveOf(stateOf(edited)), liveOf(sampler$state))
expect_identical(cutPointsOf(edited)[[1L]], three)

# ---- copy of a live sampler with a current field ---------------------------

# updateState = TRUE keeps the field the chain, so the copy's is the field
sampler <- make(updateState = TRUE)
sampler$setCutPoints(three, 1L)
invisible(sampler$run(0L, 5L))
copied <- sampler$copy()
expect_identical(copied$state, sampler$state)
expect_identical(stateOf(copied), sampler$state)
expect_identical(copied$run(0L, 5L), reloadOf(sampler)$run(0L, 5L))

# ---- dead pointer, no saved state ------------------------------------------

# setState on a reload with no field rebuilds on the record's grid, and the
# handed state goes onto it as onto a live sampler's
donor <- make()
handedState <- stateOf(donor)
sampler <- make()
sampler$setCutPoints(three, 1L)
sampler$storeState()
bare <- reloadOf(sampler)
bare$state <- NULL
expect_null(bare$setState(handedState, forceUpdate = TRUE))
expect_identical(lengths(cutPointsOf(bare)), c(3L, 100L))
# with no saved state a copy is refused as a reload is, run or not
for (runFirst in c(FALSE, TRUE)) {
  stateless <- reloadOf(make(run = runFirst))
  expect_null(stateless$state, info = runFirst)
  for (call in list(stateless$copy, function() stateless$run(0L, 1L))) {
    expect_error(call(), "without a stored state", info = runFirst)
  }
}

# ---- carried at creation ---------------------------------------------------

# a data object handed back carries the sampler's grid, weights and range to
# every creation door; a fresh object with the same rows derives its own
xWeighted <- cbind(v = rep_len(0:29, n)^2 / 100, b = x[, 2L])
sampler <- make(data = xWeighted, run = FALSE)
sampler$setCutPoints(three, 2L)
created <- sampler$data@cutPoints[[1L]]
sampler$setPredictor(2 * xWeighted[, 1L], 1L, updateCutPoints = "value")
expectRecordIsReader(sampler, "setPredictor")
# the range taken from one response and held while another is handed over
sampler$setResponse(3 * y + 10, updateScale = TRUE)
sampler$setResponse(y)
carried <- sampler$data
expect_identical(carried@cutPoints[[2L]], three)
expect_false(is.null(attr(carried@cutPoints[[1L]], "mass")))
expect_false(identical(carried@cutPoints[[1L]], created))
expect_identical(carried@response.range, range(3 * y + 10))
expect_identical(carried@y, y)
fresh <- dbarts::dbartsData(carried@x, carried@y)
spec <- dbarts::dbartsSpec(carried, control = controlWith())
doors <- list(
  new = new("dbartsSampler", sampler$control, sampler$model, carried),
  dbarts = dbarts::dbarts(carried, control = controlWith()),
  bart = dbarts::bart(
    carried,
    n.trees = 15L,
    n.samples = 2L,
    n.burn = 0L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    verbose = FALSE
  )$fit,
  dbartsSpec = new("dbartsSampler", spec$control, spec$model, spec$data)
)
for (door in names(doors)) {
  expect_identical(
    recordOf(doors[[door]])[c("cutPoints", "response.range")],
    recordOf(sampler)[c("cutPoints", "response.range")],
    info = door
  )
  expectRecordIsReader(doors[[door]], door)
}
derived <- dbarts::dbarts(fresh, control = controlWith())
expect_false(identical(derived$data@cutPoints, carried@cutPoints))
# the prior predictive draws on the sampler's grid
priorOf <- function(sampler) {
  set.seed(5L)
  dbarts::samplePriorPredictive(sampler, n.samples = 3L)
}
gridless <- carried
gridless@cutPoints <- NULL
gridless <- new("dbartsSampler", sampler$control, sampler$model, gridless)
expect_identical(gridless$data@response.range, carried@response.range)
expect_identical(priorOf(sampler), priorOf(doors$new))
expect_false(identical(priorOf(sampler), priorOf(gridless)))

# ---- xbart ignores ---------------------------------------------------------

# its folds derive their own grid and range from their training rows
crossValidated <- function(data) {
  set.seed(9L)
  dbarts::xbart(
    data,
    n.samples = 3L,
    n.burn = 2L,
    n.reps = 1L,
    n.trees = 5L,
    k = 2,
    n.threads = 1L,
    verbose = FALSE
  )
}
expect_identical(crossValidated(carried), crossValidated(fresh))

# ---- control anchor --------------------------------------------------------

# a control of several forests handed other data: the data's own anchor and
# unit, not those of the sampler the control came from
zBasis <- factor(rep_len(0:1, n))
severalOn <- function(response, control = controlWith()) {
  dbarts::dbarts(
    x,
    response,
    forests = list(
      dbarts::dbartsForests$forest(),
      dbarts::dbartsForests$forest(basis = zBasis)
    ),
    control = control
  )
}
first <- severalOn(y)
second <- severalOn(10 * y)
expect_false(identical(
  first$data@forest.calibration,
  second$data@forest.calibration
))
other <- second$data
other@forest.calibration <- NULL
other@cutPoints <- NULL
other@response.range <- NULL
handedControl <- new("dbartsSampler", first$control, first$model, other)
expect_identical(
  handedControl$data@forest.calibration,
  second$data@forest.calibration
)
expect_null(attr(first$control, "bartcore.forests")$anchor)
expect_null(attr(first$control, "bartcore.forests")$unit)
# the record holds where the response moves: a copy and a reload after a
# swap keep the anchor and unit of the response the sampler was made on
first$setResponse(10 * y)
first$storeState()
for (door in list(first$copy(), reloadOf(first))) {
  expect_identical(
    door$data@forest.calibration,
    severalOn(y)$data@forest.calibration
  )
  expect_identical(door$getLeafPrior(), first$getLeafPrior())
}

# ---- family as is ----------------------------------------------------------

# a range carried to a sampler of another family is used as it is where the
# family reads one, a count sampler reading its first number as the center
# and recording the window it holds, and is not read where it does not
gaussianRange <- make(run = FALSE)$data@response.range
counts <- rpois(n, exp(x[, 1L]))
countData <- dbarts::dbartsData(x, counts)
countData@response.range <- range(counts) + c(-1, 2)
countSampler <- dbarts::dbarts(
  countData,
  family = "nbinom",
  control = controlWith()
)
expect_identical(countSampler$data@response.range, min(counts) - 1 + c(0, 1))
binaryData <- dbarts::dbartsData(x, yBinary)
binaryData@response.range <- gaussianRange
binarySampler <- dbarts::dbarts(binaryData, control = controlWith())
expect_null(binarySampler$data@response.range)
plainBinary <- dbarts::dbarts(x, yBinary, control = controlWith())
expect_identical(binarySampler$run(0L, 5L), plainBinary$run(0L, 5L))

# ---- weights ride ----------------------------------------------------------

# a column of 30 distinct values weighs its points by their gaps, and a copy
# and a reload keep the weights from the record
sampler <- make(data = xWeighted)
mass <- attr(sampler$data@cutPoints[[1L]], "mass")
expect_identical(length(mass), 30L)
expect_identical(attr(stateOf(sampler), "cutMass")[[1L]], mass)
sampler$storeState()
for (door in list(sampler$copy(), reloadOf(sampler))) {
  expect_identical(attr(stateOf(door), "cutMass")[[1L]], mass)
  expect_identical(attr(door$data@cutPoints[[1L]], "mass"), mass)
}

# ---- misfit and order ------------------------------------------------------

# a record that fits no data is refused at creation, naming the slot and
# saying that NULL derives it; a grid out of order is sorted
recorded <- make(run = FALSE)$data
withSlot <- function(slot, value) {
  data <- recorded
  methods::slot(data, slot) <- value
  data
}
refusals <- list(
  list("cutPoints", recorded@cutPoints[1L], "one entry per predictor column"),
  list("cutPoints", list(c(0.5, 0.5), three), "repeats a point"),
  list("cutPoints", list("a", three), "finite numbers"),
  list("response.range", 1, "two finite numbers"),
  list("response.range", c(2, 1), "two finite numbers"),
  list("forest.calibration", c(anchor = -1, unit = 1), "positive")
)
for (refusal in refusals) {
  expect_error(
    dbarts::dbarts(withSlot(refusal[[1L]], refusal[[2L]])),
    paste0("'", refusal[[1L]], "' record .*", refusal[[3L]]),
    info = refusal[[3L]]
  )
  expect_error(
    dbarts::dbarts(withSlot(refusal[[1L]], refusal[[2L]])),
    "set it to NULL to derive it",
    info = refusal[[3L]]
  )
}
withMass <- three
attr(withMass, "mass") <- c(0, 1, 1, 2)
expect_error(
  dbarts::dbarts(withSlot("cutPoints", list(withMass, three))),
  "carries weights"
)
unsorted <- dbarts::dbarts(
  withSlot("cutPoints", list(rev(three), three)),
  control = controlWith()
)
expect_identical(unsorted$data@cutPoints, list(three, three))
expect_identical(cutPointsOf(unsorted), list(three, three))
# each point keeps its weight: 0.75, 0.5 and 0.25 weighing 1, 2 and 3
weighted <- rev(three)
attr(weighted, "mass") <- c(0, 1, 3, 6)
unsorted <- dbarts::dbarts(
  withSlot("cutPoints", list(weighted, three)),
  control = controlWith()
)
expect_identical(attr(unsorted$data@cutPoints[[1L]], "mass"), c(0, 3, 5, 6))
expect_identical(attr(stateOf(unsorted), "cutMass")[[1L]], c(0, 3, 5, 6))
# an ordered factor's grid is its levels', and a factor's entry is NULL
frame <- data.frame(a = x[, 1L], o = ordered(rep_len(1:4, n)))
ordinal <- dbarts::dbarts(frame, y, control = controlWith())
expect_identical(ordinal$data@cutPoints[[2L]], c(0.5, 1.5, 2.5))
moved <- ordinal$data
moved@cutPoints[[2L]] <- c(0.5, 1.5)
expect_error(
  new("dbartsSampler", ordinal$control, ordinal$model, moved),
  "an ordered factor, other cut points than those of its levels"
)
withFactor <- dbarts::dbarts(
  data.frame(a = x[, 1L], f = factor(rep_len(c("p", "q"), n))),
  y,
  control = controlWith()
)
expect_null(withFactor$data@cutPoints[[2L]])
pointed <- withFactor$data
pointed@cutPoints[[2L]] <- three
expect_error(
  new("dbartsSampler", withFactor$control, withFactor$model, pointed),
  "an unordered factor, cut points"
)
# and at setData, before anything is stored
sampler <- make()
before <- stateOf(sampler)
expect_error(
  sampler$setData(withSlot("response.range", 1)),
  "'response.range' record"
)
expect_identical(stateOf(sampler), before)
# a record carried in another column order is put in the sampler's with its
# columns and held
factorFrame <- data.frame(a = x[, 1L], f = factor(rep_len(c("p", "q"), n)))
inOrder <- dbarts::dbarts(factorFrame, y, control = controlWith())
reordered <- dbarts::dbarts(factorFrame[2:1], y, control = controlWith())
expect_null(reordered$data@cutPoints[[1L]])
inOrder$setData(reordered$data)
expectRecordIsReader(inOrder, "setData of a reordered record")
expect_null(inOrder$data@cutPoints[[2L]])

# ---- record-less object ----------------------------------------------------

# an object saved before the slots existed records nothing and derives, as a
# fresh object does
old <- recorded
attr(old, "cutPoints") <- NULL
attr(old, "response.range") <- NULL
attr(old, "forest.calibration") <- NULL
expect_false(methods::.hasSlot(old, "cutPoints"))
fromOld <- dbarts::dbarts(old, control = controlWith())
fromFresh <- dbarts::dbarts(dbarts::dbartsData(x, y), control = controlWith())
expect_identical(recordOf(fromOld), recordOf(fromFresh))
expect_identical(fromOld$run(0L, 5L), fromFresh$run(0L, 5L))

# ---- constant response -----------------------------------------------------

# several forests on a response with no spread record their anchor and unit,
# and a copy is made on them
flat <- suppressWarnings(severalOn(rep(3, n)))
expect_true(all(flat$data@forest.calibration > 0))
flatCopy <- flat$copy()
expect_identical(
  flatCopy$data@forest.calibration,
  flat$data@forest.calibration
)
expect_identical(flatCopy$getLeafPrior(), flat$getLeafPrior())
