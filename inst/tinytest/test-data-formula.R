source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that formula specification creates valid objects
trainData_df <- as.data.frame(testData)
set.seed(0)
trainData_df$weights <- runif(nrow(trainData_df))
trainData_df$offset <- rnorm(nrow(trainData_df))

modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10

expect_inherits(
  dbarts::dbartsData(modelFormula, trainData_df),
  "dbartsData"
)
# each channel named in 'data' must ARRIVE on the object, not merely be
# accepted by the call
dat <- dbarts::dbartsData(modelFormula, trainData_df, weights = weights)
expect_equal(dat@weights, trainData_df$weights)
dat <- dbarts::dbartsData(modelFormula, trainData_df, offset = offset)
expect_equal(dat@offset, trainData_df$offset)
dat <- dbarts::dbartsData(
  modelFormula,
  trainData_df,
  weights = weights,
  offset = offset
)
expect_equal(dat@weights, trainData_df$weights)
expect_equal(dat@offset, trainData_df$offset)
# and the subset applies to all three, in row order
dat <- dbarts::dbartsData(
  modelFormula,
  trainData_df,
  subset = 1:10,
  weights = weights,
  offset = offset
)
expect_equal(dat@y, trainData_df$y[1:10])
expect_equal(dat@weights, trainData_df$weights[1:10])
expect_equal(dat@offset, trainData_df$offset[1:10])
rm(dat)

testData_df <- trainData_df[1:20, ]
expect_inherits(
  dbarts::dbartsData(
    modelFormula,
    trainData_df,
    test = testData_df,
    weights = weights
  ),
  "dbartsData"
)

rm(testData_df, trainData_df)


# test that test argument creates valid objects
## test when is embedded in passed data
testData$test <- testData$x[11:20, ]
expect_inherits(
  dbarts::dbartsData(y ~ x, testData, test),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(y ~ x, testData, testData$test),
  "dbartsData"
)

## test when is in environment of formula
test <- testData$test
testData$test <- NULL
expect_inherits(
  dbarts::dbartsData(y ~ x, testData, test),
  "dbartsData"
)
expect_inherits(
  dbarts::dbartsData(y ~ x, testData, testData$x[11:20, ]),
  "dbartsData"
)
rm(test)


# test that test weights are created correctly
trainData <- as.data.frame(testData)
trainData$weights <- runif(nrow(trainData))

modelFormula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10

testData_df <- trainData[1:20, ]
data <- dbarts::dbartsData(
  modelFormula,
  trainData,
  test = testData_df,
  weights = weights
)
expect_inherits(data, "dbartsData")
expect_equal(data@weights, trainData$weights)
expect_equal(data@weights.test, testData_df$weights)
rm(data, testData_df, modelFormula, trainData)

# a one-sided formula (no response) is a legitimate dbartsData() call - the
# composed-sampler response is set later - unlike the fitting entry points
# (bart/bart/xbart), which refuse it by name instead
oneSidedData <- dbarts::dbartsData(~x, testData)
expect_inherits(oneSidedData, "dbartsData")
expect_true(all(oneSidedData@y == 0))
rm(oneSidedData)

rm(testData)

# 'data' and 'subset' are each evaluated once for a fit, at every door. A
# 'data' that draws its rows is then one draw for 'subset', written over its
# columns, and for the response, the predictors, the weights and the offset:
# every kept row has x1 > 0.5 and is whole. An `id` predictor names the row
# of `onceData` each kept row is.
set.seed(11)
onceRows <- 60L
onceData <- data.frame(
  id = as.double(seq_len(onceRows)),
  x1 = runif(onceRows),
  w = runif(onceRows, 0.5, 2),
  o = rnorm(onceRows, sd = 0.1)
)
onceData$y <- onceData$x1 + rnorm(onceRows, sd = 0.1)
dataDraws <- list()
drawData <- function(size = onceRows, replace = FALSE) {
  rows <- sample(onceRows, size, replace)
  dataDraws[[length(dataDraws) + 1L]] <<- rows
  onceData[rows, ]
}
onceControl <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 3L,
  n.samples = 2L,
  n.burn = 0L,
  updateState = FALSE,
  verbose = FALSE
)
onceDoors <- list(
  dbartsData = function(data) {
    eval(bquote(dbartsData(
      y ~ x1 + id,
      .(data),
      subset = x1 > 0.5,
      weights = w,
      offset = o
    )))
  },
  dbarts = function(data) {
    eval(bquote(dbarts(
      y ~ x1 + id,
      .(data),
      subset = x1 > 0.5,
      weights = w,
      offset = o,
      control = onceControl
    )))$data
  },
  bart = function(data) {
    eval(bquote(bart(
      y ~ x1 + id,
      .(data),
      subset = x1 > 0.5,
      weights = w,
      offset = o,
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 3L,
      n.samples = 2L,
      n.burn = 0L,
      verbose = FALSE,
      samplerOnly = TRUE
    )))$data
  }
)
onceDraws <- list(
  shuffled = quote(drawData()),
  subsample = quote(drawData(40L)),
  bootstrap = quote(drawData(replace = TRUE))
)
for (door in names(onceDoors)) {
  for (kind in names(onceDraws)) {
    info <- paste(door, kind)
    dataDraws <- list()
    built <- onceDoors[[door]](onceDraws[[kind]])
    expect_identical(length(dataDraws), 1L, info = info)
    drawn <- dataDraws[[1L]]
    rows <- as.integer(built@x[, "id"])
    expect_identical(rows, drawn[onceData$x1[drawn] > 0.5], info = info)
    expect_identical(built@y, onceData$y[rows], info = info)
    expect_identical(built@weights, onceData$w[rows], info = info)
    expect_identical(built@offset, onceData$o[rows], info = info)
  }
}
# the legacy door and the cross-validation read 'data' once too
dataDraws <- list()
invisible(bartBT(
  y ~ x1 + id,
  drawData(),
  ntree = 3L,
  ndpost = 2L,
  nskip = 0L,
  verbose = FALSE
))
expect_identical(length(dataDraws), 1L)
dataDraws <- list()
invisible(xbart(
  y ~ x1 + id,
  drawData(),
  subset = x1 > 0.5,
  n.samples = 4L,
  n.reps = 2L,
  n.burn = c(2L, 1L),
  n.test = 2,
  n.trees = 3L,
  n.threads = 1L
))
expect_identical(length(dataDraws), 1L)
# and the call a fit keeps shows what the caller wrote, not the value read
keptCall <- dbarts(y ~ x1 + id, drawData(), control = onceControl)$control@call
expect_identical(keptCall$data, quote(drawData()))
keptCall <- bart(
  y ~ x1 + id,
  drawData(),
  subset = x1 > 0.5,
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 3L,
  n.samples = 2L,
  n.burn = 0L,
  verbose = FALSE
)$call
expect_identical(keptCall$data, quote(drawData()))
expect_identical(keptCall$subset, quote(x1 > 0.5))
# 'weights' is evaluated once as well, where the data is a data frame
weightReadings <- 0L
readWeights <- function(w) {
  weightReadings <<- weightReadings + 1L
  w
}
weighted <- dbartsData(y ~ x1 + id, onceData, weights = readWeights(w))
expect_identical(weightReadings, 1L)
expect_identical(weighted@weights, onceData$w)
# With the matrix interface the predictors and the response are each
# evaluated once, in the order written: a response that names the rows the
# predictors drew is the response of the rows the fit holds, at every door.
onceX <- as.matrix(onceData[c("x1", "id")])
xDraws <- list()
drawX <- function() {
  rows <- sample(onceRows, 40L)
  xDraws[[length(xDraws) + 1L]] <<- rows
  onceX[rows, ]
}
ofDrawn <- function(values) values[xDraws[[length(xDraws)]]]
onceCounts <- rpois(onceRows, 3)
onceLevels <- factor(rep_len(c("a", "b", "c"), onceRows))
matrixDoors <- list(
  dbartsData = function() dbartsData(drawX(), ofDrawn(onceData$y)),
  dbarts = function() {
    dbarts(drawX(), ofDrawn(onceData$y), control = onceControl)$data
  },
  aft = function() {
    dbarts(
      drawX(),
      cbind(exp(ofDrawn(onceData$y)), 1),
      family = "aft",
      control = onceControl
    )$data
  },
  bart = function() {
    # written in the call: an argument handed through a helper's dots is a
    # promise already forced, and a second reading inside bart() would not show
    bart(
      drawX(),
      ofDrawn(onceData$y),
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 3L,
      n.samples = 2L,
      n.burn = 0L,
      verbose = FALSE,
      samplerOnly = TRUE
    )$data
  },
  nbinom = function() {
    # written in the call: an argument handed through a helper's dots is a
    # promise already forced, and a second reading inside bart() would not show
    bart(
      drawX(),
      ofDrawn(onceCounts),
      family = "nbinom",
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 3L,
      n.samples = 2L,
      n.burn = 0L,
      verbose = FALSE,
      samplerOnly = TRUE
    )$data
  },
  multinomial = function() {
    # written in the call: an argument handed through a helper's dots is a
    # promise already forced, and a second reading inside bart() would not show
    bart(
      drawX(),
      ofDrawn(onceLevels),
      family = "multinomial",
      n.chains = 1L,
      n.threads = 1L,
      n.trees = 3L,
      n.samples = 2L,
      n.burn = 0L,
      verbose = FALSE,
      samplerOnly = TRUE
    )$data
  },
  bartBT = function() {
    bartBT(
      drawX(),
      ofDrawn(onceData$y),
      ntree = 3L,
      ndpost = 2L,
      nskip = 0L,
      verbose = FALSE,
      sampleronly = TRUE
    )$data
  }
)
for (door in names(matrixDoors)) {
  xDraws <- list()
  built <- matrixDoors[[door]]()
  expect_identical(length(xDraws), 1L, info = door)
  rows <- as.integer(built@x[, "id"])
  expect_identical(rows, xDraws[[1L]], info = door)
  expect_identical(
    switch(
      door,
      multinomial = max.col(built@counts),
      nbinom = as.integer(built@y),
      built@y
    ),
    switch(
      door,
      multinomial = as.integer(onceLevels)[rows],
      nbinom = onceCounts[rows],
      aft = log(exp(onceData$y[rows])),
      onceData$y[rows]
    ),
    info = door
  )
}
xDraws <- list()
invisible(xbart(
  drawX(),
  ofDrawn(onceData$y),
  n.samples = 4L,
  n.reps = 2L,
  n.burn = c(2L, 1L),
  n.test = 2,
  n.trees = 3L,
  n.threads = 1L
))
expect_identical(length(xDraws), 1L)
# The values are handed on in calls that name no row: after an error inside
# a fit, the calls on the stack are no longer for 10000 rows than for 100
stackSize <- function(fit) {
  calls <- NULL
  tryCatch(
    withCallingHandlers(fit, error = function(e) calls <<- sys.calls()),
    error = function(e) NULL
  )
  sum(nchar(unlist(lapply(calls, deparse))))
}
stackSizes <- vapply(
  c(100L, 10000L),
  function(numRows) {
    big <- data.frame(x1 = runif(numRows), x2 = runif(numRows))
    big$y <- rnorm(numRows)
    bigX <- as.matrix(big[c("x1", "x2")])
    c(
      formula = stackSize(dbarts(
        y ~ x1 + nosuch,
        big,
        subset = x1 > 0.2,
        weights = x2,
        offset = x2,
        control = onceControl
      )),
      bart = stackSize(bart(y ~ x1 + nosuch, big, verbose = FALSE)),
      matrix = stackSize(dbarts(
        bigX,
        big$y,
        subset = seq_len(50L),
        weights = c(1, 2),
        control = onceControl
      )),
      bartMatrix = stackSize(bart(bigX, big$y, weights = 1:2, verbose = FALSE))
    )
  },
  numeric(4L)
)
expect_true(all(stackSizes > 0))
expect_true(all(stackSizes < 20000))
expect_true(all(abs(stackSizes[, 2L] - stackSizes[, 1L]) < 100))
rm(weightReadings, readWeights, weighted, onceX, xDraws, drawX, ofDrawn)
rm(onceCounts, onceLevels, matrixDoors, stackSize, stackSizes)
# 'subset' selects training rows and never rows of 'test'
onceTest <- onceData[1:25, ]
cutBeside <- dbartsData(
  y ~ x1 + id,
  onceData,
  test = onceTest,
  weights = w,
  subset = 25:1
)
expect_identical(as.integer(cutBeside@x[, "id"]), 25:1)
expect_identical(as.integer(cutBeside@x.test[, "id"]), 1:25)
expect_identical(cutBeside@weights.test, onceTest$w)
rm(onceRows, onceData, dataDraws, drawData, onceControl, onceDoors, onceDraws)
rm(door, kind, info, built, drawn, rows, keptCall, onceTest, cutBeside)
