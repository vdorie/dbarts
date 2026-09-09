# The na.action surface (dec-B108): the front door takes a standard
# na.action whose default drops rows with a missing RESPONSE and keeps rows
# with missing predictors, on the formula path and on the matrix (x, y)
# pair alike; na.omit, na.exclude, na.fail and na.pass keep their base
# meaning; and what an "exclude"-class record dropped pads back through
# stats::naresid, so fitted() and residuals() are at the caller's own
# length.

set.seed(311L)
n <- 120L
x1 <- runif(n)
x2 <- runif(n)
y <- 2 * x1 + x2 + rnorm(n, 0, 0.25)

missingResponse <- c(4L, 17L, 98L)
missingPredictor <- c(9L, 42L)
y[missingResponse] <- NA_real_
x1[missingPredictor] <- NA_real_

df <- data.frame(x1 = x1, x2 = x2, y = y)
xMat <- cbind(x1 = x1, x2 = x2)

quick <- list(
  n.trees = 10L,
  n.samples = 15L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
fitFormula <- function(...) {
  do.call(dbarts::bart, c(list(y ~ x1 + x2, df), quick, list(...)))
}
fitMatrix <- function(...) {
  do.call(dbarts::bart, c(list(xMat, y), quick, list(...)))
}

kept <- n - length(missingResponse)
complete <- n - length(missingResponse) - length(missingPredictor)

# --- the default, on both paths --------------------------------------------

# the package default is exported and documented, and is the formal's own
# default on all three entry points
expect_true("na.keepPredictors" %in% getNamespaceExports("dbarts"))
for (entry in list(dbarts::bart, dbarts::dbarts, dbarts::dbartsData)) {
  expect_true("na.action" %in% names(formals(entry)))
  expect_identical(
    deparse(formals(entry)[["na.action"]]),
    "dbarts::na.keepPredictors"
  )
}
# 'missing' is retired outright: no formal, and no tombstone to reach
for (entry in list(
  dbarts::bart,
  dbarts::dbarts,
  dbarts::dbartsData,
  dbarts::xbart
)) {
  expect_false("missing" %in% names(formals(entry)))
}
expect_false(
  "missing" %in%
    vapply(
      dbarts:::dbartsTombstones,
      function(e) e$name,
      character(1L)
    )
)

dataFormula <- dbarts::dbartsData(y ~ x1 + x2, df)
dataMatrix <- dbarts::dbartsData(xMat, y)
for (data in list(dataFormula, dataMatrix)) {
  expect_equal(length(data@y), kept)
  expect_equal(nrow(data@x), kept)
  # the missing predictors that survived are still there to be routed
  expect_true(anyNA(data@x[, "x1"]))
  expect_inherits(data@na.action, "exclude")
  expect_equal(
    unclass(data@na.action),
    missingResponse,
    check.attributes = FALSE
  )
}
# the two paths keep the same rows
expect_equal(unname(dataFormula@y), unname(dataMatrix@y))

# --- the base functions keep their base meaning ----------------------------

for (build in list(
  function(...) dbarts::dbartsData(y ~ x1 + x2, df, ...),
  function(...) dbarts::dbartsData(xMat, y, ...)
)) {
  expect_equal(length(build(na.action = na.omit)@y), complete)
  expect_equal(length(build(na.action = na.exclude)@y), complete)
  expect_error(build(na.action = na.fail), "missing values")
  # na.pass keeps every row, and the response completeness check then
  # refuses the missing response, exactly as before this argument existed
  expect_error(build(na.action = na.pass), "response contains missing values")
  # na.omit records what it dropped, but its record pads nothing
  expect_inherits(build(na.action = na.omit)@na.action, "omit")
  expect_inherits(build(na.action = na.exclude)@na.action, "exclude")
}

# a complete data set loses no row and records nothing, under any of them
dfComplete <- data.frame(x1 = runif(n), x2 = x2, y = runif(n))
for (action in list(na.keepPredictors, na.omit, na.exclude, na.fail, na.pass)) {
  data <- dbarts::dbartsData(y ~ x1 + x2, dfComplete, na.action = action)
  expect_equal(length(data@y), n)
  expect_null(data@na.action)
}

# --- padding ---------------------------------------------------------------

# the default and na.exclude pad the training fits back to the caller's own
# row count; na.omit does not, which is na.omit's own contract
fitDefault <- fitFormula()
expect_equal(length(fitDefault$y), kept)
expect_equal(length(fitted(fitDefault)), n)
expect_equal(length(residuals(fitDefault)), n)
expect_true(all(is.na(fitted(fitDefault)[missingResponse])))
expect_false(any(is.na(fitted(fitDefault)[-missingResponse])))
expect_equal(length(fitDefault$yhat.train.mean), n)
# the test side never lost a row, so it never pads
fitTest <- fitFormula(test = df)
expect_equal(length(fitted(fitTest, sample = "test")), n)
# an interval keeps the padded row count too
expect_equal(nrow(fitted(fitDefault, ci.level = 0.9)), n)

fitOmit <- fitFormula(na.action = na.omit)
expect_equal(length(fitted(fitOmit)), complete)
fitExclude <- fitFormula(na.action = na.exclude)
expect_equal(length(fitted(fitExclude)), n)
expect_true(all(is.na(fitted(fitExclude)[sort(c(
  missingResponse,
  missingPredictor
))])))

# the matrix path pads identically
fitMatrixDefault <- fitMatrix()
expect_equal(length(fitted(fitMatrixDefault)), n)

# --- the shape a saved-fit consumer reads ----------------------------------

# insight's case: a data frame carrying NA in the response, fit at the
# default, whose fitted() is as long as the data
insightFrame <- data.frame(y = y, x1 = x1, x2 = x2)
insightFit <- do.call(
  dbarts::bart,
  c(list(y ~ x1 + x2, insightFrame), quick)
)
expect_equal(length(fitted(insightFit)), nrow(insightFrame))
expect_equal(length(insightFit$na.action), length(missingResponse))

# --- the legacy door -------------------------------------------------------

# bartBT keeps 0.9-34's row rule, na.omit, and takes no na.action of its own
expect_false("na.action" %in% names(formals(dbarts::bartBT)))
fitLegacy <- dbarts::bartBT(
  xMat,
  y,
  ndpost = 15L,
  nskip = 5L,
  ntree = 10L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
)
expect_equal(length(fitLegacy$y), complete)
expect_equal(length(fitted(fitLegacy)), complete)

# --- na.keepPredictors itself ----------------------------------------------

# a model frame with no response loses no row, whatever is missing in it
frameNoResponse <- stats::model.frame(
  ~ x1 + x2,
  df,
  na.action = dbarts::na.keepPredictors
)
expect_equal(nrow(frameNoResponse), n)
expect_null(attr(frameNoResponse, "na.action"))
