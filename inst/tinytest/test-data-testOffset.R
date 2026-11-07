source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)
testX <- testData$X
colnames(testX) <- paste0("X.", seq_len(ncol(testX)))

# test that test offset fills in control logicals depending on specification
data <- dbarts::dbartsData(Z ~ X, testData, testX)

expect_null(data@offset)
expect_null(data@offset.test)
expect_identical(data@testUsesRegularOffset, NA)


data <- dbarts::dbartsData(Z ~ X, testData, testX, offset = 0.2)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], rep(0.2, 5))
expect_equal(data@testUsesRegularOffset, TRUE)
rm(data)


otherOffset <- 0.2 + 0.1
data <- dbarts::dbartsData(Z ~ X, testData, testX, offset = otherOffset)

expect_equal(data@offset[1:5], rep(0.3, 5))
expect_equal(data@offset.test[1:5], rep(0.3, 5))
expect_equal(data@testUsesRegularOffset, TRUE)
rm(data, otherOffset)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = NULL
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_null(data@offset.test)
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = 0.1
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], rep(0.1, 5))
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = offset + 0.1
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], rep(0.3, 5))
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = offset
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], rep(0.2, 5))
expect_equal(data@testUsesRegularOffset, TRUE)
rm(data)


set.seed(0)
otherOffset <- runif(nrow(testData$X))
data <- dbarts::dbartsData(Z ~ X, testData, testX, offset = otherOffset)

expect_equal(data@offset[1:5], otherOffset[1:5])
expect_equal(data@offset.test[1:5], otherOffset[1:5])
expect_equal(data@testUsesRegularOffset, TRUE)


expect_error(
  dbarts::dbartsData(Z ~ X, testData, testX[-1, ], offset = otherOffset),
  "'offset' was given as 'otherOffset', which cannot be evaluated on the rows of 'test'"
)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = otherOffset,
  offset.test = NULL
)

expect_equal(data@offset[1:5], otherOffset[1:5])
expect_null(data@offset.test)
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = otherOffset,
  offset.test = 0.2
)

expect_equal(data@offset[1:5], otherOffset[1:5])
expect_equal(data@offset.test[1:5], rep(0.2, 5))
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = otherOffset,
  offset.test = offset + 0.1
)

expect_equal(data@offset[1:5], otherOffset[1:5])
expect_equal(data@offset.test[1:5], otherOffset[1:5] + 0.1)
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = 0.2,
  offset.test = otherOffset
)

expect_equal(data@offset[1:5], rep(0.2, 5))
expect_equal(data@offset.test[1:5], otherOffset[1:5])
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data)


data <- dbarts::dbartsData(
  Z ~ X,
  testData,
  testX,
  offset = NULL,
  offset.test = otherOffset
)

expect_null(data@offset)
expect_equal(data@offset.test[1:5], otherOffset[1:5])
expect_equal(data@testUsesRegularOffset, FALSE)
rm(data, otherOffset)

rm(testData)

# a member named 'offset' inside the test offset expression is not the offset
# argument
set.seed(7)
nameData <- data.frame(y = rnorm(30), x = rnorm(30))
holder <- list(offset = rnorm(30), weights = runif(30) + 0.5)
rows <- 1:8
data <- dbarts::dbartsData(
  y ~ x,
  nameData,
  nameData[rows, ],
  offset = holder$offset,
  offset.test = holder$offset[rows]
)
expect_equal(data@offset.test, holder$offset[rows])
expect_equal(data@offset, holder$offset)

# and the argument itself still stands in for the name by itself
data <- dbarts::dbartsData(
  y ~ x,
  nameData,
  nameData[rows, ],
  offset = holder$offset,
  offset.test = offset[rows] + holder$offset[rows]
)
expect_equal(data@offset.test, 2 * holder$offset[rows])

data <- dbarts::dbartsData(y ~ x, nameData, weights = holder$weights)
expect_equal(data@weights, holder$weights)
rm(nameData, holder, rows, data)

# the substitution reaches every use of the name as a variable and skips only
# the member after $ or @
subTerm <- dbarts:::subTermInLanguage
repl <- quote(base)
expect_identical(
  subTerm(quote(s$offset[te]), quote(offset), repl),
  quote(s$offset[te])
)
expect_identical(subTerm(quote(s@offset), quote(offset), repl), quote(s@offset))
expect_identical(
  subTerm(quote(s@offset[te]), quote(offset), repl),
  quote(s@offset[te])
)
# the left side of $ and @ is a variable, whatever its position
expect_identical(subTerm(quote(offset$x), quote(offset), repl), quote(base$x))
expect_identical(subTerm(quote(offset@x), quote(offset), repl), quote(base@x))
expect_identical(
  subTerm(quote(offset[te]), quote(offset), repl),
  quote(base[te])
)
# the second argument of any other call (the third element) is a variable too
expect_identical(
  subTerm(quote(pmax(0, offset, 1)), quote(offset), repl),
  quote(pmax(0, base, 1))
)
expect_identical(
  subTerm(quote(f(a, offset$offset, b)), quote(offset), repl),
  quote(f(a, base$offset, b))
)

nameData <- data.frame(y = rnorm(30), x = rnorm(30))
base <- rnorm(30)
rows <- 1:8
data <- dbarts::dbartsData(
  y ~ x,
  nameData,
  nameData[rows, ],
  offset = base,
  offset.test = pmin(10, offset[rows], 20)
)
expect_equal(data@offset.test, base[rows])
rm(nameData, base, rows, data, subTerm, repl)

# a test offset given as a name bound to NULL is no test offset, as a literal
# NULL is: with a formula, through a wrapper, and in the matrix interface
set.seed(11)
nullData <- data.frame(y = rnorm(30), x = rnorm(30))
nullTest <- data.frame(x = rnorm(5))
noOffset <- NULL
data <- dbarts::dbartsData(y ~ x, nullData, nullTest, offset.test = noOffset)
expect_null(data@offset.test)
expect_identical(data@testUsesRegularOffset, FALSE)
data <- dbarts::dbartsData(
  as.matrix(nullData["x"]),
  nullData$y,
  as.matrix(nullTest),
  offset.test = noOffset
)
expect_null(data@offset.test)

# the nearest binding is the one read, a number further out being masked
o <- rnorm(5)
wrapper <- function(o) {
  dbarts::dbarts(
    y ~ x,
    nullData,
    test = nullTest,
    offset.test = o,
    control = dbarts::dbartsControl(
      verbose = FALSE,
      n.samples = 5L,
      n.burn = 5L,
      n.chains = 1L,
      n.threads = 1L
    ),
    n.trees = 5L
  )
}
expect_null(wrapper(NULL)$data@offset.test)
expect_equal(wrapper(o)$data@offset.test, o)

# and so for a count fit, with the training offset a NULL name as well
nullCounts <- data.frame(y = rpois(30, 3), x = nullData$x)
fit <- dbarts::bart(
  y ~ x,
  nullCounts,
  test = nullTest,
  family = "nbinom",
  offset = noOffset,
  offset.test = noOffset,
  n.trees = 5L,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_null(fit$offset.test)
expect_identical(dim(fit$yhat.test), c(5L, 5L))
rm(data, fit, wrapper, o, nullCounts, noOffset, nullTest, nullData)

# an offset.test forwarded through a wrapper's dots is read as written there:
# a name, an expression in 'offset', a column of 'test', through one wrapper
# or two and through dbarts(), each as the direct call reads it
set.seed(2)
dotsX <- matrix(runif(100), 50L, 2L)
dotsTestX <- dotsX[1:5, ]
dotsY <- rbinom(50L, 1L, 0.4)
dotsOffset <- rnorm(50L, sd = 0.1)
dotsTestOffset <- dotsOffset[1:5]
dotsFrame <- data.frame(
  y = dotsY,
  a = dotsX[, 1L],
  b = dotsX[, 2L],
  o = dotsOffset
)
dotsTest <- dotsFrame[1:5, ]
dotsFit <- function(...) {
  dbarts::bart(
    ...,
    n.trees = 5L,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 3L,
    verbose = FALSE
  )
}
dotsOuter <- function(...) dotsFit(...)
direct <- dotsFit(
  dotsX,
  dotsY,
  test = dotsTestX,
  offset = dotsOffset,
  offset.test = dotsTestOffset
)
expect_equal(direct$offset.test, dotsTestOffset)
viaDots <- function(fitter, ...) fitter(dotsX, dotsY, ...)
expect_identical(
  viaDots(
    dotsFit,
    test = dotsTestX,
    offset = dotsOffset,
    offset.test = dotsTestOffset
  )$yhat.test,
  direct$yhat.test
)
expect_identical(
  viaDots(
    dotsOuter,
    test = dotsTestX,
    offset = dotsOffset,
    offset.test = dotsTestOffset
  )$yhat.test,
  direct$yhat.test
)
# 'offset' in offset.test is the training offset, so the test set has as many
# rows
shifted <- dotsFit(
  dotsX,
  dotsY,
  test = dotsX,
  offset = dotsOffset,
  offset.test = offset + 1
)
expect_equal(shifted$offset.test, dotsOffset + 1)
expect_identical(
  viaDots(
    dotsOuter,
    test = dotsX,
    offset = dotsOffset,
    offset.test = offset + 1
  )$yhat.test,
  shifted$yhat.test
)
formulaDirect <- dotsFit(
  y ~ a + b,
  data = dotsFrame,
  test = dotsTest,
  offset = o,
  offset.test = o
)
expect_equal(formulaDirect$offset.test, dotsTestOffset)
expect_identical(
  dotsOuter(
    y ~ a + b,
    data = dotsFrame,
    test = dotsTest,
    offset = o,
    offset.test = o
  )$yhat.test,
  formulaDirect$yhat.test
)
# a name the wrapper's caller holds as a local of its own is read there, on
# both interfaces
heldOffset <- rep(0.3, 5L)
heldDirect <- dotsFit(
  dotsX,
  dotsY,
  test = dotsTestX,
  offset = dotsOffset,
  offset.test = heldOffset
)
heldFormulaDirect <- dotsFit(
  y ~ a + b,
  data = dotsFrame,
  test = dotsTest,
  offset = o,
  offset.test = heldOffset
)
rm(heldOffset)
callerWithLocal <- function() {
  localOffset <- rep(0.3, 5L)
  list(
    matrix = dotsFit(
      dotsX,
      dotsY,
      test = dotsTestX,
      offset = dotsOffset,
      offset.test = localOffset
    ),
    formula = dotsFit(
      y ~ a + b,
      data = dotsFrame,
      test = dotsTest,
      offset = o,
      offset.test = localOffset
    )
  )
}
heldLocally <- callerWithLocal()
expect_equal(heldLocally$matrix$offset.test, rep(0.3, 5L))
expect_identical(heldLocally$matrix$yhat.test, heldDirect$yhat.test)
expect_equal(heldLocally$formula$offset.test, rep(0.3, 5L))
expect_identical(heldLocally$formula$yhat.test, heldFormulaDirect$yhat.test)
samplerFit <- function(...) {
  dbarts::dbarts(
    ...,
    n.trees = 5L,
    control = dbarts::dbartsControl(
      n.samples = 5L,
      n.burn = 5L,
      n.chains = 1L,
      n.threads = 1L,
      verbose = FALSE
    )
  )
}
expect_equal(
  samplerFit(
    dotsX,
    dotsY,
    test = dotsTestX,
    offset = dotsOffset,
    offset.test = dotsTestOffset
  )$data@offset.test,
  dotsTestOffset
)
expect_equal(
  samplerFit(
    y ~ a + b,
    dotsFrame,
    test = dotsTest,
    offset = o,
    offset.test = o
  )$data@offset.test,
  dotsTestOffset
)
rm(
  dotsX,
  dotsTestX,
  dotsY,
  dotsOffset,
  dotsTestOffset,
  dotsFrame,
  dotsTest,
  dotsFit,
  dotsOuter,
  direct,
  viaDots,
  shifted,
  formulaDirect,
  heldDirect,
  heldFormulaDirect,
  callerWithLocal,
  heldLocally,
  samplerFit
)
