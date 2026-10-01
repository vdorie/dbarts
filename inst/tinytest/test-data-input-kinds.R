# Every entry point takes the predictor kinds a data frame takes, a sparse
# column survives a row-repeating subset, and arguments that cannot apply are
# refused or reported rather than silently dropped.

quietBart <- function(...) {
  suppressMessages(bart(
    ...,
    n.samples = 5L,
    n.burn = 5L,
    n.chains = 1L,
    n.trees = 5L,
    n.threads = 1L,
    verbose = FALSE
  ))
}

## A sparse frame column under a bootstrap subset reads the same values as its
## dense twin: repeated rows are not missing.
set.seed(1)
n <- 100L
sparse <- Matrix::sparseVector(
  x = rep(1, 20L),
  i = sort(sample(n, 20L)),
  length = n
)
sparseFrame <- data.frame(y = rnorm(n), x1 = runif(n))
sparseFrame$s <- sparse
denseFrame <- sparseFrame
denseFrame$s <- as.vector(sparse)
boot <- sample(n, n, replace = TRUE)
expect_true(anyDuplicated(boot) > 0L)
fromSparse <- dbartsData(y ~ ., sparseFrame, subset = boot)
fromDense <- dbartsData(y ~ ., denseFrame, subset = boot)
expect_false(anyNA(as.matrix(fromSparse@x)))
expect_equal(
  unname(as.matrix(fromSparse@x)[, "s"]),
  unname(fromDense@x[, "s"])
)
# a sparseFactor column likewise
factorFrame <- data.frame(y = rnorm(n))
factorFrame$f <- sparseFactor(
  factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.8, 0.1, 0.1))),
  reference = "a"
)
fromSparseFactor <- dbartsData(y ~ ., factorFrame, subset = boot)
expect_false(anyNA(as.matrix(fromSparseFactor@x)))
# with na.omit dropping a row as well
sparseFrame$x1[boot[1L]] <- NA
omitted <- dbartsData(y ~ ., sparseFrame, subset = boot, na.action = na.omit)
denseFrame$x1[boot[1L]] <- NA
omittedDense <- dbartsData(
  y ~ .,
  denseFrame,
  subset = boot,
  na.action = na.omit
)
expect_equal(
  unname(as.matrix(omitted@x)[, "s"]),
  unname(omittedDense@x[, "s"])
)

## Date, POSIXct and difftime columns are their numeric values, on both doors
## and in predict, as lm() and the indicators route read them.
set.seed(2)
timeFrame <- data.frame(
  y = rnorm(40L),
  z = as.Date("2020-01-01") + 1:40,
  p = as.POSIXct("2020-01-01", tz = "UTC") + 3600 * (1:40),
  d = as.difftime(1:40, units = "days")
)
fit <- quietBart(y ~ ., timeFrame, keepTrees = TRUE)
expect_equal(dim(predict(fit, timeFrame[1:2, ])), c(5L, 2L))
expect_equal(unname(fit$fit$data@x[, "z"]), as.double(timeFrame$z))
fit <- quietBart(timeFrame[-1L], timeFrame$y, keepTrees = TRUE)
expect_equal(dim(predict(fit, timeFrame[1:2, -1L])), c(5L, 2L))
posixlt <- list(y = rnorm(3L), w = as.POSIXlt(Sys.time() + 1:3))
expect_error(
  dbarts:::makeCategoricalModelMatrix(
    structure(list(w = posixlt$w), class = "data.frame", row.names = 1:3)
  ),
  "as.POSIXct"
)

## Integer and logical matrices are numbers, and any sparse Matrix class is
## taken, in fitting and in predict.
set.seed(3)
xInteger <- matrix(
  sample(0:5, 80L, TRUE),
  40L,
  2L,
  dimnames = list(NULL, c("a", "b"))
)
y <- rnorm(40L)
fitDouble <- quietBart(xInteger * 1.0, y, keepTrees = TRUE, seed = 1L)
fitInteger <- quietBart(xInteger, y, keepTrees = TRUE, seed = 1L)
expect_identical(fitInteger$yhat.train, fitDouble$yhat.train)
xLogical <- xInteger > 2L
fitLogical <- quietBart(xLogical, y, keepTrees = TRUE)
expect_equal(dim(predict(fitLogical, xLogical[1:3, ])), c(5L, 3L))
expect_equal(dim(predict(fitInteger, xInteger[1:3, ])), c(5L, 3L))
expect_equal(dim(quietBart(xInteger[, 1L], y)$yhat.train), c(5L, 40L))
sparseDouble <- Matrix::Matrix(xInteger * 1.0, sparse = TRUE)
fitCsc <- quietBart(sparseDouble, y, seed = 1L, sigest = 1)
fitTriplet <- quietBart(
  methods::as(sparseDouble, "TsparseMatrix"),
  y,
  seed = 1L,
  sigest = 1
)
fitRow <- quietBart(
  methods::as(sparseDouble, "RsparseMatrix"),
  y,
  seed = 1L,
  sigest = 1
)
expect_identical(fitTriplet$yhat.train, fitCsc$yhat.train)
expect_identical(fitRow$yhat.train, fitCsc$yhat.train)
expect_true(inherits(
  dbarts(Matrix::Matrix(xLogical, sparse = TRUE), y, sigest = 1),
  "dbartsSampler"
))

# a dense Matrix class is a plain matrix, and a sparseVector one sparse column
fitDense <- quietBart(Matrix::Matrix(xInteger * 1.0), y, seed = 1L)
expect_identical(fitDense$yhat.train, fitDouble$yhat.train)
sparseColumn <- methods::as(c(0, 0, 1, 0, 2)[rep(1:5, 8L)], "sparseVector")
fitVector <- quietBart(sparseColumn, y, seed = 1L, sigest = 1)
fitColumn <- quietBart(
  Matrix::Matrix(as.vector(sparseColumn), ncol = 1L, sparse = TRUE),
  y,
  seed = 1L,
  sigest = 1
)
expect_identical(fitVector$yhat.train, fitColumn$yhat.train)

## A POSIXlt column is refused by name at every door, before model.frame.
posixltFrame <- data.frame(y = rnorm(20L), x = runif(20L))
posixltFrame$w <- as.POSIXlt(
  as.POSIXct("2020-01-01", tz = "UTC") + 3600 * (1:20)
)
for (factors in c("categorical", "indicators")) {
  expect_error(quietBart(y ~ ., posixltFrame, factors = factors), "as.POSIXct")
  expect_error(
    quietBart(posixltFrame[-1L], posixltFrame$y, factors = factors),
    "as.POSIXct"
  )
  ctFrame <- posixltFrame
  ctFrame$w <- as.POSIXct(ctFrame$w)
  fit <- quietBart(y ~ ., ctFrame, factors = factors, keepTrees = TRUE)
  expect_error(predict(fit, posixltFrame[1:2, ]), "column 'w' is a POSIXlt")
}
# a POSIXlt column the formula does not use is left alone
expect_true(inherits(quietBart(y ~ x, posixltFrame), "bart"))

## A numeric column where training had a factor is refused by name, at every
## entrance a data frame reaches.
set.seed(4)
groupFrame <- data.frame(
  y = rnorm(40L),
  x1 = runif(40L),
  g = factor(rep(letters[1:4], 10L)),
  o = factor(rep(c("lo", "hi"), 20L), levels = c("lo", "hi"), ordered = TRUE)
)
fit <- quietBart(y ~ ., groupFrame, keepTrees = TRUE)
newFrame <- groupFrame[1:3, ]
newFrame$g <- 1:3
expect_error(predict(fit, newFrame), "test column 'g' is integer")
newFrame <- groupFrame[1:3, ]
newFrame$o <- c(TRUE, FALSE, TRUE)
expect_error(predict(fit, newFrame), "test column 'o' is logical")
newFrame <- groupFrame[1:3, ]
newFrame$g <- c(0, 1, 2)
expect_error(
  quietBart(y ~ ., groupFrame, test = newFrame),
  "test column 'g' is numeric"
)
expect_error(
  dbartsData(y ~ ., groupFrame, test = newFrame),
  "test column 'g' is numeric"
)
sampler <- dbarts(y ~ ., groupFrame)
expect_error(sampler$setTestPredictor(newFrame[-1L]), "test column 'g'")
# the indicators route refuses the same, on both doors
for (indicatorFit in list(
  quietBart(y ~ ., groupFrame, factors = "indicators", keepTrees = TRUE),
  quietBart(
    groupFrame[-1L],
    groupFrame$y,
    factors = "indicators",
    keepTrees = TRUE
  )
)) {
  newFrame <- groupFrame[1:3, ]
  newFrame$g <- 1:3
  expect_error(predict(indicatorFit, newFrame), "test column 'g' is integer")
}
newFrame <- groupFrame[1:3, ]
newFrame$g <- c(0, 1, 2)
# factor and character test columns still map by label
newFrame$g <- c("a", "b", "c")
expect_equal(dim(predict(fit, newFrame)), c(5L, 3L))

## A dbartsData object reports every argument that cannot reach it.
data <- dbartsData(y ~ ., groupFrame)
expect_warning(
  dbarts(data, factors = "indicators", weights = rep(2, 40L), subset = 1:10),
  "'subset', 'weights', 'factors'"
)
expect_warning(
  dbarts(data, na.action = na.omit),
  "'na.action'"
)
# a front door that stamps its own defaults reports nothing the caller did not
# write
expect_silent(suppressMessages(bart(
  data,
  n.samples = 5L,
  n.burn = 5L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)))

## n.cuts recycles a short vector and refuses a long one.
x3 <- matrix(runif(120L), 40L, 3L)
sampler <- dbarts(x3, y, control = dbartsControl(n.cuts = c(5L, 6L)))
expect_equal(sampler$data@n.cuts, c(5L, 6L, 5L))
expect_error(
  quietBart(x3[, 1:2], y, n.cuts = 5:7),
  "'n.cuts' has 3 values but the model has 2 predictor columns"
)
# the BayesTree door names its own argument
expect_error(
  suppressWarnings(bartBT(
    x3[, 1:2],
    y,
    numcut = 5:7,
    ndpost = 5L,
    nskip = 5L,
    verbose = FALSE
  )),
  "'numcut' has 3 values"
)
expect_error(
  xbart(x3[, 1:2], y, n.cuts = 5:7, n.reps = 1L, n.trees = 5L, verbose = FALSE),
  "'n.cuts' has 3 values"
)

rm(
  quietBart,
  n,
  sparse,
  sparseFrame,
  denseFrame,
  boot,
  fromSparse,
  fromDense,
  factorFrame,
  fromSparseFactor,
  omitted,
  omittedDense,
  timeFrame,
  fit,
  posixlt,
  xInteger,
  y,
  fitInteger,
  fitDouble,
  xLogical,
  fitLogical,
  sparseDouble,
  fitCsc,
  fitTriplet,
  fitRow,
  groupFrame,
  newFrame,
  sampler,
  data,
  x3
)
