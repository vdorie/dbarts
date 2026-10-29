# One refusal for test values at a factor level the training data never
# had, on every door that codes a factor against the training levels: it
# names the factor, at most five levels and a count of the rest (a factor
# with thousands of unseen levels once ran the message past R's length limit
# and was cut mid-list), and says what to do. The remedy named is the one
# each door can take: a training factor that declares the level fits and
# predicts at it where factors are kept whole, and does not under indicators,
# where a level no training row takes has no column.

set.seed(5210L)
nUL <- 80L
seenUL <- paste0("L", 1:5)
dUL <- data.frame(
  f = factor(sample(seenUL, nUL, TRUE), levels = seenUL),
  z = rnorm(nUL)
)
dUL$y <- rnorm(nUL)

# a test frame with one seen level and the given unseen ones
testUL <- function(unseen) {
  data.frame(f = factor(c("L1", unseen)), z = 0)
}
fewUL <- function(count) paste0("new", seq_len(count))

fitUL <- list(
  n.trees = 5L,
  n.samples = 6L,
  n.burn = 4L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
formulaFitUL <- function(data, ...) {
  do.call(
    bart,
    c(list(y ~ f + z, data), fitUL, list(...))
  )
}
btFitUL <- function(data, test) {
  suppressWarnings(bartBT(
    x.train = data["f"],
    y.train = data$y,
    x.test = test["f"],
    ndpost = 6L,
    nskip = 4L,
    ntree = 5L,
    nchain = 1L,
    verbose = FALSE
  ))
}

remedyWhole <- paste0(
  "refit with a training factor that declares %s, ",
  "factor(x, levels = c(levels(x), extra))"
)
oneWhole <- paste0(
  "test data factor 'f' has a level not present in the training data: ",
  "'new1'. Recode it to a training level, drop the rows that take it, or ",
  sprintf(remedyWhole, "it")
)
oneIndicators <- paste0(
  "test data factor 'f' has a level not present in the training data: ",
  "'new1'. Recode it to a training level or drop the rows that take it"
)
manyWhole <- paste0(
  "test data factor 'f' has levels not present in the training data: ",
  "'new1', 'new2', 'new3', 'new4', 'new5' and 9995 more. Recode them to ",
  "training levels, drop the rows that take them, or ",
  sprintf(remedyWhole, "them")
)
manyIndicators <- paste0(
  "test data factor 'f' has levels not present in the training data: ",
  "'new1', 'new2', 'new3', 'new4', 'new5' and 9995 more. Recode them to ",
  "training levels or drop the rows that take them"
)

messageOf <- function(expr) {
  tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
}

oneNew <- testUL(fewUL(1L))
manyNew <- testUL(fewUL(10000L))

# --- the formula door and predict, factors kept whole ---
expect_identical(
  messageOf(formulaFitUL(dUL, test = oneNew)),
  oneWhole
)
expect_identical(
  messageOf(formulaFitUL(dUL, test = manyNew)),
  manyWhole
)
fitWhole <- formulaFitUL(dUL, keepTrees = TRUE)
expect_identical(messageOf(predict(fitWhole, oneNew)), oneWhole)
expect_identical(messageOf(predict(fitWhole, manyNew)), manyWhole)
# the message is short however many levels there are
expect_true(nchar(manyWhole) < 400L)

# --- the names given: all of up to five, then a count of the rest ---
expect_identical(
  messageOf(predict(fitWhole, testUL(fewUL(2L)))),
  paste0(
    "test data factor 'f' has levels not present in the training data: ",
    "'new1' and 'new2'. Recode them to training levels, drop the rows that ",
    "take them, or ",
    sprintf(remedyWhole, "them")
  )
)
expect_true(grepl(
  "'new1', 'new2', 'new3', 'new4' and 'new5'. Recode",
  messageOf(predict(fitWhole, testUL(fewUL(5L)))),
  fixed = TRUE
))
expect_true(grepl(
  "'new1', 'new2', 'new3', 'new4', 'new5' and 1 more. Recode",
  messageOf(predict(fitWhole, testUL(fewUL(6L)))),
  fixed = TRUE
))
# a character column, and a value missing beside the unseen one
charNew <- oneNew
charNew$f <- as.character(charNew$f)
charNew$f[1L] <- NA
expect_identical(messageOf(predict(fitWhole, charNew)), oneWhole)

# --- the remedy named works where it is named: a training factor declaring
# the levels fits, and predicts at them ---
declaredUL <- dUL
declaredUL$f <- factor(as.character(dUL$f), levels = c(seenUL, fewUL(3L)))
declaredNew <- testUL(fewUL(3L))
fitDeclared <- formulaFitUL(declaredUL, test = declaredNew, keepTrees = TRUE)
expect_identical(dim(fitDeclared$yhat.test), c(6L, 4L))
expect_true(all(is.finite(predict(fitDeclared, declaredNew))))
declaredSampler <- dbarts(
  y ~ f + z,
  declaredUL,
  control = dbartsControl(
    n.trees = 5L,
    n.samples = 6L,
    n.burn = 4L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
)
expect_silent(declaredSampler$setTestPredictor(declaredNew))
# on a sampler without the levels, setTestPredictor is refused the same way
plainSampler <- dbarts(
  y ~ f + z,
  dUL,
  control = dbartsControl(
    n.trees = 5L,
    n.samples = 6L,
    n.burn = 4L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  )
)
expect_identical(
  messageOf(plainSampler$setTestPredictor(oneNew)),
  oneWhole
)
# a declared level no row takes is a level of the factor, not an unseen one
expect_identical(
  dim(predict(fitDeclared, testUL(c("L2", "new2")))),
  c(6L, 3L)
)

# --- indicators: the same declared factor is refused, so recoding or
# dropping rows is what is named, and a declared level says why ---
expect_identical(
  messageOf(formulaFitUL(dUL, test = oneNew, factors = "indicators")),
  oneIndicators
)
fitIndicators <- formulaFitUL(dUL, keepTrees = TRUE, factors = "indicators")
expect_identical(messageOf(predict(fitIndicators, oneNew)), oneIndicators)
expect_identical(messageOf(predict(fitIndicators, manyNew)), manyIndicators)
fitIndicatorsDeclared <- formulaFitUL(
  declaredUL,
  keepTrees = TRUE,
  factors = "indicators"
)
expect_identical(
  messageOf(predict(fitIndicatorsDeclared, oneNew)),
  paste0(
    oneIndicators,
    " (a declared level that no training row takes has no indicator ",
    "column)"
  )
)
# a level some training row takes predicts
expect_identical(dim(predict(fitIndicatorsDeclared, testUL("L2"))), c(6L, 2L))

# --- the BayesTree door, which is always indicators ---
expect_identical(messageOf(btFitUL(dUL, oneNew)), oneIndicators)
expect_identical(messageOf(btFitUL(dUL, manyNew)), manyIndicators)
expect_identical(
  messageOf(suppressWarnings(bart(
    x.train = dUL["f"],
    y.train = dUL$y,
    x.test = manyNew["f"],
    ndpost = 6L,
    nskip = 4L,
    ntree = 5L,
    verbose = FALSE
  ))),
  manyIndicators
)
fitBT <- suppressWarnings(bartBT(
  x.train = dUL["f"],
  y.train = dUL$y,
  ndpost = 6L,
  nskip = 4L,
  ntree = 5L,
  nchain = 1L,
  verbose = FALSE,
  keeptrees = TRUE
))
expect_identical(messageOf(predict(fitBT, oneNew["f"])), oneIndicators)
expect_identical(messageOf(predict(fitBT, manyNew["f"])), manyIndicators)
# a level the training factor declares and no row takes is refused too
expect_identical(
  messageOf(btFitUL(declaredUL, oneNew)),
  paste0(
    oneIndicators,
    " (a declared level that no training row takes has no indicator ",
    "column)"
  )
)

# --- a sparse factor column and a container, coded over the training levels ---
trainSparse <- data.frame(f = sparseFactor(dUL$f, reference = "L1"), z = dUL$z)
fitSparse <- do.call(
  bart,
  c(list(trainSparse, dUL$y), fitUL, list(keepTrees = TRUE))
)
unseenSparse <- function(count) {
  labels <- c("L1", "L2", fewUL(count))
  data.frame(
    f = sparseFactor(
      factor(labels, levels = c("L1", "L2", fewUL(count))),
      reference = "L1"
    ),
    z = 0
  )
}
expect_identical(
  messageOf(predict(fitSparse, unseenSparse(1L))),
  oneWhole
)
expect_identical(
  messageOf(predict(fitSparse, unseenSparse(8L))),
  paste0(
    "test data factor 'f' has levels not present in the training data: ",
    "'new1', 'new2', 'new3', 'new4', 'new5' and 3 more. Recode them to ",
    "training levels, drop the rows that take them, or ",
    sprintf(remedyWhole, "them")
  )
)
# an unseen reference that a row takes
takenReference <- data.frame(
  f = sparseFactor(
    "L2",
    levels = c("new1", "L1", "L2"),
    reference = "new1",
    i = 1L,
    length = 2L
  ),
  z = c(0, 0)
)
expect_identical(messageOf(predict(fitSparse, takenReference)), oneWhole)

# a container built elsewhere carries its own level table
trainMixed <- data.frame(f = dUL$f)
trainMixed$sv <- Matrix::sparseVector(
  x = runif(10L),
  i = sample(nUL, 10L),
  length = nUL
)
fitMixed <- do.call(
  bart,
  c(list(trainMixed, dUL$y), fitUL, list(keepTrees = TRUE))
)
frameMixed <- function(labels) {
  frame <- data.frame(f = factor(labels))
  frame$sv <- Matrix::sparseVector(x = 0.5, i = 2L, length = length(labels))
  dbartsData(frame, rnorm(length(labels)))@x
}
expect_identical(
  messageOf(predict(fitMixed, frameMixed(c("L1", "new1", "L2")))),
  oneWhole
)

# a sparse factor column inside the container: stored labels, and a reference
# that rows take
sparseContainerUL <- function(labels, levels, reference) {
  frame <- data.frame(
    f = sparseFactor(factor(labels, levels = levels), reference = reference)
  )
  frame$sv <- Matrix::sparseVector(x = 0.5, i = 2L, length = length(labels))
  dbartsData(frame, rnorm(length(labels)))@x
}
expect_identical(
  messageOf(predict(
    fitMixed,
    sparseContainerUL(
      c("L1", "new1", "L2", "new2"),
      c("L1", "L2", "new1", "new2"),
      "L1"
    )
  )),
  paste0(
    "test data factor 'f' has levels not present in the training data: ",
    "'new1' and 'new2'. Recode them to training levels, drop the rows that ",
    "take them, or ",
    sprintf(remedyWhole, "them")
  )
)
expect_identical(
  messageOf(predict(
    fitMixed,
    sparseContainerUL(c("L1", "new1", "L2"), c("new1", "L1", "L2"), "new1")
  )),
  oneWhole
)

# --- lists of labels elsewhere are bounded the same way ---
tooMany <- dbarts:::quotedLevelList(fewUL(7L))
expect_identical(
  tooMany,
  "'new1', 'new2', 'new3', 'new4', 'new5' and 2 more"
)
expect_identical(
  dbarts:::quotedLevelList(c("a", "b", "a")),
  "'a' and 'b'"
)
expect_error(
  pdbart(
    fitWhole,
    xind = "f",
    levs = list(f = fewUL(7L)),
    pl = FALSE
  ),
  pattern = "names levels not present in training: 'new1', 'new2', 'new3', 'new4', 'new5' and 2 more"
)

# a sampler's setPredictor and setData name the labels it never had, bounded
labelsUL <- paste0("x", 1:8)
expect_identical(
  messageOf(plainSampler$setPredictor(rep_len(labelsUL, nUL), 1L)),
  paste0(
    "column 'f' has labels 'x1', 'x2', 'x3', 'x4', 'x5' and 3 more not ",
    "among its training levels"
  )
)
swapUL <- data.frame(
  f = factor(rep_len(c(seenUL, labelsUL), nUL)),
  z = dUL$z
)
swapUL$y <- dUL$y
expect_identical(
  messageOf(plainSampler$setData(dbartsData(y ~ f + z, swapUL))),
  paste0(
    "$setData: column 'f' has labels 'x1', 'x2', 'x3', 'x4', 'x5' and 3 more ",
    "not among its training levels"
  )
)

# a hand-built indicator design, whose test factor has fewer levels than the
# training factor's, lists the test levels the same way
lengthsUL <- paste0("m", 1:8)
designUL <- data.frame(
  a = rnorm(nUL),
  f = factor(sample(lengthsUL, nUL, TRUE), levels = lengthsUL)
)
fitDesign <- suppressWarnings(bartBT(
  makeModelMatrixFromDataFrame(designUL),
  dUL$y,
  ndpost = 6L,
  nskip = 4L,
  ntree = 5L,
  verbose = FALSE,
  keeptrees = TRUE
))
expect_identical(
  messageOf(predict(
    fitDesign,
    data.frame(a = rnorm(7L), f = factor(lengthsUL[1:7]))
  )),
  paste0(
    "'test' factor 'f' does not match training's indicator columns ('test' ",
    "levels: 'm1', 'm2', 'm3', 'm4', 'm5' and 2 more); give the test factor ",
    "the training levels, factor(x, levels = ...), or fit with factors = ",
    "\"categorical\""
  )
)
