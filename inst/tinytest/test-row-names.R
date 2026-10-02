# Row names on every observation-indexed output (dec-B34): the stored draws,
# extract, fitted, residuals, predict and survivalProbabilities, on every fit
# class and both channels, across the entry paths that capture the names.

set.seed(21)
n <- 24L
trainNames <- paste0("r", seq_len(n))
x <- matrix(
  rnorm(n * 2L),
  n,
  2L,
  dimnames = list(trainNames, c("a", "b"))
)
x.test <- x[1:3, ]
testNames <- c("t1", "t2", "t3")
rownames(x.test) <- testNames
y <- x[, 1L] + rnorm(n, 0, 0.3)

# do.call, so that 'subset' arrives evaluated
quick <- function(...) {
  suppressWarnings(do.call(
    bart,
    list(
      ...,
      n.samples = 8L,
      n.burn = 4L,
      n.trees = 5L,
      n.chains = 2L,
      keepTrees = TRUE,
      verbose = FALSE
    )
  ))
}
lastNames <- function(x) {
  if (is.null(dim(x))) names(x) else dimnames(x)[[length(dim(x))]]
}

# --- class bart, matrix path: both channels and every accessor ---

fit <- quick(x, y, test = x.test)
expect_identical(fit$row.names.train, trainNames)
expect_identical(fit$row.names.test, testNames)
expect_identical(lastNames(fit$yhat.train), trainNames)
expect_identical(lastNames(fit$yhat.test), testNames)
expect_identical(names(fit$yhat.train.mean), trainNames)
expect_identical(names(fit$yhat.test.mean), testNames)
for (type in c("ev", "ppd", "bart", "loglik")) {
  expect_identical(lastNames(extract(fit, type)), trainNames, info = type)
  expect_identical(
    lastNames(extract(fit, type, combineChains = FALSE)),
    trainNames,
    info = type
  )
}
expect_identical(lastNames(extract(fit, "ev", "test")), testNames)
expect_identical(names(fitted(fit)), trainNames)
expect_identical(names(fitted(fit, sample = "test")), testNames)
expect_identical(rownames(fitted(fit, ci.level = 0.9)), trainNames)
expect_identical(names(residuals(fit)), trainNames)
for (type in c("ev", "ppd", "bart")) {
  expect_identical(lastNames(predict(fit, x.test, type)), testNames)
}
expect_identical(rownames(predict(fit, x.test, ci.level = 0.9)), testNames)
# names are the only change: the values match an unnamed comparison
expect_equal(
  unname(fitted(fit)),
  as.vector(colMeans(unname(extract(fit, "ev"))))
)

# a matrix without row names stays unnamed, as lm.fit does
unnamedX <- unname(x)
fit0 <- quick(unnamedX, y, test = unname(x.test))
expect_null(fit0$row.names.train)
expect_null(fit0$row.names.test)
expect_false("row.names.train" %in% names(fit0))
expect_null(names(fitted(fit0)))
expect_null(dimnames(extract(fit0, "ev")))
expect_null(lastNames(predict(fit0, unname(x.test))))

# --- class bart, formula and data frame paths ---

df <- data.frame(y = y, a = x[, 1L], b = x[, 2L], row.names = trainNames)
fitF <- quick(y ~ a + b, df, test = data.frame(x.test))
expect_identical(names(fitted(fitF)), trainNames)
expect_identical(names(fitted(fitF, sample = "test")), testNames)
expect_identical(
  lastNames(predict(fitF, data.frame(x.test)[2:3, ])),
  testNames[2:3]
)

# automatic data frame names count as "1".."n", as lm's
dfAuto <- df
rownames(dfAuto) <- NULL
fitAuto <- quick(y ~ a + b, dfAuto)
expect_identical(names(fitted(fitAuto)), as.character(seq_len(n)))
fitAutoX <- quick(dfAuto[, c("a", "b")], y)
expect_identical(names(fitted(fitAutoX)), as.character(seq_len(n)))

# 'subset' keeps the selected rows' names
fitSub <- quick(y ~ a + b, df, subset = 5:9)
expect_identical(names(fitted(fitSub)), trainNames[5:9])
fitSubX <- quick(x, y, subset = 5:9)
expect_identical(names(fitted(fitSubX)), trainNames[5:9])

# a factor design is a mixed container, whose names live outside it
dfFactor <- data.frame(
  a = x[, 1L],
  g = factor(rep(c("u", "v", "w"), length.out = n)),
  row.names = trainNames
)
fitFactor <- quick(dfFactor, y)
expect_identical(names(fitted(fitFactor)), trainNames)
expect_identical(
  lastNames(predict(fitFactor, dfFactor[4:5, ])),
  trainNames[4:5]
)

# a dgCMatrix keeps its row names on both paths
if (requireNamespace("Matrix", quietly = TRUE)) {
  xSparse <- Matrix::Matrix(x, sparse = TRUE)
  fitSparse <- quick(xSparse, y)
  expect_identical(names(fitted(fitSparse)), trainNames)
  expect_identical(
    lastNames(predict(fitSparse, xSparse[2:3, ])),
    trainNames[2:3]
  )
}

# na.exclude pads fitted and residuals back with the dropped rows' names, on
# the formula and matrix paths alike
yMissing <- y
yMissing[3L] <- NA
dfMissing <- df
dfMissing$y <- yMissing
for (fitE in list(
  quick(y ~ a + b, dfMissing, na.action = na.exclude),
  quick(x, yMissing, na.action = na.exclude)
)) {
  expect_identical(names(fitE$na.action), "r3")
  expect_identical(names(fitted(fitE)), trainNames)
  expect_true(is.na(fitted(fitE)[["r3"]]))
  expect_identical(names(residuals(fitE)), trainNames)
  expect_identical(lastNames(fitE$yhat.train), trainNames[-3L])
}

# an integer test matrix keeps its dimnames: its columns match by name
integerTest <- matrix(
  1:6,
  3L,
  2L,
  dimnames = list(testNames, c("b", "a"))
)
fitInt <- quick(x, y, test = integerTest)
expect_identical(lastNames(fitInt$yhat.test), testNames)
expect_identical(
  dbarts:::validateXTest(integerTest, x),
  matrix(
    c(4, 5, 6, 1, 2, 3),
    3L,
    2L,
    dimnames = list(testNames, c("a", "b"))
  )
)

# --- bartBT, which forces na.omit ---

fitBT <- suppressWarnings(bartBT(
  x,
  y,
  x.test,
  ndpost = 8L,
  nskip = 4L,
  ntree = 5L,
  verbose = FALSE
))
expect_identical(lastNames(fitBT$yhat.train), trainNames)
expect_identical(lastNames(fitBT$yhat.test), testNames)
expect_identical(names(fitted(fitBT)), trainNames)

# --- the other classes ---

category <- factor(rep(c("u", "v", "w"), length.out = n))
fitM <- quick(x, category, test = x.test, family = "multinomial")
expect_identical(dimnames(fitM$yhat.train)[[2L]], trainNames)
expect_identical(dimnames(fitM$yhat.test)[[2L]], testNames)
expect_identical(rownames(fitted(fitM)), trainNames)
expect_identical(names(fitted(fitM, "class")), trainNames)
expect_identical(rownames(residuals(fitM)), trainNames)
expect_identical(lastNames(extract(fitM, "loglik")), trainNames)
expect_identical(lastNames(extract(fitM, "ppd")), trainNames)
expect_identical(dimnames(predict(fitM, x.test))[[2L]], testNames)
expect_identical(names(predict(fitM, x.test, "class")), testNames)
expect_identical(
  dimnames(predict(fitM, x.test, ci.level = 0.9))[[1L]],
  testNames
)

fitO <- quick(
  x,
  factor(category, ordered = TRUE),
  test = x.test,
  family = "ordinal"
)
expect_identical(dimnames(fitO$yhat.train)[[2L]], trainNames)
expect_identical(lastNames(extract(fitO, "bart")), trainNames)
expect_identical(lastNames(extract(fitO, "loglik")), trainNames)
expect_identical(lastNames(extract(fitO, "bart", "test")), testNames)
expect_identical(names(fitted(fitO, "bart")), trainNames)
expect_identical(rownames(residuals(fitO)), trainNames)
expect_identical(lastNames(predict(fitO, x.test, "bart")), testNames)
expect_identical(names(predict(fitO, x.test, "class")), testNames)

counts <- rpois(n, 3)
fitN <- quick(x, counts, test = x.test, family = "nbinom")
expect_identical(names(fitted(fitN)), trainNames)
expect_identical(names(residuals(fitN)), trainNames)
for (type in c("ev", "ppd", "bart", "loglik")) {
  expect_identical(lastNames(extract(fitN, type)), trainNames, info = type)
}
for (type in c("ev", "ppd", "bart")) {
  expect_identical(lastNames(predict(fitN, x.test, type)), testNames)
}

positive <- ifelse(seq_len(n) %% 3L == 0L, 0, exp(y))
fitH <- quick(x, positive, family = "hurdle.lognormal")
expect_identical(fitH$row.names.train, trainNames)
expect_identical(names(fitted(fitH)), trainNames)
expect_identical(names(residuals(fitH)), trainNames)
expect_identical(lastNames(extract(fitH, "loglik")), trainNames)
expect_identical(lastNames(predict(fitH, x.test, "ppd")), testNames)

# --- fit-time na.exclude/na.omit parity: multinomial, ordinal and negbin ---
# (dec-B34): the four packagers store the fit's na.action, as bart does, and
# fitted/residuals pad through it the same way. A na.omit fit here also
# regression-tests a defect: an NA in x used to fail with an x/y length
# mismatch, since the derived count response was cut by 'subset' alone, not
# by the na.action's own further drop.

xMissing <- x
xMissing[3L, "a"] <- NA

fitMExclude <- quick(
  xMissing,
  category,
  family = "multinomial",
  na.action = na.exclude
)
fitMOmit <- quick(
  xMissing,
  category,
  family = "multinomial",
  na.action = na.omit
)
expect_identical(names(fitMExclude$na.action), "r3")
expect_identical(rownames(fitted(fitMExclude)), trainNames)
expect_true(all(is.na(fitted(fitMExclude)["r3", ])))
expect_identical(rownames(residuals(fitMExclude)), trainNames)
expect_identical(rownames(fitted(fitMOmit)), trainNames[-3L])
expect_identical(rownames(residuals(fitMOmit)), trainNames[-3L])
# a category offset loses the same rows
fitMOffset <- quick(
  xMissing,
  category,
  family = "multinomial",
  offset = matrix(0.1, n, 3L),
  na.action = na.exclude
)
expect_identical(rownames(fitted(fitMOffset)), trainNames)

fitOExclude <- quick(
  xMissing,
  factor(category, ordered = TRUE),
  family = "ordinal",
  na.action = na.exclude
)
fitOOmit <- quick(
  xMissing,
  factor(category, ordered = TRUE),
  family = "ordinal",
  na.action = na.omit
)
expect_identical(names(fitted(fitOExclude, "bart")), trainNames)
expect_true(is.na(fitted(fitOExclude, "bart")[["r3"]]))
expect_identical(rownames(residuals(fitOExclude)), trainNames)
expect_identical(names(fitted(fitOOmit, "bart")), trainNames[-3L])
expect_identical(rownames(residuals(fitOOmit)), trainNames[-3L])

fitNExclude <- quick(
  xMissing,
  counts,
  family = "nbinom",
  na.action = na.exclude
)
fitNOmit <- quick(xMissing, counts, family = "nbinom", na.action = na.omit)
expect_identical(names(fitted(fitNExclude)), trainNames)
expect_true(is.na(fitted(fitNExclude)[["r3"]]))
expect_identical(names(residuals(fitNExclude)), trainNames)
expect_identical(names(fitted(fitNOmit)), trainNames[-3L])
expect_identical(names(residuals(fitNOmit)), trainNames[-3L])

# hurdle stores the zero component's na.action and pads the same way.
# A hurdle fit that genuinely drops a row this way hits an unrelated,
# pre-existing routability refusal (the positive component's own 'test' is
# always the full, un-reduced design matrix, so a row na.action removes from
# training is still present, and now unroutable, in that 'test'), so this
# checks the padding machinery on a fit already trained at the reduced row
# count, its na.action attached by hand exactly as bart2Hurdle would have
# set it from a working zero component.
fitHDropped <- quick(x[-3L, ], positive[-3L], family = "hurdle.lognormal")
fitHExclude <- fitHDropped
fitHExclude$na.action <- structure(3L, class = "exclude", names = "r3")
expect_identical(names(fitted(fitHExclude)), trainNames)
expect_true(is.na(fitted(fitHExclude)[["r3"]]))
expect_identical(names(residuals(fitHExclude)), trainNames)
fitHOmit <- fitHDropped
fitHOmit$na.action <- structure(3L, class = "omit", names = "r3")
expect_identical(names(fitted(fitHOmit)), trainNames[-3L])

# --- survival: AFT and the discrete-time hazard ---

if (requireNamespace("survival", quietly = TRUE)) {
  time <- rexp(n, 1)
  status <- rep(c(1, 1, 0), length.out = n)
  response <- survival::Surv(time, status)
  fitA <- quick(x, response, family = "aft")
  expect_identical(
    lastNames(survivalProbabilities(fitA, times = c(0.5, 1))),
    trainNames
  )
  expect_identical(
    lastNames(survivalProbabilities(fitA, times = 1, newdata = x.test)),
    testNames
  )

  # defect: an NA in x under na.omit used to crash an aft fit with a status
  # length mismatch, since the status vector was cut by 'subset' alone,
  # outside dbartsData(), never by the na.action's own further drop
  xAftMissing <- x
  xAftMissing[3L, "a"] <- NA
  fitAftOmit <- quick(
    xAftMissing,
    survival::Surv(time, status),
    family = "aft",
    na.action = na.omit
  )
  expect_identical(length(fitAftOmit$y), n - 1L)
  expect_identical(names(fitted(fitAftOmit)), trainNames[-3L])

  # the person-period rows are named by make.unique over the subjects, on the
  # training and test channels alike; a subject named "r1.1" takes its own
  # name and pushes the second period of "r1" past it
  hazardX <- x[1:6, ]
  rownames(hazardX)[2L] <- "r1.1"
  hazardTime <- c(2.5, 0.5, 1.5, 2.5, 0.5, 1.5)
  hazardStatus <- c(1, 1, 0, 1, 1, 0)
  fitZ <- quick(
    hazardX,
    survival::Surv(hazardTime, hazardStatus),
    test = x.test[1:2, ],
    family = quote(hazard(breaks = c(0, 1, 2, 3)))
  )
  subjectNames <- rownames(hazardX)
  periods <- c(3L, 1L, 2L, 3L, 1L, 2L)
  expect_identical(
    fitZ$row.names.train,
    make.unique(rep(subjectNames, periods))
  )
  expect_identical(
    fitZ$row.names.train[1:4],
    c("r1", "r1.2", "r1.3", "r1.1")
  )
  expect_identical(
    fitZ$row.names.test,
    c("t1", "t2", "t1.1", "t2.1", "t1.2", "t2.2")
  )
  expect_identical(lastNames(fitZ$yhat.train), fitZ$row.names.train)
  expect_identical(lastNames(survivalProbabilities(fitZ)), c("t1", "t2"))
  expect_identical(
    lastNames(survivalProbabilities(fitZ, newdata = x.test[2:3, ])),
    c("t2", "t3")
  )
  fitZ0 <- quick(
    hazardX,
    survival::Surv(hazardTime, hazardStatus),
    family = quote(hazard(breaks = c(0, 1, 2, 3)))
  )
  expect_identical(lastNames(survivalProbabilities(fitZ0)), subjectNames)

  hazardFrame <- data.frame(
    time = hazardTime,
    status = hazardStatus,
    a = hazardX[, 1L],
    b = hazardX[, 2L],
    row.names = subjectNames
  )
  fitZF <- quick(
    survival::Surv(time, status) ~ a + b,
    hazardFrame,
    family = quote(hazard(breaks = c(0, 1, 2, 3)))
  )
  expect_identical(fitZF$row.names.train, fitZ0$row.names.train)
  expect_identical(lastNames(survivalProbabilities(fitZF)), subjectNames)

  # the matrix interface expands before the na.action runs, so an incomplete
  # subject drops as person-period rows, which the record names
  hazardXMissing <- hazardX
  hazardXMissing[4L, 2L] <- NA
  fitZE <- quick(
    hazardXMissing,
    survival::Surv(hazardTime, hazardStatus),
    na.action = na.exclude,
    family = quote(hazard(breaks = c(0, 1, 2, 3)))
  )
  allNames <- fitZ0$row.names.train
  dropped <- c("r4", "r4.1", "r4.2")
  expect_identical(names(fitZE$na.action), dropped)
  expect_identical(fitZE$row.names.train, setdiff(allNames, dropped))
  expect_identical(names(fitted(fitZE)), allNames)
  expect_identical(names(which(is.na(fitted(fitZE)))), dropped)

  # the same record without row names still pads
  hazardXUnnamed <- unname(hazardXMissing)
  fitZU <- quick(
    hazardXUnnamed,
    survival::Surv(hazardTime, hazardStatus),
    na.action = na.exclude,
    family = quote(hazard(breaks = c(0, 1, 2, 3)))
  )
  expect_identical(as.vector(fitZU$na.action), 7:9)
  expect_identical(which(is.na(fitted(fitZU))), 7:9)

  # the formula path's na.action runs BEFORE expansion, at the subject
  # level; its record is restated over the dropped subject's person-period
  # rows, so fitted() pads as it does on the matrix path
  hazardFrameMissing <- hazardFrame
  hazardFrameMissing$a[4L] <- NA
  fitZFE <- quick(
    survival::Surv(time, status) ~ a + b,
    hazardFrameMissing,
    na.action = na.exclude,
    family = quote(hazard(breaks = c(0, 1, 2, 3)))
  )
  expect_identical(names(fitZFE$na.action), dropped)
  expect_identical(fitZFE$row.names.train, setdiff(allNames, dropped))
  expect_identical(names(fitted(fitZFE)), allNames)
  # a subject with no time keeps its first-period row, as the expansion gives
  # every subject, so the record names that row and fitted() pads it
  hazardFrameMissing$time[4L] <- NA
  fitZFN <- quick(
    survival::Surv(time, status) ~ a + b,
    hazardFrameMissing,
    na.action = na.exclude,
    family = quote(hazard(breaks = c(0, 1, 2, 3)))
  )
  expect_identical(names(fitZFN$na.action), dropped[1L])
  expect_identical(fitZFN$row.names.train, setdiff(allNames, dropped))
  expect_identical(names(fitted(fitZFN)), setdiff(allNames, dropped[-1L]))
}

# --- the sampler's test setters keep the record in step with the rows ---

sampler <- dbarts(
  x,
  y,
  test = x.test,
  control = dbartsControl(n.chains = 1L, n.threads = 1L, verbose = FALSE)
)
sampler$setTestPredictor(x.test[, 1L] + 1, column = 1L)
expect_identical(sampler$data@rowNames$test, testNames)
newTest <- x[1:4, ]
rownames(newTest) <- paste0("n", 1:4)
sampler$setTestPredictor(newTest)
expect_identical(sampler$data@rowNames$test, rownames(newTest))
sampler$setTestPredictorAndOffset(x.test, NULL)
expect_identical(sampler$data@rowNames$test, testNames)
sampler$setTestPredictor(NULL)
expect_null(sampler$data@rowNames$test)
expect_identical(sampler$data@rowNames$train, trainNames)

# --- the stored arrays are named in place: extract hands them back uncopied ---

if (capabilities("profmem")) {
  fitG <- quick(x, y, combineChains = TRUE)
  stored <- fitG$yhat.train
  invisible(tracemem(stored))
  copies <- capture.output({
    evDraws <- extract(fitG, "ev")
    bartDraws <- extract(fitG, "bart")
    means <- fitted(fitG)
  })
  untracemem(stored)
  expect_identical(copies, character(0))
  expect_true(identical(evDraws, stored))
  expect_identical(names(means), trainNames)
}

# --- objects from before the names: stripped parts give unnamed output ---

stripped <- fit
stripped$row.names.train <- NULL
stripped$row.names.test <- NULL
stripped$yhat.train <- unname(stripped$yhat.train)
stripped$yhat.train.mean <- unname(stripped$yhat.train.mean)
expect_null(names(fitted(stripped)))
expect_null(dimnames(extract(stripped, "ev")))
expect_null(names(residuals(stripped)))

oldData <- dbartsData(x, y, x.test)
expect_identical(oldData@rowNames, list(train = trainNames, test = testNames))
attr(oldData, "rowNames") <- NULL
expect_false(methods::.hasSlot(oldData, "rowNames"))
expect_null(dbarts:::dataRowNames(oldData, "train"))
oldFit <- suppressWarnings(bart(
  oldData,
  n.samples = 4L,
  n.burn = 2L,
  n.trees = 5L,
  verbose = FALSE
))
expect_false("row.names.train" %in% names(oldFit))
expect_null(names(fitted(oldFit)))
expect_null(names(fitted(oldFit, sample = "test")))
