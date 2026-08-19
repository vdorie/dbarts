source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that bart2 yields similar results to bart
n.burn <- 200L
n.sims <- 400L
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = n.sims,
  nskip = n.burn,
  ntree = 50L,
  nchain = 4L,
  nthread = 1L,
  verbose = FALSE
)

bart2Fit <- dbarts::bart2(
  testData$x,
  testData$y,
  n.samples = n.sims,
  n.burn = n.burn,
  n.trees = 50L,
  n.chains = 4L,
  n.threads = 1L,
  keepTrees = FALSE,
  verbose = FALSE
)

expect_inherits(bartFit, "bart")
expect_inherits(bart2Fit, "bart")
expect_true(
  sqrt(mean((bartFit$yhat.train.mean - bart2Fit$yhat.train.mean)^2)) /
    sd(testData$y) <
    0.1
)

rm(bart2Fit, bartFit, n.sims, n.burn)

rm(testData)

# bart() gains subset, storage, and family, appended and forwarded to
# dbarts(). Each is exercised tiny.

set.seed(202)
nS10 <- 30L
xS10 <- matrix(rnorm(nS10 * 2L), nS10, dimnames = list(NULL, c("a", "b")))
yS10 <- rnorm(nS10)
yBinS10 <- rbinom(nS10, 1L, 0.5)
quickS10 <- list(
  ndpost = 5L,
  nskip = 2L,
  ntree = 3L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
)

# subset actually subsets, and reaches dbarts() the same way pre-subsetting
# x/y directly would: both build the same @x/@y, so the draws are identical
idxS10 <- c(1:10, 15:20)
subsetFit <- do.call(
  dbarts::bart,
  c(list(xS10, yS10, subset = idxS10, seed = 11L), quickS10)
)
directFit <- do.call(
  dbarts::bart,
  c(list(xS10[idxS10, ], yS10[idxS10], seed = 11L), quickS10)
)
expect_equal(length(subsetFit$yhat.train.mean), length(idxS10))
expect_identical(subsetFit$yhat.train, directFit$yhat.train)
expect_identical(subsetFit$sigma, directFit$sigma)

# storage reaches the control and changes the returned draws
doubleFit <- do.call(
  dbarts::bart,
  c(list(xS10, yS10, seed = 22L, keepsampler = TRUE), quickS10)
)
singleFit <- do.call(
  dbarts::bart,
  c(
    list(xS10, yS10, seed = 22L, storage = "single", keepsampler = TRUE),
    quickS10
  )
)
expect_equal(doubleFit$fit$control@storage, "double")
expect_equal(singleFit$fit$control@storage, "single")
expect_false(identical(doubleFit$yhat.train, singleFit$yhat.train))

# family = "logistic"/"aft" reach the right path and package as "bart"
logisticFit <- do.call(
  dbarts::bart,
  c(list(xS10, yBinS10, family = "logistic", seed = 33L), quickS10)
)
expect_inherits(logisticFit, "bart")
expect_equal(logisticFit$family, "logistic")

aftYS10 <- cbind(abs(yS10) + 0.1, rep(1L, nS10))
aftFit <- do.call(
  dbarts::bart,
  c(list(xS10, aftYS10, family = "aft", seed = 44L), quickS10)
)
expect_inherits(aftFit, "bart")
expect_equal(aftFit$family, "aft")

# the four own-class families are refused by name, pointing at bart2 - not
# match.arg's generic "should be one of" message
expect_error(
  dbarts::bart(xS10, yS10, family = "multinomial"),
  pattern = "bart2"
)
expect_error(dbarts::bart(xS10, yS10, family = "ordinal"), pattern = "bart2")
expect_error(dbarts::bart(xS10, yS10, family = "nbinom"), pattern = "bart2")
expect_error(
  dbarts::bart(xS10, yS10, family = "hurdle.lognormal"),
  pattern = "bart2"
)
expect_error(
  dbarts::bart(xS10, yS10, family = "nbinom"),
  pattern = "\"nbinom\""
)

# the appended formals break no existing bart() abbreviation
abbrevFit1 <- dbarts::bart(
  xS10,
  yS10,
  ntre = 3L,
  ndpost = 5L,
  nskip = 2L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  seed = 55L
)
expect_inherits(abbrevFit1, "bart")
abbrevFit2 <- dbarts::bart(
  xS10,
  yS10,
  ntree = 3L,
  ndpo = 5L,
  nskip = 2L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  seed = 55L
)
expect_inherits(abbrevFit2, "bart")

rm(
  nS10,
  xS10,
  yS10,
  yBinS10,
  quickS10,
  idxS10,
  subsetFit,
  directFit,
  doubleFit,
  singleFit,
  logisticFit,
  aftYS10,
  aftFit,
  abbrevFit1,
  abbrevFit2
)

# keeptrees = TRUE, keepsampler = FALSE keeps $fit anyway (keepsampler's
# default IS keeptrees; an explicit FALSE override must not lose it) -
# predict and extract("trees") both work
set.seed(303)
nKT <- 40L
xKT <- matrix(rnorm(nKT * 2L), nKT, 2L)
yKT <- xKT[, 1L] + rnorm(nKT)
fitKT <- dbarts::bart(
  xKT,
  yKT,
  ndpost = 5L,
  nskip = 4L,
  ntree = 4L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE,
  seed = 66L,
  keeptrees = TRUE,
  keepsampler = FALSE
)
expect_false(is.null(fitKT$fit))
expect_silent(predict(fitKT, xKT[1:3, ]))
expect_silent(extract(fitKT, "trees"))

rm(nKT, xKT, yKT, fitKT)

# bart()'s x/y route (factors = "indicators", no stored level table) names
# the mismatching factor and its levels rather than blaming the column
# count when a test factor's levels differ from training's
set.seed(404)
nF <- 60L
dF <- data.frame(
  x1 = rnorm(nF),
  f = factor(sample(c("a", "b", "c"), nF, TRUE))
)
dF$y <- dF$x1 *
  2 +
  ifelse(dF$f == "a", 3, ifelse(dF$f == "b", -3, 0)) +
  rnorm(nF, 0, 0.2)
trF <- dF[1:40, c("x1", "f")]
ytrF <- dF$y[1:40]
teFewer <- dF[41:60, c("x1", "f")]
teFewer <- teFewer[teFewer$f != "c", ]
teFewer$f <- droplevels(teFewer$f)
expect_error(
  dbarts::bart(
    trF,
    ytrF,
    teFewer,
    ndpost = 5L,
    nskip = 4L,
    ntree = 5L,
    nchain = 1L,
    verbose = FALSE,
    seed = 4L
  ),
  pattern = "'test' factor 'f' does not match training's indicator columns"
)
# a test factor declaring MORE levels than training's is named by the level
# tables themselves, before the drop-pattern replay indexes past its end;
# the count of undeclared levels is unbounded and does not change the refusal
teExtra <- dF[41:60, c("x1", "f")]
teExtra$f <- factor(teExtra$f, levels = c("a", "b", "c", "d"))
teMany <- dF[41:60, c("x1", "f")]
teMany$f <- factor(teMany$f, levels = c("a", "b", "c", paste0("z", 1:5000)))
expect_error(
  dbarts::bart(
    trF,
    ytrF,
    teExtra,
    ndpost = 5L,
    nskip = 4L,
    ntree = 5L,
    nchain = 1L,
    verbose = FALSE,
    seed = 4L
  ),
  pattern = "'test' factor 'f' declares 4 levels but the training design declared 3"
)
expect_error(
  dbarts::bart(
    trF,
    ytrF,
    teMany,
    ndpost = 5L,
    nskip = 4L,
    ntree = 5L,
    nchain = 1L,
    verbose = FALSE,
    seed = 4L
  ),
  pattern = "'test' factor 'f' declares 5003 levels"
)
# predict() reaches the same funnel with the fit's stored drop pattern
fitF <- dbarts::bart(
  trF,
  ytrF,
  ndpost = 5L,
  nskip = 4L,
  ntree = 5L,
  nchain = 1L,
  verbose = FALSE,
  seed = 4L,
  keeptrees = TRUE
)
expect_error(
  predict(fitF, teExtra),
  pattern = "'test' factor 'f' declares 4 levels but the training design declared 3"
)
expect_error(
  predict(fitF, teMany),
  pattern = "'test' factor 'f' declares 5003 levels"
)
expect_error(
  predict(fitF, teFewer),
  pattern = "'test' factor 'f' does not match training's indicator columns"
)
# a character test column is expanded as a factor, so its distinct values
# are compared to the training table the same way
teChar <- dF[41:60, c("x1", "f")]
teChar$f <- as.character(teChar$f)
teChar$f[1L] <- "d"
expect_error(
  predict(fitF, teChar),
  pattern = "'test' factor 'f' declares 4 levels but the training design declared 3"
)
rm(nF, dF, trF, ytrF, teFewer, teExtra, teMany, teChar, fitF)

# bart()'s missing-predictor refusal cannot advise an argument bart()
# itself rejects (missing = "incorporate", bart2()/dbarts() only); it
# points at those front doors instead
set.seed(505)
nMiss <- 40L
xMiss <- matrix(rnorm(nMiss * 2L), nMiss, 2L)
xMiss[1L, 1L] <- NA_real_
yMiss <- rnorm(nMiss)
expect_error(
  dbarts::bart(
    xMiss,
    yMiss,
    ndpost = 3L,
    nskip = 2L,
    ntree = 3L,
    verbose = FALSE
  ),
  pattern = "use bart2\\(\\) or dbarts\\(\\)"
)
expect_error(
  dbarts::bart(
    xMiss,
    yMiss,
    ndpost = 3L,
    nskip = 2L,
    ntree = 3L,
    verbose = FALSE,
    missing = "incorporate"
  ),
  pattern = "unused argument"
)
rm(nMiss, xMiss, yMiss)
