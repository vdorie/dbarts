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

# S10 (docs/plans/bart2-argument-consolidation.md 4.5, 8.7): bart() gains
# subset, storage, and family, appended and forwarded to dbarts(). Each is
# exercised tiny.

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

# the appended formals break no existing bart() abbreviation (6.8)
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
