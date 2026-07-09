source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that posterior predictive distribution samples use correct sigma
n.samples <- 7L
n.chains <- 2L
n.obs <- length(testData$y)
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  verbose = FALSE,
  ndpost = n.samples,
  nskip = 0L,
  nchain = n.chains,
  ntree = 25L,
  nthread = 1L
)
set.seed(0L)
samples.ppd <- extract(bartFit, type = "ppd")

set.seed(0L)
samples.pm <- extract(bartFit)
for (i in seq_len(n.obs)) {
  expect_equal(
    samples.pm[, i] + rnorm(n.samples * n.chains, 0, bartFit$sigma),
    samples.ppd[, i]
  )
}

set.seed(0L)
samples.ppd <- extract(bartFit, type = "ppd", combineChains = FALSE)

set.seed(0L)
samples.pm <- extract(bartFit, combineChains = FALSE)
for (i in seq_len(n.obs)) {
  expect_equal(
    samples.pm[,, i] +
      matrix(rnorm(n.samples * n.chains, 0, bartFit$sigma), nrow = n.chains),
    samples.ppd[,, i]
  )
}

rm(i, samples.pm, samples.ppd, bartFit, n.obs, n.chains, n.samples)

# extract(type = "ppd", sample = "train") must weight its noise draws by
# the fit's case weights (object$weights, not the "weigths" typo it used to
# read): a row with a large weight gets a tight PPD, a small-weight row a
# wide one, even though every row shares the same posterior sigma draws
set.seed(0L)
n <- 200L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n, 0, 0.3)
w <- rep(c(100, 0.01), each = n / 2L)
weightedFit <- dbarts::bart2(
  x,
  y,
  weights = w,
  n.samples = 200L,
  n.burn = 100L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  combineChains = TRUE,
  verbose = FALSE
)
ppd.weighted <- extract(weightedFit, type = "ppd", sample = "train")
sd.highWeight <- apply(ppd.weighted[, w > 1, drop = FALSE], 2L, sd)
sd.lowWeight <- apply(ppd.weighted[, w < 1, drop = FALSE], 2L, sd)
expect_true(all(sd.highWeight < sd.lowWeight))

rm(n, x, y, w, weightedFit, ppd.weighted, sd.highWeight, sd.lowWeight)

rm(testData)

# rbart_vi's extract/predict/fitted PPD paths must thread case weights the
# same way bart2's do (weights for sample = "train", weights.test for
# "test"): a row with a large weight gets a tight PPD, a small-weight row a
# wide one, even though every row shares the same posterior sigma draws
set.seed(0L)
n <- 200L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n, 0, 0.3)
g <- factor(rep(1:4, length.out = n))
w <- rep(c(100, 0.01), each = n / 2L)
weightedRFit <- dbarts::rbart_vi(
  x,
  y,
  group.by = g,
  weights = w,
  n.samples = 200L,
  n.burn = 100L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
ppd.weighted <- extract(weightedRFit, type = "ppd", sample = "train")
sd.highWeight <- apply(ppd.weighted[, w > 1, drop = FALSE], 2L, sd)
sd.lowWeight <- apply(ppd.weighted[, w < 1, drop = FALSE], 2L, sd)
expect_true(all(sd.highWeight < sd.lowWeight))

# fitted()'s ppd path delegates to extract() for rbart, so it inherits the
# same weighting
ppd.fitted <- fitted(weightedRFit, type = "ppd", sample = "train")
expect_true(is.numeric(ppd.fitted))

# predict.rbart has no train/test sample of its own (it always predicts on
# fresh newdata), so it takes a 'weights' argument directly instead, mirroring
# predict.bart's own convention rather than reading object$weights.test
ppd.predict <- predict(weightedRFit, x, g, type = "ppd", weights = w)
sd.highWeight.pred <- apply(ppd.predict[, w > 1, drop = FALSE], 2L, sd)
sd.lowWeight.pred <- apply(ppd.predict[, w < 1, drop = FALSE], 2L, sd)
expect_true(all(sd.highWeight.pred < sd.lowWeight.pred))

rm(
  n,
  x,
  y,
  g,
  w,
  weightedRFit,
  ppd.weighted,
  sd.highWeight,
  sd.lowWeight,
  ppd.fitted,
  ppd.predict,
  sd.highWeight.pred,
  sd.lowWeight.pred
)
