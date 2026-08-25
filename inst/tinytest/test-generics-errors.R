source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that predict fails if sampler not saved
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 20,
  nskip = 5,
  ntree = 5L,
  verbose = FALSE
)
# the pattern matters: bartFit is built with a namespace-qualified call, whose
# call[[1L]] is a `::` call rather than a bare symbol. An unanchored
# expect_error passed even while the guard died with "the condition has length
# > 1" instead of naming 'keeptrees'.
expect_error(
  predict(bartFit, testData$x, n.threads = 1L),
  pattern = "keeptrees"
)
expect_error(extract(bartFit, type = "trees"), pattern = "keeptrees")

bart2Fit <- dbarts::bart2(
  testData$x,
  testData$y,
  n.samples = 20L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = FALSE,
  keepTrainingFits = FALSE
)
expect_error(
  predict(bart2Fit, testData$x, n.threads = 1L),
  pattern = "keepTrees"
)
expect_error(extract(bart2Fit, type = "trees"), pattern = "keepTrees")
expect_error(extract(bart2Fit, sample = "train"), pattern = "keepTrainingFits")

rm(bartFit, bart2Fit)

# predict.bart/predict.rbart did not call refuseUnusedGenericArgs: a caller
# typing the sibling family's offset formal name ('offset.test', used by
# predict.bartNegbin) instead of this fit's own 'offset' had it silently
# vanish into '...' instead of applied
bart2FitKT <- dbarts::bart2(
  testData$x,
  testData$y,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_error(
  predict(bart2FitKT, testData$x, offset.test = rep(0.1, nrow(testData$x))),
  "'offset.test' is not used by predict on a bart fit",
  fixed = TRUE
)
rm(bart2FitKT)

groupBy <- rep(1:2, length.out = nrow(testData$x))
rbartFitKT <- dbarts::rbart_vi(
  testData$y ~ testData$x,
  group.by = groupBy,
  n.samples = 10L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_error(
  predict(
    rbartFitKT,
    testData$x,
    group.by = groupBy,
    offset.test = rep(0.1, nrow(testData$x))
  ),
  "'offset.test' is not used by predict on a rbart fit",
  fixed = TRUE
)
rm(rbartFitKT, groupBy)

rm(testData)
