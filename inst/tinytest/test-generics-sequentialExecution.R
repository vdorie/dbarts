source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that fixed sample mode when run sequentially gives same predictions as sequential updates mode
set.seed(0L)
pred.bart <- dbarts::bart2(
  testData$x, testData$y, testData$x,
  n.samples = 5, n.burn = 0L, n.trees = 4L,
  k = 2, n.chains = 1L, n.threads = 1L,
  keepTrees = TRUE, verbose = FALSE
)$yhat.test

set.seed(0L)
sampler <- dbarts::dbarts(
  testData$x, testData$y,
  control = dbartsControl(
    n.samples = 5, n.burn = 0L, n.trees = 4L,
    n.chains = 1L, n.threads = 1L,
    keepTrees = TRUE, updateState = FALSE
  )
)
sampler$sampleTreesFromPrior()
for (i in seq_len(5L))
  invisible(sampler$run(0L, 1L))
pred.dbarts <- sampler$predict(testData$x)

expect_equal(pred.bart, t(pred.dbarts))

rm(pred.dbarts, i, sampler, pred.bart)

# test that sequentially running samples don't overflow with fixed trees
sampler <- dbarts::dbarts(
  testData$x, testData$y,
  control = dbartsControl(
    n.samples = 5, n.burn = 0L, n.trees = 4L, n.chains = 1L, n.threads = 1L,
    keepTrees = TRUE, updateState = FALSE
  )
)
sampler$sampleTreesFromPrior()
for (i in seq_len(6L))
  invisible(sampler$run(0L, 1L))

expect_inherits(sampler, "dbartsSampler")

rm(i, sampler)

rm(testData)

