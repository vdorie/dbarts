source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that fixed sample mode when run sequentially gives same predictions as sequential updates mode
set.seed(0L)
pred.bart <- dbarts::bart2(
  testData$x,
  testData$y,
  testData$x,
  n.samples = 5,
  n.burn = 0L,
  n.trees = 4L,
  k = 2,
  n.chains = 1L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)$yhat.test

set.seed(0L)
sampler <- dbarts::dbarts(
  testData$x,
  testData$y,
  control = dbartsControl(
    n.samples = 5,
    n.burn = 0L,
    n.trees = 4L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    updateState = FALSE
  )
)
# the prior draw must reach the engine: the live trees are stumps before the
# call and carry splits after it (the yhat comparison below cannot see this -
# both sides draw their trees the same way)
treesBefore <- sampler$getTrees(current = TRUE)
expect_equal(nrow(treesBefore), 4L)
expect_equal(unique(treesBefore$var), -1L)
sampler$sampleTreesFromPrior()
treesAfter <- sampler$getTrees(current = TRUE)
expect_true(nrow(treesAfter) > nrow(treesBefore))
expect_true(any(treesAfter$var > 0L))

for (i in seq_len(5L)) {
  invisible(sampler$run(0L, 1L))
}
pred.dbarts <- sampler$predict(testData$x)

expect_equal(pred.bart, t(pred.dbarts))

rm(pred.dbarts, i, sampler, pred.bart, treesBefore, treesAfter)

# test that sequentially running samples don't overflow with fixed trees
sampler <- dbarts::dbarts(
  testData$x,
  testData$y,
  control = dbartsControl(
    n.samples = 5,
    n.burn = 0L,
    n.trees = 4L,
    n.chains = 1L,
    n.threads = 1L,
    keepTrees = TRUE,
    updateState = FALSE
  )
)
sampler$sampleTreesFromPrior()
for (i in seq_len(6L)) {
  invisible(sampler$run(0L, 1L))
}

expect_inherits(sampler, "dbartsSampler")

rm(i, sampler)

rm(testData)
