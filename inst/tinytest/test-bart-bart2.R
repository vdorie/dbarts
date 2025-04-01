source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that bart2 yields similar results to bart
n.burn <- 200L
n.sims <- 400L
bartFit <- dbarts::bart(
  testData$x, testData$y,
  ndpost = n.sims, nskip = n.burn, ntree = 50L, nchain = 4L, nthread = 1L,
  verbose = FALSE
)

bart2Fit <- dbarts::bart2(
  testData$x, testData$y,
  n.samples = n.sims, n.burn = n.burn, n.trees = 50L, n.chains = 4L, n.threads = 1L,
  keepTrees = FALSE, verbose = FALSE
)

expect_inherits(bartFit, "bart")
expect_inherits(bart2Fit, "bart")
expect_true(
  sqrt(mean((bartFit$yhat.train.mean - bart2Fit$yhat.train.mean)^2)) / sd(testData$y) < 0.1
)

rm(bart2Fit, bartFit, n.sims, n.burn)

rm(testData)

