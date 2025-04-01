source(system.file("common", "multithreadData.R", package = "dbarts"), local = TRUE)

n.sims <- 5L
n.burn <- 0L
n.tree <- 3L
weights <- c(
  rep(1, floor(.9 * nrow(testData$x))),
  rep(2, nrow(testData$x) - floor(.9 * nrow(testData$x)))
)

set.seed(99L)
bartFit <- dbarts::bart(
  testData$x, testData$y, weights = weights,
  ndpost = n.sims, nskip = n.burn, ntree = n.tree, verbose = FALSE,
  nthread = 2L
)

expect_equal(bartFit$sigest, 1.08105693951868)
expect_equal(
  bartFit$sigma,
  c(
    1.25489299085944, 1.23929900925602, 1.1457037595997,
    1.14346729349942, 1.08495804241639
  )
)
expect_equal(
  bartFit$yhat.train[n.sims,1L:5L],
  c(
    0.468714161005541, 0.468714161005541, 0.468714161005541,
    0.468714161005541, -0.166345766678715
  )
)

rm(testData)

