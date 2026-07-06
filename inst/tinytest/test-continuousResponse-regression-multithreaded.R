source(
  system.file("common", "multithreadData.R", package = "dbarts"),
  local = TRUE
)

n.sims <- 5L
n.burn <- 0L
n.tree <- 3L
weights <- c(
  rep(1, floor(.9 * nrow(testData$x))),
  rep(2, nrow(testData$x) - floor(.9 * nrow(testData$x)))
)

set.seed(99L)
bartFit <- dbarts::bart(
  testData$x,
  testData$y,
  weights = weights,
  ndpost = n.sims,
  nskip = n.burn,
  ntree = n.tree,
  verbose = FALSE,
  nthread = 2L
)

expect_equal(bartFit$sigest, 1.08105693951868)
expect_equal(
  bartFit$sigma,
  c(
    1.23724513909869,
    1.14277060463474,
    1.14262103884284,
    1.12976218154258,
    1.12824224176549
  )
)
expect_equal(
  bartFit$yhat.train[n.sims, 1L:5L],
  c(
    0.165454366580265,
    0.165454366580265,
    0.165454366580265,
    0.165454366580265,
    0.165454366580265
  )
)

rm(testData)
