source(
  system.file("common", "multithreadData.R", package = "dbarts"),
  local = TRUE
)

# test that multithreaded within chains matches single threaded
## a slice of the source set: nthread doesn't change a single-chain fit's
## code path, so size doesn't affect what's tested; keep it just weak enough
## to run quickly
testData$x <- testData$x[1:5001L, , drop = FALSE]
testData$y <- testData$y[1:5001L]
n.sims <- 15L
n.burn <- 0L
n.tree <- 3L

set.seed(99L)
singleThreadedFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = n.sims,
  nskip = n.burn,
  ntree = n.tree,
  verbose = FALSE,
  nthread = 1L
)

set.seed(99L)
multiThreadedFit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = n.sims,
  nskip = n.burn,
  ntree = n.tree,
  verbose = FALSE,
  nthread = 2L
)

expect_equal(singleThreadedFit$sigma, multiThreadedFit$sigma)
expect_equal(
  singleThreadedFit$yhat.train[n.sims, 1L:5L],
  multiThreadedFit$yhat.train[n.sims, 1L:5L]
)
expect_equal(
  singleThreadedFit$yhat.train.mean[1L:5L],
  multiThreadedFit$yhat.train.mean[1L:5L]
)
expect_null(multiThreadedFit$yhat.test)
expect_null(multiThreadedFit$yhat.test.mean)
expect_equal(
  singleThreadedFit$varcount[n.sims, ],
  multiThreadedFit$varcount[n.sims, ]
)

rm(testData)

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# test that multiple chains single threads runs correctly
set.seed(99L)
oldSeed <- .Random.seed
fit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 120L,
  nskip = 100L,
  ntree = 50L,
  nthread = 1L,
  nchain = 2L,
  combinechains = FALSE,
  verbose = FALSE
)
expect_true(any(oldSeed != .Random.seed))

rm(fit, oldSeed)

# test that multiple chains multiple threads runs correctly
set.seed(99L)
oldSeed <- .Random.seed
fit <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 120L,
  nskip = 100L,
  ntree = 50L,
  nthread = 2L,
  nchain = 2L,
  combinechains = FALSE,
  verbose = FALSE
)

expect_identical(dim(fit$yhat.train), c(2L, 120L, nrow(testData$x)))
# check that the random seeds didn't end up the same/the chains got different
# results
expect_true(
  mean(abs(fit$yhat.train[1L, , ] - fit$yhat.train[2L, , ])) > 1.0e-6
)
# chains seed from R's stream, so the fit advances it and set.seed alone
# makes multithreaded results reproducible
expect_false(identical(oldSeed, .Random.seed))
set.seed(99L)
fit2 <- dbarts::bart(
  testData$x,
  testData$y,
  ndpost = 120L,
  nskip = 100L,
  ntree = 50L,
  nthread = 2L,
  nchain = 2L,
  combinechains = FALSE,
  verbose = FALSE
)
expect_equal(fit$yhat.train, fit2$yhat.train)

rm(fit2, fit, oldSeed)

rm(testData)
