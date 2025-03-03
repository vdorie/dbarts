context("predict")

source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

test_that("combine/uncombine chains and convert from/to bart style works correctly",
{
  n.chains <- 3L
  n.samples <- 5L
  n.obs <- 7L
  
  dbartsSamples <- array(seq_len(n.chains * n.samples * n.obs),
                         c(n.obs, n.samples, n.chains),
                         dimnames = list(as.character(seq_len(n.obs)), NULL, NULL))

  bartSamples <- dbarts:::convertSamplesFromDbartsToBart(dbartsSamples, n.chains)
  
  expect_equal(dim(bartSamples), c(n.chains, n.samples, n.obs))
  expect_equal(dimnames(bartSamples)[c(3L, 2L, 1L)], dimnames(dbartsSamples))
  
  for (k in seq_len(n.chains))
    expect_equal(t(bartSamples[k,,]), dbartsSamples[,,k])
  
  bartSamples.cc <- dbarts:::convertSamplesFromDbartsToBart(dbartsSamples, n.chains, combineChains = TRUE)
  expect_equal(dim(bartSamples.cc), c(n.chains * n.samples, n.obs))
  expect_equal(colnames(bartSamples.cc), dimnames(dbartsSamples)[[1L]])
  
  for (k in seq_len(n.chains))
    expect_equal(bartSamples[k,,], bartSamples.cc[seq_len(n.samples) + (k - 1L) * n.samples,])
  
  
  expect_equal(dbarts:::uncombineChains(bartSamples.cc, n.chains), bartSamples)
  expect_equal(dbarts:::combineChains(bartSamples), bartSamples.cc)
  
  expect_equal(dbarts:::convertSamplesFromBartsToDbarts(bartSamples, n.chains),
               dbartsSamples)
  expect_equal(dbarts:::convertSamplesFromBartsToDbarts(bartSamples.cc, n.chains, uncombineChains = TRUE),
               dbartsSamples)
})

test_that("predict fails if sampler not saved", {
  bartFit <- bart(testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, verbose = FALSE)
  expect_error(predict(bartFit, testData$x, n.threads = 1L))
})

test_that("predict gives same result as x_train with linear data", {
  bartFit <- bart(testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, verbose = FALSE, keeptrees = TRUE)
  predictions <- predict(bartFit, testData$x, n.threads = 1L)
  expect_equal(predictions, bartFit$yhat.train)
  
  bartFit <- bart(testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, nchain = 4L, nthread = 1L, verbose = FALSE, keeptrees = TRUE)
  predictions <- predict(bartFit, testData$x, n.threads = 1L)
  expect_equal(predictions, bartFit$yhat.train)
})

test_that("predict gives same result when single or multi-threaded", {
  bartFit <- bart(testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, verbose = FALSE, keeptrees = TRUE)
  predictions <- predict(bartFit, testData$x, n.threads = 2L)
  expect_equal(predictions, bartFit$yhat.train)
  
  bartFit <- bart(testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L, nchain = 4L, nthread = 1L, verbose = FALSE, keeptrees = TRUE)
  predictions <- predict(bartFit, testData$x, n.threads = 2L)
  expect_equal(predictions, bartFit$yhat.train)
})

test_that("extract and fitted give correct results", {
  n.chains <- 4L
  n.samples <- 20L
  bartFit <- bart(testData$x, testData$y, testData$x[1:10,], ndpost = n.samples, nskip = 5, ntree = 5L, nchain = n.chains, verbose = FALSE)
  
  expect_equal(extract(bartFit), bartFit$yhat.train)
  expect_equal(fitted(bartFit), bartFit$yhat.train.mean)
  
  expect_equal(extract(bartFit, sample = "test"), bartFit$yhat.test)
  expect_equal(fitted(bartFit, sample = "test"), bartFit$yhat.test.mean)
  
  extracted <- extract(bartFit, combineChains = FALSE)
  for (i in seq_len(n.chains))
    expect_equal(extracted[i,,], bartFit$yhat.train[seq_len(n.samples) + (i - 1L) * n.samples,])
  
  bartFit <- bart(testData$x, testData$y, testData$x[1:10,], ndpost = n.samples, nskip = 5, ntree = 5L, nchain = n.chains, verbose = FALSE, combinechains = FALSE)
  extracted <- extract(bartFit)
  for (i in seq_len(n.chains))
    expect_equal(extracted[seq_len(n.samples) + (i - 1L) * n.samples,], bartFit$yhat.train[i,,])
})

test_that("posterior predictive distribution samples use correct sigma", {
  n.samples <- 7L
  n.chains  <- 2L
  n.obs <- length(testData$y)
  bartFit <- bart(testData$x, testData$y, verbose = FALSE,
                  ndpost = n.samples, nskip = 0L, nchain = n.chains,
                  ntree = 25L, nthread = 1L)
  set.seed(0)
  samples.ppd <- extract(bartFit, type = "ppd")
  
  set.seed(0)
  samples.pm  <- extract(bartFit)
  for (i in seq_len(n.obs))
    expect_equal(samples.pm[,i] + rnorm(n.samples * n.chains, 0, bartFit$sigma),
                 samples.ppd[,i])

  set.seed(0)
  samples.ppd <- extract(bartFit, type = "ppd", combineChains = FALSE)
  
  set.seed(0)
  samples.pm  <- extract(bartFit, combineChains = FALSE)
  for (i in seq_len(n.obs))
    expect_equal(samples.pm[,,i] + matrix(rnorm(n.samples * n.chains, 0, bartFit$sigma), nrow = n.chains),
                 samples.ppd[,,i])
})

test_that("fixed sample mode when run sequentially gives same predictions as sequential updates mode", {
  set.seed(0)
  pred.bart <- bart2(testData$x, testData$y, testData$x, n.samples = 5, n.burn = 0L, n.trees = 4L,
                     k = 2, n.chains = 1L, n.threads = 1L, keepTrees = TRUE, verbose = FALSE)$yhat.test
  
  set.seed(0)
  sampler <- dbarts(testData$x, testData$y,
                    control = dbartsControl(n.samples = 5, n.burn = 0L, n.trees = 4L,
                    n.chains = 1L, n.threads = 1L, keepTrees = TRUE, updateState = FALSE))
  sampler$sampleTreesFromPrior()
  for (i in seq_len(5L))
    invisible(sampler$run(0L, 1L))
  pred.dbarts <- sampler$predict(testData$x)
  
  expect_equal(pred.bart, t(pred.dbarts))
})

test_that("sequentially running samples don't overflow with fixed trees", {
  sampler <- dbarts(testData$x, testData$y,
                    control = dbartsControl(n.samples = 5, n.burn = 0L, n.trees = 4L, n.chains = 1L, n.threads = 1L, keepTrees = TRUE, updateState = FALSE))
  sampler$sampleTreesFromPrior()
  for (i in seq_len(6L))
    invisible(sampler$run(0L, 1L))
  
  expect_is(sampler, "dbartsSampler")
})

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

test_that("predict gives same result as x_train with binary data", {
  bartFit <- bart(y.train = testData$Z, x.train = testData$X, ndpost = 20, nskip = 5, ntree = 5L, k = 4.5, verbose = FALSE, keeptrees = TRUE)
  predictions <- predict(bartFit, testData$X, type = "bart", n.threads = 1L)
  expect_equal(predictions, bartFit$yhat.train)
  
  bartFit <- bart(y.train = testData$Z, x.train = testData$X, ndpost = 20, nskip = 5, ntree = 5L, k = 4.5, nchain = 4L, nthread = 1L, verbose = FALSE, keeptrees = TRUE)
  predictions <- predict(bartFit, testData$X, type = "bart", n.threads = 1L)
  expect_equal(predictions, bartFit$yhat.train)
})

