source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

# test that predict gives same result as x_train with linear data
bartFit <- dbarts::bart(
  testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L,
  verbose = FALSE, keeptrees = TRUE
)
predictions <- predict(bartFit, testData$x, n.threads = 1L)
expect_equal(predictions, bartFit$yhat.train)

bartFit <- dbarts::bart(
  testData$x, testData$y, ndpost = 20, nskip = 5, ntree = 5L,
  nchain = 4L, nthread = 1L,
  verbose = FALSE, keeptrees = TRUE
)
predictions <- predict(bartFit, testData$x, n.threads = 1L)
expect_equal(predictions, bartFit$yhat.train)

rm(predictions, bartFit)

# test that extract and fitted give correct results
n.chains <- 4L
n.samples <- 20L
bartFit <- dbarts::bart(
  testData$x, testData$y, testData$x[1L:10L,],
  ndpost = n.samples, nskip = 5L, ntree = 5L, nchain = n.chains,
  verbose = FALSE
)

expect_equal(extract(bartFit), bartFit$yhat.train)
expect_equal(fitted(bartFit), bartFit$yhat.train.mean)

expect_equal(extract(bartFit, sample = "test"), bartFit$yhat.test)
expect_equal(fitted(bartFit, sample = "test"), bartFit$yhat.test.mean)

extracted <- extract(bartFit, combineChains = FALSE)
for (i in seq_len(n.chains)) {
  expect_equal(
    extracted[i,,],
    bartFit$yhat.train[seq_len(n.samples) + (i - 1L) * n.samples,]
  )
}

bartFit <- dbarts::bart(
  testData$x, testData$y, testData$x[1L:10L,],
  ndpost = n.samples, nskip = 5, ntree = 5L, nchain = n.chains,
  verbose = FALSE, combinechains = FALSE
)
extracted <- extract(bartFit)
for (i in seq_len(n.chains)) {
  expect_equal(
    extracted[seq_len(n.samples) + (i - 1L) * n.samples,],
    bartFit$yhat.train[i,,]
  )
}

rm(i, extracted, bartFit, n.samples, n.chains)

rm(testData)


source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that predict gives same result as x_train with binary data
bartFit <- dbarts::bart(
  y.train = testData$Z, x.train = testData$X,
  ndpost = 20L, nskip = 5L, ntree = 5L, k = 4.5,
  verbose = FALSE, keeptrees = TRUE
)
predictions <- predict(bartFit, testData$X, type = "bart", n.threads = 1L)
expect_equal(predictions, bartFit$yhat.train)

bartFit <- dbarts::bart(
  y.train = testData$Z, x.train = testData$X,
  ndpost = 20L, nskip = 5L, ntree = 5L, k = 4.5, nchain = 4L, nthread = 1L,
  verbose = FALSE, keeptrees = TRUE
)
predictions <- predict(bartFit, testData$X, type = "bart", n.threads = 1L)
expect_equal(predictions, bartFit$yhat.train)

rm(predictions, bartFit)

rm(testData)

