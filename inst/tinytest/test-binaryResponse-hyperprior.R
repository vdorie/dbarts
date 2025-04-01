source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

# test that basic probit example with flat hyperprior superior to default
n.sims <- 200L
n.burn <- 100L

set.seed(99L)
bartFit <- dbarts::bart(
  y.train = testData$Z, x.train = testData$X,
  ntree = 50L, ndpost = n.sims, nskip = n.burn,
  verbose = FALSE
)

set.seed(99L)
bartFit.flat <- dbarts::bart(
  y.train = testData$Z, x.train = testData$X,
  ntree = 50L, ndpost = n.sims, nskip = n.burn,
  k = chi(1, Inf),
  verbose = FALSE
)

expect_true(
  cor(qnorm(testData$p), colMeans(bartFit$yhat.train))
  < cor(qnorm(testData$p), colMeans(bartFit.flat$yhat.train))
)
rm(bartFit.flat, bartFit, n.burn, n.sims)

# test_that binary model with k hyperprior is reproducible when multithreaded
fit1 <- dbarts::bart2(
  testData$X[1L:100L,], testData$Z[1L:100L],
  n.trees = 5L, n.samples = 100L, n.burn = 0L,
  n.threads = 2L, n.chains = 2L, rngSeed = 99L,
  verbose = FALSE
)
fit2 <- dbarts::bart2(
  testData$X[1L:100L,], testData$Z[1L:100L],
  n.trees = 5L, n.samples = 100L, n.burn = 0L,
  n.threads = 2L, n.chains = 2L, rngSeed = 99L,
  verbose = FALSE
)
expect_equal(fit1$yhat.train, fit2$yhat.train)
rm(fit2, fit1)

source(system.file("common", "almostLinearBinaryData.R", package = "dbarts"), local = TRUE)

fitSubset  <- 1L:100L
testSubset <- 101L:200L

fitData <- list(y = testData$y[fitSubset], x = testData$x[fitSubset,])
mu <- testData$mu[testSubset]

glmFit <- stats::glm(y ~ x, fitData, family = binomial(link = "probit"))

predictData <- list(x = testData$x[testSubset,])
mu.hat.glm <- predict(glmFit, newdata = predictData)


set.seed(99L)
bartFit <- dbarts::bart(
  testData$x[fitSubset,], testData$y[fitSubset], testData$x[testSubset,],
  binaryOffset = testData$offset,
  verbose = FALSE
)
mu.hat.bart <- colMeans(bartFit$yhat.test)

# test that binary example using close to linear function provides sensible results
expect_true(cor(mu, mu.hat.glm) < cor(mu, mu.hat.bart))
expect_true((range(mu.hat.bart) * 1.2)[1L] >= range(mu)[1L])
expect_true((range(mu.hat.bart) * 1.2)[2L] <= range(mu)[2L])

# test_that binary example using a flat prior is similar to default in tuned model
set.seed(99L)
bartFit.flat <- dbarts::bart(
  testData$x[fitSubset,], testData$y[fitSubset], testData$x[testSubset,],
  binaryOffset = testData$offset,
  verbose = FALSE, k = chi(1, Inf)
)
mu.hat.bart.flat <- colMeans(bartFit.flat$yhat.test)


expect_true(cor(mu.hat.bart, mu.hat.bart.flat) > 0.99)
expect_true(median(bartFit.flat$k) < 3)

rm(
  mu.hat.bart.flat, bartFit.flat,
  mu.hat.bart, bartFit,
  mu.hat.glm, predictData, glmFit,
  mu, fitData, testSubset, fitSubset
)

rm(testData)

