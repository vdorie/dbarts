context("bart w/binary response")

source(system.file("common", "probitData.R", package = "dbarts"), local = TRUE)

test_that("basic probit example passes regression test", {
  n.burn <- 10L
  n.sims <- 100L
  
  set.seed(99)
  bartFit <- bart(y.train = testData$Z, x.train = testData$X, ntree = 50L, ndpost = n.sims, nskip = n.burn,
                  k = 4.5, verbose = FALSE)
  
  expect_equal(bartFit$yhat.train[n.sims, 1:5], c(0.360083720993859, 0.213898385154795, 0.514888279642085, 0.402547682652614, 0.0641376173491096))
  expect_identical(bartFit$yhat.test, NULL)
  expect_equal(bartFit$varcount[n.sims,], c(19, 26, 24))
  
  expect_equal(extract(bartFit), pnorm(bartFit$yhat.train))
})

test_that("basic probit example with offset regression test", {
  n.burn <- 10L
  n.sims <- 100L
  
  set.seed(99)
  bartFit <- bart(y.train = testData$Z, x.train = testData$X, ntree = 50L, ndpost = n.sims, nskip = n.burn,
                  k = 4.5, binaryOffset = 0.1, verbose = FALSE)
  
  n.sims <- nrow(bartFit$yhat.train)
  
  expect_equal(bartFit$yhat.train[n.sims, 1:5], c(0.157043723005439, 0.649674546901119, 0.392826725618914, 0.510142912804732, -0.27263185358599))
  expect_identical(bartFit$yhat.test, NULL)
  expect_equal(bartFit$varcount[n.sims,], c(32, 21, 21))
})

test_that("basic probit example with flat hyperprior superior to default", {
  n.sims <- 200L
  n.burn <- 100L
  
  set.seed(99)
  bartFit <- bart(y.train = testData$Z, x.train = testData$X, ntree = 50L, ndpost = n.sims, nskip = n.burn,
                  verbose = FALSE)
  
  set.seed(99)
  bartFit.flat <- bart(y.train = testData$Z, x.train = testData$X, ntree = 50L, ndpost = n.sims, nskip = n.burn,
                       k = chi(1, Inf), verbose = FALSE)
  
  expect_true(cor(qnorm(testData$p), colMeans(bartFit$yhat.train)) <
              cor(qnorm(testData$p), colMeans(bartFit.flat$yhat.train)))
})

test_that("binary model with k hyperprior is reproducible when multithreaded", {
  fit1 <- bart2(testData$X[1:100,], testData$Z[1:100], n.trees = 5L,
                n.samples = 100L, n.burn = 0L, verbose = FALSE,
                n.threads = 2L, n.chains = 2L, rngSeed = 99)
  fit2 <- bart2(testData$X[1:100,], testData$Z[1:100], n.trees = 5L,
                n.samples = 100L, n.burn = 0L, verbose = FALSE,
                n.threads = 2L, n.chains = 2L, rngSeed = 99)
  expect_equal(fit1$yhat.train, fit2$yhat.train)
})

source(system.file("common", "almostLinearBinaryData.R", package = "dbarts"), local = TRUE)

fitSubset  <- 1:100
testSubset <- 101:200

fitData <- list(y = testData$y[fitSubset], x = testData$x[fitSubset,])
mu <- testData$mu[testSubset]

glmFit <- glm(y ~ x, fitData, family = binomial(link = "probit"))

predictData <- list(x = testData$x[testSubset,])
mu.hat.glm <- predict(glmFit, newdata = predictData)

set.seed(99)
bartFit <- bart(testData$x[fitSubset,], testData$y[fitSubset], testData$x[testSubset,],
                binaryOffset = testData$offset, verbose = FALSE)
mu.hat.bart <- colMeans(bartFit$yhat.test)

test_that("binary example using close to linear function provides sensible results", {
  expect_true(cor(mu, mu.hat.glm) < cor(mu, mu.hat.bart))
  expect_true((range(mu.hat.bart) * 1.2)[1] >= range(mu)[1])
  expect_true((range(mu.hat.bart) * 1.2)[2] <= range(mu)[2])
})

test_that("binary example using a flat prior is similar to default in tuned model", {
  set.seed(99)
  bartFit.flat <- bart(testData$x[fitSubset,], testData$y[fitSubset], testData$x[testSubset,],
                       binaryOffset = testData$offset, verbose = FALSE, k = chi(1, Inf))
  mu.hat.bart.flat <- colMeans(bartFit.flat$yhat.test)
  
  
  expect_true(cor(mu.hat.bart, mu.hat.bart.flat) > 0.99)
  expect_true(median(bartFit.flat$k) < 3)
})

rm(fitSubset, testSubset, fitData, mu, glmFit, predictData,
   mu.hat.glm, bartFit, mu.hat.bart, testData)

