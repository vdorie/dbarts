context("bart w/binary response")

source(system.file("common", "probitData.R", package = "dbarts"))

test_that("basic probit example passes regression test", {
  n.sims <- 1500L

  set.seed(99)
  bartFit <- bart(y.train = testData$Z, x.train = testData$X, ntree = 50L, ndpost = n.sims, k = 4.5, verbose = FALSE)

  expect_equal(bartFit$yhat.train[n.sims, 1:5], c(0.177794039813855, 0.587350502370275, 0.925468249706357, 0.369701008671282, 0.0575505774064299))
  expect_identical(bartFit$yhat.test, NULL)
  expect_equal(bartFit$varcount[n.sims,], c(14, 28, 27))
})

test_that("basic probit example with offset regression test", {
  set.seed(99)
  bartFit <- bart(y.train = testData$Z, x.train = testData$X, ntree = 50L, k = 4.5, binaryOffset = 0.1, verbose = FALSE)

  n.sims <- nrow(bartFit$yhat.train)
  
  expect_equal(bartFit$yhat.train[n.sims, 1:5], c(-0.10806177061999, 0.403378273546752, 0.917491762059623, 0.45804917334381, -0.0441185079281768))
  expect_identical(bartFit$yhat.test, NULL)
  expect_equal(bartFit$varcount[n.sims,], c(25, 21, 23))
})

# rm(testData)

source(system.file("common", "almostLinearBinaryData.R", package = "dbarts"))

test_that("binary example using close to linear function provides sensible results", {
  fitSubset  <- 1:100
  testSubset <- 101:200

  fitData <- list(y = testData$y[fitSubset], x = testData$x[fitSubset,])
  glmFit <- glm(y ~ x, fitData, family = binomial(link = "probit"))
  
  predictData <- list(x = testData$x[testSubset,])
  mu.hat.glm <- predict(glmFit, newdata = predictData)
  
  set.seed(99)
  bartFit <- bart(testData$x[fitSubset,], testData$y[fitSubset], testData$x[testSubset,],
                  binaryOffset = testData$offset, verbose = FALSE)
  
  mu.hat.bart <- colMeans(bartFit$yhat.test)
  
  mu <- testData$mu[testSubset]

  expect_true(cor(mu, mu.hat.glm) < cor(mu, mu.hat.bart))
  expect_true((range(mu.hat.bart) * 1.2)[1] >= range(mu)[1])
  expect_true((range(mu.hat.bart) * 1.2)[2] <= range(mu)[2])
})

