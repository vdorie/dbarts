context("rbart")

source(system.file("common", "friedmanData.R", package = "dbarts"))

n.g <- 5L
g <- sample(n.g, length(testData$y), replace = TRUE)
sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

testData$y <- testData$y + b[g]
testData$g <- g
testData$b <- b
rm(g, b)

test_that("rbart fails with invalid group.by", {
  expect_error(rbart_vi(y ~ x, testData, group.by = NA))
  expect_error(rbart_vi(y ~ x, testData, group.by = not_a_symbol))
  expect_error(rbart_vi(y ~ x, testData, group.by = testData$g[-1L]))
  expect_error(rbart_vi(y ~ x, testData, group.by = "not a factor"))
})

test_that("rbart finds group.by", {
  df <- as.data.frame(testData$x)
  colnames(df) <- paste0("x_", seq_len(ncol(testData$x)))
  df$y <- testData$y
  df$g <- testData$g
  expect_is(rbart_vi(y ~ . - g, df, group.by = g,
                     n.samples = 1L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                     n.trees = 25L, n.threads = 1L),
            "rbart")
  
  g <- df$g
  df$g <- NULL
  expect_is(rbart_vi(y ~ . - g, df, group.by = g,
                     n.samples = 1L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                     n.trees = 25L, n.threads = 1L),
            "rbart")
  
  y <- testData$y
  x <- testData$x
  expect_is(rbart_vi(y ~ x, group.by = g,
                     n.samples = 1L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                     n.trees = 25L, n.threads = 1L),
            "rbart")
  
  rm(y, x, g, df)
})

test_that("rbart runs example", {
  rbartFit <- rbart_vi(y ~ x, testData, group.by = g,
                       n.samples = 40L, n.burn = 10L, n.thin = 2L, n.chains = 2L,
                       n.trees = 25L, n.threads = 1L)
  expect_equal(dim(rbartFit$yhat.train), c(2L, 40L %/% 2L, length(testData$y)))
  expect_equal(length(rbartFit$yhat.train.mean), length(testData$y))
  expect_equal(dim(rbartFit$ranef), c(2L, 40L %/% 2L, length(unique(testData$g))))
  expect_equal(dim(rbartFit$first.tau), c(2L, 10L %/% 2L))
  expect_equal(dim(rbartFit$first.sigma), c(2L, 10L %/% 2L))
  expect_equal(dim(rbartFit$tau), c(2L, 40L %/% 2L))
  expect_equal(dim(rbartFit$sigma), c(2L, 40L %/% 2L))
  
  expect_true(length(unique(rbartFit$ranef)) > 1L)

  # check for one chain
  rbartFit <- rbart_vi(y ~ x, testData, group.by = g,
                       n.samples = 80L, n.burn = 20L, n.thin = 2L, n.chains = 1L,
                       n.trees = 25L, n.threads = 1L)
  expect_equal(dim(rbartFit$yhat.train), c(80L %/% 2L, length(testData$y)))
  expect_equal(length(rbartFit$yhat.train.mean), length(testData$y))
  expect_equal(dim(rbartFit$ranef), c(80L %/% 2L, length(unique(testData$g))))
  expect_equal(length(rbartFit$first.tau), 20L %/% 2L)
  expect_equal(length(rbartFit$first.sigma), 20L %/% 2L)
  expect_equal(length(rbartFit$tau), 80L %/% 2L)
  expect_equal(length(rbartFit$sigma), 80L %/% 2L)
  
  expect_true(length(unique(rbartFit$ranef)) > 1L)
})

test_that("rbart passes regression test", {
  df <- as.data.frame(testData$x)
  colnames(df) <- paste0("x_", seq_len(ncol(testData$x)))
  df$y <- testData$y
  df$g <- testData$g
  
  set.seed(99)
  rbartFit <- rbart_vi(y ~ . - g, df, group.by = g,
                       n.samples = 1L, n.burn = 5L, n.thin = 1L, n.chains = 1L,
                       n.trees = 25L, n.threads = 1L)
  
  expect_equal(as.numeric(rbartFit$ranef),
               c(0.548756620200975, 1.98489377073739, -0.123942881873723, -0.643642914323586, 3.02981874312062))
})

test_that("rbart compares favorably to lmer for nonlinear models", {
  skip_if_not_installed("lme4")
  lme4 <- asNamespace("lme4")
  
  f <- function(x) {
      10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]
  }
  
  set.seed(99)
  sigma <- 1.0
  n     <- 500
  
  x  <- matrix(runif(n * 10), n, 10)
  Ey <- f(x)
  y  <- rnorm(n, Ey, sigma)
  
  n.g <- 15
  g <- sample(n.g, length(y), replace = TRUE)
  sigma.b <- 1.5
  b <- rnorm(n.g, 0, sigma.b)
  
  y <- y + b[g]
  
  df <- as.data.frame(x)
  colnames(df) <- paste0("x_", seq_len(ncol(x)))
  df$y <- y
  df$g <- g
  
  
  rbartFit <- rbart_vi(y ~ . - g, df, group.by = g,
                       n.samples = 200L, n.burn = 100L, n.thin = 2L, n.chains = 2L,
                       n.trees = 50L, n.threads = 1L)
  ranef.rbart <- rbartFit$ranef.mean
  
  lmerFit <- lme4$lmer(y ~ . - g + (1 | g), df)
  ranef.lmer <- lme4$ranef.merMod(lmerFit)[[1L]][[1L]]
  
  expect_true(sqrt(mean((b - ranef.rbart)^2)) < sqrt(mean((b - ranef.lmer)^2)))
  
  
  f <- function(x) {
      10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]
  }
  
  rho <- 0.4
  p.y <- pnorm((Ey - mean(Ey)) / sd(Ey) + rho * .75 * b[g])
  set.seed(99)
  y <- rbinom(n, 1L, p.y)
  df <- as.data.frame(x)
  colnames(df) <- paste0("x_", seq_len(ncol(x)))
  df$y <- y
  df$g <- g
  
  rbartFit <- rbart_vi(y ~ . - g, df, group.by = g,
                       n.samples = 240L, n.burn = 120L, n.thin = 3L, n.chains = 2L,
                       n.trees = 50L, n.threads = 1L)
  ranef.rbart <- rbartFit$ranef.mean
  
  lmerFit <- lme4$glmer(y ~ . - g + (1 | g), df, family = binomial(link = probit))
  ranef.lmer <- lme4$ranef.merMod(lmerFit)[[1L]][[1L]]
  
  rbart.mu.hat <- apply(rbartFit$yhat.train, 3, mean)
  lmer.mu.hat  <- predict(lmerFit)
  expect_true(sqrt(mean((rbart.mu.hat - Ey)^2)) < sqrt(mean((lmer.mu.hat - Ey)^2)))
})

