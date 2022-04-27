context("rbart")

source(system.file("common", "friedmanData.R", package = "dbarts"), local = TRUE)

n.g <- 5L
if (getRversion() >= "3.6.0") {
  oldSampleKind <- RNGkind()[3L]
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}
g <- sample(n.g, length(testData$y), replace = TRUE)
if (getRversion() >= "3.6.0") {
  suppressWarnings(RNGkind(sample.kind = oldSampleKind))
  rm(oldSampleKind)
}

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
                     n.trees = 25L, n.threads = 1L, verbose = FALSE),
            "rbart")
  
  g <- df$g
  df$g <- NULL
  expect_is(rbart_vi(y ~ . , df, group.by = g,
                     n.samples = 1L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                     n.trees = 25L, n.threads = 1L, verbose = FALSE),
            "rbart")
  
  y <- testData$y
  x <- testData$x
  expect_is(rbart_vi(y ~ x, group.by = g,
                     n.samples = 1L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                     n.trees = 25L, n.threads = 1L, verbose = FALSE),
            "rbart")
})

test_that("works with multiple threads", {
  x <- testData$x
  y <- testData$y
  g <- factor(testData$g)
  
  set.seed(0)
  expect_is(rbart_vi(y ~ x, group.by = g,
                     n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                     n.trees = 25L, n.threads = 2L, verbose = FALSE),
            "rbart")
})

test_that("is reproducible", {
  x <- testData$x
  y <- testData$y
  g <- factor(testData$g)
  
  fit1 <- rbart_vi(y ~ x, group.by = g,
                   n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                   n.trees = 3L, n.threads = 2L, verbose = FALSE,
                   seed = 0L)
  fit2 <- rbart_vi(y ~ x, group.by = g,
                   n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                   n.trees = 3L, n.threads = 2L, verbose = FALSE,
                   seed = 0L)
  
  set.seed(0)
  seeds <- sample.int(.Machine$integer.max, 2L)
  
  set.seed(seeds[1L])
  fit3 <- rbart_vi(y ~ x, group.by = g,
                   n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                   n.trees = 3L, n.threads = 1L, verbose = FALSE)
  
  set.seed(seeds[2L])
  fit4 <- rbart_vi(y ~ x, group.by = g,
                   n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                   n.trees = 3L, n.threads = 1L, verbose = FALSE)
  
  expect_equal(fit1$yhat.train, fit2$yhat.train)
  expect_equal(fit1$ranef, fit2$ranef)
  
  
  yhat <- aperm(array(c(fit3$yhat.train, fit4$yhat.train),
                      c(dim(fit3$yhat.train), 2L)),
                c(3L, 1L, 2L))
  expect_equal(yhat, fit1$yhat.train)
  
  ranef <- aperm(array(c(fit3$ranef, fit4$ranef),
                       c(dim(fit3$ranef), 2L)),
                 c(3L, 1L:2L))
  expect_equal(as.vector(ranef), as.vector(fit1$ranef))
})

test_that("extract works at baseline", {
  x <- testData$x
  y <- testData$y
  g <- factor(testData$g)
  
  set.seed(0)
  rbartFit <- rbart_vi(y ~ x, group.by = g,
                       n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                       n.trees = 25L, n.threads = 1L, verbose = FALSE)
  expect_equal(extract(rbartFit, type = "bart", combineChains = FALSE), rbartFit$yhat.train)
  expect_equal(extract(rbartFit, type = "ranef", combineChains = FALSE), rbartFit$ranef)
  expect_equal(extract(rbartFit, type = "ev", combineChains = FALSE), rbartFit$yhat.train + unname(rbartFit$ranef[,,as.character(g)]))
  
  ppd <- extract(rbartFit, type = "ppd", combineChains = FALSE)
  sigma.hat <- apply(ppd - extract(rbartFit, type = "ev", combineChains = FALSE), c(1L, 2L), sd)
  expect_true(cor(as.vector(sigma.hat), as.vector(rbartFit$sigma)) >= 0.85) # silly test with 7 samples
  
  set.seed(0)
  rbartFit.2 <- rbart_vi(y ~ x, group.by = g,
                         n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                         n.trees = 25L, n.threads = 1L, combineChains = TRUE,
                         verbose = FALSE)
  expect_equal(extract(rbartFit,   type = "ev", combineChains = TRUE),
               extract(rbartFit.2, type = "ev", combineChains = TRUE))
  expect_equal(extract(rbartFit,   type = "ev", combineChains = FALSE),
               extract(rbartFit.2, type = "ev", combineChains = FALSE))
})

test_that("fitted works correctly", {
  x <- testData$x
  y <- testData$y
  g <- factor(testData$g)

  
  rbartFit <- rbart_vi(y ~ x, group.by = g, test = x, group.by.test = g, offset.test = 5,
                       n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                       n.trees = 25L, n.threads = 1L, verbose = FALSE)
  expect_equal(as.vector(rbartFit$yhat.train), as.vector(rbartFit$yhat.test) - 5)
  expect_equal(apply(rbartFit$yhat.train + unname(rbartFit$ranef[,,as.character(g)]), 3L, mean), fitted(rbartFit))
  expect_equal(apply(rbartFit$yhat.test  + unname(rbartFit$ranef[,,as.character(g)]), 3L, mean), fitted(rbartFit, sample = "test"))
  expect_equal(fitted(rbartFit, type = "ranef"), rbartFit$ranef.mean)
})

test_that("predict matches fitted", {
  x <- testData$x
  y <- testData$y
  g <- factor(testData$g)

  set.seed(0)
  rbartFit.0 <- rbart_vi(y ~ x, group.by = g,
                         n.samples = 14L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                         n.trees = 25L, n.threads = 1L, keepTrees = TRUE,
                         verbose = FALSE)
  expect_equal(fitted(rbartFit.0), apply(predict(rbartFit.0, x, g), 2L, mean))

  set.seed(0)
  rbartFit.0 <- rbart_vi(y ~ x, group.by = g,
                         n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                         n.trees = 25L, n.threads = 1L, keepTrees = TRUE, combineChains = FALSE,
                         verbose = FALSE)
  expect_equal(fitted(rbartFit.0), apply(predict(rbartFit.0, x, g, combineChains = FALSE), 3L, mean))
  
  set.seed(0)
  rbartFit.1 <- rbart_vi(y ~ x, group.by = g,
                       n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                       n.trees = 25L, n.threads = 1L, keepTrees = TRUE, combineChains = TRUE,
                       verbose = FALSE)
  expect_equal(fitted(rbartFit.1), apply(predict(rbartFit.1, x, g), 2L, mean))
  expect_equal(predict(rbartFit.0, x, g), predict(rbartFit.1, x, g))
  expect_equal(predict(rbartFit.0, x, g, combineChains = FALSE), predict(rbartFit.1, x, g, combineChains = FALSE))
})

test_that("works with missing levels", {
  n.train <- 80L
  x <- testData$x[seq_len(n.train),]
  y <- testData$y[seq_len(n.train)]
  g <- factor(testData$g[seq_len(n.train)])

  x.test <- testData$x[seq.int(n.train + 1L, nrow(testData$x)),]
  g.test <- factor(testData$g[seq.int(n.train + 1L, nrow(testData$x))], levels(g))
  levels(g.test)[5L] <- "6"
  
  # check that predict works when we've fit with missing levels
  rbartFit <- suppressWarnings(rbart_vi(y ~ x, group.by = g, test = x.test, group.by.test = g.test,
                                        n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                                        n.trees = 25L, n.threads = 1L, keepTrees = TRUE,
                                        verbose = FALSE))
  expect_equal(apply(predict(rbartFit, x.test, g.test), 2L, mean), fitted(rbartFit, sample = "test"))
  expect_equal(apply(predict(rbartFit, x.test, g.test, combineChains = FALSE), 3L, mean), fitted(rbartFit, sample = "test"))
  
  # check that predicts works for completely new levels
  levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
  set.seed(0)
  ranef.pred <- suppressWarnings(predict(rbartFit, x.test, g.test, type = "ranef", combineChains = FALSE))
  expect_equal(ranef.pred[,,as.character(1L:4L)], rbartFit$ranef[,,as.character(1L:4L)])
  expect_true(cor(as.numeric(rbartFit$tau), as.numeric(apply(ranef.pred[,,5L:26L], c(1L, 2L), sd))) > 0.90)
  
  # check again with combineChains as TRUE at the top level
  g.test <- droplevels(g.test)
  levels(g.test) <- c(levels(g)[-5L], "6")
  rbartFit <- suppressWarnings(rbart_vi(y ~ x, group.by = g, test = x.test, group.by.test = g.test,
                                        n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
                                        n.trees = 25L, n.threads = 1L, keepTrees = TRUE, combineChains = TRUE,
                                        verbose = FALSE))
  expect_equal(apply(predict(rbartFit, x.test, g.test), 2L, mean), fitted(rbartFit, sample = "test"))
  expect_equal(apply(predict(rbartFit, x.test, g.test, combineChains = FALSE), 3L, mean), fitted(rbartFit, sample = "test"))
  
  levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
  set.seed(0)
  ranef.pred <- suppressWarnings(predict(rbartFit, x.test, g.test, type = "ranef"))
  expect_equal(as.numeric(ranef.pred[,as.character(1L:4L)]),
               as.numeric(rbartFit$ranef[,as.character(1L:4L)]))
  expect_true(cor(as.numeric(rbartFit$tau), as.numeric(apply(ranef.pred[,5L:26L], 1L, sd))) > 0.90)
  
  # check one last time with one chain
  g.test <- droplevels(g.test)
  levels(g.test) <- c(levels(g)[-5L], "6")
  rbartFit <- suppressWarnings(rbart_vi(y ~ x, group.by = g, test = x.test, group.by.test = g.test,
                                        n.samples = 14L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
                                        n.trees = 25L, n.threads = 1L, keepTrees = TRUE,
                                        verbose = FALSE))
  expect_equal(apply(predict(rbartFit, x.test, g.test), 2L, mean), fitted(rbartFit, sample = "test"))
  levels(g.test) <- c(levels(g.test)[-5L], as.character(seq.int(7L, 28L)))
  set.seed(0)
  ranef.pred <- suppressWarnings(predict(rbartFit, x.test, g.test, type = "ranef"))
  expect_equal(as.numeric(ranef.pred[,as.character(1L:4L)]),
               as.numeric(rbartFit$ranef[,as.character(1L:4L)]))
  expect_true(cor(as.numeric(rbartFit$tau), as.numeric(apply(ranef.pred[,5L:26L], 1L, sd))) > 0.90)
})

test_that("rbart runs example", {
  rbartFit <- rbart_vi(y ~ x, testData, group.by = g,
                       n.samples = 40L, n.burn = 10L, n.thin = 2L, n.chains = 2L,
                       n.trees = 25L, n.threads = 1L,
                       verbose = FALSE)
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
                       n.trees = 25L, n.threads = 1L,
                       verbose = FALSE)
  expect_equal(dim(rbartFit$yhat.train), c(80L %/% 2L, length(testData$y)))
  expect_equal(length(rbartFit$yhat.train.mean), length(testData$y))
  expect_equal(dim(rbartFit$ranef), c(80L %/% 2L, length(unique(testData$g))))
  expect_equal(length(rbartFit$first.tau), 20L %/% 2L)
  expect_equal(length(rbartFit$first.sigma), 20L %/% 2L)
  expect_equal(length(rbartFit$tau), 80L %/% 2L)
  expect_equal(length(rbartFit$sigma), 80L %/% 2L)
  
  expect_true(length(unique(rbartFit$ranef)) > 1L)
})

test_that("rbart works with keepTrainingFits = FALSE", {
  rbartFit <- rbart_vi(y ~ x, testData, group.by = g,
                       n.samples = 2L, n.burn = 0L, n.thin = 2L, n.chains = 2L,
                       n.trees = 3L, n.threads = 1L,
                       keepTrainingFits = FALSE,
                       verbose = FALSE)
  expect_is(rbartFit, "rbart")
  expect_true(is.null(rbartFit$yhat.train))
  expect_true(is.null(rbartFit$yhat.train.mean))
})

test_that("rbart passes regression test", {
  df <- as.data.frame(testData$x)
  colnames(df) <- paste0("x_", seq_len(ncol(testData$x)))
  df$y <- testData$y
  df$g <- testData$g
  
  set.seed(99)
  rbartFit <- rbart_vi(y ~ . - g, df, group.by = g,
                       n.samples = 1L, n.burn = 5L, n.thin = 1L, n.chains = 1L,
                       n.trees = 25L, n.threads = 1L,
                       verbose = FALSE)
  
  expect_equal(as.numeric(rbartFit$ranef),
               c(1.92008577165928, 0.750130559201404, -0.837624979286864, 0.461299843412371, 3.2085237309599))
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
                       n.trees = 50L, n.threads = 1L,
                       verbose = FALSE)
  ranef.rbart <- rbartFit$ranef.mean
  
  lmerFit <- suppressWarnings(lme4$lmer(y ~ . - g + (1 | g), df))
  ranef.lmer <- lme4$ranef.merMod(lmerFit)[[1L]][[1L]]
  
  b_rmse.rbart <- sqrt(mean((b - ranef.rbart)^2))
  b_rmse.lmer <- sqrt(mean((b - ranef.lmer)^2))
  expect_true(b_rmse.rbart < b_rmse.lmer)
  
  
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
                       n.trees = 50L, n.threads = 1L,
                       verbose = FALSE)
  ranef.rbart <- rbartFit$ranef.mean
  
  glmerFit <- lme4$glmer(y ~ . - g + (1 | g), df, family = binomial(link = "probit"))
  
  rbart.mu.hat <- fitted(rbartFit)
  expect_equal(apply(extract(rbartFit), 2L, mean), rbart.mu.hat)
  glmer.mu.hat <- fitted(glmerFit, type = "response")
  
  dev.rbart <- -2 * mean(log(ifelse(df$y == 1, rbart.mu.hat, 1 - rbart.mu.hat)))
  dev.glmer <- -2 * mean(log(ifelse(df$y == 1, glmer.mu.hat, 1 - glmer.mu.hat)))
  expect_true(dev.rbart < dev.glmer)
})

