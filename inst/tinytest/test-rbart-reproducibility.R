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
rm(b, sigma.b, g, n.g)

# test that is reproducible
x <- testData$x
y <- testData$y
g <- factor(testData$g)

fit1 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 3L, n.threads = 2L, verbose = FALSE,
  seed = 0L
)
fit2 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 3L, n.threads = 2L, verbose = FALSE,
  seed = 0L
)

expect_equal(fit1$yhat.train, fit2$yhat.train)
expect_equal(fit1$ranef, fit2$ranef)
expect_equal(fit1$tau, fit2$tau)

# the chains explore independently rather than repeating one stream
expect_true(any(fit1$yhat.train[1L, , ] != fit1$yhat.train[2L, , ]))

# single-chain runs draw through R's generator and reproduce under set.seed
set.seed(0L)
fit3 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
  n.trees = 3L, n.threads = 1L, verbose = FALSE
)

set.seed(0L)
fit4 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
  n.trees = 3L, n.threads = 1L, verbose = FALSE
)

expect_equal(fit3$yhat.train, fit4$yhat.train)
expect_equal(fit3$ranef, fit4$ranef)
expect_equal(fit3$tau, fit4$tau)

rm(fit4, fit3, fit2, fit1)

rm(g, y, x)


# test that rbart passes regression test
df <- as.data.frame(testData$x)
colnames(df) <- paste0("x_", seq_len(ncol(testData$x)))
df$y <- testData$y
df$g <- testData$g

set.seed(99L)
rbartFit <- dbarts::rbart_vi(
  y ~ . - g, df, group.by = g,
  n.samples = 1L, n.burn = 5L, n.thin = 1L, n.chains = 1L,
  n.trees = 25L, n.threads = 1L,
  verbose = FALSE
)

# regenerate by replaying this file from the top: draws depend on the
# preceding fits through R's stream position and last-ulp arithmetic
expect_equal(
  as.numeric(rbartFit$ranef),
  c(
    -0.633191770812314, -0.300187735171359, -1.18793374547573,
    -0.47732899393382, 0.855568120506501
  )
)

rm(rbartFit, df)

rm(testData)

