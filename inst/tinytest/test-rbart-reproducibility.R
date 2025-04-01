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

set.seed(0L)
seeds <- sample.int(.Machine$integer.max, 2L)

set.seed(seeds[1L])
fit3 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
  n.trees = 3L, n.threads = 1L, verbose = FALSE
)

set.seed(seeds[2L])
fit4 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 5L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
  n.trees = 3L, n.threads = 1L, verbose = FALSE
)

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

rm(ranef, yhat, fit4, fit3, seeds, fit2, fit1)

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

expect_equal(
  as.numeric(rbartFit$ranef),
  c(
    1.92008577165928, 0.750130559201404, -0.837624979286864,
    0.461299843412371, 3.2085237309599
  )
)

rm(rbartFit, df)

rm(testData)

