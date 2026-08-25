source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

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

## combineChains = FALSE pinned deliberately: the per-chain independence
## check below indexes fit1$yhat.train as a raw uncombined 3-d array
fit1 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 5L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 3L,
  n.threads = 2L,
  verbose = FALSE,
  seed = 0L,
  combineChains = FALSE
)
fit2 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 5L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 3L,
  n.threads = 2L,
  verbose = FALSE,
  seed = 0L,
  combineChains = FALSE
)

expect_equal(fit1$yhat.train, fit2$yhat.train)
expect_equal(fit1$ranef, fit2$ranef)
expect_equal(fit1$tau, fit2$tau)

# the chains explore independently rather than repeating one stream
expect_true(any(fit1$yhat.train[1L, , ] != fit1$yhat.train[2L, , ]))

# single-chain runs draw through R's generator and reproduce under set.seed
set.seed(0L)
fit3 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 5L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 3L,
  n.threads = 1L,
  verbose = FALSE
)

set.seed(0L)
fit4 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 5L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 3L,
  n.threads = 1L,
  verbose = FALSE
)

expect_equal(fit3$yhat.train, fit4$yhat.train)
expect_equal(fit3$ranef, fit4$ranef)
expect_equal(fit3$tau, fit4$tau)

# a supplied seed fixes the draws on the single-chain (in-core) path too,
# independently of where R's own generator happens to stand
seededFit <- function(seed) {
  dbarts::rbart_vi(
    y ~ x,
    group.by = g,
    n.samples = 5L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 3L,
    n.threads = 1L,
    verbose = FALSE,
    seed = seed
  )
}
set.seed(1L)
fit5 <- seededFit(7L)
set.seed(2L)
fit6 <- seededFit(7L)
expect_equal(fit5$yhat.train, fit6$yhat.train)
expect_equal(fit5$ranef, fit6$ranef)
expect_equal(fit5$tau, fit6$tau)
# and the seed is what does it: a different seed moves the draws
fit7 <- seededFit(8L)
expect_false(isTRUE(all.equal(fit5$yhat.train, fit7$yhat.train)))

rm(fit7, fit6, fit5, fit4, fit3, fit2, fit1, seededFit)

rm(g, y, x)

rm(testData)
