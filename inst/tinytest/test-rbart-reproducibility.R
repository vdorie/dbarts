source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# the port's deprecation warning is pinned in test-rbart-port.R
onceState <- dbarts:::onceWarnState
onceState[["tombstone.rbart_vi"]] <- TRUE
rm(onceState)

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
  y ~ x,
  group.by = g,
  n.samples = 5L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 3L,
  n.threads = 2L,
  verbose = FALSE,
  seed = 0L
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
  seed = 0L
)

set.seed(0L)
seeds <- sample.int(.Machine$integer.max, 2L)

set.seed(seeds[1L])
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

set.seed(seeds[2L])
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

expect_equal(fit1$yhat.train, fit2$yhat.train)
expect_equal(fit1$ranef, fit2$ranef)


yhat <- aperm(
  array(c(fit3$yhat.train, fit4$yhat.train), c(dim(fit3$yhat.train), 2L)),
  c(3L, 1L, 2L)
)
expect_equal(yhat, fit1$yhat.train)

ranef <- aperm(
  array(c(fit3$ranef, fit4$ranef), c(dim(fit3$ranef), 2L)),
  c(3L, 1L:2L)
)
expect_equal(as.vector(ranef), as.vector(fit1$ranef))

rm(ranef, yhat, fit4, fit3, seeds, fit2, fit1)

# NULL reads as not given for seed and sigest, as the NA defaults do: an
# unseeded NULL fit draws exactly what an unseeded default one does
set.seed(1L)
fitDefault <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.trees = 3L,
  n.threads = 1L,
  verbose = FALSE
)
set.seed(1L)
fitNull <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 5L,
  n.burn = 0L,
  n.chains = 1L,
  n.trees = 3L,
  n.threads = 1L,
  verbose = FALSE,
  seed = NULL,
  sigest = NULL
)
expect_identical(fitNull$yhat.train, fitDefault$yhat.train)
expect_identical(fitNull$sigest, fitDefault$sigest)
rm(fitDefault, fitNull)

rm(g, y, x)


rm(testData)
