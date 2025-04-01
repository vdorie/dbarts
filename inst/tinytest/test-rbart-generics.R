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

# test that extract works at baseline
x <- testData$x
y <- testData$y
g <- factor(testData$g)

set.seed(0L)
rbartFit <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 25L, n.threads = 1L, verbose = FALSE
)
expect_equal(
  dbarts::extract(rbartFit, type = "bart", combineChains = FALSE),
  rbartFit$yhat.train
)
expect_equal(
  dbarts::extract(rbartFit, type = "ranef", combineChains = FALSE),
  rbartFit$ranef
)
expect_equal(
  dbarts::extract(rbartFit, type = "ev", combineChains = FALSE),
  rbartFit$yhat.train + unname(rbartFit$ranef[,,as.character(g)])
)

ppd <- dbarts::extract(rbartFit, type = "ppd", combineChains = FALSE)
sigma.hat <- apply(
  ppd - dbarts::extract(rbartFit, type = "ev", combineChains = FALSE),
  c(1L, 2L),
  sd
)
# silly test with 7 samples but we need something
expect_true(
  cor(as.vector(sigma.hat), as.vector(rbartFit$sigma)) >= 0.85
) 

set.seed(0L)
rbartFit.2 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 25L, n.threads = 1L, combineChains = TRUE,
  verbose = FALSE
)
expect_equal(
  dbarts::extract(rbartFit, type = "ev", combineChains = TRUE),
  dbarts::extract(rbartFit.2, type = "ev", combineChains = TRUE)
)
expect_equal(
  dbarts::extract(rbartFit, type = "ev", combineChains = FALSE),
  dbarts::extract(rbartFit.2, type = "ev", combineChains = FALSE)
)

rm(sigma.hat, ppd, rbartFit, g, y, x)



# test that fitted works correctly
x <- testData$x
y <- testData$y
g <- factor(testData$g)


rbartFit <- dbarts::rbart_vi(
  y ~ x, group.by = g, test = x, group.by.test = g, offset.test = 5,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 25L, n.threads = 1L, verbose = FALSE
)
expect_equal(
  as.vector(rbartFit$yhat.train),
  as.vector(rbartFit$yhat.test) - 5
)
expect_equal(
  apply(rbartFit$yhat.train + unname(rbartFit$ranef[,,as.character(g)]), 3L, mean),
  fitted(rbartFit)
)
expect_equal(
  apply(rbartFit$yhat.test  + unname(rbartFit$ranef[,,as.character(g)]), 3L, mean),
  fitted(rbartFit, sample = "test")
)
expect_equal(fitted(rbartFit, type = "ranef"), rbartFit$ranef.mean)

rm(rbartFit, g, y, x)


# test that predict matches fitted
x <- testData$x
y <- testData$y
g <- factor(testData$g)

set.seed(0L)
rbartFit.0 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 14L, n.burn = 0L, n.thin = 1L, n.chains = 1L,
  n.trees = 25L, n.threads = 1L, keepTrees = TRUE,
  verbose = FALSE
)
expect_equal(fitted(rbartFit.0), apply(predict(rbartFit.0, x, g), 2L, mean))

set.seed(0L)
rbartFit.0 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
  n.trees = 25L, n.threads = 1L, keepTrees = TRUE, combineChains = FALSE,
  verbose = FALSE
)
expect_equal(
  fitted(rbartFit.0),
  apply(predict(rbartFit.0, x, g, combineChains = FALSE), 3L, mean)
)

set.seed(0L)
rbartFit.1 <- dbarts::rbart_vi(
  y ~ x, group.by = g,
  n.samples = 7L, n.burn = 0L, n.thin = 1L, n.chains = 2L,
 n.trees = 25L, n.threads = 1L, keepTrees = TRUE, combineChains = TRUE,
 verbose = FALSE
)
expect_equal(
  fitted(rbartFit.1),
  apply(predict(rbartFit.1, x, g), 2L, mean)
)
expect_equal(predict(rbartFit.0, x, g), predict(rbartFit.1, x, g))
expect_equal(
  predict(rbartFit.0, x, g, combineChains = FALSE),
  predict(rbartFit.1, x, g, combineChains = FALSE)
)

rm(rbartFit.1, rbartFit.0, g, y, x)


rm(testData)

