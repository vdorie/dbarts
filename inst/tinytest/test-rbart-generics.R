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

# test that extract works at baseline
x <- testData$x
y <- testData$y
g <- factor(testData$g)
xNamed <- x
colnames(xNamed) <- paste0("x.", seq_len(ncol(x)))

set.seed(0L)
rbartFit <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE
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
  rbartFit$yhat.train + unname(rbartFit$ranef[,, as.character(g)])
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
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  combineChains = TRUE,
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
  y ~ x,
  group.by = g,
  test = xNamed,
  group.by.test = g,
  offset.test = 5,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(
  as.vector(rbartFit$yhat.train),
  as.vector(rbartFit$yhat.test) - 5
)
expect_equal(
  apply(
    rbartFit$yhat.train + unname(rbartFit$ranef[,, as.character(g)]),
    3L,
    mean
  ),
  fitted(rbartFit)
)
expect_equal(
  apply(
    rbartFit$yhat.test + unname(rbartFit$ranef[,, as.character(g)]),
    3L,
    mean
  ),
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
  y ~ x,
  group.by = g,
  n.samples = 14L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 25L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_equal(
  fitted(rbartFit.0),
  apply(predict(rbartFit.0, xNamed, g), 2L, mean)
)

set.seed(0L)
rbartFit.0 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  keepTrees = TRUE,
  combineChains = FALSE,
  verbose = FALSE
)
expect_equal(
  fitted(rbartFit.0),
  apply(predict(rbartFit.0, xNamed, g, combineChains = FALSE), 3L, mean)
)

set.seed(0L)
rbartFit.1 <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 25L,
  n.threads = 1L,
  keepTrees = TRUE,
  combineChains = TRUE,
  verbose = FALSE
)
expect_equal(
  fitted(rbartFit.1),
  apply(predict(rbartFit.1, xNamed, g), 2L, mean)
)
expect_equal(predict(rbartFit.0, xNamed, g), predict(rbartFit.1, xNamed, g))
expect_equal(
  predict(rbartFit.0, xNamed, g, combineChains = FALSE),
  predict(rbartFit.1, xNamed, g, combineChains = FALSE)
)

rm(rbartFit.1, rbartFit.0, g, y, x)

# sigma and nobs as on a bart fit: the posterior mean of the residual sd, 1 on
# a binary response, and the rows in the likelihood, those of positive weight
x <- testData$x
y <- testData$y
g <- factor(testData$g)
w <- rep(c(0, 1), c(3L, length(y) - 3L))
z <- as.numeric(y > median(y))
set.seed(0L)
weighted <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  weights = w,
  n.samples = 7L,
  n.burn = 0L,
  n.chains = 2L,
  n.trees = 10L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE
)
expect_identical(sigma(weighted), mean(weighted$sigma))
expect_identical(nobs(weighted), length(y) - 3L)
set.seed(0L)
binary <- dbarts::rbart_vi(
  z ~ x,
  group.by = g,
  weights = w,
  n.samples = 7L,
  n.burn = 0L,
  n.chains = 1L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(sigma(binary), 1)
expect_identical(nobs(binary), length(y) - 3L)
set.seed(0L)
unweighted <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 7L,
  n.burn = 0L,
  n.chains = 1L,
  n.trees = 10L,
  n.threads = 1L,
  verbose = FALSE
)
expect_identical(nobs(unweighted), length(y))
# an argument name the method does not take is warned of, in predict's words
for (generic in c("fitted", "residuals", "extract", "sigma", "nobs")) {
  expect_warning(
    match.fun(generic)(weighted, bogus = 1),
    paste0(
      "extra argument 'bogus' passed to ",
      generic,
      " on a rbart fit will be disregarded"
    ),
    fixed = TRUE,
    info = generic
  )
}
expect_warning(
  dbarts::extract(weighted, type = "trees", bogus = 1),
  "extra argument 'bogus' passed to extract on a rbart fit",
  fixed = TRUE
)
expect_error(
  dbarts::extract(weighted, type = "trees", sample = "train"),
  "'sample' is not used when type = \"trees\"",
  fixed = TRUE
)
expect_silent(dbarts::extract(weighted, type = "trees", chainNums = 1L))
rm(generic, weighted, binary, unweighted, w, z, g, y, x)


rm(testData)
