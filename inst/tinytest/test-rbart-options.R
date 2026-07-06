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

# test that rbart works with keepTrainingFits = FALSE
rbartFit <- dbarts::rbart_vi(
  y ~ x,
  testData,
  group.by = g,
  n.samples = 2L,
  n.burn = 0L,
  n.thin = 2L,
  n.chains = 2L,
  n.trees = 3L,
  n.threads = 1L,
  keepTrainingFits = FALSE,
  verbose = FALSE
)
expect_inherits(rbartFit, "rbart")
expect_true(is.null(rbartFit$yhat.train))
expect_true(is.null(rbartFit$yhat.train.mean))

rm(rbartFit)

# test that rbart works with k as a variable
# test thanks to Bruno Tancredi

k <- 1.8
expect_inherits(
  dbarts::rbart_vi(
    y ~ x,
    testData,
    group.by = g,
    n.samples = 2L,
    n.burn = 0L,
    n.thin = 2L,
    n.chains = 2L,
    n.trees = 3L,
    n.threads = 1L,
    k = k,
    keepTrainingFits = FALSE,
    verbose = FALSE
  ),
  "rbart"
)
rm(k)

# factors passes through to the data: categorical by default, indicators as
# the escape
df <- as.data.frame(testData$x)
names(df) <- paste0("x", seq_len(ncol(df)))
df$f <- factor(rep_len(c("a", "b", "c"), nrow(df)))
df$y <- testData$y + 0.5 * (df$f == "b")
df$g <- testData$g

fit.cat <- dbarts::rbart_vi(
  y ~ x1 + x2 + f,
  df,
  group.by = g,
  n.samples = 2L,
  n.burn = 0L,
  n.thin = 2L,
  n.chains = 1L,
  n.trees = 3L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(colnames(fit.cat$fit[[1L]]$data@x), c("x1", "x2", "f"))
expect_equal(fit.cat$fit[[1L]]$data@varTypes, c(0L, 0L, 1L))

fit.ind <- dbarts::rbart_vi(
  y ~ x1 + x2 + f,
  df,
  group.by = g,
  n.samples = 2L,
  n.burn = 0L,
  n.thin = 2L,
  n.chains = 1L,
  n.trees = 3L,
  n.threads = 1L,
  verbose = FALSE,
  factors = "indicators"
)
expect_true(all(fit.ind$fit[[1L]]$data@varTypes == 0L))
expect_true(ncol(fit.ind$fit[[1L]]$data@x) > 3L)

rm(fit.cat, fit.ind, df)

# dart runs inside the Gibbs sampler and reports its split probabilities
fit.dart <- dbarts::rbart_vi(
  y ~ x,
  testData,
  group.by = g,
  n.samples = 4L,
  n.burn = 2L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 3L,
  n.threads = 1L,
  verbose = FALSE,
  dart = TRUE
)
expect_equal(dim(fit.dart$varprobs), c(ncol(testData$x), 4L))
expect_equivalent(colSums(fit.dart$varprobs), rep(1, 4L))

expect_error(
  dbarts::rbart_vi(
    y ~ x,
    testData,
    group.by = g,
    n.samples = 2L,
    n.burn = 0L,
    n.thin = 1L,
    n.chains = 1L,
    n.trees = 3L,
    n.threads = 1L,
    verbose = FALSE,
    dart = TRUE,
    split.probs = c(0.5, 0.5)
  ),
  pattern = "cannot be combined"
)

rm(fit.dart)

rm(testData)
