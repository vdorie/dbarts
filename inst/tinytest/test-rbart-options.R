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

# prior.scale alone (no explicit prior =) must resolve through exact-name
# matching on the matched call, not $'s unique-partial-match onto
# prior.scale itself (which used to send 'prior' to could-not-find-function
# "prior"); called out literally, no wrapper, so the promise is the real one
expect_inherits(
  dbarts::rbart_vi(
    y ~ x,
    testData,
    group.by = g,
    prior.scale = 1.5,
    n.samples = 2L,
    n.burn = 0L,
    n.thin = 2L,
    n.chains = 1L,
    n.trees = 3L,
    n.threads = 1L,
    keepTrainingFits = FALSE,
    verbose = FALSE
  ),
  "rbart"
)

# a binary fit with no k carries the chi(1.5, 2) node hyperprior: the
# formal default is NULL, not the inert 2.0 it used to be
yBinary <- testData$y > median(testData$y)
fitBinaryK <- dbarts::rbart_vi(
  yBinary ~ x,
  testData,
  group.by = g,
  n.samples = 2L,
  n.burn = 0L,
  n.thin = 2L,
  n.chains = 1L,
  n.trees = 3L,
  n.threads = 1L,
  keepTrainingFits = FALSE,
  verbose = FALSE
)
hyperprior <- fitBinaryK$fit[[1L]]$model@node.hyperprior
expect_inherits(hyperprior, "dbartsChiHyperprior")
expect_equal(hyperprior@degreesOfFreedom, 1.5)
expect_equal(hyperprior@scale, 2)
rm(yBinary, fitBinaryK, hyperprior)

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

# the multi-chain packaging reports them as well, in the chain-major shape
# (the single-chain branch above is a different assembly)
fit.dart2 <- dbarts::rbart_vi(
  y ~ x,
  testData,
  group.by = g,
  n.samples = 4L,
  n.burn = 2L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 3L,
  n.threads = 1L,
  verbose = FALSE,
  dart = TRUE,
  combineChains = FALSE
)
varprobs2 <- fit.dart2$varprobs
expect_equal(dim(varprobs2), c(2L, 4L, ncol(testData$x)))
# guarded so a dropped channel fails this assertion too rather than aborting
# the file where a NULL meets apply()
chainSumsAreOne <- tryCatch(
  isTRUE(all.equal(
    apply(varprobs2, c(1L, 2L), sum),
    matrix(1, 2L, 4L),
    check.attributes = FALSE
  )),
  error = function(e) FALSE
)
expect_true(chainSumsAreOne)
rm(fit.dart2, varprobs2, chainSumsAreOne)

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
