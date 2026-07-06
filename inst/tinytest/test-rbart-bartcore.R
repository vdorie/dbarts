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
x <- testData$x
y <- testData$y
g <- factor(g)

# the built-in priors run in-core: one multi-chain sampler with the chain
# count on the object, tau and ranef channels filled
fit <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 6L,
  n.burn = 4L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 5L,
  n.threads = 1L,
  keepTrees = TRUE,
  verbose = FALSE,
  seed = 2L
)
expect_equal(length(fit$fit), 1L)
expect_equal(fit$fit[[1L]]$control@n.chains, 2L)
expect_equal(fit$n.chains, 2L)
expect_equal(dim(fit$tau), c(2L, 6L))
expect_equal(dim(fit$first.tau), c(2L, 4L))
expect_true(all(is.finite(fit$tau)) && all(fit$tau > 0))
expect_equal(dimnames(fit$ranef)[[3L]], levels(g))

# grouped samplers fix the response and data at creation
expect_error(fit$fit[[1L]]$setResponse(y), pattern = "grouped random effects")
expect_error(
  fit$fit[[1L]]$setData(dbarts::dbartsData(y ~ x)),
  pattern = "grouped random effects"
)

# the grouped state (b, tau) rides save/load: a re-created sampler predicts
# identically
pred <- predict(fit, x, g, combineChains = FALSE)
tempFile <- tempfile(fileext = ".rds")
saveRDS(fit, tempFile)
fit.loaded <- readRDS(tempFile)
unlink(tempFile)
expect_identical(predict(fit.loaded, x, g, combineChains = FALSE), pred)

# a custom prior function keeps the R loop and its per-chain samplers
fit.custom <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  prior = function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE),
  n.samples = 3L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(length(fit.custom$fit), 2L)
expect_true(is.null(fit.custom$n.chains))

# the gamma prior resolves by name
fit.gamma <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  prior = gamma,
  n.samples = 3L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(all(is.finite(fit.gamma$tau)))

# a custom prior passed by NAME is a symbol but not a builtin: it must
# route through the R loop, not the builtin lookup (which used to crash
# on any non-builtin symbol)
flatTauPrior <- function(x, rel.scale) dcauchy(x, 0, rel.scale * 2.5, TRUE)
fit.named <- dbarts::rbart_vi(
  y ~ x,
  group.by = g,
  prior = flatTauPrior,
  n.samples = 3L,
  n.burn = 0L,
  n.thin = 1L,
  n.chains = 1L,
  n.trees = 5L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(is.null(fit.named$n.chains))
expect_true(all(is.finite(fit.named$tau)))

rm(
  fit.gamma,
  fit.custom,
  fit.named,
  flatTauPrior,
  fit.loaded,
  pred,
  fit,
  g,
  y,
  x,
  b,
  sigma.b
)
rm(testData)
