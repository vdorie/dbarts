source(
  system.file("common", "rbartGroupData.R", package = "dbarts"),
  local = TRUE
)
x <- testData$x
y <- testData$y
g <- factor(g)

# the built-in priors run in-core: one multi-chain sampler with the chain
# count on the object, tau and ranef channels filled
# combineChains = FALSE pinned deliberately: the assertions below check the
# uncombined (n.chains x n.samples[ x n.groups]) stored shapes directly,
# now that rbart_vi's own stored-object default is combined
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
  seed = 2L,
  combineChains = FALSE
)
expect_equal(length(fit$fit), 1L)
expect_equal(fit$fit[[1L]]$control@n.chains, 2L)
expect_equal(fit$n.chains, 2L)
expect_equal(dim(fit$tau), c(2L, 6L))
expect_equal(dim(fit$first.tau), c(2L, 4L))
expect_true(all(is.finite(fit$tau)) && all(fit$tau > 0))
expect_equal(dimnames(fit$ranef)[[3L]], levels(g))

# the in-core sampler rbart_vi hands back is mutable on the response side at
# the pinned scale; re-anchoring is refused, since b and tau are held against
# the transform fixed at creation, and the whole-data conduit stays refused
# (test-grouped-swap.R carries the family-by-family cells)
expect_silent(fit$fit[[1L]]$setResponse(y))
expect_error(
  fit$fit[[1L]]$setResponse(y, updateScale = TRUE),
  pattern = "tau"
)
expect_error(
  fit$fit[[1L]]$setData(dbarts::dbartsData(y ~ x)),
  pattern = "grouped random effects"
)

# the grouped state (b, tau) rides save/load: a re-created sampler predicts
# identically
pred <- predict(fit, x, group.by = g, combineChains = FALSE)
tempFile <- tempfile(fileext = ".rds")
saveRDS(fit, tempFile)
fit.loaded <- readRDS(tempFile)
unlink(tempFile)
expect_identical(
  predict(fit.loaded, x, group.by = g, combineChains = FALSE),
  pred
)

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
expect_equal(fit.custom$n.chains, 2L)

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
expect_equal(fit.named$n.chains, 1L)
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
