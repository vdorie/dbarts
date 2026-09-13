# The 0/1 case-weight rule for the latent families. probit and ordinal carry no
# weight channel, but a weight vector of 0s and 1s does not ask for one: it
# says which rows are in the data set, which is what the active-row mask says,
# so such a vector installs as the mask at every entry point and the weights
# slot stays empty. All-ones is no mask at all, as it is no weights; any other
# value is a weighted latent likelihood and stays refused.

set.seed(20260912L)
n <- 200L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y.binary <- as.double(x[, 1L] + rnorm(n, 0, 0.3) > 0.5)
y.ordinal <- as.double(1L + (seq_len(n) %% 3L))
a <- as.double(seq_len(n) %% 4L != 1L)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  # dbartsSpec() leaves the sample count for the sampler to carry, unlike
  # dbarts(), which fills it in; every run below names its own counts anyway
  n.samples = 10L,
  updateState = FALSE,
  seed = 7L
)

probitSampler <- function(...) {
  dbarts::dbarts(x, y.binary, family = "probit", control = control, ...)
}
ordinalSampler <- function(...) {
  dbarts::dbarts(x, y.ordinal, family = "ordinal", control = control, ...)
}

# --- creation: the weights ARE the mask, bitwise -------------------------
# the weighted arm and the hand-masked arm are the same sampler, so their
# draws agree draw for draw and not merely in distribution
weighted <- probitSampler(weights = a)
expect_true(is.null(weighted$data@weights))
byMask <- probitSampler()
byMask$setActiveRows(a)
draws.weighted <- weighted$run(20L, 10L)
draws.byMask <- byMask$run(20L, 10L)
expect_identical(draws.weighted$train, draws.byMask$train)
# every row keeps its fitted value, the inactive ones included: the mask takes
# a row out of the likelihood, not out of the reported channels
expect_identical(nrow(draws.weighted$train), n)
expect_true(all(is.finite(draws.weighted$train)))
rm(weighted, byMask, draws.weighted, draws.byMask)

# ordinal reaches the same channel on the same terms
weighted.ord <- ordinalSampler(weights = a)
expect_true(is.null(weighted.ord$data@weights))
byMask.ord <- ordinalSampler()
byMask.ord$setActiveRows(a)
expect_identical(
  weighted.ord$run(20L, 10L)$train,
  byMask.ord$run(20L, 10L)$train
)
rm(weighted.ord, byMask.ord)

# --- all ones is still no weights at all, bitwise ------------------------
ones <- probitSampler(weights = rep(1, n))
plain <- probitSampler()
expect_identical(ones$run(20L, 10L)$train, plain$run(20L, 10L)$train)
ones.ord <- ordinalSampler(weights = rep(1, n))
plain.ord <- ordinalSampler()
expect_identical(ones.ord$run(20L, 10L)$train, plain.ord$run(20L, 10L)$train)
rm(ones, plain, ones.ord, plain.ord)

# --- any other value keeps the refusal ------------------------------------
for (bad in list(rep(2, n), rep(0.5, n), replace(a, 1L, 2))) {
  expect_error(
    probitSampler(weights = bad),
    "probit models do not support weights other than 0 and 1"
  )
  expect_error(
    ordinalSampler(weights = bad),
    "ordinal models do not support weights other than 0 and 1"
  )
}
rm(bad)

# --- the mutators ---------------------------------------------------------
# $setWeights installs the mask rather than refusing, and leaves the data
# object's weights slot empty
bySetWeights <- probitSampler()
bySetWeights$setWeights(a)
expect_true(is.null(bySetWeights$data@weights))
byMask <- probitSampler()
byMask$setActiveRows(a)
expect_identical(
  bySetWeights$run(20L, 10L)$train,
  byMask$run(20L, 10L)$train
)
rm(bySetWeights, byMask)

bySetWeights.ord <- ordinalSampler()
bySetWeights.ord$setWeights(a)
byMask.ord <- ordinalSampler()
byMask.ord$setActiveRows(a)
expect_identical(
  bySetWeights.ord$run(20L, 10L)$train,
  byMask.ord$run(20L, 10L)$train
)
rm(bySetWeights.ord, byMask.ord)

expect_error(
  probitSampler()$setWeights(rep(2, n)),
  "probit models do not support case weights other than 0 and 1"
)
expect_error(
  ordinalSampler()$setWeights(rep(0.5, n)),
  "ordinal models do not support case weights other than 0 and 1"
)

# $setData carries the same rule on the whole-data conduit: the replacement's
# weights install as the mask the replacement's own n sizes
bySetData <- probitSampler()
bySetData$setData(dbarts::dbartsData(x, y.binary, weights = a))
expect_true(is.null(bySetData$data@weights))
byMask <- probitSampler()
byMask$setData(dbarts::dbartsData(x, y.binary))
byMask$setActiveRows(a)
expect_identical(
  bySetData$run(20L, 10L)$train,
  byMask$run(20L, 10L)$train
)
expect_error(
  bySetData$setData(dbarts::dbartsData(x, y.binary, weights = rep(2, n))),
  "probit models do not support case weights other than 0 and 1"
)
rm(bySetData, byMask)

# --- the specification surface --------------------------------------------
# dbartsSpec() resolves the same rule but builds no sampler, so it hands the
# mask back on its own element for the caller to install
spec <- dbarts::dbartsSpec(
  suppressWarnings(dbarts::dbartsData(x, y.binary, weights = a)),
  control = control,
  family = "probit"
)
expect_identical(spec$active, a)
expect_true(is.null(spec$data@weights))
specSampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
specSampler$setActiveRows(spec$active)
byMask <- probitSampler()
byMask$setActiveRows(a)
expect_identical(
  specSampler$run(20L, 10L)$train,
  byMask$run(20L, 10L)$train
)
rm(spec, specSampler, byMask)

# --- the front doors ------------------------------------------------------
# bart() and bartBT() reach the rule through dbarts(); a masked fit reports
# every row, and survives a save/load round trip well enough to predict
fit <- suppressWarnings(dbarts::bart(
  x,
  y.binary,
  weights = a,
  n.samples = 5L,
  n.burn = 5L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 11L,
  verbose = FALSE,
  keepTrees = TRUE
))
expect_identical(ncol(fit$yhat.train), n)
expect_true(all(is.finite(fit$yhat.train)))
predicted <- predict(fit, x[1:5, , drop = FALSE])
# the store-then-save flow every kept sampler takes; the mask itself is not
# saved state, and predict reads the stored forests rather than the training
# rows, so a reloaded masked fit predicts what the live one does
fit$fit$storeState()
reloaded <- unserialize(serialize(fit, NULL))
expect_identical(predict(reloaded, x[1:5, , drop = FALSE]), predicted)
expect_identical(ncol(predicted), 5L)
expect_true(all(is.finite(predicted)))

expect_error(
  dbarts::bart(
    x,
    y.binary,
    weights = rep(0.5, n),
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  "probit models do not support weights other than 0 and 1"
)
expect_error(
  dbarts::bartBT(
    x,
    y.binary,
    weights = rep(2, n),
    ndpost = 5L,
    nskip = 5L,
    ntree = 10L,
    nchain = 1L,
    nthread = 1L,
    verbose = FALSE
  ),
  "probit models do not support weights other than 0 and 1"
)
expect_silent(suppressWarnings(dbarts::bartBT(
  x,
  y.binary,
  weights = a,
  ndpost = 5L,
  nskip = 5L,
  ntree = 10L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
)))

# cross-validation partitions the rows itself and has no sampler to install a
# mask on, so it refuses the vector every other entry point takes
expect_error(
  dbarts::xbart(
    x,
    y.binary,
    weights = a,
    family = "probit",
    n.reps = 1L,
    n.threads = 1L
  ),
  "xbart does not accept weights of 0 and 1"
)

rm(
  n,
  x,
  y.binary,
  y.ordinal,
  a,
  control,
  probitSampler,
  ordinalSampler,
  fit,
  reloaded,
  predicted
)
