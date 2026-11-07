# The 0/1 case-weight rule for the latent families. probit and ordinal have no
# weighted latent likelihood, but a weight vector of 0s and 1s does not ask for
# one: it says which rows are in the data set. Such a vector stays in the data
# object's weights slot at every entry point, and the engine composes it with
# the active-row mask, so it draws what the same mask draws while activeRows
# holds only what setActiveRows installed. All-ones is no weights at all; any
# other value is a weighted latent likelihood and stays refused.

set.seed(20260912L)
n <- 200L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y.binary <- as.double(x[, 1L] + rnorm(n, 0, 0.3) > 0.5)
y.ordinal <- as.double(1L + (seq_len(n) %% 3L))
a <- as.double(seq_len(n) %% 4L != 1L)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  # dbartsSpec() leaves the sample count for the sampler to carry, unlike
  # dbarts(), which fills it in; every run below names its own counts anyway
  n.samples = 10L,
  updateState = FALSE,
  seed = 7L
)

probitSampler <- function(...) {
  dbarts::dbarts(
    x,
    y.binary,
    family = "probit",
    control = control,
    ...,
    n.trees = 20L
  )
}
countWarnings <- function(expr) {
  count <- 0L
  withCallingHandlers(
    expr,
    warning = function(w) {
      count <<- count + 1L
      invokeRestart("muffleWarning")
    }
  )
  count
}
ordinalSampler <- function(...) {
  dbarts::dbarts(
    x,
    y.ordinal,
    family = "ordinal",
    control = control,
    ...,
    n.trees = 20L
  )
}

# --- creation: the weights draw as the mask, bitwise ----------------------
# the weighted arm and the hand-masked arm are the same sampler, so their
# draws agree draw for draw and not merely in distribution
weighted <- probitSampler(weights = a)
expect_identical(weighted$data@weights, a)
expect_null(weighted$activeRows)
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
expect_identical(weighted.ord$data@weights, a)
expect_null(weighted.ord$activeRows)
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
# a vector carrying zeros reports them before the refusal; the zeros are
# part of what that case checks
expectRefusal <- function(expr, pattern, zeros) {
  if (zeros) {
    expect_warning(
      expect_error(expr, pattern),
      "rows with 'weights' of 0 are left out of the likelihood",
      strict = TRUE
    )
  } else {
    expect_error(expr, pattern)
  }
}
for (bad in list(rep(2, n), rep(0.5, n), replace(a, 1L, 2))) {
  expectRefusal(
    probitSampler(weights = bad),
    "probit models do not support weights other than 0 and 1",
    any(bad == 0)
  )
  expectRefusal(
    ordinalSampler(weights = bad),
    "ordinal models do not support weights other than 0 and 1",
    any(bad == 0)
  )
}
rm(bad, expectRefusal)

# --- the mutators ---------------------------------------------------------
# $setWeights takes the 0/1 vector into the data object's weights slot rather
# than refusing, and installs no mask
bySetWeights <- probitSampler()
bySetWeights$setWeights(a)
expect_identical(bySetWeights$data@weights, a)
expect_null(bySetWeights$activeRows)
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
# weights ride its data object, at the replacement's own n
bySetData <- probitSampler()
bySetData$setData(dbarts::dbartsData(x, y.binary, weights = a))
expect_identical(bySetData$data@weights, a)
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
# dbartsSpec() resolves the same rule and keeps the weights on its data
# object, so a sampler built from the spec draws as the mask with nothing to
# install
spec <- dbarts::dbartsSpec(
  dbarts::dbartsData(x, y.binary, weights = a),
  control = control,
  family = "probit",
  n.trees = 20L
)
expect_null(spec$active)
expect_false("active" %in% names(spec))
expect_identical(spec$data@weights, a)
specSampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
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
fit <- dbarts::bart(
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
)
expect_identical(ncol(fit$yhat.train), n)
expect_true(all(is.finite(fit$yhat.train)))
predicted <- predict(fit, x[1:5, , drop = FALSE])
# the store-then-save flow every kept sampler takes; predict reads the stored
# forests rather than the training rows, so a reloaded masked fit predicts
# what the live one does
fit$fit$storeState()
reloaded <- unserialize(serialize(fit, NULL))
expect_identical(predict(reloaded, x[1:5, , drop = FALSE]), predicted)
expect_identical(ncol(predicted), 5L)
expect_true(all(is.finite(predicted)))

# the fit records the rows its weights leave out, and the log-likelihood
# channel reports NaN at a masked row - the engine's own convention for a row
# that is not in the model - rather than the finite value its fit would give
expect_identical(fit$active, a)
loglik <- dbarts::extract(fit, "loglik")
expect_identical(ncol(loglik), n)
expect_true(all(is.nan(loglik[, a == 0])))
expect_true(all(is.finite(loglik[, a == 1])))

# the run itself must carry the mask, not merely the sampler bart() built:
# substituting arbitrary labels at the INACTIVE rows leaves every active row's
# draw bitwise. An ordinal fit reaches its engine by a second route - the run
# adopts an engine built separately from the data object, which carries no
# weights at all - so it is pinned here too.
frontDoor <- function(y, family, weights) {
  dbarts::bart(
    x,
    y,
    family = family,
    weights = weights,
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 10L,
    n.chains = 1L,
    n.threads = 1L,
    seed = 11L,
    verbose = FALSE
  )$yhat.train
}
expect_identical(
  frontDoor(y.binary, "probit", a)[, a == 1],
  frontDoor(ifelse(a == 0, 1 - y.binary, y.binary), "probit", a)[, a == 1]
)
expect_false(identical(
  frontDoor(y.binary, "probit", a),
  frontDoor(y.binary, "probit", rep(1, n))
))
y.ordinal.flipped <- y.ordinal
y.ordinal.flipped[a == 0] <- (y.ordinal.flipped[a == 0] %% 3) + 1
expect_identical(
  frontDoor(y.ordinal, "ordinal", a)[, a == 1, , drop = FALSE],
  frontDoor(y.ordinal.flipped, "ordinal", a)[, a == 1, , drop = FALSE]
)
expect_false(identical(
  frontDoor(y.ordinal, "ordinal", a),
  frontDoor(y.ordinal, "ordinal", rep(1, n))
))
rm(frontDoor, y.ordinal.flipped)

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
expect_silent(dbarts::bartBT(
  x,
  y.binary,
  weights = a,
  ndpost = 5L,
  nskip = 5L,
  ntree = 10L,
  nchain = 1L,
  nthread = 1L,
  verbose = FALSE
))

# cross-validation takes the vector every other entry point takes, each fold
# carrying its rows' weights (test-xbart-fold-oracle.R pins the folds)
expect_silent(dbarts::xbart(
  x,
  y.binary,
  weights = a,
  family = "probit",
  n.samples = 5L,
  n.burn = c(5L, 5L),
  n.trees = 10L,
  n.reps = 1L,
  n.threads = 1L
))

# --- the weights survive the sampler's re-creation ------------------------
# The rows out are not saved state, so a sampler whose external pointer has
# died across a save and load is re-created from its data object, which
# carries the 0/1 weights. Were they lost there, every row of weight 0 would
# silently rejoin the likelihood. An inactive row's latent is not drawn, so a
# sweep that leaves the inactive latents exactly as they were is the
# membership still in force.
inactiveFrozen <- function(sampler) {
  before <- sampler$getLatents()
  invisible(sampler$run(0L, 1L))
  identical(before[a == 0], sampler$getLatents()[a == 0])
}
saved <- probitSampler(weights = a)
invisible(saved$run(20L, 5L))
saved$storeState()
kept <- saved$state
expect_null(saved$activeRows)
expect_true(inactiveFrozen(saved))

restored <- unserialize(serialize(saved, NULL))
expect_false(.Call(dbarts:::C_dbarts_bartcore_isValidPointer, restored$pointer))
expect_identical(restored$data@weights, a)
expect_true(inactiveFrozen(restored))

# $setState re-creates the same way, and $copy builds a second engine
restated <- unserialize(serialize(saved, NULL))
restated$setState(kept)
expect_true(inactiveFrozen(restated))
expect_true(inactiveFrozen(saved$copy()))

# $setData onto weightless data takes every row in: nothing of the old
# weights is re-applied at its next re-creation
swapped <- probitSampler(weights = a)
swapped$setData(dbarts::dbartsData(x, y.binary))
expect_null(swapped$activeRows)
expect_null(swapped$data@weights)
rm(swapped)

# the detector is not vacuous: weights of 1 unfreeze those latents
unmasked <- probitSampler(weights = a)
invisible(unmasked$run(20L, 5L))
unmasked$setWeights(rep(1, n))
expect_null(unmasked$activeRows)
expect_false(inactiveFrozen(unmasked))
rm(saved, kept, restored, restated, unmasked, inactiveFrozen)

# --- the zero-weight warning ----------------------------------------------
# a vector of nothing but 0s and 1s states which rows are in the data set, so
# the zeros are not an inert value to warn about - they are the whole point,
# and on these families they draw as the mask
expect_identical(
  countWarnings(dbarts::dbartsData(x, y.binary, weights = a)),
  0L
)
expect_identical(countWarnings(probitSampler(weights = a)), 0L)
expect_identical(countWarnings(ordinalSampler(weights = a)), 0L)
# a real weight vector with an inert zero among it still warns
expect_warning(
  dbarts::dbartsData(x, y.binary, weights = replace(runif(n, 0.5, 1.5), 1L, 0)),
  "rows with 'weights' of 0 are left out of the likelihood but stay in the fit; use 'subset' to drop them",
  fixed = TRUE
)

# an ordinal bart fit records the rows out too, and its log-likelihood is NaN
# at those rows, as the probit fit's is
fit.ordinal <- dbarts::bart(
  x,
  factor(y.ordinal, ordered = TRUE),
  weights = a,
  n.trees = 10L,
  n.samples = 10L,
  n.burn = 10L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_equal(as.double(fit.ordinal$active), a)
loglik.ordinal <- dbarts::extract(fit.ordinal, "loglik")
expect_true(all(is.nan(loglik.ordinal[, a == 0])))
expect_true(all(is.finite(loglik.ordinal[, a == 1])))
rm(fit.ordinal, loglik.ordinal)

rm(
  n,
  x,
  y.binary,
  y.ordinal,
  a,
  control,
  probitSampler,
  ordinalSampler,
  countWarnings,
  fit,
  reloaded,
  predicted,
  loglik
)
