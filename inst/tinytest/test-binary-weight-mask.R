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
  dbarts::dbartsData(x, y.binary, weights = a),
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

# the fit carries the mask its weights installed, and the log-likelihood
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

# --- the mask survives the sampler's re-creation --------------------------
# The mask is not saved state, so a sampler whose external pointer has died
# across a save and load is re-created from its data object - which for these
# families carries no weights at all. Were the mask not mirrored and
# re-applied, every masked row would silently rejoin the likelihood there. An
# inactive row's latent is not drawn, so a sweep that leaves the inactive
# latents exactly as they were is the mask still in force.
inactiveFrozen <- function(sampler) {
  before <- sampler$getLatents()
  invisible(sampler$run(0L, 1L))
  identical(before[a == 0], sampler$getLatents()[a == 0])
}
saved <- probitSampler(weights = a)
invisible(saved$run(20L, 5L))
saved$storeState()
kept <- saved$state
expect_identical(saved$activeRows, a)
expect_true(inactiveFrozen(saved))

restored <- unserialize(serialize(saved, NULL))
expect_false(.Call(dbarts:::C_dbarts_bartcore_isValidPointer, restored$pointer))
expect_identical(restored$activeRows, a)
expect_true(inactiveFrozen(restored))

# $setState re-creates the same way, and $copy builds a second engine
restated <- unserialize(serialize(saved, NULL))
restated$setState(kept)
expect_true(inactiveFrozen(restated))
expect_true(inactiveFrozen(saved$copy()))

# $setData clears the mask, and the record of it with it: a sampler swapped
# onto weightless data must not have one re-applied at its next re-creation
swapped <- probitSampler(weights = a)
swapped$setData(dbarts::dbartsData(x, y.binary))
expect_null(swapped$activeRows)
rm(swapped)

# the detector is not vacuous: clearing the mask unfreezes those latents
unmasked <- probitSampler(weights = a)
invisible(unmasked$run(20L, 5L))
unmasked$setActiveRows(NULL)
expect_null(unmasked$activeRows)
expect_false(inactiveFrozen(unmasked))
rm(saved, kept, restored, restated, unmasked, inactiveFrozen)

# --- the zero-weight warning ----------------------------------------------
# a vector of nothing but 0s and 1s states which rows are in the data set, so
# the zeros are not an inert value to warn about - they are the whole point,
# and on these families they are the mask outright
expect_identical(
  countWarnings(dbarts::dbartsData(x, y.binary, weights = a)),
  0L
)
expect_identical(countWarnings(probitSampler(weights = a)), 0L)
expect_identical(countWarnings(ordinalSampler(weights = a)), 0L)
# a real weight vector with an inert zero among it still warns
expect_warning(
  dbarts::dbartsData(x, y.binary, weights = replace(runif(n, 0.5, 1.5), 1L, 0)),
  "'weights' of 0 will be ignored"
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
  countWarnings,
  fit,
  reloaded,
  predicted,
  loglik
)
