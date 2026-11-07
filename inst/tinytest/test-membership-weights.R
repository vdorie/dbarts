# Weights of 0 and 1 as membership: on probit and ordinal, a 0/1 vector rides
# the data object and the engine composes it with the active-row mask, a row
# in when its weight is 1 and the mask keeps it; on logistic a 0 is the mask
# at that row; on nbinom 0/1 are counts, a 0 out under the same composition.
# Each is pinned draw for draw against the mask set by setActiveRows.

set.seed(20261107L)
n <- 120L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y.binary <- as.double(x[, 1L] + rnorm(n, 0, 0.3) > 0.5)
y.ordinal <- as.double(1L + (seq_len(n) %% 3L))
y.count <- as.double(rpois(n, exp(1 + x[, 1L])))
a <- as.double(seq_len(n) %% 4L != 1L)
b <- as.double(seq_len(n) %% 5L != 2L)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 10L,
  updateState = FALSE,
  seed = 7L
)
sampler <- function(y, family, ...) {
  dbarts::dbarts(x, y, family = family, control = control, ..., n.trees = 20L)
}
trainDraws <- function(s) s$run(20L, 10L)$train
# the zero-weight warning a vector of counts with zeros among them raises,
# pinned once below and muffled where it is not the point
muffleZeroWeights <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (
        grepl("rows with 'weights' of 0", conditionMessage(w), fixed = TRUE)
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

for (family in c("probit", "ordinal")) {
  y <- if (family == "probit") y.binary else y.ordinal

  # --- setData(sampler$data) keeps the rows of weight 0 out ---------------
  # the replacement's weights ride its data object, so the swapped sampler
  # draws what one swapped onto weightless data and masked by hand draws
  swapped <- sampler(y, family, weights = a)
  swapped$setData(swapped$data)
  expect_identical(swapped$data@weights, a, info = family)
  expect_null(swapped$activeRows, info = family)
  byMask <- sampler(y, family)
  byMask$setData(dbarts::dbartsData(x, y))
  byMask$setActiveRows(a)
  expect_identical(trainDraws(swapped), trainDraws(byMask), info = family)
  rm(swapped, byMask)

  # --- the mask composes with the weights in either order ------------------
  composed <- sampler(y, family)
  composed$setActiveRows(a * b)
  reference <- trainDraws(composed)
  weightsFirst <- sampler(y, family, weights = a)
  weightsFirst$setActiveRows(b)
  expect_identical(weightsFirst$activeRows, b, info = family)
  expect_identical(trainDraws(weightsFirst), reference, info = family)
  maskFirst <- sampler(y, family)
  maskFirst$setActiveRows(b)
  maskFirst$setWeights(a)
  expect_identical(trainDraws(maskFirst), reference, info = family)
  rm(composed, weightsFirst, maskFirst, reference)

  # --- a composition that reactivates a row draws nothing there ------------
  # weights of 1 bring the rows of weight 0 back as clearing the mask does,
  # leaving every latent as it was, and the sweeps after draw alike
  byWeights <- sampler(y, family, weights = a)
  byMask <- sampler(y, family)
  byMask$setActiveRows(a)
  invisible(byWeights$run(10L, 5L))
  invisible(byMask$run(10L, 5L))
  before <- byWeights$getLatents()
  byWeights$setWeights(rep(1, n))
  byMask$setActiveRows(NULL)
  expect_identical(byWeights$getLatents(), before, info = family)
  expect_identical(trainDraws(byWeights), trainDraws(byMask), info = family)
  rm(byWeights, byMask, before)
}
rm(family, y)

# --- logistic: a 0 is the mask at that row --------------------------------
# the counts elsewhere as they are; the row of weight 0 consumes no variate,
# so the fit draws what the fit with that row masked draws
counts <- as.double(1L + (seq_len(n) %% 3L))
zeroed <- replace(counts, a == 0, 0)
# among other counts the zeros are reported, as on every family
expect_warning(
  dbarts::dbartsData(x, y.binary, weights = zeroed),
  "rows with 'weights' of 0 are left out of the likelihood",
  strict = TRUE
)
masked <- sampler(y.binary, "logistic", weights = counts)
masked$setActiveRows(a)
zero <- muffleZeroWeights(sampler(y.binary, "logistic", weights = zeroed))
expect_identical(zero$data@weights, zeroed)
expect_identical(trainDraws(zero), trainDraws(masked))
rm(zero, masked)

# the same at setWeights, bitwise at the rows in: a row out keeps the omega
# it had, and its fitted value, which no row in reads, follows its working
# response, stated at count 1 there and at its count under the mask
bySetWeights <- sampler(y.binary, "logistic", weights = counts)
bySetWeights$setWeights(zeroed)
masked <- sampler(y.binary, "logistic", weights = counts)
masked$setActiveRows(a)
draws.bySetWeights <- trainDraws(bySetWeights)
draws.masked <- trainDraws(masked)
expect_identical(draws.bySetWeights[a == 1, ], draws.masked[a == 1, ])
expect_equal(draws.bySetWeights, draws.masked, tolerance = 1e-12)
rm(bySetWeights, masked, draws.bySetWeights, draws.masked)

# and the log-likelihood channel is NaN at the row

fit <- muffleZeroWeights(dbarts::bart(
  x,
  y.binary,
  family = "logistic",
  weights = zeroed,
  n.samples = 5L,
  n.burn = 5L,
  n.trees = 10L,
  n.chains = 1L,
  n.threads = 1L,
  seed = 11L,
  verbose = FALSE
))
loglik <- dbarts::extract(fit, "loglik")
expect_true(all(is.nan(loglik[, zeroed == 0])))
expect_true(all(is.finite(loglik[, zeroed > 0])))
rm(fit, loglik)

# a fractional count is refused in the family's words, a negative one by the
# data object first
expect_error(
  sampler(y.binary, "logistic", weights = replace(counts, 1L, 0.5)),
  "logistic weights are observation counts and must be non-negative whole numbers"
)
expect_error(
  sampler(y.binary, "logistic", weights = replace(counts, 1L, -1)),
  "'weights' must all be non-negative"
)
rm(counts, zeroed)

# --- nbinom: 0/1 weights are counts, bitwise the mask ----------------------
# a sampler created from 0/1 weights draws what one created without them and
# masked at the same rows, the center taken again from the rows in, draws
nbinomPair <- function(y, weights) {
  spec <- dbarts::dbartsSpec(
    dbarts::dbartsData(x, y, weights = weights),
    control = control,
    family = "nbinom",
    n.trees = 20L
  )
  weighted <- new("dbartsSampler", spec$control, spec$model, spec$data)
  plain <- dbarts::dbartsSpec(
    dbarts::dbartsData(x, y),
    control = control,
    family = "nbinom",
    n.trees = 20L
  )
  masked <- new("dbartsSampler", plain$control, plain$model, plain$data)
  masked$setActiveRows(weights, updateScale = TRUE)
  list(weighted = weighted, masked = masked)
}
pair <- nbinomPair(y.count, a)
expect_identical(pair$weighted$data@weights, a)
expect_null(pair$weighted$activeRows)
expect_identical(trainDraws(pair$weighted), trainDraws(pair$masked))
rm(pair)

# where the rows in hold one count the center is read from every row on both
# routes, a row of weight 0 as one of weight 1
y.flat <- replace(y.count, a == 1, 3)
expect_warning(
  pair <- nbinomPair(y.flat, a),
  "fewer than two distinct values of the response are in the likelihood"
)
expect_identical(trainDraws(pair$weighted), trainDraws(pair$masked))
rm(pair, y.flat, nbinomPair)

rm(
  n,
  x,
  y.binary,
  y.ordinal,
  y.count,
  a,
  b,
  control,
  sampler,
  trainDraws,
  muffleZeroWeights
)
