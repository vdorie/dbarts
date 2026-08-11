source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

# test that dbarts sampler settors raise errors
train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
test <- data.frame(x = testData$x, z = 1 - testData$z)

control <- dbarts::dbartsControl(
  n.burn = 0L,
  n.samples = 1L,
  n.thin = 5L,
  n.chains = 1L,
  n.threads = 1L,
  updateState = FALSE,
  verbose = FALSE
)
sampler <- dbarts::dbarts(y ~ x + z, train, test, control = control)

expect_error(
  sampler$setControl("not-a-control"),
  "'control' must inherit from dbartsControl"
)
expect_error(
  sampler$setResponse(numeric(0L)),
  paste0("y must be of length equal to ", nrow(train))
)

expect_error(
  sampler$setOffset(numeric(0L)),
  "length of replacement offset is not equal to number of observations"
)
expect_error(
  sampler$setPredictor(numeric(0L), 1L),
  "length of new x does not match y"
)
expect_error(
  sampler$setPredictor(testData$z, 3L),
  "column '3' is out of range"
)
expect_error(
  sampler$setTestPredictor(numeric(0L), 1L),
  "length of new x does not match old x.test"
)
expect_error(
  sampler$setTestPredictor(numeric(0L)),
  "number of columns in 'test' must be equal to that of 'x'"
)
expect_error(
  sampler$setTestPredictor(testData$z, 3L),
  "column '3' is out of range"
)

n <- length(testData$y)
expect_error(
  sampler$setPredictor(matrix(numeric(n * 3L), n)),
  paste0("dimension of x must be equal to ", ncol(sampler$data@x))
)
expect_error(
  sampler$setPredictor(matrix(numeric((n - 1L) * 2L), n - 1L)),
  paste0("dimension of x must be equal to ", nrow(sampler$data@x))
)
expect_error(
  sampler$setPredictor(matrix(numeric(n * 2L), n), 1L),
  "number of columns of new x does not match length of columns to replace"
)

# setResponse/setWeights/setOffset/setSigma reject the same class of bad
# input dbartsData's construction/validity reject, and do so before any
# sampler state is touched (the rejected value never reaches data@... or
# the engine)

yBad <- testData$y
yBad[1L] <- NA
yBefore <- sampler$data@y
expect_error(sampler$setResponse(yBad), "response contains missing values")
expect_identical(sampler$data@y, yBefore)

weightsBefore <- sampler$data@weights
expect_error(
  sampler$setWeights(rep(1, n - 1L)),
  "'weights' must have length equal to that of 'y'"
)
expect_identical(sampler$data@weights, weightsBefore)

weightsBad <- rep(1, n)
weightsBad[1L] <- NA
expect_error(sampler$setWeights(weightsBad), "'weights' cannot be NA")
expect_identical(sampler$data@weights, weightsBefore)

weightsBad <- rep(1, n)
weightsBad[1L] <- -1
expect_error(
  sampler$setWeights(weightsBad),
  "'weights' must all be non-negative"
)
expect_identical(sampler$data@weights, weightsBefore)

offsetBad <- rep(0, n)
offsetBad[1L] <- NA
offsetBefore <- sampler$data@offset
expect_error(sampler$setOffset(offsetBad), "'offset' contains missing values")
expect_identical(sampler$data@offset, offsetBefore)

sigmasBefore <- sampler$getSigmas()
expect_error(sampler$setSigma(c(1, 2)), "'sigma' must be of length 1")
expect_identical(sampler$getSigmas(), sigmasBefore)
expect_error(sampler$setSigma(NA_real_), "'sigma' must be positive")
expect_identical(sampler$getSigmas(), sigmasBefore)
expect_error(sampler$setSigma(-1), "'sigma' must be positive")
expect_identical(sampler$getSigmas(), sigmasBefore)

# a well-formed value is installed on the gaussian sampler the rejections above
# left untouched
sampler$setSigma(2)
expect_equal(sampler$getSigmas(), 2)

# setSigma is refused where the residual standard deviation is not a free
# parameter: the binary/latent families pin it at 1 by the model definition,
# and a heteroscedastic sampler's variance forest owns the residual scale. In
# both cases the value would be installed with the redraw that corrects it
# gated off, silently rescaling every leaf posterior precision. The predicate
# is the FAMILY, not "sigma is fixed": a gaussian resid.prior = fixed()
# sampler pins sigma too, and driving it per sweep is the supported outer-Gibbs
# conditioning idiom, asserted positively below.
xErr <- cbind(x = testData$x)
yBinary <- as.double(testData$z.0)
sampler.probit <- dbarts::dbarts(
  xErr,
  yBinary,
  family = "probit",
  control = control
)
expect_error(
  sampler.probit$setSigma(2),
  "response family fixes the residual standard deviation"
)
sampler.logistic <- dbarts::dbarts(
  xErr,
  yBinary,
  family = "logistic",
  control = control
)
expect_error(
  sampler.logistic$setSigma(2),
  "response family fixes the residual standard deviation"
)

sampler.variance <- dbarts::dbarts(
  y ~ x + z,
  train,
  control = control,
  variance = TRUE,
  n.trees.variance = 10L
)
expect_error(
  sampler.variance$setSigma(2),
  "variance forest owns the residual scale"
)

# the two permitted fixed-sigma cases: gaussian + resid.prior = fixed() (the
# outer-Gibbs conduit stan4bart's mvbart drives) and aft, whose sigma is drawn
# conjugately like gaussian's
sampler.fixed <- dbarts::dbarts(
  y ~ x + z,
  train,
  resid.prior = fixed(1),
  control = control
)
sampler.fixed$setSigma(0.4)
expect_equal(sampler.fixed$getSigmas(), 0.4)

sampler.aft <- dbarts::dbarts(
  xErr,
  cbind(time = exp(scale(testData$y)[, 1L]), status = testData$z.0),
  family = "aft",
  control = control
)
sampler.aft$setSigma(0.7)
expect_equal(sampler.aft$getSigmas(), 0.7)

# The transactional and per-observation predictor entries are refused on a
# heteroscedastic sampler: the variance forest holds its trees outside the
# forest vector revalidateAllChains loops, so an accepted change would leave
# s^2(x) routed by the predictors the forest was built with. The forced paths,
# setCutPoints and setData re-route it (setData resizes its n-sized buffers
# first), and are covered by the mutation test file.
xVariance <- as.matrix(train[, c("x", "z")])
xReplacement <- xVariance
xReplacement[, 1L] <- rev(xVariance[, 1L])
varianceRefusal <- "variance forest keeps routing observations"
# a single column without forceUpdate is the transactional path (a whole
# matrix defaults to the forced one, which is routed and tested elsewhere)
expect_error(
  sampler.variance$setPredictor(xReplacement[, 1L], 1L),
  varianceRefusal
)
expect_error(
  sampler.variance$setPredictor(
    xReplacement[, 1L],
    1L,
    forceUpdate = "partial"
  ),
  varianceRefusal
)
expect_error(
  dbarts::updatePredictorPerObservationJointly(
    list(sampler.variance),
    xReplacement[, 1L],
    "x"
  ),
  varianceRefusal
)

# resid.prior calibrates the variance forest's scale leaf rather than a sigma
# that is not a parameter here, so setModel recalibrates the leaf from the
# incoming triple instead of refusing it (S5 replacing S1's refusal). A prior
# this tight swamps the residuals: every leaf factor is drawn at essentially
# its prior value, so the surface collapses onto sigest^2 * scale wherever the
# data would have put it (it sat 20% below that pin beforehand).
varianceSurface <- sampler.variance$run(0L, 5L)$variance
expect_true(sd(varianceSurface) / mean(varianceSurface) > 0.5)
modelWithNewPrior <- sampler.variance$model
modelWithNewPrior@resid.prior@df <- 2000
modelWithNewPrior@resid.prior@quantile <- 0.001
sampler.variance$setModel(modelWithNewPrior)
recalibrated <- sampler.variance$run(0L, 5L)$variance
pinned <- attr(sampler.variance$data, "sigma")^2 * qchisq(0.999, 2000) / 2000
expect_true(sd(recalibrated) / mean(recalibrated) < 0.1)
expect_true(abs(mean(recalibrated) / pinned - 1) < 0.05)
expect_true(abs(mean(varianceSurface) / pinned - 1) > 0.1)

# and sigma stays pinned across setModel: under a variance forest it is not a
# parameter, so neither branch of the model's sigma clause may move it
pinnedSigma <- sampler.variance$getSigmas()
sampler.variance$setModel(sampler.variance$model)
expect_true(all(sampler.variance$run(0L, 5L)$sigma == pinnedSigma))
expect_equal(sampler.variance$getSigmas(), pinnedSigma)
# non-vacuity: the same recipe moves sigma where it IS a parameter
sampler.homoscedastic <- dbarts::dbarts(y ~ x + z, train, control = control)
homoscedasticSigma <- sampler.homoscedastic$getSigmas()
sampler.homoscedastic$setModel(sampler.homoscedastic$model)
expect_true(
  any(sampler.homoscedastic$run(0L, 5L)$sigma != homoscedasticSigma)
)

# a warm start carries no variance trees, so a donor and destination that
# disagree on the variance forest are refused in both directions
donor.variance <- dbarts::dbarts(
  y ~ x + z,
  train,
  control = control,
  variance = TRUE,
  n.trees.variance = 10L
)
invisible(donor.variance$run(0L, 1L))
donor.homoscedastic <- dbarts::dbarts(y ~ x + z, train, control = control)
invisible(donor.homoscedastic$run(0L, 1L))
expect_error(
  sampler.homoscedastic$installTrees(donor.variance),
  "shape-compatible"
)
expect_error(
  sampler.variance$installTrees(donor.homoscedastic),
  "shape-compatible"
)
# two heteroscedastic samplers agree on the shape, and the donor's variance
# trees now ride the install (S5); the state-level gate is in tests/cpp
expect_silent(sampler.variance$installTrees(donor.variance))

# The response-side conduits carry the same scale pin the transactional
# predictor paths carry: a variance forest's scale leaf is calibrated once,
# against the response transform in force at creation, and updateScale = TRUE
# re-anchors that transform under the calibration - the fit then runs away
# while getSigmas(), which reads the pinned sigma rather than the forest, shows
# nothing. updateScale = FALSE pins the transform and is the supported swap.
varianceScaleRefusal <- "variance forest is calibrated against the response"
expect_error(
  sampler.variance$setResponse(train$y, updateScale = TRUE),
  varianceScaleRefusal
)
expect_error(
  sampler.variance$setOffset(rep(0.5, n), updateScale = TRUE),
  varianceScaleRefusal
)
expect_silent(sampler.variance$setResponse(train$y))
expect_silent(sampler.variance$setOffset(rep(0.5, n)))

# A binary sampler's latents are drawn against y directly, so an out-of-support
# swap is silently garbage rather than an error; mutation now accepts exactly
# what creation accepts. The refusal names the family, not the length.
binaryRefusal <- "must be coded 0 or 1"
expect_error(sampler.probit$setResponse(rnorm(length(yBinary))), binaryRefusal)
expect_error(sampler.logistic$setResponse(yBinary + 0.5), binaryRefusal)
# a well-formed swap still lands
sampler.probit$setResponse(1 - yBinary)
expect_identical(sampler.probit$data@y, 1 - yBinary)

# the weight-change refusal names the family that refuses: the single message
# it used to carry told an aft, ordinal or nbinom caller about "a binary
# response" they never asked for
expect_error(
  sampler.aft$setWeights(rep(1, length(sampler.aft$data@y))),
  "aft \\(survival\\) models do not support case weights"
)

rm(
  binaryRefusal,
  varianceScaleRefusal,
  xVariance,
  xReplacement,
  varianceRefusal,
  varianceSurface,
  modelWithNewPrior,
  recalibrated,
  pinned,
  pinnedSigma,
  homoscedasticSigma,
  sampler.homoscedastic,
  donor.variance,
  donor.homoscedastic,
  n,
  sampler,
  control,
  test,
  train,
  yBad,
  yBefore,
  weightsBefore,
  weightsBad,
  offsetBad,
  offsetBefore,
  sigmasBefore,
  xErr,
  yBinary,
  sampler.probit,
  sampler.logistic,
  sampler.variance,
  sampler.fixed,
  sampler.aft
)

rm(testData)
