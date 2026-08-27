source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

# The Student-t residual surface: the
# resid.dist = student() constructor, its validation, the gaussian-only family
# refusal, and that the resolved degrees of freedom ride the model's resid.df
# attribute the C bridge reads (absent = gaussian, 0 = estimate, positive =
# fixed). A tiny smoke fit in both df modes, predict/fitted shapes, and a
# serialize+restore round-trip.

x <- testData$x
y <- testData$y

# --- constructor construction and validation -------------------------------

gaussianDist <- dbarts:::gaussian()
studentEstimate <- dbarts:::student()
studentFixed <- dbarts:::student(df = 4)

expect_true(is(gaussianDist, "dbartsResidDist"))
expect_true(is(studentEstimate, "dbartsStudentDist"))
expect_true(is.na(studentEstimate@df)) # NULL df means estimate
expect_equal(studentFixed@df, 4)

# df must be NULL (estimate) or a single positive finite number
expect_error(dbarts:::student(df = -1), "positive finite")
expect_error(dbarts:::student(df = 0), "positive finite")
expect_error(dbarts:::student(df = Inf), "positive finite")
expect_error(dbarts:::student(df = NA_real_), "positive finite")
expect_error(dbarts:::student(df = c(2, 4)), "single positive finite")
expect_error(dbarts:::student(df = "4"), "single positive finite")

# --- the resid.df attribute the bridge reads -------------------------------

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 30L,
  updateState = FALSE
)

samplerGaussian <- dbarts::dbarts(x, y, control = control)
expect_null(attr(samplerGaussian$model, "resid.df"))

samplerFixed <- dbarts::dbarts(
  x,
  y,
  control = control,
  resid.dist = student(df = 4)
)
expect_equal(attr(samplerFixed$model, "resid.df"), 4)

samplerEstimate <- dbarts::dbarts(
  x,
  y,
  control = control,
  resid.dist = student()
)
expect_equal(attr(samplerEstimate$model, "resid.df"), 0) # 0 signals estimate

# --- family refusal (R-side backstop for the C cross-check) ----------------

yBinary <- as.numeric(y > median(y))
expect_error(
  dbarts::dbarts(
    x,
    yBinary,
    control = control,
    family = "probit",
    resid.dist = student(df = 4)
  ),
  "student residuals require a continuous gaussian response"
)
# gaussian() on a binary probit fit is fine (the default error law), and
# attaches no resid.df attribute
samplerBinary <- dbarts::dbarts(
  x,
  yBinary,
  control = control,
  family = "probit",
  resid.dist = gaussian()
)
expect_null(attr(samplerBinary$model, "resid.df"))

# --- smoke fit in both df modes, predict/fitted shapes ---------------------

set.seed(7L)
fitFixed <- dbarts::bart2(
  x,
  y,
  n.trees = 25L,
  n.samples = 40L,
  n.burn = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE,
  resid.dist = student(df = 4)
)
expect_equal(dim(fitFixed$yhat.train), c(40L, length(y)))
predsFixed <- predict(fitFixed, x)
expect_equal(dim(predsFixed), c(40L, length(y)))
expect_equal(length(fitted(fitFixed)), length(y))

set.seed(7L)
fitEstimate <- dbarts::bart2(
  x,
  y,
  n.trees = 25L,
  n.samples = 40L,
  n.burn = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE,
  resid.dist = student()
)
expect_equal(dim(fitEstimate$yhat.train), c(40L, length(y)))
expect_true(all(is.finite(fitEstimate$sigma)))

# --- serialize + restore round-trip (the storeState ritual) ----------------

roundTripRobustErrors <- function(object) {
  tempFile <- tempfile()
  saveRDS(object, file = tempFile)
  on.exit(unlink(tempFile))
  readRDS(tempFile)
}

fitFixed$fit$storeState()
restored <- roundTripRobustErrors(fitFixed)
# the resolved df survives serialization on the model object
expect_equal(attr(restored$fit$model, "resid.df"), 4)
# and the restored fit reproduces its pre-serialization predictions
expect_equal(predict(restored, x), predsFixed)

rm(testData)
