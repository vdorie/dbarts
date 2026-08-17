source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

# sums of squared residuals descale to the original response scale, against
# the current (last recorded) fits
set.seed(0L)
sampler <- dbarts::dbarts(
  y ~ x + z,
  data.frame(y = testData$y, x = testData$x, z = testData$z),
  control = dbarts::dbartsControl(
    n.threads = 1L,
    n.chains = 1L,
    updateState = FALSE
  )
)
samples <- sampler$run(20L, 1L)
expect_equal(
  sampler$getSumsOfSquaredResiduals(),
  sum((testData$y - samples$train[, 1L])^2)
)

# On a K-forest the fitted value is the COMBINED location
# sum_f m_f(i) f_f(x_i), not forest 0's own total, and this getter reported the
# latter. storeSample routes the recorded train channel through that same
# combination and runs last in the sweep, so the recorded draw is the oracle -
# no reconstruction is possible from the public getters anyway, since none of
# them exposes a basis. Both K = 2 and K = 3 with a multi-column basis, where
# the per-observation multiplier is not a single scalar per forest.
kForestX <- cbind(testData$x, testData$z)
kForestControl <- dbartsControl(
  n.chains = 2L,
  n.threads = 1L,
  n.trees = 25L,
  updateState = FALSE,
  seed = 12L
)
zBasis <- cbind(1 - testData$z, testData$z)
for (bases in list(
  list(NULL, zBasis),
  list(NULL, zBasis, cbind(1, testData$x / 100))
)) {
  kForest <- dbarts(
    dbartsData(kForestX, testData$y, bases = bases),
    control = kForestControl
  )
  kResult <- kForest$run(20L, 1L)
  expect_equal(
    kForest$getSumsOfSquaredResiduals(),
    colSums((testData$y - kResult$train[, 1L, ])^2)
  )
}

rm(
  bases,
  kForest,
  kForestControl,
  kForestX,
  kResult,
  samples,
  sampler,
  testData,
  zBasis
)
