source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

# sums of squared residuals descale to the original response scale, against
# the current (last recorded) fits
set.seed(0L)
sampler <- dbarts::dbarts(
  y ~ x + z, data.frame(y = testData$y, x = testData$x, z = testData$z),
  control = dbarts::dbartsControl(n.threads = 1L, n.chains = 1L,
                                  updateState = FALSE)
)
samples <- sampler$run(20L, 1L)
expect_equal(sampler$getSumsOfSquaredResiduals(),
             sum((testData$y - samples$train[, 1L])^2))

rm(samples, sampler, testData)
