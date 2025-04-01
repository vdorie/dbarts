source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

# test that sampling from prior works correctly
train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
test <- data.frame(x = testData$x, z = 1 - testData$z)

set.seed(0L)
sampler <- dbarts::dbarts(
  y ~ x + z, train, test,
  control = dbarts::dbartsControl(n.threads = 1L, n.chains = 1L)
)

sampler$sampleTreesFromPrior()
sampler$sampleNodeParametersFromPrior()

trees <- sampler$getTrees()

observed <- sd(trees$value[trees$var == -1L])
expected <- sampler$model@node.scale / (sampler$model@node.hyperprior@k * sqrt(sampler$control@n.trees))
expect_true(abs(observed - expected) < 1.0e-3)

rm(expected, observed, trees, sampler, test, train)

rm(testData)

