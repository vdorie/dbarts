source(system.file("common", "hillData.R", package = "dbarts"), local = TRUE)

# test that sampling from prior works correctly
train <- data.frame(y = testData$y, x = testData$x, z = testData$z)
test <- data.frame(x = testData$x, z = 1 - testData$z)

set.seed(0L)
sampler <- dbarts::dbarts(
  y ~ x + z,
  train,
  test,
  control = dbarts::dbartsControl(n.threads = 1L, n.chains = 1L)
)

sampler$sampleTreesFromPrior()
sampler$sampleLeafParametersFromPrior()

trees <- sampler$getTrees()

leaves <- trees$value[trees$var == -1L]
observed <- sd(leaves)
expected <- sampler$model@leaf.scale /
  (sampler$model@leaf.hyperprior@k * sqrt(sampler$control@n.trees))
# within three standard errors of a normal sample's standard deviation; the
# number of leaves follows the prior over trees, which the grid sets
expect_true(abs(observed - expected) < 3 * expected / sqrt(2 * length(leaves)))

rm(expected, observed, leaves, trees, sampler, test, train)

rm(testData)
