# A sampleTreesFromPrior issued in the middle of a chain must return the forest
# to the zero-fit state a freshly built sampler carries. The falsifier is the
# two fit channels of one recorded sweep over the SAME design: the training
# channel is the incrementally rolled totalFits, while the test channel is
# rebuilt from the trees alone every recorded sweep, so a pre-reset sum left
# behind in totalFits separates them permanently. Each arm carries its
# no-reset control, which pins the tolerance at the last ulp.

source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 25L,
  updateState = FALSE,
  keepTrees = FALSE
)

makeSampler <- function(seed) {
  set.seed(seed)
  dbarts(x, y, test = x, control = control)
}

# the channels of a plain continued sweep agree, so the fixture itself is not
# what any arm below is measuring
sampler <- makeSampler(31L)
sampler$run(50L, 0L)
samples <- sampler$run(0L, 1L)
expect_equal(as.numeric(samples$train), as.numeric(samples$test))

# a mid-chain prior reset keeps them agreeing
sampler <- makeSampler(31L)
sampler$run(50L, 0L)
sampler$sampleTreesFromPrior()
samples <- sampler$run(0L, 1L)
expect_equal(as.numeric(samples$train), as.numeric(samples$test))

# growFromRoot consumes the post-reset state too: it rolls each tree's residual
# before rebuilding that tree's leaf map, so it reads the same cached fits
sampler <- makeSampler(17L)
sampler$run(50L, 0L)
sampler$growFromRoot(2L)
samples <- sampler$run(0L, 1L)
expect_equal(as.numeric(samples$train), as.numeric(samples$test))

sampler <- makeSampler(17L)
sampler$run(50L, 0L)
sampler$sampleTreesFromPrior()
sampler$growFromRoot(2L)
samples <- sampler$run(0L, 1L)
expect_equal(as.numeric(samples$train), as.numeric(samples$test))

rm(sampler, samples, makeSampler, x, y, testData)

# the latent families reach the same defect through their own working response,
# and separate by more than the gaussian arm does
source(
  system.file("common", "probitData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$X
z <- testData$Z

makeSampler <- function(seed) {
  set.seed(seed)
  dbarts(x, z, test = x, control = control)
}

sampler <- makeSampler(5L)
sampler$run(50L, 0L)
samples <- sampler$run(0L, 1L)
expect_equal(as.numeric(samples$train), as.numeric(samples$test))

sampler <- makeSampler(5L)
sampler$run(50L, 0L)
sampler$sampleTreesFromPrior()
samples <- sampler$run(0L, 1L)
expect_equal(as.numeric(samples$train), as.numeric(samples$test))

rm(sampler, samples, makeSampler, control, x, z, testData)
