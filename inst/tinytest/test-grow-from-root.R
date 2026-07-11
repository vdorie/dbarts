source(
  system.file("common", "friedmanData.R", package = "dbarts"),
  local = TRUE
)

x <- testData$x
y <- testData$y

makeSampler <- function(n.trees = 25L) {
  dbarts::dbarts(
    x,
    y,
    control = dbarts::dbartsControl(
      n.trees = n.trees,
      n.chains = 1L,
      updateState = FALSE,
      keepTrees = FALSE
    )
  )
}

## grow-from-root in place yields a well-formed forest with finite predictions
sampler <- makeSampler()
sampler$growFromRoot(3L)
grownFit <- as.vector(sampler$predict(x))
expect_true(all(is.finite(grownFit)))

## the exact sampler continues from the grown forest with finite fits
res <- sampler$run(0L, 10L)
expect_true(all(is.finite(res$train)))

## seeded reproducibility: same seed -> identical grown fit
set.seed(42L)
s1 <- makeSampler()
s1$growFromRoot(2L)
f1 <- as.vector(s1$predict(x))
set.seed(42L)
s2 <- makeSampler()
s2$growFromRoot(2L)
f2 <- as.vector(s2$predict(x))
expect_equal(f1, f2)

## a different seed grows a different forest (grow consumes RNG)
set.seed(7L)
s3 <- makeSampler()
s3$growFromRoot(2L)
f3 <- as.vector(s3$predict(x))
expect_false(isTRUE(all.equal(f1, f3)))

## cross-sampler workflow: donor$growFromRoot -> target$installTrees round trip
## reproduces the grown forest with no bespoke install code
donor <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.trees = 8L,
    n.chains = 1L,
    updateState = FALSE,
    keepTrees = FALSE
  )
)
donor$growFromRoot(2L)
target <- dbarts::dbarts(
  x,
  y,
  control = dbarts::dbartsControl(
    n.trees = 8L,
    n.chains = 1L,
    updateState = FALSE,
    keepTrees = FALSE
  )
)
target$installTrees(donor)
expect_equal(as.vector(target$predict(x)), as.vector(donor$predict(x)))

## n.sweeps must be a single positive integer
expect_error(makeSampler()$growFromRoot(0L), "positive integer")
expect_error(makeSampler()$growFromRoot(-2L), "positive integer")

## the linear-leaf model refuses grow-from-root with an informative error
linSampler <- dbarts::dbarts(
  x,
  y,
  node.prior = linear(1L),
  control = dbarts::dbartsControl(
    n.trees = 8L,
    n.chains = 1L,
    updateState = FALSE
  )
)
expect_error(linSampler$growFromRoot(2L), "constant-leaf")
