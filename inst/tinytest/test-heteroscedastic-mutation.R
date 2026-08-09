# A heteroscedastic sampler's variance forest lives outside the forest vector
# the mutation helpers loop, so a predictor or cut-grid change used to leave
# s^2(x) routing observations by the predictors the forest was built with.
# The forced predictor paths and setCutPoints now re-route it
# (docs/plans/variance-forest-mutation-routing.md, slice S3).

set.seed(29, sample.kind = "Rejection")

n <- 200L
x <- cbind(x1 = runif(n), x2 = runif(n), x3 = runif(n))
y <- 2 * x[, 1L] + ifelse(x[, 2L] < 0.5, 0.3, 1.5) * rnorm(n)
control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.trees = 25L,
  n.burn = 0L,
  verbose = FALSE
)
buildVarianceSampler <- function(predictors, response) {
  dbarts::dbarts(
    predictors,
    response,
    control = control,
    variance = TRUE,
    n.trees.variance = 10L
  )
}
# distinct reported values, quantized well below any routing difference: the
# reported channels carry per-observation rounding of order 1e-16, which is not
# what these bounds are counting
numDistinct <- function(values) length(unique(signif(as.vector(values), 8)))

# ---- forced predictor swap: an all-identical design admits one value ----
# every row carries the same predictors, so every tree of either forest routes
# every row to the same leaf and both surfaces are constant. The bound is 1
# because the design has one distinct row, not because of any tolerance.
sampler <- buildVarianceSampler(x, y)
invisible(sampler$run(50L, 0L))
identicalRows <- matrix(rep(c(0.4, 0.55, 0.7), each = n), n, 3L)
sampler$setPredictor(identicalRows, forceUpdate = TRUE)
swapped <- sampler$run(0L, 1L)
expect_equal(numDistinct(swapped$variance), 1L)
# non-vacuity: the mean channel, which has always been re-routed, reaches the
# same bound on the same mutated sampler
expect_equal(numDistinct(swapped$train), 1L)
# and the bound is achievable, not merely small: a fresh fit on the
# post-mutation design reports one value too
fresh <- buildVarianceSampler(identicalRows, y)
expect_equal(numDistinct(fresh$run(50L, 1L)$variance), 1L)

# ---- setCutPoints: a coarsened grid caps the surface combinatorially ----
# 2 cuts on each of 3 columns quantizes the design into at most 3^3 = 27 cells,
# and s^2(x_i) is a function of row i's cell alone. A stale partition reports a
# value per training row.
coarse <- buildVarianceSampler(x, y)
invisible(coarse$run(50L, 0L))
coarse$setCutPoints(rep(list(c(1 / 3, 2 / 3)), 3L), 1:3)
quantized <- coarse$run(0L, 1L)
expect_true(numDistinct(quantized$variance) <= 27L)
expect_true(numDistinct(quantized$train) <= 27L)

# every serialized variance threshold must be a member of the new grid: the
# rebuild matches each flattened cut value against the current cut points
# exactly, so a state carrying a stale threshold cannot restore
coarse$storeState()
expect_silent(coarse$setState(coarse$state))
expect_true(numDistinct(coarse$run(0L, 1L)$variance) <= 27L)

# ---- statistical agreement with a from-scratch fit (post-mutation-assertions
# ---- .md): a repair that re-routes but scatters the wrong factors fails here
set.seed(37, sample.kind = "Rejection")
nAgree <- 300L
xOld <- cbind(runif(nAgree), runif(nAgree))
xNew <- cbind(runif(nAgree), runif(nAgree))
yAgree <- 2 * xNew[, 1L] + ifelse(xNew[, 2L] < 0.5, 0.4, 1.6) * rnorm(nAgree)
mutated <- buildVarianceSampler(xOld, yAgree)
invisible(mutated$run(50L, 0L))
mutated$setPredictor(xNew, forceUpdate = TRUE)
sMutated <- sqrt(apply(mutated$run(100L, 100L)$variance, 1L, mean))
scratch <- buildVarianceSampler(xNew, yAgree)
sScratch <- sqrt(apply(scratch$run(150L, 100L)$variance, 1L, mean))
expect_true(cor(sMutated, sScratch) > 0.9)
expect_true(abs(mean(sMutated) / mean(sScratch) - 1) < 0.15)

rm(
  n,
  x,
  y,
  control,
  buildVarianceSampler,
  numDistinct,
  sampler,
  identicalRows,
  swapped,
  fresh,
  coarse,
  quantized,
  nAgree,
  xOld,
  xNew,
  yAgree,
  mutated,
  sMutated,
  scratch,
  sScratch
)
