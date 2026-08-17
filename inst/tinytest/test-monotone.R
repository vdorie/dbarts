# Per-variable monotone (mBART) constraints (docs/design/monotone.md): the
# `monotone` surface, its refusals, and family/prior forcing, plus a recovery
# smoke test. The exact-posterior gate lives in benchmarks/R/monotone-reference.R.

set.seed(22L)

# ---- surface: direction resolution and forcing ----

monotoneOf <- function(...) {
  sampler <- dbarts::dbarts(
    ...,
    control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L),
    resid.prior = fixed(1)
  )
  list(
    directions = attr(sampler$model, "monotone"),
    p.birth_death = sampler$model@p.birth_death,
    p.swap = sampler$model@p.swap,
    p.change = sampler$model@p.change,
    k = sampler$model@node.hyperprior
  )
}

n <- 60L
x <- matrix(runif(n * 3L), n, 3L, dimnames = list(NULL, c("a", "b", "c")))
y <- x[, 1L] - x[, 3L] + rnorm(n, 0, 0.2)

# glyphs, words, and integers all resolve to {-1, 0, +1} by column name
expect_equal(
  monotoneOf(x, y, monotone = c(a = "+", c = "-"))$directions,
  c(1L, 0L, -1L)
)
expect_equal(
  monotoneOf(x, y, monotone = c(a = "increasing", c = "decreasing"))$directions,
  c(1L, 0L, -1L)
)
expect_equal(
  monotoneOf(x, y, monotone = c(a = 1L, c = -1L))$directions,
  c(1L, 0L, -1L)
)
# an unnamed length-p vector is positional
expect_equal(
  monotoneOf(x, y, monotone = c(1L, 0L, -1L))$directions,
  c(1L, 0L, -1L)
)

# a monotone fit is forced to birth/death-only, fixed k = 2
forced <- monotoneOf(x, y, monotone = c(a = "+"))
expect_equal(forced$p.birth_death, 1)
expect_equal(forced$p.swap, 0)
expect_equal(forced$p.change, 0)
expect_inherits(forced$k, "dbartsFixedHyperprior")
expect_equal(forced$k@k, 2)

# no constraint leaves the model unmarked and the default proposals
plain <- monotoneOf(x, y)
expect_null(plain$directions)
expect_equal(plain$p.birth_death, 0.5)

# an all-zero spec is treated as no constraint
expect_null(monotoneOf(x, y, monotone = c(a = 0L, b = 0L, c = 0L))$directions)

# ---- refusals ----

# a categorical predictor cannot carry a direction
xf <- data.frame(a = runif(n), g = factor(sample(letters[1:3], n, TRUE)))
yf <- xf$a + rnorm(n, 0, 0.2)
expect_error(
  dbarts::dbarts(yf ~ ., data = xf, monotone = c(g = "+")),
  "categorical"
)

# an unrecognized name is an error
expect_error(
  monotoneOf(x, y, monotone = c(nosuchcolumn = "+")),
  "unrecognized"
)

# an explicit k hyperprior conflicts with the fixed-k monotone rule
expect_error(
  dbarts::dbarts(
    x,
    y,
    monotone = c(a = "+"),
    node.prior = normal(chi(1.5, 2)),
    control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L)
  ),
  "monotone"
)

# an explicit non-default proposal.probs conflicts with birth/death-only
expect_error(
  dbarts::dbarts(
    x,
    y,
    monotone = c(a = "+"),
    proposal.probs = c(
      birth_death = 0.6,
      swap = 0.1,
      change = 0.3,
      birth = 0.5
    ),
    control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L)
  ),
  "proposal.probs"
)

# an explicit proposal.probs = NULL is treated as absent: it succeeds
# under monotone and is bitwise identical to the defaulted call
ctrlD4 <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  n.samples = 5L,
  n.burn = 5L,
  seed = 55L
)
defaultedD4 <- dbarts::dbarts(x, y, monotone = c(a = "+"), control = ctrlD4)
defaultedD4$sampleTreesFromPrior()
samplesDefaultedD4 <- defaultedD4$run(0L, 5L)

explicitNullD4 <- dbarts::dbarts(
  x,
  y,
  monotone = c(a = "+"),
  proposal.probs = NULL,
  control = ctrlD4
)
explicitNullD4$sampleTreesFromPrior()
samplesExplicitNullD4 <- explicitNullD4$run(0L, 5L)

expect_identical(samplesDefaultedD4$train, samplesExplicitNullD4$train)
expect_identical(samplesDefaultedD4$sigma, samplesExplicitNullD4$sigma)

rm(
  ctrlD4,
  defaultedD4,
  samplesDefaultedD4,
  explicitNullD4,
  samplesExplicitNullD4
)

rm(monotoneOf, forced, plain, x, y, xf, yf)

# ---- recovery: a monotone truth is recovered with tighter intervals ----

set.seed(101L)
nRec <- 200L
xRec <- matrix(sort(runif(nRec)), nRec, 1L, dimnames = list(NULL, "x1"))
truth <- 1.5 * xRec[, 1L]
yRec <- truth + rnorm(nRec, 0, 0.3)

fitArgs <- list(
  n.trees = 50L,
  n.chains = 1L,
  n.threads = 1L,
  n.samples = 500L,
  n.burn = 300L,
  keepTrees = FALSE,
  verbose = FALSE
)
fitMono <- do.call(
  dbarts::bart2,
  c(list(xRec, yRec, monotone = c(x1 = "+")), fitArgs)
)
fitFree <- do.call(dbarts::bart2, c(list(xRec, yRec), fitArgs))

# the constrained posterior-mean fit is monotone in x (x is sorted)
expect_true(all(diff(fitMono$yhat.train.mean) > -1e-8))

# it recovers the truth
expect_true(
  sqrt(mean((fitMono$yhat.train.mean - truth)^2)) / sd(truth) < 0.25
)

# and its pointwise posterior intervals are, on average, tighter than the
# unconstrained fit's (the constraint pools information across the ordering)
sdMono <- mean(apply(fitMono$yhat.train, 2L, sd))
sdFree <- mean(apply(fitFree$yhat.train, 2L, sd))
expect_true(sdMono < sdFree)

# ---- a non-monotone truth flattens under the constraint, does not crash ----

set.seed(202L)
yBump <- sin(2 * pi * xRec[, 1L]) + rnorm(nRec, 0, 0.3)
fitBump <- do.call(
  dbarts::bart2,
  c(list(xRec, yBump, monotone = c(x1 = "+")), fitArgs)
)
fitBumpFree <- do.call(dbarts::bart2, c(list(xRec, yBump), fitArgs))
# the fit stays monotone (the true rise-then-fall is flattened, not recovered)
expect_true(all(diff(fitBump$yhat.train.mean) > -1e-8))
expect_true(all(is.finite(fitBump$yhat.train.mean)))
# the flattening bias: the unconstrained fit recovers the rise-and-fall (large
# range); the monotone fit cannot chase the descent, so its range is smaller
expect_true(
  diff(range(fitBump$yhat.train.mean)) <
    diff(range(fitBumpFree$yhat.train.mean))
)

rm(
  nRec,
  xRec,
  truth,
  yRec,
  fitArgs,
  fitMono,
  fitFree,
  sdMono,
  sdFree,
  yBump,
  fitBump,
  fitBumpFree
)

# ---- the prior draw is constrained too ----

# samplePriorPredictive installs its draws as the sampler's live leaf values,
# so under a constraint they must come from the truncated prior; an
# unconstrained draw would both misstate the prior and strand the sampler in an
# infeasible state.

set.seed(303L)
nPri <- 80L
xPri <- matrix(sort(runif(nPri)), nPri, 1L, dimnames = list(NULL, "x1"))
samplerPri <- dbarts::dbarts(
  xPri,
  xPri[, 1L] + rnorm(nPri, 0, 0.3),
  monotone = c(x1 = "+"),
  control = dbarts::dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 25L)
)

priorEv <- samplePriorPredictive(samplerPri, n.samples = 100L, type = "ev")
expect_equal(dim(priorEv), c(100L, nPri))
expect_true(all(is.finite(priorEv)))

# x is sorted, so each individual draw - not merely their average - is
# non-decreasing across the constrained column
expect_true(all(apply(priorEv, 1L, function(f) all(diff(f) > -1e-8))))

# the prior draw is not the constant zero forest it would collapse to if the
# rejection simply gave up
expect_true(mean(apply(priorEv, 1L, function(f) diff(range(f)))) > 0.1)

# and a sampler left in that state runs on
samplerPri$sampleTreesFromPrior()
samplerPri$sampleNodeParametersFromPrior()
samplesPri <- samplerPri$run(10L, 10L)
expect_true(all(is.finite(samplesPri$train)))

rm(nPri, xPri, samplerPri, priorEv, samplesPri)
