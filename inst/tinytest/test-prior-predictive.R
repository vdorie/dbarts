# samplePriorPredictive (issue #31): draws from the BART prior on a private
# freshly built sampler, so the caller's own sampler is left exactly as it
# was.

library(dbarts, quietly = TRUE)

set.seed(0)
n <- 30L
p <- 2L
x <- matrix(runif(n * p), n, p)
y <- x[, 1L] - x[, 2L] + rnorm(n, 0, 0.3)
z <- rbinom(n, 1L, pnorm(x[, 1L] - x[, 2L]))

control.plain <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 15L
)

sampler.gauss <- dbarts(y ~ x, control = control.plain)
sampler.binary <- dbarts(z ~ x, control = control.plain)
colnames(x) <- colnames(sampler.gauss$data@x)

x.test <- matrix(runif(6L * p), 6L, p)
colnames(x.test) <- colnames(sampler.gauss$data@x)

# (a) shapes: default x.test (training predictors) and an explicit x.test,
# both response families, both types
ev.gauss <- samplePriorPredictive(sampler.gauss, n.samples = 25L, type = "ev")
expect_equal(dim(ev.gauss), c(25L, n))
ppd.gauss <-
  samplePriorPredictive(sampler.gauss, n.samples = 25L, type = "ppd")
expect_equal(dim(ppd.gauss), c(25L, n))

ev.gauss.test <- samplePriorPredictive(
  sampler.gauss,
  x.test = x.test,
  n.samples = 10L,
  type = "ev"
)
expect_equal(dim(ev.gauss.test), c(10L, nrow(x.test)))

ev.binary <-
  samplePriorPredictive(sampler.binary, n.samples = 25L, type = "ev")
expect_equal(dim(ev.binary), c(25L, n))
expect_true(all(ev.binary >= 0 & ev.binary <= 1))

ppd.binary <-
  samplePriorPredictive(sampler.binary, n.samples = 25L, type = "ppd")
expect_equal(dim(ppd.binary), c(25L, n))
expect_true(all(ppd.binary %in% c(0, 1)))

ppd.binary.test <- samplePriorPredictive(
  sampler.binary,
  x.test = x.test,
  n.samples = 8L,
  type = "ppd"
)
expect_equal(dim(ppd.binary.test), c(8L, nrow(x.test)))

# (b) the caller's sampler is left untouched: fits off its live trees, and
# its residual sd, are bit-identical before and after
fitsBefore <- sampler.gauss$predict(x)
sigmasBefore <- sampler.gauss$getSigmas()
invisible(samplePriorPredictive(sampler.gauss, n.samples = 25L, type = "ppd"))
expect_identical(sampler.gauss$predict(x), fitsBefore)
expect_identical(sampler.gauss$getSigmas(), sigmasBefore)

# (c) seeding semantics. Successive calls are independent by default: each
# call constructs its private sampler fresh, seeding the engine RNG from
# R's stream when the control's seed is NA
indep1 <- samplePriorPredictive(sampler.gauss, n.samples = 20L, type = "ev")
indep2 <- samplePriorPredictive(sampler.gauss, n.samples = 20L, type = "ev")
expect_false(identical(indep1, indep2))

# set.seed() reproduces the whole draw, engine and observation layer alike
set.seed(11)
repeat1 <- samplePriorPredictive(sampler.gauss, n.samples = 20L, type = "ppd")
set.seed(11)
repeat2 <- samplePriorPredictive(sampler.gauss, n.samples = 20L, type = "ppd")
expect_identical(repeat1, repeat2)

set.seed(12)
repeatBinary1 <-
  samplePriorPredictive(sampler.binary, n.samples = 20L, type = "ppd")
set.seed(12)
repeatBinary2 <-
  samplePriorPredictive(sampler.binary, n.samples = 20L, type = "ppd")
expect_identical(repeatBinary1, repeatBinary2)

# a control-fixed seed pins the engine stream instead, so two ev calls
# repeat identically without any set.seed
control.pinned <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 15L,
  seed = 99L
)
sampler.pinned <- dbarts(y ~ x, control = control.pinned)
pinned1 <- samplePriorPredictive(sampler.pinned, n.samples = 20L, type = "ev")
pinned2 <- samplePriorPredictive(sampler.pinned, n.samples = 20L, type = "ev")
expect_identical(pinned1, pinned2)

# (d) statistical sanity: ppd adds observation noise on top of ev, so its
# aggregate variance is at least as large (gaussian). The pinned sampler
# pairs the forests across the two calls, making the comparison structural
# rather than luck across independent forest sets; binary ev draws staying
# valid probabilities is checked above
ppd.pinned <-
  samplePriorPredictive(sampler.pinned, n.samples = 20L, type = "ppd")
expect_true(var(as.vector(ppd.pinned)) >= var(as.vector(pinned1)))

# (e) validation errors
expect_error(samplePriorPredictive(1:5), pattern = "dbartsSampler")
expect_error(
  samplePriorPredictive(sampler.gauss, type = "bogus"),
  pattern = "should be one of"
)

# (f) a heteroscedastic sampler has no scalar sigma to add as observation
# noise - s(x) comes from a prior draw of the variance forest, which the ppd
# path does not make - so "ppd" is refused. "ev" never reaches that draw and
# is a legitimate mean-surface prior draw, so it keeps working.
sampler.variance <- dbarts(
  y ~ x,
  control = control.plain,
  variance = varianceForest(n.trees = 5L)
)
expect_error(
  samplePriorPredictive(sampler.variance, n.samples = 5L, type = "ppd"),
  pattern = "heteroscedastic sampler"
)
ev.variance <- samplePriorPredictive(
  sampler.variance,
  n.samples = 5L,
  type = "ev"
)
expect_equal(dim(ev.variance), c(5L, n))
expect_true(all(is.finite(ev.variance)))

rm(
  sampler.variance,
  ev.variance,
  n,
  p,
  x,
  y,
  z,
  control.plain,
  control.pinned,
  sampler.gauss,
  sampler.binary,
  sampler.pinned,
  x.test,
  ev.gauss,
  ppd.gauss,
  ev.gauss.test,
  ev.binary,
  ppd.binary,
  ppd.binary.test,
  fitsBefore,
  sigmasBefore,
  indep1,
  indep2,
  repeat1,
  repeat2,
  repeatBinary1,
  repeatBinary2,
  pinned1,
  pinned2,
  ppd.pinned
)
