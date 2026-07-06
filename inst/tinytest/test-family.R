# the family argument: auto dispatch, forced gaussian on 0/1 responses, and
# the public logistic family

set.seed(11)
n <- 300L
x <- matrix(runif(n * 3L), n)
f <- 3 * x[, 1L] - 1.5
y.binary <- rbinom(n, 1L, plogis(f))
y.continuous <- f + rnorm(n, 0, 0.5)

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE
)

# auto preserves the existing dispatch
sampler.auto <- dbarts(y.binary ~ x, control = control)
expect_equal(sampler.auto$control@family, "probit")
expect_true(sampler.auto$control@binary)
expect_equal(
  dbarts(y.continuous ~ x, control = control)$control@family,
  "gaussian"
)

# gaussian on a 0/1 response is allowed and fits a continuous model
sampler.gauss <- dbarts(y.binary ~ x, family = "gaussian", control = control)
expect_false(sampler.gauss$control@binary)
samples.gauss <- sampler.gauss$run(50L, 50L)
expect_true(
  all(is.finite(samples.gauss$sigma)) &&
    length(unique(samples.gauss$sigma)) > 1L
)

# binary families need a 0/1 response
expect_error(
  dbarts(y.continuous ~ x, family = "probit", control = control),
  pattern = "requires a response coded 0/1"
)

control.bc <- dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = 50L)
sampler.logit <- dbarts(y.binary ~ x, family = "logistic", control = control.bc)
expect_equal(sampler.logit$control@family, "logistic")
expect_inherits(sampler.logit$model@node.hyperprior, "dbartsChiHyperprior")
expect_equal(sampler.logit$model@node.scale, pi * sqrt(3))

# fits live on the latent logit scale and recover the signal
samples.logit <- sampler.logit$run(300L, 300L)
p.hat <- rowMeans(plogis(samples.logit$train))
expect_true(cor(p.hat, plogis(f)) > 0.8)
expect_true(mean(p.hat[y.binary == 1L]) > mean(p.hat[y.binary == 0L]))

# the family survives save/load: the pointer is recreated from the stored
# state with the same response model (a single chain draws through R's
# generator, so bitwise continuation needs the seed reset)
serialized <- tempfile(fileext = ".rds")
saveRDS(sampler.logit, serialized)
set.seed(101)
result.a <- sampler.logit$run(0L, 5L, updateState = FALSE)
sampler.loaded <- readRDS(serialized)
expect_equal(sampler.loaded$control@family, "logistic")
expect_true(all(sampler.loaded$getLatents() > 0)) # omega, not probit z
set.seed(101)
result.b <- sampler.loaded$run(0L, 5L, updateState = FALSE)
expect_identical(result.a, result.b)
unlink(serialized)

# setControl cannot silently change the family
newControl <- sampler.logit$control
sampler.logit$setControl(newControl)
expect_equal(sampler.logit$control@family, "logistic")

# bart2 forwards the argument
fit.bart2 <- bart2(
  y.binary ~ x,
  family = "gaussian",
  n.samples = 30L,
  n.burn = 30L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)
expect_true(!is.null(fit.bart2$sigma))

# probit refuses weights at the R layer through the wrapper (the bridge keeps
# the same refusal as a backstop for direct-API consumers): a weighted probit
# has no tractable latent-variable form. Logistic weights are covered in
# test-weighted-logistic.R
expect_error(
  bart2(
    y.binary ~ x,
    weights = runif(n, 0.5, 1.5),
    n.samples = 5L,
    n.burn = 5L,
    n.trees = 25L,
    n.chains = 1L,
    n.threads = 1L,
    verbose = FALSE
  ),
  pattern = "probit models do not support weights"
)

# the wrappers record the family and transform through its link
fit.probit <- bart2(
  y.binary ~ x,
  n.samples = 40L,
  n.burn = 40L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_equal(fit.probit$family, "probit")
expect_equal(extract(fit.probit, "ev"), pnorm(extract(fit.probit, "bart")))

fit.logit <- bart2(
  y.binary ~ x,
  family = "logistic",
  n.samples = 40L,
  n.burn = 40L,
  n.trees = 25L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE,
  keepTrees = TRUE
)
expect_equal(fit.logit$family, "logistic")
latents <- extract(fit.logit, "bart")
expect_equal(extract(fit.logit, "ev"), plogis(latents))
expect_equal(
  predict(fit.logit, x, type = "ev"),
  plogis(predict(fit.logit, x, type = "bart"))
)
expect_equal(
  fitted(fit.logit),
  apply(plogis(latents), length(dim(latents)), mean)
)

# fits saved before the family element existed fall back to probit
expect_equal(dbarts:::probabilityFromLatents(0.5, list()), pnorm(0.5))
