# the offset path of the binary families: logistic and probit fits where the
# offset carries real signal, both in-core and through setOffset. The BART
# latent absorbs a fixed offset, so the probability is plogis/pnorm of the
# latent PLUS the offset, and a large offset shift moves the mean fitted
# probability.

set.seed(11)
n <- 200L
x <- matrix(runif(n * 2L), n)
f <- 1.5 * x[, 1L] - 0.5
o <- 2 * x[, 2L] - 1 # offset carries signal in column 2
z.logit <- rbinom(n, 1L, plogis(f + o))

control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 10L,
  updateState = FALSE
)

# logistic + offset fits on the latent logit scale
sampler <- dbarts(
  z.logit ~ x,
  family = "logistic",
  offset = o,
  control = control
)
expect_equal(sampler$control@family, "logistic")
expect_equal(sampler$data@offset, o)

samples <- sampler$run(200L, 200L)
expect_true(all(is.finite(samples$train)))
latents <- sampler$getLatents()
expect_true(all(is.finite(latents)) && all(latents > 0)) # omega, not probit z

# the fitted probability combines the latent with the offset and tracks truth
p.hat <- rowMeans(plogis(samples$train + o))
expect_true(all(is.finite(p.hat)))
expect_true(cor(p.hat, plogis(f + o)) > 0.7)

# predictions respond to the offset: a large positive shift drives the mean
# fitted probability toward 1, a large negative shift toward 0
sampler$setOffset(o + 5)
expect_equal(sampler$data@offset, o + 5)
samples.up <- sampler$run(200L, 200L)
expect_true(all(is.finite(samples.up$train)))
p.hat.up <- rowMeans(plogis(samples.up$train + o + 5))
expect_true(mean(p.hat.up) > 0.9)

sampler$setOffset(o - 5)
samples.dn <- sampler$run(200L, 200L)
expect_true(all(is.finite(samples.dn$train)))
p.hat.dn <- rowMeans(plogis(samples.dn$train + o - 5))
expect_true(mean(p.hat.dn) < 0.1)

# probit + offset for contrast: signed-z latents, not omega, and the same
# offset-driven probability shift through pnorm
set.seed(22)
z.probit <- rbinom(n, 1L, pnorm(f + o))
sampler.probit <- dbarts(
  z.probit ~ x,
  family = "probit",
  offset = o,
  control = control
)
expect_equal(sampler.probit$control@family, "probit")
samples.probit <- sampler.probit$run(200L, 200L)
expect_true(all(is.finite(samples.probit$train)))
p.hat.probit <- rowMeans(pnorm(samples.probit$train + o))
expect_true(cor(p.hat.probit, pnorm(f + o)) > 0.7)

sampler.probit$setOffset(o + 5)
samples.probit.up <- sampler.probit$run(200L, 200L)
expect_true(all(is.finite(samples.probit.up$train)))
expect_true(mean(rowMeans(pnorm(samples.probit.up$train + o + 5))) > 0.9)

rm(
  sampler,
  sampler.probit,
  samples,
  samples.up,
  samples.dn,
  samples.probit,
  samples.probit.up,
  p.hat,
  p.hat.up,
  p.hat.dn,
  p.hat.probit,
  latents,
  control,
  x,
  f,
  o,
  z.logit,
  z.probit,
  n
)
