# $setActiveRows redraws the latent of every row it switches from inactive to
# active, against the current fit and before it returns, so a larger sampler
# that redraws the mask from the fit - the latents integrated out - draws the
# mask and those latents jointly. Nothing else moves: a row that stays active,
# stays inactive or is switched off keeps its latent, and a mask that
# reactivates no row draws nothing, so reinstalling the mask in force leaves
# the stream where it was.

set.seed(20261005L)
n <- 120L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y.binary <- as.double(x[, 1L] + rnorm(n, 0, 0.3) > 0.5)
y.ordinal <- as.double(1L + (seq_len(n) %% 3L))
y.counts <- as.double(rpois(n, exp(1 + x[, 1L])))
# rows 1, 5, 9, ... start inactive; the second mask switches every other one
# of those back in and switches rows 2, 6, 10, ... off
first <- as.double(seq_len(n) %% 4L != 1L)
second <- as.double(seq_len(n) %% 8L != 1L & seq_len(n) %% 4L != 2L)
reactivated <- first == 0 & second == 1
expect_true(sum(reactivated) > 10L)

control <- dbarts::dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 20L,
  n.samples = 5L,
  updateState = FALSE,
  seed = 13L
)
samplers <- list(
  probit = function() {
    dbarts::dbarts(x, y.binary, family = "probit", control = control)
  },
  ordinal = function() {
    dbarts::dbarts(x, y.ordinal, family = "ordinal", control = control)
  },
  logistic = function() {
    dbarts::dbarts(x, y.binary, family = "logistic", control = control)
  },
  nbinom = function() {
    dbarts::dbarts(x, y.counts, family = "nbinom", control = control)
  }
)

for (family in names(samplers)) {
  sampler <- samplers[[family]]()
  twin <- samplers[[family]]()
  invisible(sampler$run(20L, 1L))
  invisible(twin$run(20L, 1L))
  sampler$setActiveRows(first)
  twin$setActiveRows(first)
  invisible(sampler$run(10L, 1L))
  invisible(twin$run(10L, 1L))

  # the mask in force, installed again: no latent moves and no variate is
  # drawn, so the twin that skips the call draws the same sweeps
  held <- sampler$getLatents()
  sampler$setActiveRows(first)
  expect_identical(sampler$getLatents(), held, info = family)
  expect_identical(
    sampler$run(0L, 3L)$train,
    twin$run(0L, 3L)$train,
    info = family
  )

  # a mask that reactivates rows: exactly those rows' latents are redrawn
  held <- sampler$getLatents()
  sampler$setActiveRows(second)
  moved <- sampler$getLatents() != held
  expect_true(all(moved[reactivated]), info = family)
  expect_false(any(moved[!reactivated]), info = family)
  expect_true(all(is.finite(sampler$getLatents())), info = family)

  # clearing the mask reactivates every row still out
  held <- sampler$getLatents()
  sampler$setActiveRows(NULL)
  moved <- sampler$getLatents() != held
  expect_true(all(moved[second == 0]), info = family)
  expect_false(any(moved[second == 1]), info = family)
  expect_true(all(is.finite(sampler$run(0L, 3L)$train)), info = family)
}

# a probit latent redrawn on reactivation is a draw given the row's response:
# its sign is the response's
probit <- samplers$probit()
invisible(probit$run(20L, 1L))
probit$setActiveRows(first)
invisible(probit$run(10L, 1L))
probit$setActiveRows(NULL)
expect_true(all((probit$getLatents() > 0) == (y.binary == 1)))

# 0/1 case weights on probit are the mask, so $setWeights redraws the same
# latents $setActiveRows does
byWeights <- samplers$probit()
byMask <- samplers$probit()
for (mask in list(first, second, rep(1, n))) {
  invisible(byWeights$run(5L, 1L))
  invisible(byMask$run(5L, 1L))
  byWeights$setWeights(mask)
  byMask$setActiveRows(mask)
  expect_identical(byWeights$getLatents(), byMask$getLatents())
}
expect_identical(byWeights$run(0L, 3L)$train, byMask$run(0L, 3L)$train)

# a gaussian sampler holds no latent, so a mask that comes and goes draws
# nothing there: it is still setWeights(a) and back
masked <- dbarts::dbarts(x, x[, 1L] + y.binary, control = control)
weighted <- dbarts::dbarts(x, x[, 1L] + y.binary, control = control)
invisible(masked$run(10L, 1L))
invisible(weighted$run(10L, 1L))
masked$setActiveRows(first)
weighted$setWeights(first)
invisible(masked$run(5L, 1L))
invisible(weighted$run(5L, 1L))
masked$setActiveRows(second)
weighted$setWeights(second)
expect_identical(masked$run(0L, 3L)$train, weighted$run(0L, 3L)$train)

rm(
  n,
  x,
  y.binary,
  y.ordinal,
  y.counts,
  first,
  second,
  reactivated,
  control,
  samplers,
  family,
  sampler,
  twin,
  held,
  moved,
  probit,
  byWeights,
  byMask,
  mask,
  masked,
  weighted
)
