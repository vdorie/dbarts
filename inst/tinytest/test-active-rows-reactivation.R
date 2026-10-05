# $setActiveRows redraws the latent of every row it switches from inactive to
# active, against the current fit and before it returns, so a larger sampler
# that redraws the mask from the fit - the latents integrated out - draws the
# mask and those latents jointly. Nothing else moves: a row that stays active,
# stays inactive or is switched off keeps its latent, and a mask that
# reactivates no row draws nothing, so reinstalling the mask in force leaves
# the stream where it was. A Student-t sampler does the same for a row whose
# case weight leaves zero through $setWeights.

set.seed(20261005L)
n <- 120L
x <- matrix(runif(n * 2L), n, 2L, dimnames = list(NULL, c("x1", "x2")))
y.binary <- as.double(x[, 1L] + rnorm(n, 0, 0.3) > 0.5)
y.ordinal <- as.double(1L + (seq_len(n) %% 3L))
y.counts <- as.double(rpois(n, exp(1 + x[, 1L])))
nu <- 4
y.heavy <- 2 * x[, 1L] + 0.5 * rt(n, nu)
w.positive <- 0.5 + runif(n)
# rows 1, 5, 9, ... start out; the second pattern brings every other one of
# those back in and takes rows 2, 6, 10, ... out
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
byMask <- function(sampler, rows) sampler$setActiveRows(rows)
# zeros in the case weights, the rows that are in at their own weights
byWeights <- function(sampler, rows) {
  sampler$setWeights(if (is.null(rows)) w.positive else rows * w.positive)
}
student <- function() {
  dbarts::dbarts(
    x,
    y.heavy,
    weights = w.positive,
    family = dbarts:::student(df = nu),
    control = control
  )
}
arms <- list(
  probit = list(
    make = function() {
      dbarts::dbarts(x, y.binary, family = "probit", control = control)
    },
    install = byMask
  ),
  ordinal = list(
    make = function() {
      dbarts::dbarts(x, y.ordinal, family = "ordinal", control = control)
    },
    install = byMask
  ),
  logistic = list(
    make = function() {
      dbarts::dbarts(x, y.binary, family = "logistic", control = control)
    },
    install = byMask
  ),
  nbinom = list(
    make = function() {
      dbarts::dbarts(x, y.counts, family = "nbinom", control = control)
    },
    install = byMask
  ),
  student = list(make = student, install = byMask),
  student.weights = list(make = student, install = byWeights)
)

for (arm in names(arms)) {
  install <- arms[[arm]]$install
  sampler <- arms[[arm]]$make()
  twin <- arms[[arm]]$make()
  invisible(sampler$run(20L, 1L))
  invisible(twin$run(20L, 1L))
  install(sampler, first)
  install(twin, first)
  invisible(sampler$run(10L, 1L))
  invisible(twin$run(10L, 1L))

  # the pattern in force, installed again: no latent moves and no variate is
  # drawn, so the twin that skips the call draws the same sweeps
  held <- sampler$getLatents()
  install(sampler, first)
  expect_identical(sampler$getLatents(), held, info = arm)
  expect_identical(
    sampler$run(0L, 3L)$train,
    twin$run(0L, 3L)$train,
    info = arm
  )

  # a pattern that brings rows back: exactly those rows' latents are redrawn
  held <- sampler$getLatents()
  install(sampler, second)
  moved <- sampler$getLatents() != held
  expect_true(all(moved[reactivated]), info = arm)
  expect_false(any(moved[!reactivated]), info = arm)
  expect_true(all(is.finite(sampler$getLatents())), info = arm)

  # lifting it brings back every row still out
  held <- sampler$getLatents()
  install(sampler, NULL)
  moved <- sampler$getLatents() != held
  expect_true(all(moved[second == 0]), info = arm)
  expect_false(any(moved[second == 1]), info = arm)
  expect_true(all(is.finite(sampler$run(0L, 3L)$train)), info = arm)
}

# 0/1 case weights on probit are the mask, so $setWeights redraws the same
# latents $setActiveRows does
weighted <- arms$probit$make()
masked <- arms$probit$make()
for (rows in list(first, second, rep(1, n))) {
  invisible(weighted$run(5L, 1L))
  invisible(masked$run(5L, 1L))
  weighted$setWeights(rows)
  masked$setActiveRows(rows)
  expect_identical(weighted$getLatents(), masked$getLatents())
}
expect_identical(weighted$run(0L, 3L)$train, masked$run(0L, 3L)$train)

# The redrawn Student-t scale is a draw from its conditional at the fit and
# residual scale in force when the weight comes back: lambda_i is gamma with
# shape (nu + 1) / 2 and rate (nu + w_i (y_i - f_i)^2 / sigma^2) / 2. Forty
# rounds of thirty rows; a scale kept from while the row was out sits at a
# mean probability transform of 0.61.
sampler <- student()
out <- which(first == 0)
transformed <- numeric(0)
for (round in seq_len(40L)) {
  sampler$setWeights(first * w.positive)
  draw <- sampler$run(3L, 1L)
  sampler$setWeights(w.positive)
  rate <- 0.5 *
    (nu +
      w.positive[out] *
        (y.heavy[out] - draw$train[out, 1L])^2 /
        draw$sigma[1L]^2)
  transformed <- c(
    transformed,
    pgamma(sampler$getLatents()[out], shape = 0.5 * (nu + 1), rate = rate)
  )
  invisible(sampler$run(2L, 1L))
}
expect_true(abs(mean(transformed) - 0.5) < 0.03)
expect_true(ks.test(transformed, "punif")$p.value > 0.01)

rm(
  n,
  x,
  y.binary,
  y.ordinal,
  y.counts,
  nu,
  y.heavy,
  w.positive,
  first,
  second,
  reactivated,
  control,
  byMask,
  byWeights,
  student,
  arms,
  arm,
  install,
  sampler,
  twin,
  held,
  moved,
  weighted,
  masked,
  rows,
  out,
  transformed,
  round,
  draw,
  rate
)
