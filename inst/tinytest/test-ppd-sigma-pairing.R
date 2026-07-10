# sampleFromPPD's gaussian noise pairs each posterior draw with its own
# sigma draw. For a multi-chain fit, the noise is drawn in object$sigma's
# own chain-fastest order and, when the caller wants combined output,
# reshaped (not redrawn) with the package's own combineChains() - the same
# helper the stored draws themselves go through - so a combined and a split
# ppd draw from the same seed agree bit-for-bit after accounting for row
# order. See docs/plans/ppd-sigma-pairing.md.

# 1. layout invariance: a combined ppd draw equals the split ppd draw
# reshaped by combineChains(), same seed, same fit
set.seed(1, sample.kind = "Rejection")
n <- 80L
x <- matrix(runif(n * 2L), n, 2L)
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.5)

fit.mc <- bart2(
  y ~ x,
  n.samples = 20L,
  n.burn = 15L,
  n.trees = 20L,
  n.chains = 3L,
  n.threads = 1L,
  verbose = FALSE
)

set.seed(0L)
ppd.split <- extract(fit.mc, type = "ppd", combineChains = FALSE)
set.seed(0L)
ppd.comb <- extract(fit.mc, type = "ppd")

expect_equal(dim(ppd.split), c(3L, 20L, n))
expect_equal(dim(ppd.comb), c(60L, n))
expect_identical(dbarts:::combineChains(ppd.split), ppd.comb)

rm(fit.mc, ppd.split, ppd.comb, x, y, n)

# 2. within-draw coupling: a short, thinly-informed multi-chain fit whose
# sigma varies noticeably draw to draw. The combined ppd row for draw s must
# use draw s's own sigma - checked by bit-exact recomputation against
# ev + rnorm(..., sd) drawn in sigma's native order and reshaped to match,
# exactly as sampleFromPPD does. Under the old code (plain recycling of the
# chain-interleaved sigma against chain-blocked ev rows) this comparison
# fails for any chain past the first.
set.seed(2, sample.kind = "Rejection")
n <- 15L
x <- matrix(runif(n * 2L), n, 2L)
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.5)

fit.vol <- bart2(
  y ~ x,
  n.samples = 12L,
  n.burn = 5L,
  n.trees = 5L,
  n.chains = 4L,
  n.threads = 1L,
  verbose = FALSE
)

set.seed(4L)
ppd.vol <- extract(fit.vol, type = "ppd")
ev.vol <- extract(fit.vol, type = "ev")

n.chains <- 4L
n.samples <- 12L
set.seed(4L)
noise.vol <- rnorm(
  n * n.chains * n.samples,
  0,
  rep_len(fit.vol$sigma, n * n.chains * n.samples)
)
noise.vol <- dbarts:::combineChains(array(noise.vol, c(n.chains, n.samples, n)))
expect_identical(ev.vol + noise.vol, ppd.vol)

rm(fit.vol, ppd.vol, ev.vol, noise.vol, n.chains, n.samples, x, y, n)

# 3. single chain: unaffected by the fix, matches ev + rnorm(sigma) exactly
# (no chain dimension to reorder)
set.seed(3, sample.kind = "Rejection")
n <- 40L
x <- matrix(runif(n * 2L), n, 2L)
y <- 2 * x[, 1L] - x[, 2L] + rnorm(n, 0, 0.5)

fit.sc <- bart2(
  y ~ x,
  n.samples = 25L,
  n.burn = 15L,
  n.trees = 20L,
  n.chains = 1L,
  n.threads = 1L,
  verbose = FALSE
)

set.seed(5L)
ppd.sc <- extract(fit.sc, type = "ppd")
ev.sc <- extract(fit.sc, type = "ev")

set.seed(5L)
expect_identical(
  ev.sc + rnorm(length(ev.sc), 0, rep_len(fit.sc$sigma, length(ev.sc))),
  ppd.sc
)

rm(fit.sc, ppd.sc, ev.sc, x, y, n)

# 4. weighted gaussian, multi-chain combined: the noise sd is sigma / sqrt(w)
# with weights aligned per observation, on top of the same reshape fix
set.seed(4, sample.kind = "Rejection")
n <- 60L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n, 0, 0.4)
w <- rep_len(c(1, 4, 9), n)

fit.w <- bart2(
  y ~ x,
  weights = w,
  n.samples = 15L,
  n.burn = 10L,
  n.trees = 15L,
  n.chains = 3L,
  n.threads = 1L,
  verbose = FALSE
)

set.seed(6L)
ppd.w <- extract(fit.w, type = "ppd")
ev.w <- extract(fit.w, type = "ev")

n.chains <- 3L
n.samples <- 15L
set.seed(6L)
sd.w <- rep_len(fit.w$sigma, n * n.chains * n.samples) *
  rep(sqrt(1 / w), each = n.chains * n.samples)
noise.w <- rnorm(n * n.chains * n.samples, 0, sd.w)
noise.w <- dbarts:::combineChains(array(noise.w, c(n.chains, n.samples, n)))
expect_identical(ev.w + noise.w, ppd.w)

rm(fit.w, ppd.w, ev.w, sd.w, noise.w, n.chains, n.samples, x, y, w, n)

# 5. rbart_vi: multi-chain ppd shape and the same layout invariance
set.seed(5, sample.kind = "Rejection")
n <- 60L
x <- matrix(runif(n * 2L), n, 2L)
g <- factor(rep_len(1:4, n))
y <- x[, 1L] + rnorm(4L, 0, 1)[as.integer(g)] + rnorm(n, 0, 0.4)

fit.r <- rbart_vi(
  y ~ x,
  group.by = g,
  n.samples = 10L,
  n.burn = 5L,
  n.thin = 1L,
  n.chains = 2L,
  n.trees = 15L,
  n.threads = 1L,
  verbose = FALSE
)

set.seed(0L)
ppd.r.split <- extract(fit.r, type = "ppd", combineChains = FALSE)
set.seed(0L)
ppd.r.comb <- extract(fit.r, type = "ppd")

expect_equal(dim(ppd.r.split), c(2L, 10L, n))
expect_equal(dim(ppd.r.comb), c(20L, n))
expect_identical(dbarts:::combineChains(ppd.r.split), ppd.r.comb)

rm(fit.r, ppd.r.split, ppd.r.comb, x, y, g, n)
