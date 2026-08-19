# sampleFromPPD's gaussian noise pairs each posterior draw with its own
# sigma draw. For a multi-chain fit, the noise is drawn in object$sigma's
# own chain-fastest order and, when the caller wants combined output,
# reshaped (not redrawn) with the package's own combineChains() - the same
# helper the stored draws themselves go through - so a combined and a split
# ppd draw from the same seed agree bit-for-bit after accounting for row
# order.

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
# fit.vol$sigma is stored combined (chain-major); normalized to the split
# (n.chains x n.samples) matrix, its as.vector() is chain-fastest - the
# order sampleFromPPD draws noise in before re-combining to match a
# combined ev, exactly as below
sigma.split <- dbarts:::uncombineChains(fit.vol$sigma, n.chains)
set.seed(4L)
noise.vol <- rnorm(
  n * n.chains * n.samples,
  0,
  rep_len(as.vector(sigma.split), n * n.chains * n.samples)
)
noise.vol <- dbarts:::combineChains(array(noise.vol, c(n.chains, n.samples, n)))
expect_identical(ev.vol + noise.vol, ppd.vol)

rm(
  fit.vol,
  ppd.vol,
  ev.vol,
  noise.vol,
  sigma.split,
  n.chains,
  n.samples,
  x,
  y,
  n
)

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
sigma.split.w <- dbarts:::uncombineChains(fit.w$sigma, n.chains)
set.seed(6L)
sd.w <- rep_len(as.vector(sigma.split.w), n * n.chains * n.samples) *
  rep(sqrt(1 / w), each = n.chains * n.samples)
noise.w <- rnorm(n * n.chains * n.samples, 0, sd.w)
noise.w <- dbarts:::combineChains(array(noise.w, c(n.chains, n.samples, n)))
expect_identical(ev.w + noise.w, ppd.w)

rm(
  fit.w,
  ppd.w,
  ev.w,
  sd.w,
  noise.w,
  sigma.split.w,
  n.chains,
  n.samples,
  x,
  y,
  w,
  n
)

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

# 6. binary layout invariance: the binary branches draw rbinom against the
# split-layout probabilities and reshape the outcome (the draw depends on ev,
# so it cannot be reshaped after the fact), so a combined and a split ppd
# draw from the same seed agree bit-for-bit. Under the pre-fix code the
# combined draw consumed the RNG stream in the chain-blocked row order and
# this comparison failed for every chain past the first.
set.seed(6, sample.kind = "Rejection")
n <- 60L
x <- matrix(runif(n * 2L), n, 2L)
y.b <- rbinom(n, 1L, pnorm(0.8 * x[, 1L] - 0.4))

# probit
fit.pb <- bart2(
  y.b ~ x,
  n.samples = 20L,
  n.burn = 20L,
  n.trees = 15L,
  n.chains = 3L,
  n.threads = 1L,
  verbose = FALSE
)
set.seed(0L)
ppd.pb.split <- extract(fit.pb, type = "ppd", combineChains = FALSE)
set.seed(0L)
ppd.pb.comb <- extract(fit.pb, type = "ppd")
expect_equal(dim(ppd.pb.split), c(3L, 20L, n))
expect_equal(dim(ppd.pb.comb), c(60L, n))
expect_identical(dbarts:::combineChains(ppd.pb.split), ppd.pb.comb)

# weighted logistic: draws are binomial(w, p) counts, so also test that the
# reshape composes with the per-observation weight recycling
w <- rep_len(c(1L, 3L, 5L), n)
fit.lg <- bart2(
  y.b ~ x,
  weights = w,
  family = "logistic",
  n.samples = 20L,
  n.burn = 20L,
  n.trees = 15L,
  n.chains = 3L,
  n.threads = 1L,
  verbose = FALSE
)
set.seed(0L)
ppd.lg.split <- extract(fit.lg, type = "ppd", combineChains = FALSE)
set.seed(0L)
ppd.lg.comb <- extract(fit.lg, type = "ppd")
expect_identical(dbarts:::combineChains(ppd.lg.split), ppd.lg.comb)
# counts respect the per-observation weight ceiling in both layouts
expect_true(all(ppd.lg.comb <= rep(w, each = 60L)))

rm(
  fit.pb,
  fit.lg,
  ppd.pb.split,
  ppd.pb.comb,
  ppd.lg.split,
  ppd.lg.comb,
  w,
  x,
  y.b,
  n
)

# 7. chain-major order pin: this is the one assertion in this file that
# checks the combined sigma vector against an independent ground truth
# (a same-seed uncombined fit) rather than against combineChains()/
# uncombineChains() themselves, which would pass under any layout those two
# helpers happen to agree on. Chain 1's whole run must come first.
set.seed(7, sample.kind = "Rejection")
n <- 30L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] + rnorm(n, 0, 0.5)

n.chains <- 3L
n.samples <- 8L
fit.comb <- bart2(
  y ~ x,
  n.samples = n.samples,
  n.burn = 5L,
  n.trees = 8L,
  n.chains = n.chains,
  n.threads = 1L,
  combineChains = TRUE,
  seed = 123L,
  verbose = FALSE
)
fit.split <- bart2(
  y ~ x,
  n.samples = n.samples,
  n.burn = 5L,
  n.trees = 8L,
  n.chains = n.chains,
  n.threads = 1L,
  combineChains = FALSE,
  seed = 123L,
  verbose = FALSE
)
expect_identical(fit.comb$sigma[seq_len(n.samples)], fit.split$sigma[1L, ])

rm(fit.comb, fit.split, x, y, n, n.chains, n.samples)
