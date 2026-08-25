# AFT log-normal survival family on the bartcore engine (src/bartcore/,
# docs/design/survival.md). Exercised through the internal bartcore surface,
# with the per-observation status on the control's bartcore.survival
# attribute (as the public survival surface sets it). The exact-posterior
# gate lives in benchmarks/R/aft-exact.R.

source(
  system.file("common", "bartcoreHandle.R", package = "dbarts"),
  local = TRUE
)

set.seed(21L)
n <- 200L
p <- 3L
x <- matrix(runif(n * p), n, p)
f <- 1.5 * sin(pi * x[, 1L]) + x[, 2L] - 0.5 * x[, 3L]
sigma.true <- 0.5
log.t <- f + sigma.true * rnorm(n)

# a seeded control makes each chain's Mersenne twister deterministic, so the
# reduction below can compare bitwise
control <- dbartsControl(
  n.chains = 1L,
  n.threads = 1L,
  n.trees = 50L,
  updateState = FALSE,
  seed = 271L
)

aftSampler <- function(y, status, weights = NULL, groups = NULL) {
  sampler <- dbarts(x, y, weights = weights, control = control)
  ctrl <- sampler$control
  attr(ctrl, "bartcore.survival") <- as.numeric(status)
  if (!is.null(groups)) {
    attr(ctrl, "bartcore.groups") <- list(
      indices = as.integer(groups),
      n.groups = nlevels(factor(groups)),
      prior = "cauchy",
      rel.scale = sd(y),
      n.steps = 1L
    )
  }
  sampler$control <- ctrl
  dbarts:::bartcoreSampler(sampler, family = "aft")
}

# ---- reduction: all-uncensored aft == gaussian on log T, bitwise ----

samp.g <- dbarts(x, log.t, control = control)
bc.g <- dbarts:::bartcoreSampler(samp.g)
res.g <- bartcoreRun(bc.g, 100L, 200L)

bc.a <- aftSampler(log.t, rep(1, n)) # every observation an event
res.a <- bartcoreRun(bc.a, 100L, 200L)

expect_identical(res.g$train, res.a$train)
expect_identical(res.g$sigma, res.a$sigma)

# ---- recovery under censoring, and the naive downward bias corrected ----

recover <- function(censor.rate) {
  set.seed(11L)
  # censor by an independent time chosen to hit the target rate in expectation
  cens.time <- f +
    quantile(sigma.true * rnorm(2000L), 1 - censor.rate) +
    sigma.true * rnorm(n)
  status <- as.numeric(log.t <= cens.time)
  obs.log.t <- ifelse(status == 1, log.t, cens.time)

  bc <- aftSampler(obs.log.t, status)
  res <- bartcoreRun(bc, 200L, 400L)
  fit.aft <- rowMeans(res$train)

  # ignoring the censoring underestimates the mean log-time
  samp.naive <- dbarts(x, obs.log.t, control = control)
  res.naive <- bartcoreRun(
    dbarts:::bartcoreSampler(samp.naive),
    200L,
    400L
  )
  list(
    rate = mean(status == 0),
    cor = cor(fit.aft, f),
    sigma = mean(res$sigma),
    mean.aft = mean(fit.aft),
    mean.naive = mean(rowMeans(res.naive$train)),
    lat = bartcoreGetLatents(bc),
    status = status,
    obs = obs.log.t
  )
}

for (rate in c(0.2, 0.5)) {
  r <- recover(rate)
  # fit still tracks the signal
  expect_true(r$cor > 0.8)
  # sigma stays in a sane band around the truth (loose at test scale)
  expect_true(r$sigma > 0.3 && r$sigma < 0.8)
  # AFT recovers a higher mean log-time than the censoring-ignoring fit,
  # correcting its downward bias (the model extrapolates the censored tail)
  expect_true(r$mean.aft > r$mean.naive)
  # censored latents sit at or above their observed censoring time; events
  # keep their observed log event time exactly
  cens <- r$status == 0
  expect_true(all(r$lat[cens] >= r$obs[cens] - 1e-8))
  expect_equal(r$lat[!cens], r$obs[!cens])
}

# ---- setResponse under censoring redraws the latents, keeps status ----

set.seed(3L)
cens.time <- f + 0.3 + sigma.true * rnorm(n)
status <- as.numeric(log.t <= cens.time)
obs.log.t <- ifelse(status == 1, log.t, cens.time)
bc.mut <- aftSampler(obs.log.t, status)
invisible(bartcoreRun(bc.mut, 100L, 1L))
# shift every log-time up by 1; the fit should move up with it
bartcoreSetResponse(bc.mut, obs.log.t + 1)
res.mut <- bartcoreRun(bc.mut, 20L, 20L)
expect_equal(dim(res.mut$train), c(n, 20L))
lat.mut <- bartcoreGetLatents(bc.mut)
expect_true(all(lat.mut[status == 0] >= obs.log.t[status == 0] + 1 - 1e-8))
expect_equal(lat.mut[status == 1], obs.log.t[status == 1] + 1)

# ---- grouped composition smoke: AFT + random intercepts (riAFTBART) ----

groups <- rep(1:8, length.out = n)
group.shift <- rnorm(8L, 0, 0.5)[groups]
log.t.g <- f + group.shift + sigma.true * rnorm(n)
cens.g <- f + 0.4 + sigma.true * rnorm(n)
status.g <- as.numeric(log.t.g <= cens.g)
obs.g <- ifelse(status.g == 1, log.t.g, cens.g)
bc.grouped <- aftSampler(obs.g, status.g, groups = groups)
res.grouped <- bartcoreRun(bc.grouped, 100L, 100L)
expect_equal(dim(res.grouped$train), c(n, 100L))
expect_true(all(is.finite(res.grouped$train)))
expect_true(!is.null(res.grouped$ranef))
expect_equal(dim(res.grouped$ranef), c(8L, 100L))

# ---- refusals: weights and post-creation setData on an AFT sampler ----

expect_error(
  aftSampler(log.t, rep(1, n), weights = runif(n) + 0.5),
  "weight"
)
expect_error(
  bartcoreSetData(bc.a, samp.g$data),
  "aft"
)

# ---- public surface: Surv / two-column ingestion, predict, and the
# ---- survivalProbabilities generic ----

set.seed(8L)
cens <- f + 0.4 + sigma.true * rnorm(n)
status.s <- as.numeric(log.t <= cens)
time.s <- exp(ifelse(status.s == 1, log.t, cens)) # observed time (not logged)

# two-column (time, status) matrix with family = "aft"
fit.2col <- bart2(
  x,
  cbind(time.s, status.s),
  family = "aft",
  n.trees = 50L,
  n.burn = 100L,
  n.samples = 200L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L,
  keepTrees = TRUE
)
expect_identical(fit.2col[["family"]], "aft")
# aft carries sigma and returns the linear predictor E[log T | x] (no
# probability transform), so it tracks the signal on the log scale
expect_false(is.null(fit.2col[["sigma"]]))
expect_true(cor(fitted(fit.2col), f) > 0.8)

# a Surv-like object (built without importing survival) auto-dispatches to
# aft and gives an identical fit
surv <- structure(
  cbind(time = time.s, status = status.s),
  class = "Surv",
  type = "right"
)
fit.surv <- bart2(
  x,
  surv,
  n.trees = 50L,
  n.burn = 100L,
  n.samples = 200L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L,
  keepTrees = TRUE
)
expect_identical(fit.surv[["family"]], "aft")
expect_equal(fitted(fit.2col), fitted(fit.surv))

# predict returns log-scale linear-predictor draws; median time is exp of it
x.new <- matrix(runif(5L * p), 5L, p)
pr <- predict(fit.2col, x.new)
expect_equal(ncol(pr), 5L)

# survivalProbabilities returns DRAWS (draws x times x observations); the
# posterior-mean curve is monotone decreasing in [0, 1]
times <- c(0.5, 1, 2, 4)
sp <- survivalProbabilities(fit.2col, times, newdata = x.new)
n.draws <- nrow(fit.2col$yhat.train)
expect_equal(dim(sp), c(n.draws, length(times), 5L))
expect_true(all(sp >= 0 & sp <= 1))
sp.mean <- apply(sp, c(2L, 3L), mean)
expect_true(all(apply(sp.mean, 2L, function(curve) all(diff(curve) <= 1e-8))))
# every individual draw's curve is monotone too (each is an exact normal tail)
expect_true(all(sp[, -1L, ] <= sp[, -length(times), ] + 1e-8))

# refusals through the public surface
expect_error(dbarts(x, log.t, family = "aft"), "two-column|Surv")
expect_error(
  bart2(x, cbind(c(-1, time.s[-1]), status.s), family = "aft", n.chains = 1L),
  "positive"
)
# the training-fit path (no newdata) spans the training observations
sp.train <- survivalProbabilities(fit.2col, times)
expect_equal(dim(sp.train), c(n.draws, length(times), n))

# an explicitly conflicting family with a Surv response errors instead of
# silently becoming aft
expect_error(
  bart2(x, surv, family = "gaussian", n.chains = 1L, verbose = FALSE),
  "aft"
)
expect_error(
  dbarts(x, surv, family = "probit"),
  "aft"
)

# a factor status: survival::Surv codes it as multi-state ("mright"), which
# is detected with a hint; a data.frame factor status hints the same way
surv.mright <- structure(
  cbind(time = time.s, status = status.s + 1),
  class = "Surv",
  type = "mright"
)
expect_error(dbarts(x, surv.mright, family = "aft"), "factor")
expect_error(
  dbarts(
    x,
    data.frame(time = time.s, status = factor(status.s)),
    family = "aft"
  ),
  "factor"
)

# the formula interface is documented out of scope for aft: every call shape
# points at the matrix interface instead of failing hostilely downstream
surv.df <- data.frame(t = time.s, s = status.s, x1 = x[, 1L])
# a plain formula with family = "aft"
expect_error(dbarts(t ~ x1, surv.df, family = "aft"), "matrix interface")
# a Surv-like response through the formula path with family "auto" (guarded
# in dbartsData, before any arithmetic can trip survival's Ops.Surv)
expect_error(
  dbarts(y ~ x1, data = list(y = surv, x1 = x[, 1L])),
  "matrix interface"
)
# as a survreg user would type them, with the real survival package
if (requireNamespace("survival", quietly = TRUE)) {
  expect_error(
    dbarts(survival::Surv(t, s) ~ x1, surv.df, family = "aft"),
    "matrix interface"
  )
  expect_error(
    dbarts(survival::Surv(t, s) ~ x1, surv.df),
    "matrix interface"
  )
  expect_error(
    bart2(survival::Surv(t, s) ~ x1, surv.df, verbose = FALSE),
    "matrix interface"
  )
}

# the rbart method refuses: rbart_vi cannot fit an aft model
rbart.stub <- structure(list(), class = "rbart")
expect_error(survivalProbabilities(rbart.stub, times = 1), "rbart")

# non-aft fits are refused by the bart method
fit.gauss <- bart2(
  x,
  log.t,
  n.trees = 25L,
  n.burn = 50L,
  n.samples = 50L,
  n.chains = 1L,
  verbose = FALSE,
  seed = 7L
)
expect_error(survivalProbabilities(fit.gauss, times = 1), "aft")

# multi-chain conventions: combineChains collapses the chain margin
fit.chains <- bart2(
  x,
  cbind(time.s, status.s),
  family = "aft",
  n.trees = 25L,
  n.burn = 50L,
  n.samples = 50L,
  n.chains = 2L,
  n.threads = 1L,
  verbose = FALSE,
  seed = 7L
)
sp.comb <- survivalProbabilities(fit.chains, times)
expect_equal(dim(sp.comb), c(100L, length(times), n))
sp.unc <- survivalProbabilities(fit.chains, times, combineChains = FALSE)
expect_equal(dim(sp.unc), c(2L, 50L, length(times), n))
# the combined result is the uncombined one with chains stacked sample-major
expect_equal(sp.comb[1:50, , ], sp.unc[1L, , , ])
expect_equal(sp.comb[51:100, , ], sp.unc[2L, , , ])

# ground truth: a fit packaged uncombined (same seed, so identical draws)
# carries yhat.train and sigma with explicit chain margins; the probability
# at any (chain, sample, time, obs) is the exact normal upper tail, which
# pins the sigma-to-draw alignment
fit.chains2 <- bart2(
  x,
  cbind(time.s, status.s),
  family = "aft",
  n.trees = 25L,
  n.burn = 50L,
  n.samples = 50L,
  n.chains = 2L,
  n.threads = 1L,
  combineChains = FALSE,
  verbose = FALSE,
  seed = 7L
)
sp.unc2 <- survivalProbabilities(fit.chains2, times, combineChains = FALSE)
expect_equal(sp.unc, sp.unc2)
expect_equal(
  sp.unc2[2L, 17L, 3L, 5L],
  pnorm(
    (log(times[3L]) - fit.chains2$yhat.train[2L, 17L, 5L]) /
      fit.chains2$sigma[2L, 17L],
    lower.tail = FALSE
  )
)
